"""Joint Anderson Acceleration (DIIS) and self-consistent mean-field solver."""

from typing import cast

import numpy as np

from kitaev_mf.domain.lattice import UnitCell
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.ports import LinalgBackendPort
from kitaev_mf.domain.solver import SolverConfig, SolverMetrics
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.correlators import (
    compute_correlators_from_bloch,
    compute_gauge_residuals,
)
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian


class JointAndersonSolver:
    """Joint DIIS solver operating on concatenated state vectors X."""

    def __init__(self, history_depth: int = 6, reg: float = 1e-12) -> None:
        self.history_depth = history_depth
        self.reg = reg
        self.reset()

    def reset(self) -> None:
        """Clear historical subspace buffers."""
        self._history_x: list[np.ndarray] = []
        self._history_r: list[np.ndarray] = []

    def step(self, x: np.ndarray, f_x: np.ndarray) -> np.ndarray:
        """Perform DIIS extrapolation from current state x and fixed-point output f_x."""
        residual = f_x - x

        self._history_x.append(np.asarray(x, dtype=np.float64).copy())
        self._history_r.append(np.asarray(residual, dtype=np.float64).copy())

        if len(self._history_x) > self.history_depth:
            self._history_x.pop(0)
            self._history_r.pop(0)

        m = len(self._history_x)
        if m == 1:
            return cast(np.ndarray, np.asarray(f_x, dtype=np.float64).copy())

        # Build Gram matrix Omega_{pq} = <R_p, R_q> + reg * delta_{pq} * <R_p, R_p>
        omega = np.zeros((m, m), dtype=np.float64)
        for p in range(m):
            for q in range(p, m):
                dot_val = float(np.dot(self._history_r[p], self._history_r[q]))
                omega[p, q] = dot_val
                omega[q, p] = dot_val

        diag_scale = np.diag(omega)
        for p in range(m):
            omega[p, p] += self.reg * max(diag_scale[p], 1e-14)

        # Assemble (m + 1) x (m + 1) KKT system: [[Omega, -1], [-1^T, 0]] [c, nu]^T = [0, -1]^T
        kkt = np.zeros((m + 1, m + 1), dtype=np.float64)
        kkt[:m, :m] = omega
        kkt[:m, m] = -1.0
        kkt[m, :m] = -1.0

        rhs = np.zeros(m + 1, dtype=np.float64)
        rhs[m] = -1.0

        try:
            sol = np.linalg.solve(kkt, rhs)
            weights = sol[:m]
        except np.linalg.LinAlgError:
            # Fallback to last fixed-point iteration if KKT matrix is singular
            return cast(np.ndarray, np.asarray(f_x, dtype=np.float64).copy())

        # Construct accelerated state: X_next = sum_j c_j * (X_j + R_j) = sum_j c_j * f(X_j)
        x_next = np.zeros_like(x, dtype=np.float64)
        for j in range(m):
            x_next += weights[j] * (self._history_x[j] + self._history_r[j])

        return cast(np.ndarray, x_next)


def solve_self_consistent_mean_field(
    cell: UnitCell,
    params: ModelParameters,
    initial_state: MeanFieldState,
    config: SolverConfig,
    backend: LinalgBackendPort | None = None,
) -> tuple[MeanFieldState, SolverMetrics]:
    """Execute the self-consistent mean-field loop until convergence or max_iterations."""
    k_grid = cell.make_k_grid(*config.k_grid_shape)
    curr_state = initial_state
    solver = JointAndersonSolver(
        history_depth=config.diis_history_depth,
        reg=config.diis_reg,
    )

    residual_history: list[float] = []
    converged = False
    diff_u = 1.0
    max_g = 1.0
    step_count = 0

    for _ in range(config.max_iterations):
        step_count += 1
        # 1. Build Bloch Hamiltonian at current mean-field state
        h_batch = build_bloch_hamiltonian(k_grid, cell, params, curr_state, backend)

        # 2. Diagonalize and integrate reciprocal spectral projector over BZ
        u_new, v_new = compute_correlators_from_bloch(h_batch, cell, k_grid, backend)

        # 3. Compute gauge residuals: g_s^gamma = tr(G^gamma V_s)
        g_res = compute_gauge_residuals(v_new)

        # 4. Target multiplier update: dual gradient ascent lambda_target = lambda - alpha_lambda * g
        lambda_target = curr_state.lambdas - config.lambda_step * g_res

        # 5. Convergence check
        diff_u = float(np.max(np.abs(u_new - curr_state.U)))
        max_g = float(np.max(np.abs(g_res)))
        total_res = max(diff_u, max_g)
        residual_history.append(total_res)

        if diff_u < config.tol_correlator and max_g < config.tol_gauge:
            converged = True
            curr_state = MeanFieldState(U=u_new, V=v_new, lambdas=lambda_target)
            break

        # 6. Base Picard mixing target
        u_target = (1.0 - config.mixing_alpha) * curr_state.U + config.mixing_alpha * u_new
        v_target = (1.0 - config.mixing_alpha) * curr_state.V + config.mixing_alpha * v_new

        # 7. Concatenate into joint DIIS state vectors
        x_curr = np.concatenate(
            [curr_state.U.ravel(), curr_state.V.ravel(), curr_state.lambdas.ravel()]
        )
        f_target = np.concatenate([u_target.ravel(), v_target.ravel(), lambda_target.ravel()])

        # 8. DIIS extrapolation
        x_next = solver.step(x_curr, f_target)

        # 9. Unpack into new MeanFieldState
        len_u = curr_state.U.size
        len_v = curr_state.V.size

        u_next = x_next[:len_u].reshape(curr_state.U.shape)
        v_raw = x_next[len_u : len_u + len_v].reshape(curr_state.V.shape)
        v_next = 0.5 * (v_raw - np.swapaxes(v_raw, -1, -2))
        lambda_next = x_next[len_u + len_v :].reshape(curr_state.lambdas.shape)

        curr_state = MeanFieldState(U=u_next, V=v_next, lambdas=lambda_next)

    metrics = SolverMetrics(
        converged=converged,
        iterations=step_count,
        final_correlator_residual=diff_u,
        final_gauge_residual=max_g,
        residual_history=tuple(residual_history),
    )

    return curr_state, metrics
