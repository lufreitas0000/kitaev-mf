"""Self-consistent mean-field solver and fixed-point algorithm orchestrator."""

from __future__ import annotations

import numpy as np
import torch

from kitaev_mf.domain.lattice import UnitCell
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.ports import LinalgBackendPort
from kitaev_mf.domain.solver import SolverConfig, SolverMetrics
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.correlators import (
    compute_correlators_from_bloch,
    compute_gauge_residuals,
    compute_spectral_projector,
)
from kitaev_mf.services.fixed_point import (
    BroydenSolver,
    DampedPicard,
    FixedPointAccelerator,
    JFNKSolver,
    JointAndersonSolver,
    PlainFixedPoint,
)
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian

__all__ = [
    "JointAndersonSolver",
    "compute_commutator_stationarity",
    "solve_self_consistent_mean_field",
]


def _create_accelerator(config: SolverConfig) -> FixedPointAccelerator:
    """Instantiate accelerator strategy according to configuration."""
    if config.method == "fp":
        return PlainFixedPoint()
    elif config.method == "picard":
        return DampedPicard(alpha=config.mixing_alpha)
    elif config.method == "diis":
        return JointAndersonSolver(
            history_depth=config.diis_history_depth,
            reg=config.diis_reg,
            beta=config.mixing_alpha,
        )
    elif config.method == "broyden":
        return BroydenSolver(alpha=config.broyden_alpha)
    elif config.method == "jfnk":
        return JFNKSolver(
            gmres_tol=config.jfnk_gmres_tol,
            gmres_maxiter=config.jfnk_gmres_maxiter,
            eps_fd=config.jfnk_eps,
        )
    else:
        raise ValueError(f"Unknown solver method: {config.method}")


def compute_commutator_stationarity(
    h_batch: torch.Tensor,
    p_batch: torch.Tensor,
) -> float:
    """Evaluate BZ-averaged Frobenius norm of the commutator ||[H(k), P(k)]||_F.

    Parameters:
        h_batch: Complex Bloch Hamiltonian tensor of shape (N_k, dim, dim).
        p_batch: Complex spectral projector tensor of shape (N_k, dim, dim).

    Returns:
        Scalar average commutator norm over all k-points.
    """
    comm = torch.matmul(h_batch, p_batch) - torch.matmul(p_batch, h_batch)
    frob_per_k = torch.linalg.norm(comm, dim=(-2, -1))
    return float(torch.mean(frob_per_k).item())


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
    accelerator = _create_accelerator(config)

    symmetrize_output = config.skew_mode != "raw"
    retract_step = config.skew_mode != "raw"

    residual_history: list[float] = []
    contraction_factors: list[float] = []
    comm_residuals: list[float] = []
    converged = False
    diff_u = 1.0
    diff_v = 1.0
    max_g = 1.0
    step_count = 0
    total_map_evals = 0

    prev_res_norm: float | None = None

    def evaluate_phi(
        state_in: MeanFieldState,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, torch.Tensor]:
        """Evaluate one pass of the self-consistent map Phi: state -> (U', V', lambda', P_k)."""
        h_batch = build_bloch_hamiltonian(k_grid, cell, params, state_in, backend)
        p_k = compute_spectral_projector(h_batch, backend)
        u_out, v_out = compute_correlators_from_bloch(
            h_batch,
            cell,
            k_grid,
            backend=backend,
            symmetrize_onsite=symmetrize_output,
        )
        g_res = compute_gauge_residuals(v_out)
        lam_out = state_in.lambdas - config.lambda_step * g_res
        return u_out, v_out, lam_out, p_k

    for _ in range(config.max_iterations):
        step_count += 1
        total_map_evals += 1

        # 1. Evaluate self-consistent mapping Phi
        u_new, v_new, lambda_target, p_batch = evaluate_phi(curr_state)

        # 2. Gauge residuals on integrated output
        g_res = compute_gauge_residuals(v_new)

        # 3. Commutator stationarity: commute P_batch with Hamiltonian rebuilt from (u_new, v_new)
        state_integrated = MeanFieldState(U=u_new, V=v_new, lambdas=lambda_target)
        h_rebuilt = build_bloch_hamiltonian(k_grid, cell, params, state_integrated, backend)
        eps_comm = compute_commutator_stationarity(h_rebuilt, p_batch)
        comm_residuals.append(eps_comm)

        # 4. Residual norm evaluation
        diff_u = float(np.max(np.abs(u_new - curr_state.U)))
        diff_v = float(np.max(np.abs(v_new - curr_state.V)))
        max_g = float(np.max(np.abs(g_res)))
        total_res = max(diff_u, diff_v, max_g)
        residual_history.append(total_res)

        # Contraction factor mu_k = ||R_{k+1}|| / ||R_k||
        if prev_res_norm is not None and prev_res_norm > 1e-15:
            contraction_factors.append(total_res / prev_res_norm)
        prev_res_norm = total_res

        # 5. Convergence check
        if (
            diff_u < config.tol_correlator
            and diff_v < config.tol_correlator
            and max_g < config.tol_gauge
        ):
            converged = True
            curr_state = state_integrated
            break

        # 6. Assemble state vectors
        x_curr = np.concatenate(
            [curr_state.U.ravel(), curr_state.V.ravel(), curr_state.lambdas.ravel()]
        )
        phi_vec = np.concatenate([u_new.ravel(), v_new.ravel(), lambda_target.ravel()])

        # 7. Residual callback for matrix-free directional derivatives (used by JFNK)
        def residual_callback(
            x_trial: np.ndarray,
            state_ref: MeanFieldState = curr_state,
        ) -> np.ndarray:
            len_u = state_ref.U.size
            len_v = state_ref.V.size
            u_t = x_trial[:len_u].reshape(state_ref.U.shape)
            v_t = x_trial[len_u : len_u + len_v].reshape(state_ref.V.shape)
            if retract_step:
                v_t = 0.5 * (v_t - np.swapaxes(v_t, -1, -2))
            lam_t = x_trial[len_u + len_v :].reshape(state_ref.lambdas.shape)
            state_t = MeanFieldState(U=u_t, V=v_t, lambdas=lam_t)
            u_p, v_p, lam_p, _ = evaluate_phi(state_t)
            phi_t = np.concatenate([u_p.ravel(), v_p.ravel(), lam_p.ravel()])
            diff_vec = phi_t - x_trial
            return np.asarray(diff_vec, dtype=np.float64)

        # 8. Acceleration proposal
        x_next, extra_evals = accelerator.step(x_curr, phi_vec, residual_fn=residual_callback)
        total_map_evals += extra_evals

        # 9. Unpack into new MeanFieldState
        len_u = curr_state.U.size
        len_v = curr_state.V.size

        u_next = x_next[:len_u].reshape(curr_state.U.shape)
        v_raw = x_next[len_u : len_u + len_v].reshape(curr_state.V.shape)
        if retract_step:
            v_next = 0.5 * (v_raw - np.swapaxes(v_raw, -1, -2))
        else:
            v_next = v_raw
        lambda_next = x_next[len_u + len_v :].reshape(curr_state.lambdas.shape)

        curr_state = MeanFieldState(U=u_next, V=v_next, lambdas=lambda_next)

    metrics = SolverMetrics(
        converged=converged,
        iterations=step_count,
        final_correlator_residual=diff_u,
        final_onsite_residual=diff_v,
        final_gauge_residual=max_g,
        residual_history=tuple(residual_history),
        contraction_factors=tuple(contraction_factors),
        comm_residuals=tuple(comm_residuals),
        map_evaluations=total_map_evals,
    )

    return curr_state, metrics
