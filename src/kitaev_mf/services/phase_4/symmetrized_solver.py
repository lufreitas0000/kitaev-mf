"""Ultra-low field symmetrized fixed-point solver on the reduced D_sym manifold.

Operates directly on the 7-parameter C3- and sublattice-symmetric manifold:
    x = (lambda, v0, v1, u0, uz, uperp, ucross)
eliminating numerical symmetry breaking and accelerating convergence down
to ultra-low fields h in [1e-6, 1e-2].
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from kitaev_mf.domain.lattice import UnitCell, make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverMetrics
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.fixed_point import JointAndersonSolver
from kitaev_mf.services.phase_4.symmetry_reduction import (
    SymmetricParameters,
    extract_symmetric_parameters,
    pack_symmetric_state,
)
from kitaev_mf.services.correlators import (
    compute_correlators_from_bloch,
    compute_gauge_residuals,
)


@dataclass(frozen=True)
class SymmetrizedSolverResult:
    """Output bundle of the symmetrized fixed-point solve."""

    sym_params: SymmetricParameters
    state: MeanFieldState
    metrics: SolverMetrics
    dirac_gap_mf: float
    dirac_gap_canonical: float
    magnetization: float


def solve_symmetrized_fixed_point(
    params: ModelParameters,
    cell: UnitCell | None = None,
    init_sym: SymmetricParameters | None = None,
    nk: int = 24,
    max_iter: int = 60,
    tol: float = 1e-9,
    diis_depth: int = 6,
    diis_beta: float = 0.5,
    lambda_step: float = 0.2,
) -> SymmetrizedSolverResult:
    """Solve the mean-field fixed point directly in the symmetric D_sym subspace using DIIS acceleration.

    Args:
        params: Model exchange couplings and magnetic field h.
        cell: Honeycomb unit cell (defaults to 2-site primitive cell).
        init_sym: Initial symmetric guess (defaults to pure Kitaev state).
        nk: BZ k-mesh dimension (nk x nk).
        max_iter: Maximum fixed-point iterations.
        tol: Convergence tolerance for parameter updates and gauge residual.
        diis_depth: History depth for Joint Anderson Acceleration.
        diis_beta: Relaxation parameter for Anderson acceleration.
        lambda_step: Step size for KKT gauge Lagrange multiplier update.

    Returns:
        SymmetrizedSolverResult with converged parameters and physical observables.
    """
    if cell is None:
        cell = make_honeycomb_2site()

    backend = TorchBackend()
    k_grid = cell.make_k_grid(nk, nk)

    if init_sym is None:
        init_sym = SymmetricParameters.pure_kitaev()

    current_sym = init_sym
    accelerator = JointAndersonSolver(history_depth=diis_depth, beta=diis_beta, reg=1e-12)
    residual_history: list[float] = []
    converged = False

    for it in range(1, max_iter + 1):
        # 1. Assemble symmetric state
        state = pack_symmetric_state(current_sym, cell)

        # 2. Build Bloch Hamiltonian across the k-grid
        h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, backend)

        # 3. Compute BZ-integrated correlators
        u_new, v_new = compute_correlators_from_bloch(
            h_batch, cell, k_grid, backend=backend, symmetrize_onsite=True
        )
        gauge_res = compute_gauge_residuals(v_new)
        max_gauge_res = float(np.max(np.abs(gauge_res)))

        # 4. Project updated correlators onto symmetric parameters
        temp_state = MeanFieldState(U=u_new, V=v_new, lambdas=state.lambdas)
        extracted_sym = extract_symmetric_parameters(temp_state)

        # 5. Multiplier gradient step
        g_mean = float(np.mean(gauge_res))
        lam_target = current_sym.lambda_val - lambda_step * g_mean

        # 6. DIIS extrapolation on 7D parameter vector
        vec_curr = current_sym.to_vector()
        vec_target = np.array(
            [
                lam_target,
                extracted_sym.v0,
                extracted_sym.v1,
                extracted_sym.u0,
                extracted_sym.uz,
                extracted_sym.uperp,
                extracted_sym.ucross,
            ],
            dtype=np.float64,
        )

        vec_next, _ = accelerator.step(vec_curr, vec_target)
        diff = float(np.max(np.abs(vec_target - vec_curr)))
        residual_history.append(diff)

        current_sym = SymmetricParameters.from_vector(vec_next)

        if diff < tol and max_gauge_res < tol:
            converged = True
            break

    # Final converged state
    final_state = pack_symmetric_state(current_sym, cell)

    # Compute Dirac gap at K = (1/3, 2/3)
    k_k = np.array([[1.0 / 3.0, 2.0 / 3.0]])
    h_bloch = build_bloch_hamiltonian(k_k, cell, params, final_state, backend)
    if hasattr(h_bloch, "detach"):
        h_k = h_bloch.detach().cpu().numpy()[0]
    elif hasattr(h_bloch, "numpy"):
        h_k = h_bloch.numpy()[0]
    else:
        h_k = np.asarray(h_bloch[0])
    evals = np.linalg.eigvalsh(h_k)
    gap_mf = float(evals[4] - evals[3])
    gap_canonical = 4.0 * gap_mf
    magnetization = 2.0 * current_sym.v0

    metrics = SolverMetrics(
        converged=converged,
        iterations=it,
        final_correlator_residual=diff,
        final_onsite_residual=float(abs(current_sym.v0 - extracted_sym.v0)),
        final_gauge_residual=max_gauge_res,
        residual_history=residual_history,
        contraction_factors=[],
        comm_residuals=[],
    )

    return SymmetrizedSolverResult(
        sym_params=current_sym,
        state=final_state,
        metrics=metrics,
        dirac_gap_mf=gap_mf,
        dirac_gap_canonical=gap_canonical,
        magnetization=magnetization,
    )
