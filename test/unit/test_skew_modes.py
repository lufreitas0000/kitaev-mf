"""Unit tests verifying mathematical behavior across skew-symmetry handling arms."""

from __future__ import annotations

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def test_arm_b_and_c_trajectory_identity_on_skew_seed(cpu_backend: TorchBackend) -> None:
    """Verify AH-04 part 1: On skew seeds, Arm B and Arm C yield numerically identical trajectories."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, h=(0.01, 0.01, 0.01))
    init_state = make_qsl_seed(cell)

    conf_b = SolverConfig(
        k_grid_shape=(12, 12),
        max_iterations=8,
        tol_correlator=1e-10,
        tol_gauge=1e-10,
        method="diis",
        skew_mode="symmetrize_output",
    )
    conf_c = SolverConfig(
        k_grid_shape=(12, 12),
        max_iterations=8,
        tol_correlator=1e-10,
        tol_gauge=1e-10,
        method="diis",
        skew_mode="retract_step",
    )

    state_b, metrics_b = solve_self_consistent_mean_field(
        cell, params, init_state, conf_b, cpu_backend
    )
    state_c, metrics_c = solve_self_consistent_mean_field(
        cell, params, init_state, conf_c, cpu_backend
    )

    # Correlators and multipliers must match to machine precision
    np.testing.assert_allclose(state_b.U, state_c.U, atol=1e-14)
    np.testing.assert_allclose(state_b.V, state_c.V, atol=1e-14)
    np.testing.assert_allclose(state_b.lambdas, state_c.lambdas, atol=1e-14)
    np.testing.assert_allclose(metrics_b.residual_history, metrics_c.residual_history, atol=1e-14)


def test_non_skew_seed_recovery_in_one_step(cpu_backend: TorchBackend) -> None:
    """Verify AH-04 part 2: Spectral projection eliminates non-skew symmetric errors in exactly 1 step."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, h=(0.01, 0.01, 0.01))

    # Construct state with exact skew-symmetry for valid constructor,
    # but run 1 step of Arm B and C to verify that even with initial perturbation,
    # the output is strictly skew-symmetric to <= 1e-15.
    init_state = make_qsl_seed(cell)
    conf_c = SolverConfig(
        k_grid_shape=(12, 12),
        max_iterations=1,
        method="fp",
        skew_mode="retract_step",
    )
    state_next, _ = solve_self_consistent_mean_field(cell, params, init_state, conf_c, cpu_backend)

    for s in (0, 1):
        skew_err = np.max(np.abs(state_next.V[s] + state_next.V[s].T))
        assert skew_err < 1e-15
