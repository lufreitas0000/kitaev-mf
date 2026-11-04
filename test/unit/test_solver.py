"""Test Joint Anderson Acceleration (DIIS) and self-consistent mean-field solver."""

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.solver import (
    JointAndersonSolver,
    solve_self_consistent_mean_field,
)


def test_joint_anderson_solver_linear_system() -> None:
    """Verify JointAndersonSolver achieves rapid residual norm reduction on a coupled saddle system."""
    dim = 8
    rng = np.random.default_rng(42)

    k_primal = rng.standard_normal((4, 4))
    k_primal = k_primal.T @ k_primal + np.eye(4)
    k_dual = rng.standard_normal((4, 4))
    k_dual = k_dual.T @ k_dual + np.eye(4)
    c_block = 0.3 * rng.standard_normal((4, 4))

    mat_saddle = np.block([[k_primal, c_block], [-c_block.T, k_dual]])
    step_size = 0.8 / float(np.linalg.norm(mat_saddle, 2))
    x_target = rng.standard_normal(dim)
    b_vec = -mat_saddle @ x_target

    def f_map(x_vec: np.ndarray) -> np.ndarray:
        return np.asarray(x_vec - step_size * (mat_saddle @ x_vec + b_vec), dtype=np.float64)

    solver = JointAndersonSolver(history_depth=4, reg=1e-12)
    x_curr = np.zeros(dim)

    for _ in range(60):
        f_val = f_map(x_curr)
        residual_norm = np.linalg.norm(f_val - x_curr)
        if residual_norm < 1e-10:
            break
        x_curr = solver.step(x_curr, f_val)

    assert np.linalg.norm(x_curr - x_target) < 1e-8


def test_pure_kitaev_one_step_convergence(cpu_backend: TorchBackend) -> None:
    """Verify that pure Kitaev state converges within 3 iterations with zero gauge violation."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    # Seed close to analytical fixed point
    init_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5248645)
    config = SolverConfig(
        k_grid_shape=(24, 24),
        max_iterations=10,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        mixing_alpha=0.8,
    )

    final_state, metrics = solve_self_consistent_mean_field(
        cell, params, init_state, config, cpu_backend
    )

    assert metrics.converged is True
    assert metrics.iterations <= 3
    assert metrics.final_correlator_residual < 1e-8
    assert metrics.final_onsite_residual < 1e-12
    assert metrics.final_gauge_residual < 1e-12
    assert np.allclose(final_state.lambdas, 0.0, atol=1e-12)


def test_finite_perturbation_convergence(cpu_backend: TorchBackend) -> None:
    """Verify convergence under finite Heisenberg J and Zeeman field with dynamic multiplier shift."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, J=0.04, h=(0.02, 0.01, 0.03))
    init_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5)
    config = SolverConfig(
        k_grid_shape=(16, 16),
        max_iterations=60,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        mixing_alpha=0.5,
        lambda_step=0.2,
        diis_history_depth=6,
    )

    final_state, metrics = solve_self_consistent_mean_field(
        cell, params, init_state, config, cpu_backend
    )

    assert metrics.converged is True
    assert metrics.final_correlator_residual < 1e-8
    assert metrics.final_onsite_residual < 1e-8
    assert metrics.final_gauge_residual < 1e-8
    assert np.linalg.norm(final_state.lambdas) > 1e-4
