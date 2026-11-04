"""TDD contract test suite for production energetics and spectral consistency (M1.4.4)."""

import numpy as np

from kitaev_mf.domain.lattice import (
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.energetics import (
    check_energy_consistency,
    compute_gaussian_energy,
)
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def test_energy_breakdown_pure_kitaev_normalization() -> None:
    """Verify exact algebraic normalization of pure-Kitaev ground-state energy."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    u_matter_exact = -0.5248657458506
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=u_matter_exact)

    breakdown = compute_gaussian_energy(state, cell, params)
    assert np.isclose(breakdown.per_bond, -0.5248657458506, atol=1e-12)
    assert np.isclose(breakdown.per_cell, -1.5745972375518, atol=1e-12)
    assert np.isclose(breakdown.per_site, -0.7872986187759, atol=1e-12)
    assert np.isclose(breakdown.E_field, 0.0, atol=1e-14)
    assert np.isclose(breakdown.E_multiplier, 0.0, atol=1e-14)


def test_energy_consistency_pure_kitaev(cpu_backend: TorchBackend) -> None:
    """Verify |(E_quad,0 - C) - E_total^G| < 1e-10 in the pure Kitaev ground state."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    init_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5)

    config = SolverConfig(
        k_grid_shape=(24, 24),
        max_iterations=10,
        tol_correlator=1e-10,
        tol_gauge=1e-10,
        mixing_alpha=0.8,
    )
    state, metrics = solve_self_consistent_mean_field(cell, params, init_state, config, cpu_backend)
    assert metrics.converged is True

    k_grid = cell.make_k_grid(24, 24)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)

    residual = check_energy_consistency(state, cell, params, h_batch, cpu_backend, tolerance=1e-10)
    assert residual < 1e-10


def test_energy_consistency_under_perturbations(cpu_backend: TorchBackend) -> None:
    """Verify spectral and Wick energy consistency under finite J, field h, and multipliers."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, J=0.04, h=(0.02, 0.01, 0.03))
    init_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5)

    config = SolverConfig(
        k_grid_shape=(16, 16),
        max_iterations=50,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        mixing_alpha=0.5,
        lambda_step=0.2,
        diis_history_depth=6,
    )
    state, metrics = solve_self_consistent_mean_field(cell, params, init_state, config, cpu_backend)
    assert metrics.converged is True

    k_grid = cell.make_k_grid(16, 16)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)

    residual = check_energy_consistency(state, cell, params, h_batch, cpu_backend, tolerance=1e-8)
    assert residual < 1e-8


def test_supercell_energy_density_folding_parity() -> None:
    """Verify energy per physical site matches between 2-site and 4-site unit cells."""
    cell_2 = make_honeycomb_2site()
    cell_4 = make_honeycomb_4site_rectangular()
    params = ModelParameters(K=-1.0, Gamma=0.08, h=(0.02, 0.02, 0.02))

    state_2 = MeanFieldState.from_pure_kitaev(cell_2, u_gauge=1.0, u_matter=-0.3)
    state_4 = MeanFieldState.from_pure_kitaev(cell_4, u_gauge=1.0, u_matter=-0.3)

    energy_2 = compute_gaussian_energy(state_2, cell_2, params)
    energy_4 = compute_gaussian_energy(state_4, cell_4, params)

    assert np.isclose(energy_2.per_site, energy_4.per_site, atol=1e-12)
