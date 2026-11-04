"""End-to-end integration tests for self-consistent MMFT solver, observables, and energetics."""

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
from kitaev_mf.services.observables import (
    compute_plaquette_fluxes,
    compute_sublattice_magnetizations,
)
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def test_system_pure_kitaev_convergence_and_invariants(
    cpu_backend: TorchBackend,
) -> None:
    """Verify solver convergence, vanishing magnetization, W_p=+1, and spectral consistency."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    init_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5)

    config = SolverConfig(
        k_grid_shape=(18, 18),
        max_iterations=15,
        tol_correlator=1e-10,
        tol_gauge=1e-10,
        mixing_alpha=0.8,
    )
    state, metrics = solve_self_consistent_mean_field(cell, params, init_state, config, cpu_backend)
    assert metrics.converged is True
    assert metrics.final_gauge_residual < 1e-12
    assert metrics.final_onsite_residual < 1e-12

    sigma_mag, spin_mag = compute_sublattice_magnetizations(state, cell)
    assert np.allclose(sigma_mag, 0.0, atol=1e-12)
    assert np.allclose(spin_mag, 0.0, atol=1e-12)

    fluxes = compute_plaquette_fluxes(state, cell)
    assert np.allclose(fluxes, 1.0, atol=1e-12)

    k_grid = cell.make_k_grid(18, 18)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)
    residual = check_energy_consistency(state, cell, params, h_batch, cpu_backend, tolerance=1e-10)
    assert residual < 1e-10


def test_system_field_polarized_response_and_energy_consistency(
    cpu_backend: TorchBackend,
) -> None:
    """Verify finite-field polarization alignment and spectral-versus-Wick energy equality."""
    cell = make_honeycomb_2site()
    h_field = (0.03, 0.02, 0.04)
    params = ModelParameters(K=-1.0, J=0.03, h=h_field)
    init_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5)

    config = SolverConfig(
        k_grid_shape=(12, 12),
        max_iterations=45,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        mixing_alpha=0.5,
        lambda_step=0.2,
        diis_history_depth=6,
    )
    state, metrics = solve_self_consistent_mean_field(cell, params, init_state, config, cpu_backend)
    assert metrics.converged is True

    sigma_mag, spin_mag = compute_sublattice_magnetizations(state, cell)
    for s in range(cell.num_sublattices):
        for a in range(3):
            assert sigma_mag[s, a] * h_field[a] > 0.0
            assert np.isclose(spin_mag[s, a], 0.5 * sigma_mag[s, a], atol=1e-14)

    k_grid = cell.make_k_grid(12, 12)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)
    residual = check_energy_consistency(state, cell, params, h_batch, cpu_backend, tolerance=1e-8)
    assert residual < 1e-8


def test_system_2site_vs_4site_scale_parity(
    cpu_backend: TorchBackend,
) -> None:
    """Verify energy per physical site and flux parity match between 2-site and 4-site cells."""
    cell_2 = make_honeycomb_2site()
    cell_4 = make_honeycomb_4site_rectangular()
    params = ModelParameters(K=-1.0)

    config_2 = SolverConfig(
        k_grid_shape=(24, 24),
        max_iterations=10,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        mixing_alpha=0.8,
    )
    config_4 = SolverConfig(
        k_grid_shape=(12, 12),
        max_iterations=10,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        mixing_alpha=0.8,
    )

    init_2 = MeanFieldState.from_pure_kitaev(cell_2, u_gauge=1.0, u_matter=-0.5)
    init_4 = MeanFieldState.from_pure_kitaev(cell_4, u_gauge=1.0, u_matter=-0.5)

    state_2, metrics_2 = solve_self_consistent_mean_field(
        cell_2, params, init_2, config_2, cpu_backend
    )
    state_4, metrics_4 = solve_self_consistent_mean_field(
        cell_4, params, init_4, config_4, cpu_backend
    )
    assert metrics_2.converged is True
    assert metrics_4.converged is True

    energy_2 = compute_gaussian_energy(state_2, cell_2, params)
    energy_4 = compute_gaussian_energy(state_4, cell_4, params)

    assert np.isclose(energy_2.per_site, energy_4.per_site, atol=1e-3)

    flux_2 = compute_plaquette_fluxes(state_2, cell_2)
    flux_4 = compute_plaquette_fluxes(state_4, cell_4)
    assert np.allclose(flux_2, 1.0, atol=1e-12)
    assert np.allclose(flux_4, 1.0, atol=1e-12)

    sigma_2, _ = compute_sublattice_magnetizations(state_2, cell_2)
    sigma_4, _ = compute_sublattice_magnetizations(state_4, cell_4)
    assert np.allclose(sigma_2, 0.0, atol=1e-12)
    assert np.allclose(sigma_4, 0.0, atol=1e-12)
