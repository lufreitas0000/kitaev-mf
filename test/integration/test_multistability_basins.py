"""Integration tests verifying weak-field multistability basin separation and symmetry."""

from __future__ import annotations

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.observables import (
    compute_plaquette_fluxes,
    compute_sublattice_magnetizations,
)
from kitaev_mf.services.seeds import make_polarized_seed, make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def test_weak_field_basin_separation_and_symmetry(cpu_backend: TorchBackend) -> None:
    """Verify PH-02 and PH-03: QSL and Polarized basins remain separated with exact C3 isotropy."""
    cell = make_honeycomb_2site()
    h_val = 0.02
    h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
    params = ModelParameters(K=-1.0, h=h_vec)

    config = SolverConfig(
        k_grid_shape=(16, 16),
        max_iterations=20,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        method="diis",
    )

    # 1. QSL Basin
    qsl_seed = make_qsl_seed(cell)
    state_qsl, metrics_qsl = solve_self_consistent_mean_field(
        cell, params, qsl_seed, config, cpu_backend
    )
    assert metrics_qsl.converged is True

    # 2. Polarized Basin
    pol_seed = make_polarized_seed(cell, axis="111")
    state_pol, metrics_pol = solve_self_consistent_mean_field(
        cell, params, pol_seed, config, cpu_backend
    )
    assert metrics_pol.converged is True

    # Energy ordering (PH-02)
    e_qsl = compute_gaussian_energy(state_qsl, cell, params).per_site
    e_pol = compute_gaussian_energy(state_pol, cell, params).per_site
    assert e_qsl < e_pol - 0.2, f"Expected QSL ground state, got e_qsl={e_qsl}, e_pol={e_pol}"

    # Wilson loop flux separation (PH-02)
    w_qsl = compute_plaquette_fluxes(state_qsl, cell)[0]
    w_pol = compute_plaquette_fluxes(state_pol, cell)[0]
    assert w_qsl > 0.99
    assert abs(w_pol) < 0.05

    # Magnetization isotropy and sublattice symmetry (PH-03)
    sig_qsl, _ = compute_sublattice_magnetizations(state_qsl, cell)
    np.testing.assert_allclose(sig_qsl[0], sig_qsl[1], atol=1e-12)
    np.testing.assert_allclose(sig_qsl[0, 0], sig_qsl[0, 1], atol=1e-12)
    np.testing.assert_allclose(sig_qsl[0, 1], sig_qsl[0, 2], atol=1e-12)


def test_high_field_classical_convergence(cpu_backend: TorchBackend) -> None:
    """Verify PH-06: at h=5.0, state converges to fully polarized along [111]."""
    cell = make_honeycomb_2site()
    h_val = 5.0
    h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
    params = ModelParameters(K=-1.0, h=h_vec)

    config = SolverConfig(
        k_grid_shape=(16, 16),
        max_iterations=40,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        method="diis",
    )

    pol_seed = make_polarized_seed(cell, axis="111")
    state_pol, metrics_pol = solve_self_consistent_mean_field(
        cell, params, pol_seed, config, cpu_backend
    )
    assert metrics_pol.converged is True

    # Energy per site should approach classical E = -h - 0.5 = -5.5
    e_val = compute_gaussian_energy(state_pol, cell, params).per_site
    e_classical = -h_val - 0.5
    np.testing.assert_allclose(e_val, e_classical, atol=0.01)

    # Wilson loop flux should be zero
    w_val = compute_plaquette_fluxes(state_pol, cell)[0]
    assert abs(w_val) < 1e-4

    # Magnetization should be saturated along [111] with norm ~ 1.0
    sig, _ = compute_sublattice_magnetizations(state_pol, cell)
    m_mag = np.linalg.norm(sig[0])
    np.testing.assert_allclose(m_mag, 1.0, atol=1e-3)
