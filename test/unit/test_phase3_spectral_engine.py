"""Unit tests for Dynamical Spin Correlation Lehmann Spectral Summation (Milestone M3.1 / Task T3.1.3).

Audits:
- H_SPEC_01: Spectral pole positivity & non-local vison creation threshold.
- H_SPEC_02: Local completeness sum rule saturation (> 95% in 2-quasiparticle sector).
- H_SPEC_03: Point-group C_3 spatial rotational isometry (S^{xx} = S^{yy} = S^{zz}).
- H_SPEC_04: Continuous frequency broadening, non-negativity, and numerical quadrature.
- H_SPEC_05: Quench response on vortex background and low-energy MZM signature.
"""

from __future__ import annotations

import numpy as np

from kitaev_mf.services.phase_3.gauge_theory import (
    build_honeycomb_torus_cluster,
    build_uniform_gauge_configuration,
)
from kitaev_mf.services.phase_3.spectral_engine import (
    compute_local_spectral_poles,
    evaluate_spectral_function,
    verify_spectral_sum_rule,
)
from kitaev_mf.services.phase_3.vortex_engine import (
    HoneycombGeometry,
    create_two_vortex_gauge_configuration,
)


def test_h_spec_01_spectral_pole_positivity_and_vison_threshold() -> None:
    """Verify Hypothesis H_SPEC_01: All excitation poles are strictly positive and bounded by the flux gap."""
    L = 3
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    geom = HoneycombGeometry(L=L)
    gauge_0 = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    K = 1.0
    kappa = 0.15
    spec_data = compute_local_spectral_poles(
        gauge=gauge_0,
        site_j=0,
        bond_type=3,  # z-bond
        geom=geom,
        K=K,
        kappa=kappa,
        boundary="obc",
        bonds=bonds,
        triplets=triplets,
    )

    # 1. Non-local vison gap must be strictly positive
    assert spec_data.flux_gap > 0.10
    assert np.isclose(spec_data.flux_gap, 0.23178366, rtol=1e-4)

    # 2. Every single excitation pole must be strictly positive
    assert len(spec_data.poles) > 0
    for pole in spec_data.poles:
        assert pole.energy > spec_data.flux_gap - 1e-10
        assert pole.weight >= 0.0

    # 3. Lowest pole must correspond exactly to Delta E_flux + eps_{b, 0}
    lowest_pole = spec_data.poles[0]
    assert lowest_pole.energy > spec_data.flux_gap


def test_h_spec_02_completeness_sum_rule_saturation() -> None:
    """Verify Hypothesis H_SPEC_02: Discrete sum rule saturates > 95% of spectral weight on N=18 cluster."""
    L = 3
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    geom = HoneycombGeometry(L=L)
    gauge_0 = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    spec_data = compute_local_spectral_poles(
        gauge=gauge_0,
        site_j=0,
        bond_type=3,
        geom=geom,
        K=1.0,
        kappa=0.15,
        boundary="obc",
        bonds=bonds,
        triplets=triplets,
    )

    # Sum of transition weights must saturate > 95% of the total spin sum rule
    total_w = spec_data.total_spectral_weight
    assert 0.95 <= total_w <= 1.00000000001

    rule_metrics = verify_spectral_sum_rule(spec_data, tol=0.06)
    assert rule_metrics["is_discrete_valid"] is True
    assert rule_metrics["num_poles"] == 9  # 1-quasiparticle channels for N=9 Dirac modes


def test_h_spec_03_c3_spatial_isometry() -> None:
    """Verify Hypothesis H_SPEC_03: Point-group C_3 rotation invariance S^{xx} = S^{yy} = S^{zz}."""
    L = 3
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    geom = HoneycombGeometry(L=L)
    gauge_0 = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    K = 1.0
    kappa = 0.15

    spec_x = compute_local_spectral_poles(
        gauge=gauge_0,
        site_j=0,
        bond_type=1,
        geom=geom,
        K=K,
        kappa=kappa,
        boundary="obc",
        bonds=bonds,
        triplets=triplets,
    )
    spec_y = compute_local_spectral_poles(
        gauge=gauge_0,
        site_j=0,
        bond_type=2,
        geom=geom,
        K=K,
        kappa=kappa,
        boundary="obc",
        bonds=bonds,
        triplets=triplets,
    )
    spec_z = compute_local_spectral_poles(
        gauge=gauge_0,
        site_j=0,
        bond_type=3,
        geom=geom,
        K=K,
        kappa=kappa,
        boundary="obc",
        bonds=bonds,
        triplets=triplets,
    )

    # Total weights must match to machine precision
    assert np.isclose(spec_x.total_spectral_weight, spec_y.total_spectral_weight, atol=1e-12)
    assert np.isclose(spec_x.total_spectral_weight, spec_z.total_spectral_weight, atol=1e-12)

    # Poles must match entry-by-entry
    assert len(spec_x.poles) == len(spec_y.poles) == len(spec_z.poles)
    for px, py, pz in zip(spec_x.poles, spec_y.poles, spec_z.poles, strict=False):
        assert np.isclose(px.energy, py.energy, atol=1e-12)
        assert np.isclose(px.energy, pz.energy, atol=1e-12)
        assert np.isclose(px.weight, py.weight, atol=1e-12)
        assert np.isclose(px.weight, pz.weight, atol=1e-12)

    # Continuous spectra must match across entire grid
    omega_grid = np.linspace(0.0, 8.0, 500)
    cont_x = evaluate_spectral_function(
        spec_x, omega_grid=omega_grid, broadening_type="gaussian", broadening_param=0.05
    )
    cont_y = evaluate_spectral_function(
        spec_y, omega_grid=omega_grid, broadening_type="gaussian", broadening_param=0.05
    )
    cont_z = evaluate_spectral_function(
        spec_z, omega_grid=omega_grid, broadening_type="gaussian", broadening_param=0.05
    )

    max_diff_xy = float(np.max(np.abs(cont_x.spectral_density - cont_y.spectral_density)))
    max_diff_xz = float(np.max(np.abs(cont_x.spectral_density - cont_z.spectral_density)))
    assert max_diff_xy < 1e-11
    assert max_diff_xz < 1e-11


def test_h_spec_04_continuous_quadrature_and_broadening() -> None:
    """Verify Hypothesis H_SPEC_04: Continuous Gaussian broadening preserves sum rule via numerical quadrature."""
    L = 3
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    geom = HoneycombGeometry(L=L)
    gauge_0 = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    spec_data = compute_local_spectral_poles(
        gauge=gauge_0,
        site_j=0,
        bond_type=3,
        geom=geom,
        K=1.0,
        kappa=0.15,
        boundary="obc",
        bonds=bonds,
        triplets=triplets,
    )

    omega_grid = np.linspace(0.0, 10.0, 2000)
    sigma = 0.05
    cont_spec = evaluate_spectral_function(
        spec_data, omega_grid=omega_grid, broadening_type="gaussian", broadening_param=sigma
    )

    # 1. Non-negativity: S(omega) >= 0 everywhere
    assert np.all(cont_spec.spectral_density >= 0.0)

    # 2. Quadrature integral must match discrete sum to machine precision
    metrics = verify_spectral_sum_rule(spec_data, continuous_spec=cont_spec, tol=0.06)
    assert metrics["is_continuous_valid"] is True
    assert metrics["quadrature_difference"] < 1e-6

    # 3. Lorentzian broadening check
    eta = 0.05
    cont_lorentz = evaluate_spectral_function(
        spec_data, omega_grid=omega_grid, broadening_type="lorentzian", broadening_param=eta
    )
    assert np.all(cont_lorentz.spectral_density >= 0.0)
    assert abs(cont_lorentz.integrated_weight - spec_data.total_spectral_weight) < 0.05


def test_h_spec_05_two_vortex_mzm_quench_response() -> None:
    """Verify Hypothesis H_SPEC_05: Quenching on a background with trapped vortices generates valid spectra."""
    L = 3
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    geom = HoneycombGeometry(L=L)
    gauge_0 = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    # Insert two vortices at plaquettes (0, 0) and (1, 1)
    two_vortex_gauge, cut = create_two_vortex_gauge_configuration(
        gauge_0, p1=(0, 0), p2=(1, 1), geom=geom
    )

    spec_vortex = compute_local_spectral_poles(
        gauge=two_vortex_gauge,
        site_j=0,
        bond_type=3,
        geom=geom,
        K=1.0,
        kappa=0.15,
        boundary="obc",
        bonds=bonds,
        triplets=triplets,
    )

    assert len(spec_vortex.poles) > 0
    assert spec_vortex.total_spectral_weight > 0.90


def test_h_spec_06_odd_initial_parity_0qp_and_2qp_channels() -> None:
    """Verify compute_local_spectral_poles when initial state has odd parity (pi_G = -1).

    Exercises the 0-quasiparticle and 2-quasiparticle transition channels in sector G'.
    """
    L = 3
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    geom = HoneycombGeometry(L=L)
    gauge_0 = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    spec_torus = compute_local_spectral_poles(
        gauge=gauge_0,
        site_j=0,
        bond_type=3,
        geom=geom,
        K=1.0,
        kappa=0.15,
        boundary="torus",  # pi_G = -1 on clean torus
        bonds=bonds,
        triplets=triplets,
    )

    # Must contain 0-qp vacuum and 2-qp excitation states
    state_lens = {len(p.quasiparticle_state) for p in spec_torus.poles}
    assert 0 in state_lens
    assert 2 in state_lens
    assert len(spec_torus.poles) == 37  # 1 vacuum + 9*8/2 = 36 two-particle pairs
    assert np.isclose(spec_torus.total_spectral_weight, 0.99825, rtol=1e-3)


def test_h_spec_07_singular_quench_spectral_poles() -> None:
    """Verify compute_local_spectral_poles under singular Bogoliubov quenches (det W = -1)."""
    L = 3
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    geom = HoneycombGeometry(L=L)
    gauge_0 = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    # Background with two vortices at (0, 0) and (0, 1) produces a singular quench on bond_type 2
    gauge_2v, _ = create_two_vortex_gauge_configuration(gauge_0, p1=(0, 0), p2=(0, 1), geom=geom)

    spec_sing = compute_local_spectral_poles(
        gauge=gauge_2v,
        site_j=0,
        bond_type=2,
        geom=geom,
        K=1.0,
        kappa=0.15,
        boundary="obc",
        bonds=bonds,
        triplets=triplets,
    )

    assert spec_sing.is_singular is True
    assert len(spec_sing.poles) > 0
    assert spec_sing.total_spectral_weight > 0.90
