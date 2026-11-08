"""Unit tests and scientific contracts for Task T3.4.1 (4-Majorana Mean-Field Dynamical Spin Response).

Verifies contracts:
- H_MF_01: 4-Majorana Spin Operator & Exact Sum Rule (sum w_mn + |<sigma>|^2 = 1.0 to < 1e-12)
- H_MF_02: Convolution Bubble Spectrum & Quadrature (Numerical quadrature integral matches discrete sum to < 0.01)
- H_MF_03: Gap Comparison: Mean-Field vs Exact Quench (omega_onset^MF > 2.0 |K| >> Delta E_2v)
- H_MF_04: C3 Isotropy & Site Invariance (Identical spectra across x, y, z and sites j to < 1e-10)
"""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.phase_3.mean_field_spectral import (
    compare_mean_field_vs_exact_quench_thresholds,
    compute_mean_field_spectral_poles,
    evaluate_mean_field_continuous_spectrum,
)
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


@pytest.fixture(scope="module")
def qsl_cluster_and_state() -> tuple[RealSpaceTorusCluster, ModelParameters, MeanFieldState]:
    """Shared L=4 honeycomb torus cluster with isotropic FM Kitaev QSL state."""
    cell = make_honeycomb_2site()
    cluster = RealSpaceTorusCluster(L=4)
    params = ModelParameters(K=-1.0)
    state = make_qsl_seed(cell)
    return cluster, params, state


def test_h_mf_01_sum_rule_analytical_saturation(
    qsl_cluster_and_state: tuple[RealSpaceTorusCluster, ModelParameters, MeanFieldState],
) -> None:
    """H_MF_01: Verify analytical saturation of the 4-Majorana spin sum rule to machine precision."""
    cluster, params, state = qsl_cluster_and_state

    for gamma in (1, 2, 3):
        mf_data = compute_mean_field_spectral_poles(
            cluster, params, state, site_j=0, spin_component=gamma
        )

        # 1. Total spectral weight + ground state expectation squared must equal 1.0
        tot_sum = mf_data.total_spectral_weight + mf_data.ground_state_magnetization
        assert np.isclose(tot_sum, 1.0, atol=1e-12), (
            f"Sum rule violated for gamma={gamma}: total sum = {tot_sum}"
        )

        # 2. In unpolarized QSL, ground state magnetization is zero
        assert mf_data.ground_state_magnetization < 1e-12

        # 3. All poles have non-negative weights
        for p in mf_data.poles:
            assert p.weight >= 0.0
            assert p.energy > 0.0


def test_h_mf_02_convolution_bubble_continuous_quadrature(
    qsl_cluster_and_state: tuple[RealSpaceTorusCluster, ModelParameters, MeanFieldState],
) -> None:
    """H_MF_02: Verify continuous spectral broadening and numerical quadrature integral."""
    cluster, params, state = qsl_cluster_and_state
    mf_data = compute_mean_field_spectral_poles(cluster, params, state, site_j=0, spin_component=3)

    omega_grid = np.linspace(0.0, 10.0, 1000)
    cont_spec = evaluate_mean_field_continuous_spectrum(mf_data, omega_grid=omega_grid, sigma=0.06)

    # Continuous integral int (d omega / 2 pi) S(omega) matches discrete sum to < 0.01
    quadrature_diff = abs(cont_spec.integrated_weight - mf_data.total_spectral_weight)
    assert quadrature_diff < 0.01, f"Quadrature difference {quadrature_diff} exceeds tolerance 0.01"


def test_h_mf_03_gap_demarcation_mf_vs_exact_quench(
    qsl_cluster_and_state: tuple[RealSpaceTorusCluster, ModelParameters, MeanFieldState],
) -> None:
    """H_MF_03: Contrast mean-field two-particle onset against exact non-perturbative vison gap."""
    cluster, params, state = qsl_cluster_and_state
    mf_data = compute_mean_field_spectral_poles(cluster, params, state, site_j=0, spin_component=3)

    exact_vison_threshold = 0.4564  # Two-vison creation onset from Track A

    comp = compare_mean_field_vs_exact_quench_thresholds(
        exact_vison_gap=exact_vison_threshold, mf_data=mf_data
    )

    # 1. Mean-field threshold is shifted up into the two-quasiparticle continuum (> 2.0 |K|)
    assert comp["mean_field_threshold"] > 2.0, (
        f"Expected MF threshold > 2.0, got {comp['mean_field_threshold']}"
    )

    # 2. Substantial positive threshold shift due to absence of non-perturbative flux loops
    assert comp["threshold_shift"] > 1.5
    assert comp["relative_threshold_shift"] > 3.0


def test_h_mf_04_c3_isotropy_and_site_invariance(
    qsl_cluster_and_state: tuple[RealSpaceTorusCluster, ModelParameters, MeanFieldState],
) -> None:
    """H_MF_04: Verify exact C3 spin isotropy and site translational invariance on torus."""
    cluster, params, state = qsl_cluster_and_state

    # 1. C3 spin isotropy: gamma in {1, 2, 3} (x, y, z)
    weights_gamma = []
    min_energies_gamma = []
    for gamma in (1, 2, 3):
        data = compute_mean_field_spectral_poles(
            cluster, params, state, site_j=0, spin_component=gamma
        )
        weights_gamma.append(data.total_spectral_weight)
        min_energies_gamma.append(data.min_excitation_energy)

    assert np.allclose(weights_gamma, 1.0, atol=1e-12)
    assert np.isclose(min_energies_gamma[0], min_energies_gamma[1], atol=1e-10)
    assert np.isclose(min_energies_gamma[1], min_energies_gamma[2], atol=1e-10)

    # 2. Site translational invariance: sites 0, 1, 4, 7
    site_weights = []
    for s in (0, 1, 4, 7):
        data = compute_mean_field_spectral_poles(cluster, params, state, site_j=s, spin_component=3)
        site_weights.append(data.total_spectral_weight)

    assert np.allclose(site_weights, 1.0, atol=1e-12)
