"""Unit tests for Two-Vortex MZM Spectroscopy, Real-Space Local Density of States, and arXiv:2209.12930 Reproduction (Task T3.2.1).

Audits:
- H_MZM_01: Midgap MZM in-gap isolation and spectral gap demarcation.
- H_MZM_02: Reproduction of Figure 3 (3-tier STS spectra: Core vs Ring vs Bulk).
- H_MZM_03: Reproduction of Figure 4 (Spatial exponential decay and fitted localization length xi_M).
- H_MZM_04: Local completeness sum-rule preservation across defect backgrounds.
"""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_mf.services.phase_3.mzm_spectroscopy import (
    build_two_vortex_mzm_system,
    reproduce_figure_3_spectrogram,
    reproduce_figure_4_spatial_decay,
)


def test_h_mzm_01_midgap_mode_isolation_and_spectral_gap() -> None:
    """Verify Hypothesis H_MZM_01: Two-vortex background traps two isolated mid-gap MZM states."""
    L = 6
    K = 1.0
    kappa = 0.15
    system = build_two_vortex_mzm_system(L=L, p1=(1, 1), p2=(4, 4), K=K, kappa=kappa)

    # 1. Vortex separation must be significant (R >= 4a)
    assert system.vortex_distance >= 4.0

    # 2. In-gap MZM mode energies
    eps_0, eps_1 = system.mzm_energies
    assert eps_0 < 0.45
    assert eps_1 < 0.45
    assert eps_0 > 0.0
    assert eps_1 >= eps_0

    # 3. Bulk band gap and clean isolation ratio
    bulk_gap = system.bulk_gap_estimate  # 6 * sqrt(3) * kappa ~ 1.5588
    assert np.isclose(bulk_gap, 6.0 * np.sqrt(3.0) * kappa, rtol=1e-5)

    eps_2 = system.diag.energies[2]
    assert eps_2 > 1.0  # Bulk continuum onset
    assert system.in_gap_ratio < 0.45  # eps_1 / eps_2 << 1

    # 4. Total MZM density normalization across lattice
    assert np.isclose(float(np.sum(system.mzm_density)), 1.0, atol=1e-12)


def test_h_mzm_02_figure_3_spectrogram_reproduction() -> None:
    """Verify Hypothesis H_MZM_02: Reproduction of arXiv:2209.12930 Figure 3 (Core vs Ring vs Bulk)."""
    L = 6
    system = build_two_vortex_mzm_system(L=L, p1=(1, 1), p2=(4, 4), K=1.0, kappa=0.15)
    fig3 = reproduce_figure_3_spectrogram(system, sigma=0.06)

    # 1. Core site response: dominant MZM zero-bias peak
    assert fig3.core_lowest_energy <= 0.25
    idx_subgap = np.where(fig3.omega <= 0.30)[0]
    idx_continuum = np.where((fig3.omega > 0.40) & (fig3.omega <= 1.50))[0]

    max_core_zero_bias = float(np.max(fig3.core_spectrum.spectral_density[idx_subgap]))
    max_core_continuum = float(np.max(fig3.core_spectrum.spectral_density[idx_continuum]))
    assert max_core_zero_bias > 10.0
    assert max_core_zero_bias > 2.0 * max_core_continuum

    # 2. Ring site response: suppressed MZM zero-bias peak
    max_ring_zero_bias = float(np.max(fig3.ring_spectrum.spectral_density[idx_subgap]))
    assert max_ring_zero_bias < 1.0
    assert max_ring_zero_bias < 0.10 * max_core_zero_bias

    # 3. Bulk site response: zero subgap weight, onset at two-vison threshold
    assert fig3.bulk_lowest_energy >= 0.40
    max_bulk_zero_bias = float(np.max(fig3.bulk_spectrum.spectral_density[idx_subgap]))
    assert max_bulk_zero_bias < 0.50

    max_bulk_continuum = float(np.max(fig3.bulk_spectrum.spectral_density[idx_continuum]))
    assert max_bulk_continuum > 8.0


def test_h_mzm_03_figure_4_spatial_decay_and_localization_length() -> None:
    """Verify Hypothesis H_MZM_03: Reproduction of arXiv:2209.12930 Figure 4 (Spatial decay and xi_M)."""
    L = 6
    system = build_two_vortex_mzm_system(L=L, p1=(1, 1), p2=(4, 4), K=1.0, kappa=0.15)
    fig4 = reproduce_figure_4_spatial_decay(system, omega_cutoff=0.35, max_fit_distance=3.5)

    # 1. Low-energy weight at vortex core vs bulk
    assert fig4.low_energy_weights[0] > 0.25  # Core r=0
    # Bulk weight at r >= 3.0 drops to zero or near-zero
    bulk_idxs = np.where(fig4.radial_distances >= 3.0)[0]
    assert np.any(fig4.low_energy_weights[bulk_idxs] == 0.0)

    # 2. MZM wavefunction density monotonically decreases over short distance
    assert fig4.mzm_densities[0] > fig4.mzm_densities[1]
    assert fig4.mzm_densities[1] > fig4.mzm_densities[2]

    # 3. Fitted localization length xi_M in physical range [1.5, 4.5] a
    assert 1.5 <= fig4.fitted_xi_M <= 4.5
    assert fig4.fit_slope < 0.0
    assert fig4.fit_r_squared >= 0.85


def test_h_mzm_04_sum_rule_preservation_across_vortex_background() -> None:
    """Verify Hypothesis H_MZM_04: Local completeness sum rule holds across Core, Ring, and Bulk sites."""
    L = 6
    system = build_two_vortex_mzm_system(L=L, p1=(1, 1), p2=(4, 4), K=1.0, kappa=0.15)
    fig3 = reproduce_figure_3_spectrogram(system, sigma=0.06)

    # Integrated weights must all saturate >= 90%
    w_core = float(sum(p.weight for p in fig3.core_poles))
    w_ring = float(sum(p.weight for p in fig3.ring_poles))
    w_bulk = float(sum(p.weight for p in fig3.bulk_poles))

    assert w_core >= 0.90
    assert w_ring >= 0.90
    assert w_bulk >= 0.90

    # Numerical continuum quadratures match discrete sums to within convolution boundary limits
    assert abs(fig3.core_spectrum.integrated_weight - w_core) < 0.05
    assert abs(fig3.ring_spectrum.integrated_weight - w_ring) < 0.05
    assert abs(fig3.bulk_spectrum.integrated_weight - w_bulk) < 0.05
