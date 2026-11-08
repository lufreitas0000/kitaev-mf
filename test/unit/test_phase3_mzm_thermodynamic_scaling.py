"""Unit tests for Phase 3 Task T3.W1.2: MZM Thermodynamic Scaling & Sharpness.

Validates Hypothesis Contracts:
- H_MZM_SCALE_01: Exponential suppression of MZM splitting with distance R.
- H_MZM_SCALE_02: Convergence of zero-bias MZM peak energy as L -> infty (omega_MZM <= 0.05 for L >= 16).
- H_MZM_SCALE_03: Saturated zero-bias sharpness ratio S_subgap >= 0.15 across all L in [6, 20].
- H_MZM_SCALE_04: Bulk-site subgap vanishing and gap protection (omega_bulk >= 0.20 |K| for L >= 12).
"""

from kitaev_mf.services.phase_3.mzm_thermodynamic_scaling import (
    compute_thermodynamic_zero_bias_sharpness,
    measure_mzm_hybridization_scaling,
)


def test_h_mzm_scale_01_exponential_suppression_commensurate() -> None:
    """Verify H_MZM_SCALE_01: MZM splitting decays exponentially across commensurate sizes 3 | L."""
    data, fit = measure_mzm_hybridization_scaling(L_values=[6, 12, 18, 24])

    assert fit.r_squared >= 0.85, f"Expected R^2 >= 0.85, got {fit.r_squared:.4f}"
    assert 1.5 <= fit.xi_M <= 4.5, f"Expected xi_M in [1.5, 4.5]a, got {fit.xi_M:.3f}"
    assert data[-1].delta_eps_0 < 0.005, (
        f"Expected splitting at L=24 to be < 0.005, got {data[-1].delta_eps_0}"
    )


def test_h_mzm_scale_02_zero_bias_convergence() -> None:
    """Verify H_MZM_SCALE_02: zero-bias MZM peak energy converges toward zero with L."""
    results = compute_thermodynamic_zero_bias_sharpness(L_values=[6, 8, 12, 16, 20])
    omega_by_L = {r.L: r.omega_MZM for r in results}

    # As L increases from 6 to 20, zero-bias peak converges
    assert omega_by_L[6] >= 0.10, f"L=6 should have hybridized peak: {omega_by_L[6]}"
    assert omega_by_L[16] <= 0.05, f"L=16 peak should be <= 0.05: {omega_by_L[16]}"
    assert omega_by_L[20] <= 0.01, f"L=20 peak should be <= 0.01: {omega_by_L[20]}"


def test_h_mzm_scale_03_zero_bias_sharpness_saturation() -> None:
    """Verify H_MZM_SCALE_03: zero-bias sharpness ratio S_subgap >= 0.15 for all L."""
    results = compute_thermodynamic_zero_bias_sharpness(L_values=[6, 8, 12, 16, 20])
    for r in results:
        assert r.sharpness_ratio >= 0.15, (
            f"L={r.L}: Expected sharpness >= 0.15, got {r.sharpness_ratio:.4f}"
        )
        assert r.total_spectral_weight >= 0.70, (
            f"L={r.L}: Expected total weight >= 0.70, got {r.total_spectral_weight:.4f}"
        )


def test_h_mzm_scale_04_bulk_site_gap_protection() -> None:
    """Verify H_MZM_SCALE_04: bulk-site quench has zero subgap weight and protected gap >= 0.20."""
    results = compute_thermodynamic_zero_bias_sharpness(L_values=[6, 8, 12, 16, 20])
    for r in results:
        # Bulk minimum excitation frequency must respect the two-vison threshold >= 0.20 |K|
        assert r.bulk_min_energy >= 0.20, (
            f"L={r.L}: Expected bulk minimum excitation >= 0.20 |K|, got {r.bulk_min_energy:.4f}"
        )
        # Bulk subgap weight below omega_c = 0.12 |K| must strictly vanish (< 1e-5)
        assert r.bulk_subgap_weight < 1e-5, (
            f"L={r.L}: Expected bulk subgap weight < 1e-5, got {r.bulk_subgap_weight:.4e}"
        )
