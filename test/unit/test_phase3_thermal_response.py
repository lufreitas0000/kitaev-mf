"""Unit tests and scientific contracts for Task T3.3.2 (Finite-Temperature Response & Detailed Balance).

Verifies contracts:
- H_TEMP_01: KMS Detailed Balance Exactness (|S(-omega)/S(omega) - exp(-beta*omega)| < 1e-10)
- H_TEMP_02: Thermal MZM Peak Suppression & Broadening (Monotonic decrease, > 50% drop from T=0.02 to 0.20)
- H_TEMP_03: Thermal Subgap Infill (Bulk subgap weight increases monotonically with T)
- H_TEMP_04: Total Sum Rule Conservation (Two-sided integral conserved across T)
"""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_mf.services.phase_3.mzm_spectroscopy import (
    TwoVortexMZMSystem,
    build_two_vortex_mzm_system,
)
from kitaev_mf.services.phase_3.spectral_engine import (
    compute_local_spectral_poles,
    evaluate_thermal_spectral_function,
    verify_kms_detailed_balance,
)


@pytest.fixture(scope="module")
def shared_vortex_system() -> TwoVortexMZMSystem:
    """Shared L=6 two-vortex cluster for thermal tests."""
    return build_two_vortex_mzm_system(L=6, p1=(1, 1), p2=(4, 4), K=1.0, kappa=0.15)


def test_h_temp_01_kms_detailed_balance(shared_vortex_system: TwoVortexMZMSystem) -> None:
    """H_TEMP_01: Verify exact KMS detailed balance relation S(-omega) = exp(-beta * omega) S(omega)."""
    sys = shared_vortex_system
    p1_bonds = sys.geom.plaquette_bonds(1, 1)
    core_site = p1_bonds[0][0]

    spec_core = compute_local_spectral_poles(
        gauge=sys.gauge,
        site_j=core_site,
        bond_type=3,
        geom=sys.geom,
        K=1.0,
        kappa=0.15,
        boundary="obc",
        bonds=sys.bonds,
        triplets=sys.triplets,
        diag_initial=sys.diag,
    )

    omega_grid = np.linspace(-4.0, 4.0, 801)

    for T in (0.05, 0.10, 0.20):
        spec_T = evaluate_thermal_spectral_function(
            spec_core,
            omega_grid=omega_grid,
            temperature=T,
            broadening_param=0.06,
        )

        res = verify_kms_detailed_balance(
            spec_T,
            temperature=T,
            omega_test_points=[0.10, 0.20, 0.30, 0.50, 0.80],
            tol=1e-10,
        )

        assert res["is_valid"], (
            f"KMS detailed balance violated at T={T}: max error = {res['max_relative_error']}"
        )
        assert res["max_relative_error"] < 1e-10


def test_h_temp_02_thermal_mzm_peak_suppression(
    shared_vortex_system: TwoVortexMZMSystem,
) -> None:
    """H_TEMP_02: Verify monotonic decrease of MZM zero-bias peak height with increasing temperature."""
    sys = shared_vortex_system
    p1_bonds = sys.geom.plaquette_bonds(1, 1)
    core_site = p1_bonds[0][0]

    spec_core = compute_local_spectral_poles(
        gauge=sys.gauge,
        site_j=core_site,
        bond_type=3,
        geom=sys.geom,
        K=1.0,
        kappa=0.15,
        boundary="obc",
        bonds=sys.bonds,
        triplets=sys.triplets,
        diag_initial=sys.diag,
    )

    omega_grid = np.linspace(-4.0, 4.0, 801)
    temperatures = [0.02, 0.05, 0.10, 0.20]
    peak_heights: list[float] = []

    for T in temperatures:
        spec_T = evaluate_thermal_spectral_function(
            spec_core, omega_grid=omega_grid, temperature=T, broadening_param=0.06
        )
        # Low positive bias region around MZM peak
        idx_peak = np.where((omega_grid >= 0.05) & (omega_grid <= 0.25))[0]
        peak = float(np.max(spec_T.spectral_density[idx_peak]))
        peak_heights.append(peak)

    # 1. Monotonic decrease with T
    for i in range(len(peak_heights) - 1):
        assert peak_heights[i] > peak_heights[i + 1], (
            f"Peak height must decrease monotonically: {peak_heights}"
        )

    # 2. Peak amplitude drops by more than 50% between T=0.02 and T=0.20
    suppression_ratio = peak_heights[-1] / peak_heights[0]
    assert suppression_ratio < 0.50, (
        f"Expected suppression ratio < 0.50, got {suppression_ratio:.3f}"
    )


def test_h_temp_03_thermal_subgap_infill(
    shared_vortex_system: TwoVortexMZMSystem,
) -> None:
    """H_TEMP_03: Verify that bulk subgap spectral weight increases monotonically with temperature."""
    sys = shared_vortex_system
    geom = sys.geom
    best_dist = -1.0
    bulk_site = 0
    for m2 in range(geom.L):
        for m1 in range(geom.L):
            d = min(geom.torus_distance((m1, m2), sys.p1), geom.torus_distance((m1, m2), sys.p2))
            if d > best_dist:
                best_dist = d
                bulk_site = geom.plaquette_bonds(m1, m2)[0][0]

    spec_bulk = compute_local_spectral_poles(
        gauge=sys.gauge,
        site_j=bulk_site,
        bond_type=3,
        geom=geom,
        K=1.0,
        kappa=0.15,
        boundary="obc",
        bonds=sys.bonds,
        triplets=sys.triplets,
        diag_initial=sys.diag,
    )

    omega_grid = np.linspace(-4.0, 4.0, 801)
    domega = omega_grid[1] - omega_grid[0]
    subgap_idx = np.where((omega_grid >= 0.0) & (omega_grid <= 0.30))[0]

    temperatures = [0.00, 0.02, 0.05, 0.10, 0.20]
    subgap_weights: list[float] = []

    for T in temperatures:
        spec_T = evaluate_thermal_spectral_function(
            spec_bulk, omega_grid=omega_grid, temperature=T, broadening_param=0.06
        )
        w_subgap = float(np.sum(spec_T.spectral_density[subgap_idx]) * domega / (2.0 * np.pi))
        subgap_weights.append(w_subgap)

    # Monotonic increase in bulk subgap weight
    for i in range(len(subgap_weights) - 1):
        assert subgap_weights[i] < subgap_weights[i + 1], (
            f"Bulk subgap weight must increase monotonically: {subgap_weights}"
        )

    # Subgap weight at T=0.20 exceeds T=0.00 by at least 10x
    ratio = subgap_weights[-1] / max(subgap_weights[0], 1e-12)
    assert ratio > 10.0, f"Expected subgap infill ratio > 10x, got {ratio:.1f}x"


def test_h_temp_04_sum_rule_conservation(
    shared_vortex_system: TwoVortexMZMSystem,
) -> None:
    """H_TEMP_04: Verify conservation of two-sided integrated spectral weight across temperatures."""
    sys = shared_vortex_system
    p1_bonds = sys.geom.plaquette_bonds(1, 1)
    core_site = p1_bonds[0][0]

    spec_core = compute_local_spectral_poles(
        gauge=sys.gauge,
        site_j=core_site,
        bond_type=3,
        geom=sys.geom,
        K=1.0,
        kappa=0.15,
        boundary="obc",
        bonds=sys.bonds,
        triplets=sys.triplets,
        diag_initial=sys.diag,
    )

    # Broad frequency grid covering all poles up to 8.0 |K|
    omega_grid = np.linspace(-8.0, 8.0, 1601)

    # Integrated weight at T=0
    spec_0 = evaluate_thermal_spectral_function(
        spec_core, omega_grid=omega_grid, temperature=0.0, broadening_param=0.06
    )
    w0 = spec_0.integrated_weight

    for T in (0.05, 0.10, 0.20):
        spec_T = evaluate_thermal_spectral_function(
            spec_core, omega_grid=omega_grid, temperature=T, broadening_param=0.06
        )
        wT = spec_T.integrated_weight
        # Quadrature difference should be less than 5% across temperatures
        rel_diff = abs(wT - w0) / w0
        assert rel_diff < 0.05, f"Sum rule violated at T={T}: rel_diff = {rel_diff:.4f}"
