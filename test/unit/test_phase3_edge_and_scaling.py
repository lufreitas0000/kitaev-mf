"""Unit tests and scientific contracts for Task T3.3.1 (Cylinder Geometry, Edge Modes & Scaling).

Verifies contracts:
- H_EDGE_01: Cylinder Geometry & Skew-Symmetry (Exact skew-symmetry norm error < 1e-15, link truncation)
- H_EDGE_02: Boundary Localization of In-Gap States (Chiral edge mode localization > 90% on outer layers)
- H_EDGE_03: Edge vs. Bulk Local Dynamic Response (Edge low-bias weight > 5.0 vs Bulk subgap silence < 0.05)
- H_EDGE_04: Exponential Decoupling of MZM Doublet (Monotonic decrease in delta_eps across L in {4, 6, 8})
"""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_mf.services.phase_3.edge_and_scaling import (
    analyze_cylinder_edge_modes,
    build_honeycomb_cylinder_system,
    compute_cylinder_edge_vs_bulk_spectra,
    sweep_finite_size_mzm_scaling,
)
from kitaev_mf.services.phase_3.gauge_theory import (
    assemble_matter_skew_matrix,
    build_honeycomb_cylinder_cluster,
    build_uniform_gauge_configuration,
)


def test_h_edge_01_cylinder_geometry_and_skew_symmetry():
    """H_EDGE_01: Verify cylinder cluster generation, bond counts, and skew-symmetry to machine precision."""
    Lx, Ly = 6, 8
    num_sites, bonds, triplets = build_honeycomb_cylinder_cluster(Lx=Lx, Ly=Ly)

    # 1. Total sites N = 2 * Lx * Ly
    assert num_sites == 2 * Lx * Ly == 96

    # 2. Total bonds: on torus it would be 3 * Lx * Ly; on cylinder with open m2,
    # the x-bonds crossing m2=0 to m2=-1 are omitted, so exactly Lx bonds are cut.
    assert len(bonds) == 3 * Lx * Ly - Lx

    # 3. Verify matter skew matrix assembly and exact anti-symmetry
    gauge = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)
    A = assemble_matter_skew_matrix(gauge, K=1.0, kappa=0.15, triplets=triplets)

    skew_error = float(np.linalg.norm(A + A.T))
    assert skew_error < 1e-15, f"Skew-symmetry violated: ||A + A^T|| = {skew_error}"

    # Verify hermiticity of i*A
    H_matter = 1j * A
    herm_error = float(np.linalg.norm(H_matter - H_matter.conj().T))
    assert herm_error < 1e-15, f"Hermiticity of i*A violated: {herm_error}"


def test_h_edge_02_chiral_edge_mode_boundary_localization():
    """H_EDGE_02: Verify emergence of in-gap chiral edge states localized on cylinder boundaries."""
    Lx, Ly = 6, 12
    system = build_honeycomb_cylinder_system(Lx=Lx, Ly=Ly, K=1.0, kappa=0.15)
    modes = analyze_cylinder_edge_modes(system, max_modes=4)

    # 1. In-gap mode existence below bulk gap
    bulk_gap = system.bulk_gap_estimate  # ~1.5588
    assert modes["in_gap_count"] >= 4, f"Expected >= 4 in-gap modes, got {modes['in_gap_count']}"

    # 2. Lowest excitation energy is nearly gapless chiral edge branch
    eps_0 = modes["energies"][0]
    assert eps_0 < 1e-4, f"Lowest chiral edge mode should be nearly gapless, got {eps_0}"

    # 3. Wavefunction boundary localization on outer layers
    boundary_fraction_0 = modes["boundary_fractions"][0]
    assert boundary_fraction_0 > 0.90, (
        f"Lowest edge mode must have > 90% boundary localization, got {boundary_fraction_0:.4f}"
    )


def test_h_edge_03_edge_vs_bulk_dynamic_spin_response():
    """H_EDGE_03: Verify edge STS low-bias conductance vs bulk subgap silence."""
    Lx, Ly = 4, 6
    system = build_honeycomb_cylinder_system(Lx=Lx, Ly=Ly, K=1.0, kappa=0.15)
    omega_grid = np.linspace(0.0, 3.0, 150)
    cont_edge, cont_bulk = compute_cylinder_edge_vs_bulk_spectra(
        system, omega_grid=omega_grid, sigma=0.06
    )

    # Frequency idx around omega = 0.10 |K| (subgap below bulk vison threshold ~0.45)
    idx_low = np.where((omega_grid >= 0.05) & (omega_grid <= 0.15))[0]

    edge_low_weight = float(np.max(cont_edge.spectral_density[idx_low]))
    bulk_low_weight = float(np.max(cont_bulk.spectral_density[idx_low]))

    # Edge site must have strong low-bias spectral weight from chiral edge cotunneling
    assert edge_low_weight > 5.0, f"Expected edge low-bias weight > 5.0, got {edge_low_weight}"

    # Bulk site on the same cylinder must have strictly suppressed subgap weight
    assert bulk_low_weight < 0.05, f"Expected bulk subgap weight < 0.05, got {bulk_low_weight}"

    # Contrast ratio between edge and bulk at low bias
    contrast = edge_low_weight / max(bulk_low_weight, 1e-12)
    assert contrast > 50.0, f"Expected edge/bulk contrast > 50x, got {contrast:.1f}x"


def test_h_edge_04_finite_size_mzm_exponential_decoupling():
    """H_EDGE_04: Verify monotonic decrease of MZM splitting delta_eps as cluster size L increases."""
    records = sweep_finite_size_mzm_scaling(L_values=[4, 6, 8], K=1.0, kappa=0.15)

    assert len(records) == 3
    rec_4, rec_6, rec_8 = records

    # 1. Distances scale with L
    assert rec_4["vortex_distance"] < rec_6["vortex_distance"] < rec_8["vortex_distance"]

    # 2. MZM splitting delta_eps decreases monotonically toward 0
    delta_4 = rec_4["delta_eps"]
    delta_6 = rec_6["delta_eps"]
    delta_8 = rec_8["delta_eps"]

    assert delta_4 > delta_6 > delta_8, (
        f"MZM splitting must decrease monotonically: L=4 ({delta_4:.4e}), L=6 ({delta_6:.4e}), L=8 ({delta_8:.4e})"
    )

    # 3. For L=8, the splitting is strongly suppressed
    assert delta_8 < 0.01, f"L=8 MZM splitting should be < 0.01, got {delta_8:.4e}"

    # 4. In-gap ratio improves with size
    assert rec_8["in_gap_ratio"] <= rec_6["in_gap_ratio"] < rec_4["in_gap_ratio"]
