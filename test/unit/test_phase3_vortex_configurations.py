"""Unit test suite for Phase 3 Task T3.1.1: Static Flux Configurations & Two-Vortex State Generation.

Audits:
- Hypothesis H_VOR_01: Dual branch cut exactness and two-vortex isolation.
- Hypothesis H_VOR_02: MZM sub-gap isolation in the chiral phase.
- Hypothesis H_VOR_03: Vortex translation selection rule under local spin quench.
- Hypothesis H_VOR_04: Additional vortex-pair creation on unshared bulk bonds.
"""

import numpy as np
import pytest

from kitaev_mf.services.phase_3.gauge_theory import (
    assemble_matter_skew_matrix,
    build_honeycomb_torus_cluster,
    build_uniform_gauge_configuration,
    diagonalize_matter_hamiltonian,
)
from kitaev_mf.services.phase_3.vortex_engine import (
    HoneycombGeometry,
    apply_local_spin_quench,
    categorize_spin_quench,
    create_two_vortex_gauge_configuration,
    evaluate_all_plaquette_fluxes,
    identify_vortices,
)


def test_h_vor_01_dual_branch_cut_exactness() -> None:
    """Verify Hypothesis H_VOR_01: Dual branch cut flips flux on target plaquettes exclusively."""
    L = 4
    geom = HoneycombGeometry(L=L)
    num_sites, bonds, _ = build_honeycomb_torus_cluster(L=L)
    base_gauge = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    # 1. Base gauge is vortex-free (W_p = +1 everywhere)
    init_fluxes = evaluate_all_plaquette_fluxes(base_gauge, geom)
    assert len(init_fluxes) == L * L
    assert all(w == 1 for w in init_fluxes.values())
    assert identify_vortices(base_gauge, geom) == []

    # 2. Insert two vortices at p1=(0, 0) and p2=(2, 2)
    p1 = (0, 0)
    p2 = (2, 2)
    vortex_gauge, cut_bonds = create_two_vortex_gauge_configuration(base_gauge, p1, p2, geom)

    assert len(cut_bonds) >= 1
    # Check that cut bonds are all distinct
    assert len(cut_bonds) == len(set(cut_bonds))

    # 3. Verify exactly two vortices exist at p1 and p2
    vortex_fluxes = evaluate_all_plaquette_fluxes(vortex_gauge, geom)
    vortices = identify_vortices(vortex_gauge, geom)
    assert set(vortices) == {p1, p2}
    assert vortex_fluxes[p1] == -1
    assert vortex_fluxes[p2] == -1

    # Remaining L^2 - 2 plaquettes must have W_p = +1
    clean_count = sum(1 for p, w in vortex_fluxes.items() if w == 1)
    assert clean_count == L * L - 2


def test_h_vor_02_mzm_in_gap_isolation() -> None:
    """Verify Hypothesis H_VOR_02: Isolated two-vortex background produces sub-gap MZM mode."""
    L = 4
    geom = HoneycombGeometry(L=L)
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    base_gauge = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    kappa = 0.15
    K = 1.0

    # 1. Vacuum spectrum (vortex-free)
    A_clean = assemble_matter_skew_matrix(base_gauge, K=K, kappa=kappa, triplets=triplets)
    diag_clean = diagonalize_matter_hamiltonian(A_clean, gauge=base_gauge, boundary="torus", L=L)
    lowest_clean_mode = diag_clean.energies[0]

    # 2. Two-vortex spectrum with vortices at (0, 0) and (2, 2)
    p1 = (0, 0)
    p2 = (2, 2)
    vortex_gauge, _ = create_two_vortex_gauge_configuration(base_gauge, p1, p2, geom)
    A_vor = assemble_matter_skew_matrix(vortex_gauge, K=K, kappa=kappa, triplets=triplets)
    diag_vor = diagonalize_matter_hamiltonian(A_vor, gauge=vortex_gauge, boundary="torus", L=L)

    lowest_vor_mode = diag_vor.energies[0]

    # Verify that the two-vortex state binds an isolated in-gap sub-gap mode:
    # On L=4 (which misses Dirac points), clean lowest mode is ~ 1.400 due to finite-size quantization.
    # On L=3 or L=6, clean lowest mode is exactly 3*sqrt(3)*kappa ~ 0.7794.
    assert lowest_vor_mode < 0.35
    assert lowest_clean_mode > 1.30
    subgap_ratio = lowest_vor_mode / lowest_clean_mode
    assert subgap_ratio < 0.25, f"Subgap ratio {subgap_ratio} is not isolated below bulk gap"

    # Ensure the sub-gap mode is isolated from the rest of the vortex spectrum
    assert diag_vor.energies[1] - lowest_vor_mode > 0.5


def test_h_vor_03_vortex_translation_vs_pair_creation() -> None:
    """Verify Hypotheses H_VOR_03 & H_VOR_04: Local spin quench translates vortex or creates pair."""
    L = 4
    geom = HoneycombGeometry(L=L)
    num_sites, bonds, _ = build_honeycomb_torus_cluster(L=L)
    base_gauge = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    # Reference state: two vortices at (0, 0) and (2, 2)
    vortex_gauge, _ = create_two_vortex_gauge_configuration(base_gauge, (0, 0), (2, 2), geom)
    assert set(identify_vortices(vortex_gauge, geom)) == {(0, 0), (2, 2)}

    # Case A: Quench bond adjacent to plaquette (0, 0)
    # Site s_A(0, 0) with z-bond (bond_type=3) directly bounds plaquette (0, 0)
    site_near = geom.site_index(0, 0, 0)
    gauge_translated, flipped_bond_A = apply_local_spin_quench(
        vortex_gauge, site_j=site_near, bond_type=3, geom=geom
    )

    analysis_A = categorize_spin_quench(vortex_gauge, gauge_translated, geom)
    assert analysis_A["transition_type"] == "vortex_translation"
    assert analysis_A["initial_vortex_count"] == 2
    assert analysis_A["quenched_vortex_count"] == 2
    # The vortex at (0, 0) moved to adjacent plaquette sharing the z-bond
    assert (0, 0) in analysis_A["displaced_vortices"][0]

    # Case B: Quench distant bulk bond away from (0, 0) and (2, 2)
    # E.g., cell (0, 2) with z-bond
    site_distant = geom.site_index(0, 2, 0)
    gauge_pair_created, flipped_bond_B = apply_local_spin_quench(
        vortex_gauge, site_j=site_distant, bond_type=3, geom=geom
    )

    analysis_B = categorize_spin_quench(vortex_gauge, gauge_pair_created, geom)
    assert analysis_B["transition_type"] == "pair_creation"
    assert analysis_B["initial_vortex_count"] == 2
    assert analysis_B["quenched_vortex_count"] == 4
    # Both original vortices remain present
    assert (0, 0) in analysis_B["quenched_vortices"]
    assert (2, 2) in analysis_B["quenched_vortices"]


def test_h_vor_04_torus_distance_and_geometry() -> None:
    """Verify HoneycombGeometry torus distance symmetry and coordinate bijection."""
    geom = HoneycombGeometry(L=4)

    # 1. Coordinate roundtrip
    for s in range(geom.num_sites):
        m1, m2, sub = geom.site_coordinates(s)
        assert geom.site_index(m1, m2, sub) == s

    # 2. Torus distance symmetry and identity
    p1 = (1, 1)
    p2 = (3, 2)
    assert geom.torus_distance(p1, p1) == 0.0
    d12 = geom.torus_distance(p1, p2)
    d21 = geom.torus_distance(p2, p1)
    assert np.isclose(d12, d21)
    assert d12 > 0.0

    # 3. Periodic equivalence including arbitrary periodic multiples
    p_img = (p1[0] + geom.L, p1[1])
    assert np.isclose(geom.torus_distance(p1, p_img), 0.0)
    p_img_large = (p1[0] + 5 * geom.L, p1[1] - 3 * geom.L)
    assert np.isclose(geom.torus_distance(p1, p_img_large), 0.0)
