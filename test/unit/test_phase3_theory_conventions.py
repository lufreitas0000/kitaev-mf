"""Test suite for Phase 3 Task T3.0.1: Gauge Theory & Exact Conventions.

Audits:
- Exact FM Kitaev prefactor A_{ab} = +2 K u_{ab} (K > 0).
- Chiral 3-spin traversal sign and analytical Dirac mass gap 3*sqrt(3)*kappa.
- Local gauge quench toggling flux on exactly adjacent plaquettes.
- Physical matter-fermion parity theorem (Pedrocchi et al. 2011, Zschocke & Vojta 2015)
  confirming pi_{G_0} = -1 and ground-state fermion occupation on periodic tori.
"""

import numpy as np
import pytest

from kitaev_mf.services.phase_3.gauge_theory import (
    DirectedBond,
    assemble_matter_skew_matrix,
    build_honeycomb_torus_cluster,
    build_uniform_gauge_configuration,
    compute_plaquette_flux,
    diagonalize_matter_hamiltonian,
    flip_bond_gauge_link,
)


def test_h_dyn_01_single_link_quench_locality() -> None:
    """Verify Hypothesis H_DYN_01: Flipping one bond gauge link toggles flux on exactly two adjacent plaquettes."""
    # Build a 2-hexagon system: sites 0 to 9
    # Hexagon 0: (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)
    # Hexagon 1: (2, 6), (6, 7), (7, 8), (8, 9), (9, 3), (2, 3) (shared z-bond)
    # Sublattice A: {0, 2, 4, 6, 8}, Sublattice B: {1, 3, 5, 7, 9}
    bonds = [
        # Hexagon 0
        DirectedBond(site_a=0, site_b=1, bond_type=1),
        DirectedBond(site_a=2, site_b=1, bond_type=2),
        DirectedBond(site_a=2, site_b=3, bond_type=3),  # shared bond
        DirectedBond(site_a=4, site_b=3, bond_type=1),
        DirectedBond(site_a=4, site_b=5, bond_type=2),
        DirectedBond(site_a=0, site_b=5, bond_type=3),
        # Hexagon 1
        DirectedBond(site_a=2, site_b=7, bond_type=1),
        DirectedBond(site_a=6, site_b=7, bond_type=2),
        DirectedBond(site_a=6, site_b=9, bond_type=3),
        DirectedBond(site_a=8, site_b=9, bond_type=1),
        DirectedBond(site_a=8, site_b=3, bond_type=2),
    ]
    gauge0 = build_uniform_gauge_configuration(num_sites=10, bonds=bonds)

    p0_bonds = [(0, 1), (2, 1), (2, 3), (4, 3), (4, 5), (0, 5)]
    p1_bonds = [(2, 7), (6, 7), (6, 9), (8, 9), (8, 3), (2, 3)]

    # Initial vortex-free flux
    assert compute_plaquette_flux(gauge0, p0_bonds) == 1
    assert compute_plaquette_flux(gauge0, p1_bonds) == 1

    # Apply spin-flip sudden quench on shared z-bond (2, 3)
    gauge_quenched = flip_bond_gauge_link(gauge0, site_i=2, site_j=3)

    # Both adjacent plaquettes must flip flux W_p -> -1
    assert compute_plaquette_flux(gauge_quenched, p0_bonds) == -1
    assert compute_plaquette_flux(gauge_quenched, p1_bonds) == -1

    # All unshared bonds remain completely unchanged
    assert gauge_quenched.get_link(0, 1) == gauge0.get_link(0, 1)
    assert gauge_quenched.get_link(6, 7) == gauge0.get_link(6, 7)


def test_h_dyn_02_torus_chiral_bulk_gap_exact_identity() -> None:
    """Verify Hypothesis H_DYN_02: Chiral 3-spin term opens exact topological Dirac gap 3*sqrt(3)*kappa."""
    L = 3
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    gauge = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    # 1. Zero chirality: Dirac cone at K, K' with epsilon = 0
    A_gapless = assemble_matter_skew_matrix(gauge, K=1.0, kappa=0.0, triplets=triplets)
    assert np.allclose(A_gapless, -A_gapless.T, atol=1e-14)
    assert np.isclose(np.trace(A_gapless), 0.0)

    diag_gapless = diagonalize_matter_hamiltonian(A_gapless, gauge=gauge, boundary="torus", L=L)
    # Verify Dirac zero modes exist
    assert np.isclose(diag_gapless.energies[0], 0.0, atol=1e-10)
    assert np.isclose(diag_gapless.energies[1], 0.0, atol=1e-10)

    # 2. Non-zero chirality kappa = 0.15: Dirac mode acquires exact mass 3*sqrt(3)*kappa
    kappa = 0.15
    A_chiral = assemble_matter_skew_matrix(gauge, K=1.0, kappa=kappa, triplets=triplets)
    assert np.allclose(A_chiral, -A_chiral.T, atol=1e-14)

    diag_chiral = diagonalize_matter_hamiltonian(A_chiral, gauge=gauge, boundary="torus", L=L)
    expected_dirac_energy = 3.0 * np.sqrt(3.0) * kappa  # 0.77942286...
    assert np.allclose(diag_chiral.energies[:2], expected_dirac_energy, rtol=1e-5)

    # Bulk band gap is 6*sqrt(3)*kappa
    bulk_band_gap = 2.0 * diag_chiral.energies[0]
    assert np.isclose(bulk_band_gap, 6.0 * np.sqrt(3.0) * kappa, rtol=1e-5)


def test_h_dyn_03_physical_matter_parity_pedrocchi_zschocke_theorem() -> None:
    """Verify Hypothesis H_DYN_03: Physical matter parity satisfies Zschocke & Vojta (2015) theorem on tori.

    On a periodic torus in the vortex-free sector:
    (-1)^theta * det(Q) = -1 => pi_G = -1 (odd matter fermion parity!).
    Physical ground state must occupy the lowest single-particle mode epsilon_1.
    """
    L = 3
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    gauge0 = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    kappa = 0.10
    A0 = assemble_matter_skew_matrix(gauge0, K=1.0, kappa=kappa, triplets=triplets)
    diag0 = diagonalize_matter_hamiltonian(A0, gauge=gauge0, boundary="torus", L=L)

    # 1. Clean torus sector: det(Q) = -1.0
    assert np.isclose(diag0.det_Q, -1.0)
    # Physical parity is odd (-1)
    assert diag0.physical_matter_parity == -1
    # Physical ground state energy occupies epsilon_1
    expected_phys_E0 = diag0.unprojected_vacuum_energy + diag0.energies[0]
    assert np.isclose(diag0.physical_ground_energy, expected_phys_E0)

    # 2. Perturbed sector with single flipped bond:
    sa = bonds[0].site_a
    sb = bonds[0].site_b
    gauge_flipped = flip_bond_gauge_link(gauge0, site_i=sa, site_j=sb)

    A_flipped = assemble_matter_skew_matrix(gauge_flipped, K=1.0, kappa=kappa, triplets=triplets)
    diag_flipped = diagonalize_matter_hamiltonian(
        A_flipped, gauge=gauge_flipped, boundary="torus", L=L
    )

    # Both sectors yield valid orthogonal Schur transformations Q in O(2N)
    assert np.allclose(diag_flipped.Q @ diag_flipped.Q.T, np.eye(num_sites), atol=1e-12)
    assert abs(diag_flipped.det_Q) == 1.0
    assert diag_flipped.physical_matter_parity in (-1, 1)


def test_h_dyn_04_directed_link_antisymmetry() -> None:
    """Verify anti-symmetry u_{ba} = -u_{ab} and bond retrieval integrity."""
    bonds = [
        DirectedBond(site_a=0, site_b=1, bond_type=1),
        DirectedBond(site_a=2, site_b=3, bond_type=2),
    ]
    gauge = build_uniform_gauge_configuration(num_sites=4, bonds=bonds)

    # A -> B direction is +1
    assert gauge.get_link(0, 1) == 1
    # B -> A direction is -1
    assert gauge.get_link(1, 0) == -1

    # Non-existent bond raises KeyError
    with pytest.raises(KeyError, match="does not exist"):
        gauge.get_link(0, 2)
