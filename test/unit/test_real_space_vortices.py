"""Unit tests for real-space vortex geometry, branch cuts, and Wilson loops (HYP-W1-02)."""

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


def test_flux_free_vacuum_uniformity() -> None:
    """HYP-W1-02: Verify that in the unperturbed vacuum, W_p = +1.0 across all L^2 plaquettes."""
    L = 4
    cluster = RealSpaceTorusCluster(L=L)
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    state = make_qsl_seed(cell)

    A = cluster.assemble_skew_hamiltonian(params, state)
    omega, _ = cluster.solve_covariance_complex_eigh(A)
    fluxes = cluster.compute_plaquette_fluxes(omega)

    assert len(fluxes) == L * L
    for p, w in fluxes.items():
        assert np.isclose(w, 1.0, atol=1e-12), f"Plaquette {p} has non-unit flux W_p = {w}"


def test_single_bond_flip_adjacent_vortex_pair() -> None:
    """HYP-W1-02: Flipping one bond must invert exactly two adjacent plaquettes (W = -1.0)."""
    L = 6
    cluster = RealSpaceTorusCluster(L=L)
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    state = make_qsl_seed(cell)

    # Flip bond (0, 0, 0): z-bond in cell (0, 0)
    # Connects plaquette (0, 0) and plaquette (1, L-1)
    flipped = [(0, 0, 0)]
    A = cluster.assemble_skew_hamiltonian(params, state, flipped_bonds=flipped)
    omega, _ = cluster.solve_covariance_complex_eigh(A)
    fluxes = cluster.compute_plaquette_fluxes(omega)

    defects = [p for p, w in fluxes.items() if w < 0.0]
    assert len(defects) == 2, f"Expected 2 defect plaquettes, got {len(defects)}: {defects}"
    expected_defects = {(0, 0), (1, L - 1)}
    assert set(defects) == expected_defects

    for p in expected_defects:
        assert np.isclose(fluxes[p], -1.0, atol=1e-12)

    for p, w in fluxes.items():
        if p not in expected_defects:
            assert np.isclose(w, 1.0, atol=1e-12)


def test_dual_graph_shortest_path_separated_pair() -> None:
    """HYP-W1-02: Shortest-path cut between distant plaquettes yields strictly two vortex endpoints."""
    L = 8
    cluster = RealSpaceTorusCluster(L=L)
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    state = make_qsl_seed(cell)

    p1 = (0, 0)
    p2 = (4, 4)
    cut = cluster.generate_vortex_pair_cut(p1, p2)
    assert len(cut) > 0

    A = cluster.assemble_skew_hamiltonian(params, state, flipped_bonds=cut)
    omega, _ = cluster.solve_covariance_complex_eigh(A)
    fluxes = cluster.compute_plaquette_fluxes(omega)

    defects = [p for p, w in fluxes.items() if w < 0.0]
    assert set(defects) == {p1, p2}, f"Expected defects at {p1} and {p2}, got {defects}"

    assert np.isclose(fluxes[p1], -1.0, atol=1e-12)
    assert np.isclose(fluxes[p2], -1.0, atol=1e-12)

    vacuum_fluxes = [w for p, w in fluxes.items() if p not in {p1, p2}]
    assert np.allclose(vacuum_fluxes, 1.0, atol=1e-12)


def test_two_vison_pair_creation_energy_and_in_gap_mode() -> None:
    """HYP-W2-01 & HYP-W2-02: Vison pair costs positive energy and introduces sub-gap localized mode.

    Note: On torus geometries with 3 | L, discrete Dirac resonance zero modes occur in the
    flux-free sector, which require physical parity projection. On non-resonant geometries
    (3 does not divide L, e.g. L=4, 8), the Gaussian vacuum is gapped and defect creation energy
    is strictly positive.
    """
    L = 8
    cluster = RealSpaceTorusCluster(L=L)
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    state = make_qsl_seed(cell)

    p1 = (0, 0)
    p2 = (4, 4)
    cut = cluster.generate_vortex_pair_cut(p1, p2)

    A_vac = cluster.assemble_skew_hamiltonian(params, state)
    A_vor = cluster.assemble_skew_hamiltonian(params, state, flipped_bonds=cut)

    omega_vac, gs_e_vac = cluster.solve_covariance_complex_eigh(A_vac)
    omega_vor, gs_e_vor = cluster.solve_covariance_complex_eigh(A_vor)

    # Positive creation energy
    delta_E = (gs_e_vor - gs_e_vac) / 2.0
    assert delta_E > 0.0, f"Expected positive defect creation energy, got {delta_E}"

    # Spectral comparison: in-gap bound mode
    w_vac, _ = np.linalg.eigh(1j * A_vac.astype(np.complex128))
    w_vor, _ = np.linalg.eigh(1j * A_vor.astype(np.complex128))

    pos_vac = np.sort(w_vac[w_vac > 1e-10])
    pos_vor = np.sort(w_vor[w_vor > 1e-10])

    bulk_gap = pos_vac[0]
    vortex_mode = pos_vor[0]

    assert vortex_mode < bulk_gap, f"Vortex mode {vortex_mode} not below bulk gap {bulk_gap}"
