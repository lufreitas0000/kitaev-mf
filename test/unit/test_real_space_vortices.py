"""Unit tests for real-space vortex geometry, branch cuts, and Wilson loops (HYP-W1-02)."""

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import InhomogeneousMeanFieldState
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


def test_unpinned_vortex_pair_self_consistent_trapping() -> None:
    """K-020, K-021, K-022: Unpinned vortex seed stably converges to exact pi-fluxes at h=0."""
    L = 4
    cluster = RealSpaceTorusCluster(L=L)
    params = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))
    vac_state = InhomogeneousMeanFieldState.from_pure_kitaev(num_cells=cluster.num_cells)

    p1 = (0, 0)
    p2 = (2, 2)
    vor_seed, cut = cluster.create_unpinned_vortex_pair_seed(vac_state, p1, p2)
    assert len(cut) > 0

    # Solve unpinned: flipped_bonds=None!
    omega_vor, final_state, met = cluster.solve_inhomogeneous_purification(
        params=params,
        initial_state=vor_seed,
        max_iter=15,
        tol=1e-7,
        accelerator="diis",
        flipped_bonds=None,
    )
    assert met["converged"] is True
    assert met["iterations"] <= 6

    # Verify exact topological flux trapping
    fluxes = cluster.compute_plaquette_fluxes(omega_vor)
    defects = [p for p, w in fluxes.items() if w < 0.0]
    assert set(defects) == {p1, p2}, f"Expected defects at {p1} and {p2}, got {defects}"
    assert np.isclose(fluxes[p1], -1.0, atol=1e-12)
    assert np.isclose(fluxes[p2], -1.0, atol=1e-12)

    # Vacuum plaquettes must be strictly +1.0
    for p, w in fluxes.items():
        if p not in {p1, p2}:
            assert np.isclose(w, 1.0, atol=1e-12)

    # Inverted gauge link along the branch cut remains -1.0
    m1, m2, b_idx = cut[0]
    c = m1 + L * m2
    gamma = [3, 1, 2][b_idx]
    assert np.isclose(final_state.U[c, b_idx, gamma, gamma], -1.0, atol=1e-12)


def test_unpinned_vortex_pair_weak_field_trapping() -> None:
    """K-022: Unpinned two-vortex pair remains topologically trapped under finite field h || [111]."""
    L = 4
    cluster = RealSpaceTorusCluster(L=L)
    params = ModelParameters(K=-1.0, h=(0.02, 0.02, 0.02))
    vac_state = InhomogeneousMeanFieldState.from_pure_kitaev(num_cells=cluster.num_cells)

    p1 = (0, 0)
    p2 = (2, 2)
    vor_seed, _ = cluster.create_unpinned_vortex_pair_seed(vac_state, p1, p2)

    # Solve unpinned with moderate iterations
    omega_vor, final_state, met = cluster.solve_inhomogeneous_purification(
        params=params,
        initial_state=vor_seed,
        max_iter=20,
        tol=1e-4,
        accelerator="diis",
        mixing_alpha=0.4,
        lambda_step=0.15,
        flipped_bonds=None,
    )

    # Check that topological fluxes remain negative exclusively at p1 and p2
    fluxes = cluster.compute_plaquette_fluxes(omega_vor)
    defects = [p for p, w in fluxes.items() if w < 0.0]
    assert set(defects) == {p1, p2}, (
        f"Defects should remain strictly at {p1} and {p2}, got {defects}"
    )
    assert fluxes[p1] < -0.95
    assert fluxes[p2] < -0.95


def test_unpinned_vortex_healing_under_high_field() -> None:
    """K-022: Under strong field h > h_c2, unpinned vortex branch cut heals back to uniform vacuum."""
    L = 4
    cluster = RealSpaceTorusCluster(L=L)
    # High polarized field h = 0.50 >> h_c2+ = 0.40
    params = ModelParameters(K=-1.0, h=(0.50, 0.50, 0.50))
    vac_state = InhomogeneousMeanFieldState.from_pure_kitaev(num_cells=cluster.num_cells)

    p1 = (0, 0)
    p2 = (1, 1)
    vor_seed, _ = cluster.create_unpinned_vortex_pair_seed(vac_state, p1, p2)

    # Under high field, the system collapses to the polarized state
    omega, final_state, met = cluster.solve_inhomogeneous_purification(
        params=params,
        initial_state=vor_seed,
        max_iter=20,
        tol=1e-5,
        accelerator="diis",
        flipped_bonds=None,
    )

    # Gauge links along the cut heal (no longer trapped at -1.0)
    fluxes = cluster.compute_plaquette_fluxes(omega)
    # At high field, all plaquette Wilson loops heal to positive or near-vacuum values
    defects = [p for p, w in fluxes.items() if w < -0.5]
    assert len(defects) == 0, (
        f"Expected vortices to heal under high field, but found defects: {defects}"
    )


def test_inhomogeneous_site_vacancy_local_relaxation() -> None:
    """K-021: Non-magnetic site vacancy relaxes without translational averaging."""
    L = 4
    cluster = RealSpaceTorusCluster(L=L)
    params = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))
    initial_state = InhomogeneousMeanFieldState.from_pure_kitaev(num_cells=cluster.num_cells)

    vacancy_site = 1  # site B in cell (0, 0)
    vacancies = {vacancy_site}

    omega, final_state, met = cluster.solve_inhomogeneous_purification(
        params=params,
        initial_state=initial_state,
        max_iter=15,
        tol=1e-5,
        accelerator="diis",
        site_vacancies=vacancies,
    )
    assert met["converged"] is True

    # The Hamiltonian block for the vacancy must be zeroed
    A = cluster.assemble_skew_hamiltonian(params, final_state, site_vacancies=vacancies)
    sl = cluster.majorana_slice(vacancy_site)
    assert np.allclose(A[sl, :], 0.0, atol=1e-15)
    assert np.allclose(A[:, sl], 0.0, atol=1e-15)


def test_plaquette_boundary_sites_topology() -> None:
    """K-020, K-021: Verify 6-site hexagonal perimeter and 2-site shared interface."""
    L = 6
    cluster = RealSpaceTorusCluster(L=L)

    p0 = (0, 0)
    sites0 = cluster.compute_plaquette_boundary_sites(p0)
    assert len(sites0) == 6

    # 3 A-sites (sublattice 0) and 3 B-sites (sublattice 1)
    sublattices = [s % 2 for s in sites0]
    assert sublattices.count(0) == 3
    assert sublattices.count(1) == 3

    # Neighboring plaquette (0, 1) shares exactly the z-link (2 sites: 1 A, 1 B)
    p_nbr = (0, 1)
    sites_nbr = cluster.compute_plaquette_boundary_sites(p_nbr)
    shared = sites0 & sites_nbr
    assert len(shared) == 2
    assert {s % 2 for s in shared} == {0, 1}


def test_vortex_mode_core_localization_weight() -> None:
    """K-022, K-023: Sub-gap vortex bound mode exhibits strong core localization (W_core / W_unif >= 2.5x)."""
    L = 8
    cluster = RealSpaceTorusCluster(L=L)
    params = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))
    state = InhomogeneousMeanFieldState.from_pure_kitaev(num_cells=cluster.num_cells)

    p1 = (0, 0)
    p2 = (4, 4)
    cut = cluster.generate_vortex_pair_cut(p1, p2)
    A = cluster.assemble_skew_hamiltonian(params, state, flipped_bonds=cut)

    w, v = np.linalg.eigh(1j * A.astype(np.complex128))
    pos_mask = w > 1e-10
    pos_w = w[pos_mask]
    pos_v = v[:, pos_mask]

    # Lowest excitation is the in-gap vortex core mode
    idx0 = int(np.argmin(pos_w))
    v0 = pos_v[:, idx0]
    w_core, w_uniform, ratio = cluster.compute_vortex_mode_core_weight(v0, p1, p2)

    assert np.isclose(w_uniform, 12.0 / (2.0 * L * L), atol=1e-12)
    assert ratio >= 2.5, f"Expected vortex core enrichment >= 2.5x, got {ratio:.2f}x"
    assert w_core > 0.30, f"Expected >= 30% wavefunction weight on 12 core sites, got {w_core:.4f}"

    # Delocalized bulk mode in the continuum
    idx_bulk = len(pos_w) // 2
    v_bulk = pos_v[:, idx_bulk]
    w_core_bulk, _, ratio_bulk = cluster.compute_vortex_mode_core_weight(v_bulk, p1, p2)
    assert ratio_bulk < 1.5, f"Expected delocalized bulk ratio < 1.5x, got {ratio_bulk:.2f}x"


def test_finite_size_vs_chiral_gap_demarcation() -> None:
    """[B2], K-023: Demarcate finite-size quantization gap Delta_FS from topological chiral Dirac mass.

    On resonant tori (3 | L, e.g. L=6), the Dirac point lies on the momentum grid,
    yielding exact Dirac zero modes Delta_FS = 0 at h = 0.
    On non-resonant tori (3 does not divide L, e.g. L=4, 8), the grid misses the Dirac
    point, producing an apparent finite-size quantization gap Delta_FS > 0.
    """
    # 1. Resonant geometry L=6: Delta_FS = 0 exactly
    c6 = RealSpaceTorusCluster(L=6)
    cell = make_honeycomb_2site()
    state_vac = make_qsl_seed(cell)
    params_0 = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))

    A6_0 = c6.assemble_skew_hamiltonian(params_0, state_vac)
    w6_0 = np.linalg.eigvalsh(1j * A6_0.astype(np.complex128))
    # Count exact Dirac zero modes (|w| < 1e-12)
    zero_modes_6 = np.sum(np.abs(w6_0) < 1e-12)
    assert zero_modes_6 > 0, "Resonant L=6 torus must exhibit exact Dirac zero modes at h=0"

    # 2. Non-resonant geometries: Delta_FS > 0 purely from momentum quantization
    c4 = RealSpaceTorusCluster(L=4)
    A4_0 = c4.assemble_skew_hamiltonian(params_0, state_vac)
    w4_0 = np.linalg.eigvalsh(1j * A4_0.astype(np.complex128))
    delta_fs_4 = float(np.min(np.abs(w4_0)))
    assert np.isclose(delta_fs_4, 1.049729, atol=1e-5), (
        f"L=4 finite-size gap unexpected: {delta_fs_4}"
    )

    c8 = RealSpaceTorusCluster(L=8)
    A8_0 = c8.assemble_skew_hamiltonian(params_0, state_vac)
    w8_0 = np.linalg.eigvalsh(1j * A8_0.astype(np.complex128))
    delta_fs_8 = float(np.min(np.abs(w8_0)))
    assert np.isclose(delta_fs_8, 0.828427, atol=1e-5), (
        f"L=8 finite-size gap unexpected: {delta_fs_8}"
    )
    # As L increases from 4 to 8, Delta_FS scales towards 0 as O(1/L)
    assert delta_fs_8 < delta_fs_4


def test_exact_12_majorana_pfaffian_vacuum_and_vortex_at_h0() -> None:
    """Contract H2: At h=0, exact 12-Majorana Pfaffian equals factored link product W_p = +-1.0."""
    L = 4
    cluster = RealSpaceTorusCluster(L=L)
    cell = make_honeycomb_2site()
    params0 = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))
    state = make_qsl_seed(cell)
    p0 = (0, 0)

    # 1. Vacuum state
    A0 = cluster.assemble_skew_hamiltonian(params0, state)
    omega0, _ = cluster.solve_covariance_complex_eigh(A0)
    w_fact_vac = cluster.compute_plaquette_fluxes(omega0)[p0]
    w_pfaff_vac = cluster.compute_plaquette_pfaffian_flux(omega0, p0)

    assert np.isclose(w_fact_vac, 1.0, atol=1e-12)
    assert np.isclose(w_pfaff_vac, 1.0, atol=1e-12)
    assert np.isclose(w_fact_vac, w_pfaff_vac, atol=1e-12)

    # 2. Vortex pair state (flipped bond (0, 0, 0) flips plaquettes (0, 0) and (1, L-1))
    cut = [(0, 0, 0)]
    A_vor = cluster.assemble_skew_hamiltonian(params0, state, flipped_bonds=cut)
    omega_vor, _ = cluster.solve_covariance_complex_eigh(A_vor)

    w_fact_vor = cluster.compute_plaquette_fluxes(omega_vor)[p0]
    w_pfaff_vor = cluster.compute_plaquette_pfaffian_flux(omega_vor, p0)

    assert np.isclose(w_fact_vor, -1.0, atol=1e-12)
    assert np.isclose(w_pfaff_vor, -1.0, atol=1e-12)
    assert np.isclose(w_fact_vor, w_pfaff_vor, atol=1e-12)


def test_exact_pfaffian_vs_factored_flux_quadratic_scaling() -> None:
    """Contract H2: At h > 0, discrepancy |W_p^fact - Pf(Omega|_partial p)| scales quadratically as O(h^2)."""
    L = 4
    cluster = RealSpaceTorusCluster(L=L)
    cell = make_honeycomb_2site()
    seed = make_qsl_seed(cell)
    inhom_seed = InhomogeneousMeanFieldState.from_homogeneous(seed, cluster.num_cells)
    p0 = (0, 0)

    h_vals = [0.02, 0.04]
    diffs: list[float] = []

    for h_val in h_vals:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)
        omega, _, met = cluster.solve_inhomogeneous_purification(
            params=params,
            initial_state=inhom_seed,
            max_iter=30,
            tol=1e-6,
        )
        assert met["converged"]

        w_fact = cluster.compute_plaquette_fluxes(omega)[p0]
        w_pfaff = cluster.compute_plaquette_pfaffian_flux(omega, p0)
        diff = float(abs(w_fact - w_pfaff))
        diffs.append(diff)
        # Quantum gauge fluctuations soften flux degradation: Pf > W_fact
        assert w_pfaff > w_fact

    # Doubling field from 0.02 to 0.04 must multiply discrepancy by 2^2 = 4.0
    scaling_ratio = diffs[1] / diffs[0]
    assert np.isclose(scaling_ratio, 4.0, atol=0.05), (
        f"Expected quadratic scaling ratio ~4.0, got {scaling_ratio:.4f}"
    )
