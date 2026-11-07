"""Unit tests for finite-torus geometry, matter-parity selection rules, and physical state projectors.

Validates Backlog Items B-004, B-010, and B-011.
"""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_mf.services.exact_matter import (
    build_antisymmetric_majorana_matrix_A,
    build_matter_matrix_T,
    build_reference_links,
    compute_matter_vacuum_parity,
    compute_physical_ground_state_energy,
    compute_torus_gauge_sign,
    compute_torus_plaquette_fluxes,
    compute_torus_wilson_loops,
    coord_to_idx,
    is_physical_state,
)


def compute_explicit_clifford_permutation_sign(N1: int, N2: int) -> int:
    """Compute the exact Clifford permutation sign between gauge Majoranas and U bond product."""
    N = N1 * N2
    # B ordering: (0, j, a) for A, then (1, k, a) for B
    maj_list_B = []
    for j in range(N):
        for a in range(3):
            maj_list_B.append((0, j, a))
    for k in range(N):
        for a in range(3):
            maj_list_B.append((1, k, a))
    maj_to_idx = {m: idx for idx, m in enumerate(maj_list_B)}

    def coord(j: int) -> tuple[int, int]:
        return (j % N1, j // N1)

    def P(j: int, a: int) -> int:
        m, n = coord(j)
        if a == 0:
            return coord_to_idx(m + 1, n, N1, N2)
        elif a == 1:
            return coord_to_idx(m, n + 1, N1, N2)
        else:
            return coord_to_idx(m, n, N1, N2)

    maj_list_U = []
    for a in range(3):
        for j in range(N):
            maj_list_U.append((0, j, a))
            maj_list_U.append((1, P(j, a), a))

    perm = [maj_to_idx[m] for m in maj_list_U]
    # Parity of permutation:
    inv_count = 0
    n_elem = len(perm)
    for i in range(n_elem):
        for j in range(i + 1, n_elem):
            if perm[i] > perm[j]:
                inv_count += 1

    sign_U = (-1) ** inv_count
    # Collecting the phase from Clifford factoring:
    phase = sign_U * ((-1) ** (N + N * (N - 1) // 2))
    return int(phase)


@pytest.mark.parametrize(
    ("N1", "N2"),
    [
        (2, 2),
        (2, 3),
        (3, 3),
        (2, 4),
        (3, 4),
        (4, 4),
    ],
)
def test_torus_geometry_and_plaquette_flux_identities(N1: int, N2: int) -> None:
    """Verify exact combinatorial flux identities from notes/05_finite_torus_conventions.md."""
    # 1. Reference gauge fluxes
    u_ref = build_reference_links(N1, N2)
    W_ref = compute_torus_plaquette_fluxes(u_ref, N1, N2)
    assert np.all(W_ref == 1), "All plaquettes in reference gauge must have W_r = +1"
    assert int(np.prod(W_ref)) == 1

    # 2. Random gauge: prod_r W_r must be identically +1
    rng = np.random.default_rng(42 + N1 * 10 + N2)
    u_rand = rng.choice([-1, 1], size=(N1, N2, 3))
    W_rand = compute_torus_plaquette_fluxes(u_rand, N1, N2)
    assert int(np.prod(W_rand)) == 1, "Product of all plaquette fluxes must be identically +1"

    # 3. Single bond flip must flip exactly two adjacent plaquettes
    u_flip = u_ref.copy()
    u_flip[0, 0, 0] = -u_flip[0, 0, 0]  # flip x-bond at (0, 0)
    W_flip = compute_torus_plaquette_fluxes(u_flip, N1, N2)
    flipped_indices = np.argwhere(W_flip == -1)
    assert len(flipped_indices) == 2, (
        f"Single bond flip must create exactly 2 vortices, got {len(flipped_indices)}"
    )

    # 4. Gauge transformation invariance
    # Local site signs eta_A(r), eta_B(r)
    eta_A = rng.choice([-1, 1], size=(N1, N2))
    eta_B = rng.choice([-1, 1], size=(N1, N2))
    u_gauge = np.zeros_like(u_rand)
    for m in range(N1):
        for n in range(N2):
            # x-bond: A_(m, n) -> B_(m+1, n)
            u_gauge[m, n, 0] = eta_A[m, n] * u_rand[m, n, 0] * eta_B[(m + 1) % N1, n]
            # y-bond: A_(m, n) -> B_(m, n+1)
            u_gauge[m, n, 1] = eta_A[m, n] * u_rand[m, n, 1] * eta_B[m, (n + 1) % N2]
            # z-bond: A_(m, n) -> B_(m, n)
            u_gauge[m, n, 2] = eta_A[m, n] * u_rand[m, n, 2] * eta_B[m, n]

    W_gauge = compute_torus_plaquette_fluxes(u_gauge, N1, N2)
    assert np.array_equal(W_gauge, W_rand), (
        "Plaquette fluxes must be strictly invariant under local gauge transformations"
    )

    # Primitive Wilson loops must also be invariant under local gauge transformations
    w1_rand, w2_rand = compute_torus_wilson_loops(u_rand, N1, N2)
    w1_gauge, w2_gauge = compute_torus_wilson_loops(u_gauge, N1, N2)
    assert (w1_rand, w2_rand) == (w1_gauge, w2_gauge), (
        "Wilson loops must be gauge-invariant under local transformations"
    )


@pytest.mark.parametrize(
    ("N1", "N2"),
    [
        (2, 2),
        (2, 3),
        (3, 2),
        (3, 3),
        (2, 4),
        (4, 2),
        (3, 4),
        (4, 3),
        (4, 4),
    ],
)
def test_exact_clifford_permutation_sign_theorem(N1: int, N2: int) -> None:
    """Prove that sigma_geom = (-1)**(N1 + N2) matches explicit Clifford permutation sign."""
    sigma_analytical = (-1) ** (N1 + N2)
    sigma_clifford = compute_explicit_clifford_permutation_sign(N1, N2)
    assert sigma_analytical == sigma_clifford, (
        f"Analytical sign (-1)^({N1}+{N2})={sigma_analytical} differs from "
        f"Clifford permutation sign {sigma_clifford}"
    )


def test_gaussian_vacuum_parity_against_fock_oracle() -> None:
    """Verify pi_M(|0_M; u>) = sign(det(T)) against exact Jordan-Wigner many-body Fock diagonalizer."""
    # Test on N = 3 complex fermions (dimension 2^3 = 8)
    N = 3
    dim = 1 << N

    # Construct Jordan-Wigner annihilation operators f_m
    f = []
    for m in range(N):
        fm = np.zeros((dim, dim), dtype=np.float64)
        for state in range(dim):
            if (state >> m) & 1:
                state_new = state ^ (1 << m)
                sign = (-1) ** bin(state & ((1 << m) - 1)).count("1")
                fm[state_new, state] = float(sign)
        f.append(fm)

    # Parity matrix: prod_m (1 - 2 f_m^\dagger f_m) = diag((-1)**popcount)
    pi_M_mat = np.diag([float((-1) ** bin(state).count("1")) for state in range(dim)])

    rng = np.random.default_rng(101)
    for _ in range(5):
        # Generate random real orthogonal U, V
        U, _ = np.linalg.qr(rng.standard_normal((N, N)))
        V, _ = np.linalg.qr(rng.standard_normal((N, N)))
        s = np.array([1.0, 2.0, 3.0])
        T = U @ np.diag(s) @ V.T

        X = 0.5 * (U + V)
        Y = 0.5 * (U - V)

        # Bogoliubov annihilation operators: alpha_n = sum_m [ X_mn f_m + Y_mn f_m^\dagger ]
        alpha = []
        for n in range(N):
            an = np.zeros((dim, dim), dtype=np.float64)
            for m in range(N):
                an += X[m, n] * f[m] + Y[m, n] * f[m].T
            alpha.append(an)

        # The Gaussian vacuum |vac> is the null space of all alpha_n
        stacked = np.vstack(alpha)
        _, _, vh = np.linalg.svd(stacked)
        vac = vh[-1]
        vac = vac / np.linalg.norm(vac)

        # Expectation of fermion parity in the Gaussian vacuum:
        parity_exp = float(vac.T @ pi_M_mat @ vac)
        assert np.isclose(abs(parity_exp), 1.0, atol=1e-10)

        # Compare with compute_matter_vacuum_parity:
        vac_parity = compute_matter_vacuum_parity(T)
        assert vac_parity == int(np.round(parity_exp)), (
            f"compute_matter_vacuum_parity ({vac_parity}) must match Fock expectation ({parity_exp})"
        )


def test_matter_parity_selection_on_torus_clusters() -> None:
    """Verify physical state projection and ground state energetics on torus clusters."""
    # Cluster 1: 4x4 torus (N=16 unit cells, N1+N2=8 even)
    N1, N2 = 4, 4
    u_ref = build_reference_links(N1, N2)
    T_ref = build_matter_matrix_T(u_ref, N1, N2, K=1.0)
    A_ref = build_antisymmetric_majorana_matrix_A(T_ref)
    assert np.allclose(A_ref, -A_ref.T), "Majorana matrix A must be antisymmetric"

    vac_parity = compute_matter_vacuum_parity(T_ref)
    gauge_sign = compute_torus_gauge_sign(N1, N2, u=u_ref)

    # In 4x4 periodic reference gauge, sign(det T) = -1 while gauge_sign = +1.
    # Therefore, the unprojected Gaussian vacuum has odd parity and is UNPHYSICAL (Pedrocchi et al. 2011).
    assert gauge_sign == 1
    assert vac_parity == -1

    # Verify physical ground state energy is shifted by lowest single-particle mode s_1
    E_phys, vac_is_phys, E_vac = compute_physical_ground_state_energy(T_ref, N1, N2, u=u_ref)
    assert vac_is_phys is False
    s_sorted = np.sort(np.linalg.svd(T_ref, compute_uv=False))
    assert np.isclose(E_phys, E_vac + s_sorted[0])

    # Zero quasiparticles is unphysical, 1 quasiparticle is physical:
    occ_zero = np.zeros(N1 * N2, dtype=np.int64)
    occ_one = np.zeros(N1 * N2, dtype=np.int64)
    occ_one[0] = 1
    assert is_physical_state(T_ref, occ_zero, N1, N2, u=u_ref) is False
    assert is_physical_state(T_ref, occ_one, N1, N2, u=u_ref) is True

    # Under antiperiodic twist (-1, 1), the Gaussian vacuum becomes physical:
    u_twisted = build_reference_links(N1, N2, twists=(-1, 1))
    T_twisted = build_matter_matrix_T(u_twisted, N1, N2, K=1.0)
    vac_parity_tw = compute_matter_vacuum_parity(T_twisted)
    gauge_sign_tw = compute_torus_gauge_sign(N1, N2, u=u_twisted, twists=(-1, 1))
    assert vac_parity_tw == 1
    assert gauge_sign_tw == 1
    E_phys_tw, vac_is_phys_tw, E_vac_tw = compute_physical_ground_state_energy(
        T_twisted, N1, N2, u=u_twisted, twists=(-1, 1)
    )
    assert vac_is_phys_tw is True
    assert np.isclose(E_phys_tw, E_vac_tw)
    assert is_physical_state(T_twisted, occ_zero, N1, N2, u=u_twisted, twists=(-1, 1)) is True

    # Cluster 2: 2x3 torus (N=6 unit cells, N1+N2=5 odd)
    # Here sigma_geom = (-1)^5 = -1. In reference gauge (u=+1), gauge sign is -1.
    N1, N2 = 2, 3
    u_ref_23 = build_reference_links(N1, N2)
    T_ref_23 = build_matter_matrix_T(u_ref_23, N1, N2, K=1.0)
    vac_parity_23 = compute_matter_vacuum_parity(T_ref_23)
    gauge_sign_23 = compute_torus_gauge_sign(N1, N2, u=u_ref_23)

    assert gauge_sign_23 == -1
    # det(T_ref_23) > 0 so vac_parity is +1 != -1.
    assert vac_parity_23 == 1

    E_phys_23, vac_is_phys_23, E_vac_23 = compute_physical_ground_state_energy(
        T_ref_23, N1, N2, u=u_ref_23
    )
    assert vac_is_phys_23 is False
    # Ground state must be shifted by the lowest quasiparticle energy s_1
    s_23 = np.sort(np.linalg.svd(T_ref_23, compute_uv=False))
    assert np.isclose(E_phys_23, E_vac_23 + s_23[0])

    # In this sector, 1 quasiparticle state is physical, zero quasiparticles is unphysical:
    occ_zero_23 = np.zeros(N1 * N2, dtype=np.int64)
    occ_one_23 = np.zeros(N1 * N2, dtype=np.int64)
    occ_one_23[0] = 1
    assert is_physical_state(T_ref_23, occ_zero_23, N1, N2, u=u_ref_23) is False
    assert is_physical_state(T_ref_23, occ_one_23, N1, N2, u=u_ref_23) is True
