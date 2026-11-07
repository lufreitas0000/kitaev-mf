"""Exact matter solver, torus geometry, and physical matter-parity selection rules (Track 4 / B-010, B-011).

Implements exact combinatorial gauge constraints and fermion parity projection on finite
periodic honeycomb tori according to Pedrocchi et al. (2011) and Zschocke & Vojta (2015).
"""

from __future__ import annotations

from collections import deque

import numpy as np
import scipy.integrate as integrate


def coord_to_idx(m: int, n: int, N1: int, N2: int) -> int:
    """Map 2D cell coordinate (m, n) modulo (N1, N2) to a 1D unit-cell index 0 <= j < N1*N2."""
    return (n % N2) * N1 + (m % N1)


def idx_to_coord(idx: int, N1: int, N2: int) -> tuple[int, int]:
    """Map 1D unit-cell index to 2D cell coordinate (m, n)."""
    return idx % N1, (idx // N1) % N2


def build_reference_links(
    N1: int,
    N2: int,
    twists: tuple[int, int] = (1, 1),
) -> np.ndarray:
    """Construct canonical link array u[m, n, alpha] in {+1, -1} with optional boundary twists.

    Ordering of alpha: 0 -> 'x', 1 -> 'y', 2 -> 'z'.
    Twist tau_1 is placed on x-bonds across the m = N1 - 1 seam.
    Twist tau_2 is placed on y-bonds across the n = N2 - 1 seam.
    """
    if twists[0] not in (-1, 1) or twists[1] not in (-1, 1):
        raise ValueError(f"Twists must be in {{+1, -1}}, got {twists}")

    u = np.ones((N1, N2, 3), dtype=np.int64)
    if twists[0] == -1:
        u[N1 - 1, :, 0] = -1
    if twists[1] == -1:
        u[:, N2 - 1, 1] = -1
    return u


def compute_torus_plaquette_fluxes(u: np.ndarray, N1: int, N2: int) -> np.ndarray:
    """Compute all N1*N2 plaquette fluxes W_r according to notes/05_finite_torus_conventions.md.

    Cycle based at A_r:
        A_r -> B_{r+n1} (x) -> A_{r+n1-n2} (y) -> B_{r+n1-n2} (z)
        -> A_{r-n2} (x) -> B_r (y) -> A_r (z).
    Stored link product:
        W_r = u^x_r * u^y_{r+n1-n2} * u^z_{r+n1-n2} * u^x_{r-n2} * u^y_{r-n2} * u^z_r.
    """
    W = np.zeros((N1, N2), dtype=np.int64)
    for m in range(N1):
        for n in range(N2):
            ux1 = u[m, n, 0]
            uy1 = u[(m + 1) % N1, (n - 1) % N2, 1]
            uz1 = u[(m + 1) % N1, (n - 1) % N2, 2]
            ux2 = u[m, (n - 1) % N2, 0]
            uy2 = u[m, (n - 1) % N2, 1]
            uz2 = u[m, n, 2]
            W[m, n] = int(ux1 * uy1 * uz1 * ux2 * uy2 * uz2)
    return W


def compute_torus_wilson_loops(u: np.ndarray, N1: int, N2: int) -> tuple[int, int]:
    """Compute the two primitive noncontractible Wilson loops W1 and W2.

    W1 = prod_{m=0}^{N1-1} u^z_{(m, 0)} u^x_{(m-1, 0)}
    W2 = prod_{n=0}^{N2-1} u^z_{(0, n)} u^y_{(0, n-1)}
    """
    w1 = 1
    for m in range(N1):
        w1 *= int(u[m, 0, 2] * u[(m - 1) % N1, 0, 0])

    w2 = 1
    for n in range(N2):
        w2 *= int(u[0, n, 2] * u[0, (n - 1) % N2, 1])

    return w1, w2


def compute_torus_gauge_sign(
    N1: int,
    N2: int,
    u: np.ndarray | None = None,
    twists: tuple[int, int] = (1, 1),
) -> int:
    """Compute the physical parity eigenvalue sigma_geom * prod_e u_e.

    sigma_geom = (-1)**(N1 + N2).
    """
    if u is None:
        u = build_reference_links(N1, N2, twists=twists)

    prod_u = int(np.prod(u))
    sigma_geom = (-1) ** (N1 + N2)
    return int(sigma_geom * prod_u)


def build_matter_matrix_T(
    u: np.ndarray,
    N1: int,
    N2: int,
    K: float = 1.0,
) -> np.ndarray:
    """Construct N x N bipartite matter coupling matrix T_{r, r'} (N = N1 * N2).

    T_{r, r'} = 2K [ u^z_r delta_{r', r} + u^x_r delta_{r', r+n1} + u^y_r delta_{r', r+n2} ]
    """
    N = N1 * N2
    T = np.zeros((N, N), dtype=np.float64)

    for m in range(N1):
        for n in range(N2):
            idx_r = coord_to_idx(m, n, N1, N2)

            # z bond: target is r' = (m, n)
            idx_z = idx_r
            T[idx_r, idx_z] += 2.0 * K * float(u[m, n, 2])

            # x bond: target is r' = (m+1, n)
            idx_x = coord_to_idx(m + 1, n, N1, N2)
            T[idx_r, idx_x] += 2.0 * K * float(u[m, n, 0])

            # y bond: target is r' = (m, n+1)
            idx_y = coord_to_idx(m, n + 1, N1, N2)
            T[idx_r, idx_y] += 2.0 * K * float(u[m, n, 1])

    return T


def build_antisymmetric_majorana_matrix_A(T: np.ndarray) -> np.ndarray:
    """Assemble the 2N x 2N antisymmetric Majorana matrix A = [[0, T], [-T^T, 0]]."""
    N = T.shape[0]
    A = np.zeros((2 * N, 2 * N), dtype=np.float64)
    A[:N, N:] = T
    A[N:, :N] = -T.T
    return A


def compute_matter_vacuum_parity(T: np.ndarray, tol: float = 1e-10) -> int:
    """Compute the unprojected Gaussian matter vacuum fermion parity pi_M(|0_M; u>).

    By the Thouless theorem, pi_M(|0_M; u>) = det(U) * det(V) = sign(det(T)).
    Returns 0 if T possesses a zero singular value (s_1 < tol).
    """
    s = np.linalg.svd(T, compute_uv=False)
    if float(np.min(s)) < tol:
        return 0

    sign, _ = np.linalg.slogdet(T)
    return int(sign)


def compute_matter_state_parity(
    T: np.ndarray,
    occupations: np.ndarray,
    tol: float = 1e-10,
) -> int:
    """Compute fermion parity of a state with mode occupations: pi_M * (-1)**(sum n_k)."""
    vac_parity = compute_matter_vacuum_parity(T, tol=tol)
    if vac_parity == 0:
        return 0
    n_tot = int(np.sum(occupations))
    return int(vac_parity * ((-1) ** n_tot))


def is_physical_state(
    T: np.ndarray,
    occupations: np.ndarray,
    N1: int,
    N2: int,
    u: np.ndarray | None = None,
    twists: tuple[int, int] = (1, 1),
    tol: float = 1e-10,
) -> bool:
    """Determine whether a given quasiparticle Fock state satisfies the physical parity rule."""
    state_parity = compute_matter_state_parity(T, occupations, tol=tol)
    if state_parity == 0:
        # Zero-mode subspace contains both even and odd physical states
        return True

    gauge_sign = compute_torus_gauge_sign(N1, N2, u=u, twists=twists)
    return state_parity == gauge_sign


def compute_physical_ground_state_energy(
    T: np.ndarray,
    N1: int,
    N2: int,
    u: np.ndarray | None = None,
    twists: tuple[int, int] = (1, 1),
    tol: float = 1e-10,
) -> tuple[float, bool, float]:
    """Compute the true physical ground state energy in the given gauge sector.

    Returns:
        (E_phys, vacuum_is_physical, E_vac)
        where:
            E_vac = -0.5 * sum_k s_k(T)
            If vacuum is physical or s_1 < tol, E_phys = E_vac.
            If vacuum is unphysical and s_1 >= tol, E_phys = E_vac + s_1.
    """
    s = np.linalg.svd(T, compute_uv=False)
    s_sorted = np.sort(s)
    E_vac = float(-0.5 * np.sum(s_sorted))

    vac_parity = compute_matter_vacuum_parity(T, tol=tol)
    gauge_sign = compute_torus_gauge_sign(N1, N2, u=u, twists=twists)

    if vac_parity == 0 or s_sorted[0] < tol:
        # Zero mode present: exact 2-fold degeneracy between even/odd Fock states
        return E_vac, True, E_vac

    if vac_parity == gauge_sign:
        return E_vac, True, E_vac
    else:
        # Vacuum has wrong parity; project onto 1-particle state
        E_phys = float(E_vac + s_sorted[0])
        return E_phys, False, E_vac


def compute_thermodynamic_vacuum_quadrature(
    K: float = 1.0,
    epsabs: float = 1e-11,
) -> tuple[float, float]:
    """Compute the thermodynamic vacuum energy and Fock width density by 2D BZ quadrature.

    Returns:
        (e_vac, w_fock)
        where:
            e_vac = -K / (4 pi^2) int int sqrt(3 + 2 cos t1 + 2 cos t2 + 2 cos(t1 - t2)) dt1 dt2
            w_fock = -2 * e_vac
    """

    def integrand(t1: float, t2: float) -> float:
        val = 3.0 + 2.0 * np.cos(t1) + 2.0 * np.cos(t2) + 2.0 * np.cos(t1 - t2)
        return float(np.sqrt(max(val, 0.0)))

    res, _ = integrate.dblquad(
        integrand,
        -np.pi,
        np.pi,
        -np.pi,
        np.pi,
        epsabs=epsabs,
    )
    e_vac = float(-res * K / (4.0 * np.pi**2))
    w_fock = float(-2.0 * e_vac)
    return e_vac, w_fock


def get_plaquette_neighbors(
    m: int,
    n: int,
    N1: int,
    N2: int,
) -> list[tuple[int, int, int, int, int]]:
    """Return the 6 dual triangular lattice neighbors of plaquette (m, n).

    Each neighbor is represented as (nxt_m, nxt_n, bond_m, bond_n, alpha).
    Flipping the bond (bond_m, bond_n, alpha) flips plaquettes (m, n) and (nxt_m, nxt_n).
    """
    neighbors: list[tuple[int, int, int, int, int]] = []
    # 1. (m, n+1) via x-bond at (m, n)
    neighbors.append(((m) % N1, (n + 1) % N2, m, n, 0))
    # 2. (m, n-1) via x-bond at (m, n-1)
    neighbors.append(((m) % N1, (n - 1) % N2, m, (n - 1) % N2, 0))
    # 3. (m-1, n+1) via z-bond at (m, n)
    neighbors.append(((m - 1) % N1, (n + 1) % N2, m, n, 2))
    # 4. (m+1, n-1) via z-bond at (m+1, n-1)
    neighbors.append(((m + 1) % N1, (n - 1) % N2, (m + 1) % N1, (n - 1) % N2, 2))
    # 5. (m+1, n) via y-bond at (m+1, n-1)
    neighbors.append(((m + 1) % N1, (n) % N2, (m + 1) % N1, (n - 1) % N2, 1))
    # 6. (m-1, n) via y-bond at (m, n-1)
    neighbors.append(((m - 1) % N1, (n) % N2, m, (n - 1) % N2, 1))
    return neighbors


def find_dual_path(
    p1: tuple[int, int],
    p2: tuple[int, int],
    N1: int,
    N2: int,
) -> list[tuple[int, int, int]]:
    """Find a shortest sequence of bonds to flip along the dual lattice between plaquettes p1 and p2."""
    if p1 == p2:
        return []

    queue: deque[tuple[int, int]] = deque([p1])
    parent: dict[tuple[int, int], tuple[int, int] | None] = {p1: None}
    bond_used: dict[tuple[int, int], tuple[int, int, int]] = {}

    while queue:
        curr = queue.popleft()
        if curr == p2:
            break
        for nxt_m, nxt_n, bm, bn, alpha in get_plaquette_neighbors(curr[0], curr[1], N1, N2):
            nxt = (nxt_m, nxt_n)
            if nxt not in parent:
                parent[nxt] = curr
                bond_used[nxt] = (bm, bn, alpha)
                queue.append(nxt)

    curr = p2
    bonds: list[tuple[int, int, int]] = []
    while curr != p1:
        bonds.append(bond_used[curr])
        p = parent[curr]
        assert p is not None
        curr = p

    return bonds


def create_vison_pair_links(
    N1: int,
    N2: int,
    p1: tuple[int, int],
    p2: tuple[int, int],
    twists: tuple[int, int] = (1, 1),
) -> tuple[np.ndarray, list[tuple[int, int, int]]]:
    """Construct link array creating an isolated vison pair at plaquettes p1 and p2."""
    u = build_reference_links(N1, N2, twists=twists)
    bonds = find_dual_path(p1, p2, N1, N2)
    for bm, bn, alpha in bonds:
        u[bm, bn, alpha] *= -1
    return u, bonds


def evaluate_vison_pair_energy(
    N1: int,
    N2: int,
    p1: tuple[int, int],
    p2: tuple[int, int],
    K: float = 1.0,
    twists: tuple[int, int] = (1, 1),
) -> dict[str, float]:
    """Compute physical and vacuum vison pair creation energy relative to the flux-free sector."""
    # Reference sector
    u_0 = build_reference_links(N1, N2, twists=twists)
    T_0 = build_matter_matrix_T(u_0, N1, N2, K=K)
    E_phys_0, is_phys_0, E_vac_0 = compute_physical_ground_state_energy(
        T_0, N1, N2, u=u_0, twists=twists
    )

    # Vison pair sector
    u_pair, _ = create_vison_pair_links(N1, N2, p1, p2, twists=twists)
    T_pair = build_matter_matrix_T(u_pair, N1, N2, K=K)
    E_phys_pair, is_phys_pair, E_vac_pair = compute_physical_ground_state_energy(
        T_pair, N1, N2, u=u_pair, twists=twists
    )

    return {
        "dE_phys": float(E_phys_pair - E_phys_0),
        "dE_vac": float(E_vac_pair - E_vac_0),
        "E_phys_pair": float(E_phys_pair),
        "E_vac_pair": float(E_vac_pair),
        "E_phys_0": float(E_phys_0),
        "E_vac_0": float(E_vac_0),
        "is_phys_pair": 1.0 if is_phys_pair else 0.0,
        "is_phys_0": 1.0 if is_phys_0 else 0.0,
    }


def extrapolate_vison_gap_infinite_separation(
    L_list: list[int],
    K: float = 1.0,
) -> dict[str, float]:
    """Perform 1/L linear extrapolation of maximum-separation vison pair energy to thermodynamic limit.

    For each torus size L, computes the true global physical ground state across all 4 boundary twist
    sectors for both the flux-free and the 2-vison sector, yielding the physical pair creation energy
    dE_pair(L) = E_pair(L) - E_0(L).
    """
    inv_L = [1.0 / L for L in L_list]
    dE_list: list[float] = []

    for L in L_list:
        # 1. Global flux-free ground state over all twists
        E0_global = 1e9
        for tw in [(1, 1), (-1, 1), (1, -1), (-1, -1)]:
            u = build_reference_links(L, L, twists=tw)
            T = build_matter_matrix_T(u, L, L, K=K)
            E_phys, _, _ = compute_physical_ground_state_energy(T, L, L, u=u, twists=tw)
            if E_phys < E0_global:
                E0_global = E_phys

        # 2. Global pair ground state over all twists
        p1 = (0, 0)
        p2 = (L // 2, L // 2)
        E_pair_global = 1e9
        for tw in [(1, 1), (-1, 1), (1, -1), (-1, -1)]:
            u_pair, _ = create_vison_pair_links(L, L, p1, p2, twists=tw)
            T_pair = build_matter_matrix_T(u_pair, L, L, K=K)
            E_phys, _, _ = compute_physical_ground_state_energy(T_pair, L, L, u=u_pair, twists=tw)
            if E_phys < E_pair_global:
                E_pair_global = E_phys

        dE_list.append(float(E_pair_global - E0_global))

    # Linear fit: dE = dE_inf + c / L
    poly = np.polyfit(inv_L, dE_list, 1)
    dE_inf = float(poly[1])
    slope = float(poly[0])

    return {
        "dE_inf": dE_inf,
        "slope": slope,
        "mean_dE": float(np.mean(dE_list)),
    }
