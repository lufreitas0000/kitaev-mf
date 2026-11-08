"""Finite chiral-matter subspace transport; not the four-flavor finite-field MF model.

H=(i/4)c^T A c on a periodic honeycomb graph with NN K and prescribed
NNN kappa. Transport is defined only for an isolated central rank-two spectral
subspace with nonsingular successive overlaps and exact site-gauge closure.
The returned one-particle holonomy is not a projected many-body braid matrix.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ChiralVortexTransportResult:
    """Discrete geometric holonomy and diagnostics of an isolated central subspace."""

    enclosing: bool
    loop_bonds: tuple[tuple[int, int, int], ...]
    enclosed_sites: frozenset[int]
    holonomy_matrix: np.ndarray
    eigenvalue_phases: np.ndarray
    berry_phase: float
    gauge_closure_error: float
    min_subgap_splitting: float
    min_bulk_gap: float
    min_subspace_gap: float
    max_subgap_splitting: float
    min_overlap_singular_value: float


def build_dual_graph(
    L: int,
) -> tuple[
    dict[tuple[int, int], list[tuple[tuple[int, int], tuple[int, int, int]]]],
    dict[tuple[int, int, int], list[tuple[int, int]]],
]:
    """Build dual triangular lattice adjacency and bond-to-plaquette incidence map."""
    bond_to_plaquettes: dict[tuple[int, int, int], list[tuple[int, int]]] = {}
    for m2 in range(L):
        for m1 in range(L):
            p = (m1, m2)
            bonds = [
                (m1, m2, 0),
                (m1, (m2 + 1) % L, 1),
                (m1, (m2 + 1) % L, 2),
                ((m1 - 1) % L, (m2 + 1) % L, 0),
                ((m1 - 1) % L, (m2 + 1) % L, 1),
                (m1, m2, 2),
            ]
            for b in bonds:
                bond_to_plaquettes.setdefault(b, []).append(p)

    adj: dict[tuple[int, int], list[tuple[tuple[int, int], tuple[int, int, int]]]] = {
        (m1, m2): [] for m2 in range(L) for m1 in range(L)
    }
    for b, pls in bond_to_plaquettes.items():
        pA, pB = pls[0], pls[1]
        adj[pA].append((pB, b))
        adj[pB].append((pA, b))
    return adj, bond_to_plaquettes


def generate_vortex_pair_cut(
    L: int,
    p1: tuple[int, int],
    p2: tuple[int, int],
) -> list[tuple[int, int, int]]:
    """Shortest dual-lattice branch cut connecting plaquette p1 to plaquette p2."""
    if p1 == p2:
        return []
    adj, _ = build_dual_graph(L)
    visited = {p1}
    queue: deque[tuple[tuple[int, int], list[tuple[int, int, int]]]] = deque([(p1, [])])
    while queue:
        curr, path = queue.popleft()
        if curr == p2:
            return path
        for nbr, b in adj[curr]:
            if nbr not in visited:
                visited.add(nbr)
                queue.append((nbr, path + [b]))
    return []


def bond_endpoints(L: int, m1: int, m2: int, b_idx: int) -> tuple[int, int]:
    """Return physical site indices (s_A, s_B) connected by bond (m1, m2, b_idx)."""
    disp = [(0, 0), (0, -1), (-1, 0)]
    s_A = 2 * ((m1 % L) + L * (m2 % L))
    d1, d2 = disp[b_idx]
    s_B = 2 * (((m1 + d1) % L) + L * ((m2 + d2) % L)) + 1
    return s_A, s_B


def plaquette_sites(L: int, m1: int, m2: int) -> set[int]:
    """Return the 6 vertex site indices surrounding hexagonal plaquette (m1, m2)."""
    bonds = [
        (m1, m2, 0),
        (m1, (m2 + 1) % L, 1),
        (m1, (m2 + 1) % L, 2),
        ((m1 - 1) % L, (m2 + 1) % L, 0),
        ((m1 - 1) % L, (m2 + 1) % L, 1),
        (m1, m2, 2),
    ]
    pts: set[int] = set()
    for b in bonds:
        s_A, s_B = bond_endpoints(L, *b)
        pts.add(s_A)
        pts.add(s_B)
    return pts


def build_homology_braiding_loop(
    L: int,
    center_plaquette: tuple[int, int],
    radius: int = 1,
) -> tuple[list[tuple[int, int, int]], frozenset[int], tuple[int, int]]:
    """Construct a closed dual-lattice bond loop delta(R) surrounding center_plaquette.

    Returns:
        Tuple of (ordered_loop_bonds, enclosed_site_set_R, start_plaquette_on_loop).
    """
    if L < 4 or radius < 0 or 2 * radius + 3 > L:
        raise ValueError("Loop region must fit in an L >= 4 torus")
    _, bond_to_plaq = build_dual_graph(L)
    c1, c2 = center_plaquette
    r_sites: set[int] = set()
    for dm1 in range(-radius, radius + 1):
        for dm2 in range(-radius, radius + 1):
            r_sites |= plaquette_sites(L, (c1 + dm1) % L, (c2 + dm2) % L)

    delta_r: list[tuple[int, int, int]] = []
    for b in bond_to_plaq:
        s_A, s_B = bond_endpoints(L, *b)
        if (s_A in r_sites) ^ (s_B in r_sites):
            delta_r.append(b)

    adj_loop: dict[tuple[int, int], list[tuple[tuple[int, int], tuple[int, int, int]]]] = {}
    for b in delta_r:
        pA, pB = bond_to_plaq[b]
        adj_loop.setdefault(pA, []).append((pB, b))
        adj_loop.setdefault(pB, []).append((pA, b))

    if not adj_loop or any(len(edges) != 2 for edges in adj_loop.values()):
        raise ValueError("Region boundary is not one simple dual cycle")
    p_start = sorted(adj_loop.keys())[0]
    ordered_bonds: list[tuple[int, int, int]] = []
    curr_p = p_start
    visited_bonds: set[tuple[int, int, int]] = set()
    while len(ordered_bonds) < len(delta_r):
        for nxt_p, b in adj_loop[curr_p]:
            if b not in visited_bonds:
                visited_bonds.add(b)
                ordered_bonds.append(b)
                curr_p = nxt_p
                break
        else:
            raise ValueError("Region boundary has disconnected cycles")
    if curr_p != p_start:
        raise ValueError("Region boundary is not closed")

    return ordered_bonds, frozenset(r_sites), p_start


def assemble_chiral_bdg_generator(
    L: int,
    K: float = 1.0,
    kappa: float = 0.12,
    bond_signs: dict[tuple[int, int, int], float] | None = None,
) -> np.ndarray:
    """Assemble gauge-covariant real skew-symmetric BdG generator A({u}) with Peierls links."""
    n_sites = 2 * L * L
    A = np.zeros((n_sites, n_sites), dtype=np.float64)
    signs = bond_signs if bond_signs is not None else {}

    def get_u(m1: int, m2: int, b_idx: int) -> float:
        return float(signs.get((m1 % L, m2 % L, b_idx), 1.0))

    disp = [(0, 0), (0, -1), (-1, 0)]
    for m2 in range(L):
        for m1 in range(L):
            s_A = 2 * (m1 + L * m2)
            neighbors_B: list[tuple[int, int, float]] = []
            for b_idx, (d1, d2) in enumerate(disp):
                n1 = (m1 + d1) % L
                n2 = (m2 + d2) % L
                s_B = 2 * (n1 + L * n2) + 1
                u_b = get_u(m1, m2, b_idx)
                A[s_A, s_B] += -2.0 * K * u_b
                A[s_B, s_A] -= -2.0 * K * u_b
                neighbors_B.append((b_idx, s_B, u_b))

            for idx1, idx2 in ((0, 1), (1, 2), (2, 0)):
                _, s_B1, u1 = neighbors_B[idx1]
                _, s_B2, u2 = neighbors_B[idx2]
                A[s_B1, s_B2] += -2.0 * kappa * (u1 * u2)
                A[s_B2, s_B1] -= -2.0 * kappa * (u1 * u2)

            neighbors_A = [
                (0, 2 * ((m1 % L) + L * (m2 % L)), get_u(m1, m2, 0)),
                (1, 2 * ((m1 % L) + L * ((m2 + 1) % L)), get_u(m1, m2 + 1, 1)),
                (2, 2 * (((m1 + 1) % L) + L * (m2 % L)), get_u(m1 + 1, m2, 2)),
            ]
            for idx1, idx2 in ((0, 1), (1, 2), (2, 0)):
                _, s_A1, u1 = neighbors_A[idx1]
                _, s_A2, u2 = neighbors_A[idx2]
                A[s_A1, s_A2] += -2.0 * kappa * (u1 * u2)
                A[s_A2, s_A1] -= -2.0 * kappa * (u1 * u2)

    return A


def _ordered_plaquettes(
    L: int, bonds: list[tuple[int, int, int]], start: tuple[int, int]
) -> list[tuple[int, int]]:
    """Validate a simple closed dual cycle without assuming its edge ordering."""
    _, incidence = build_dual_graph(L)
    if not bonds or len(set(bonds)) != len(bonds):
        raise ValueError("Loop must contain distinct bonds")
    path = [start]
    for bond in bonds:
        pair = incidence.get(bond, [])
        if len(pair) != 2 or path[-1] not in pair:
            raise ValueError("Bonds do not form an ordered dual path")
        path.append(pair[1] if pair[0] == path[-1] else pair[0])
    if path[-1] != start or len(set(path[:-1])) != len(path) - 1:
        raise ValueError("Loop must be a simple closed dual cycle")
    return path


def compute_vortex_braiding_holonomy(
    L: int = 10,
    p_fixed: tuple[int, int] = (4, 4),
    loop_center: tuple[int, int] = (4, 4),
    loop_radius: int = 1,
    K: float = 1.0,
    kappa: float = 0.12,
    n_interp_steps: int = 6,
    custom_loop_bonds: Sequence[tuple[int, int, int]] | None = None,
    custom_enclosed_sites: frozenset[int] | None = None,
    min_overlap: float = 0.1,
) -> ChiralVortexTransportResult:
    """Transport a rank-two spectral subspace through a closed bond interpolation.

    ``min_overlap`` is an explicit discretization admission threshold: refine
    the path if it fails. It is not a universal physical tolerance. Isolation
    is checked relative to matrix scale at every sample; no continuous-path
    isolation theorem or physical parity projection is asserted.
    """
    if L < 4 or n_interp_steps < 2:
        raise ValueError("Require L >= 4 and at least two interpolation steps")
    if not np.isfinite(K) or K <= 0 or not np.isfinite(kappa) or kappa == 0:
        raise ValueError("Require finite K > 0 and finite nonzero kappa")
    if not 0 < min_overlap < 1:
        raise ValueError("min_overlap must lie strictly between zero and one")
    for point in (p_fixed, loop_center):
        if any(v < 0 or v >= L for v in point):
            raise ValueError("Plaquette coordinates must be in [0, L)")
    if (custom_loop_bonds is None) != (custom_enclosed_sites is None):
        raise ValueError("Supply both a custom loop and its enclosed sites")
    if custom_loop_bonds is None:
        ordered_bonds, r_sites, p_start = build_homology_braiding_loop(
            L, loop_center, radius=loop_radius
        )
    else:
        ordered_bonds = list(custom_loop_bonds)
        r_sites = frozenset(custom_enclosed_sites or ())
        if not ordered_bonds or any(s < 0 or s >= 2 * L * L for s in r_sites):
            raise ValueError("Invalid custom loop or enclosed site")
        _, incidence = build_dual_graph(L)
        first = incidence.get(ordered_bonds[0], [])
        if not first:
            raise ValueError("Invalid loop bond")
        # Either orientation may be supplied; exactly one must trace the cycle.
        valid = []
        for candidate in first:
            try:
                _ordered_plaquettes(L, ordered_bonds, candidate)
                valid.append(candidate)
            except ValueError:
                continue
        if not valid:
            raise ValueError("Custom bonds do not form an ordered closed loop")
        p_start = valid[0]
    path = _ordered_plaquettes(L, ordered_bonds, p_start)
    if p_fixed in path:
        raise ValueError("Stationary vortex must not lie on the transport path")
    _, incidence = build_dual_graph(L)
    boundary = {
        b
        for b in incidence
        if (bond_endpoints(L, *b)[0] in r_sites) != (bond_endpoints(L, *b)[1] in r_sites)
    }
    if set(ordered_bonds) != boundary:
        raise ValueError("Loop must equal the site-region gauge boundary")
    cut_init = generate_vortex_pair_cut(L, p_fixed, p_start)
    curr_signs: dict[tuple[int, int, int], float] = dict.fromkeys(cut_init, -1.0)
    A_init = assemble_chiral_bdg_generator(L, K=K, kappa=kappa, bond_signs=curr_signs)
    mid = A_init.shape[0] // 2

    def frame(A: np.ndarray) -> tuple[np.ndarray, float, float, float]:
        values, vectors = np.linalg.eigh(1j * A)
        splitting = float(np.max(np.abs(values[mid - 1 : mid + 1])))
        bulk = float(min(abs(values[mid - 2]), abs(values[mid + 1])))
        separation = bulk - splitting
        scale = float(np.max(np.abs(values)))
        if separation <= 1e-10 * scale:
            raise ValueError("Central rank-two subspace is not spectrally isolated")
        return vectors[:, mid - 1 : mid + 1], splitting, bulk, separation

    V_init, initial_split, initial_bulk, initial_gap = frame(A_init)
    V_curr = V_init.copy()
    min_subgap, max_subgap = initial_split, initial_split
    min_bulk, min_gap, min_singular = initial_bulk, initial_gap, 1.0
    for bond in ordered_bonds:
        u_start = curr_signs.get(bond, 1.0)
        for step in range(1, n_interp_steps + 1):
            curr_signs[bond] = float(u_start * np.cos(np.pi * step / n_interp_steps))
            A_step = assemble_chiral_bdg_generator(L, K=K, kappa=kappa, bond_signs=curr_signs)
            V_next, splitting, bulk, gap = frame(A_step)
            min_subgap, max_subgap = min(min_subgap, splitting), max(max_subgap, splitting)
            min_bulk, min_gap = min(min_bulk, bulk), min(min_gap, gap)
            left, singular, right = np.linalg.svd(V_next.conj().T @ V_curr)
            min_singular = min(min_singular, float(singular[-1]))
            if singular[-1] < min_overlap:
                raise ValueError("Transport overlap is singular or unresolved; refine the path")
            V_curr = V_next @ (left @ right)
        curr_signs[bond] = -u_start
    A_final = assemble_chiral_bdg_generator(L, K=K, kappa=kappa, bond_signs=curr_signs)
    gauge = np.array([-1.0 if s in r_sites else 1.0 for s in range(2 * L * L)])
    gauge_error = float(np.linalg.norm(gauge[:, None] * A_final * gauge[None, :] - A_init))
    if gauge_error > 1e-12 * np.linalg.norm(A_init):
        raise ValueError("Final Hamiltonian does not close under the site gauge")
    raw = V_init.conj().T @ (gauge[:, None] * V_curr)
    if np.linalg.norm(raw.conj().T @ raw - np.eye(2)) > 1e-8:
        raise ValueError("Final transported subspace does not close")
    left, _, right = np.linalg.svd(raw)
    unitary = left @ right
    phases = np.angle(np.linalg.eigvals(unitary))
    return ChiralVortexTransportResult(
        enclosing=plaquette_sites(L, *p_fixed).issubset(r_sites),
        loop_bonds=tuple(ordered_bonds),
        enclosed_sites=r_sites,
        holonomy_matrix=unitary,
        eigenvalue_phases=phases,
        berry_phase=float(np.mean(np.abs(phases))),
        gauge_closure_error=gauge_error,
        min_subgap_splitting=min_subgap,
        min_bulk_gap=min_bulk,
        min_subspace_gap=min_gap,
        max_subgap_splitting=max_subgap,
        min_overlap_singular_value=min_singular,
    )
