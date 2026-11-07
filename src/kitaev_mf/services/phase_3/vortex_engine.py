"""Real-space static flux configurations, dual branch cuts, and two-vortex states.

Implements the geometric and gauge-theoretic engine for:
- Generating minimal-length branch cuts on the dual triangular lattice.
- Creating isolated two-vortex backgrounds G(p1, p2) with W_p1 = -1, W_p2 = -1.
- Characterizing Majorana zero mode (MZM) sub-gap state emergence.
- Evaluating local spin-quench transitions: vortex translation vs. pair creation.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from kitaev_mf.services.phase_3.gauge_theory import (
    DirectedBond,
    GaugeConfiguration,
    compute_plaquette_flux,
    flip_bond_gauge_link,
)


@dataclass(frozen=True)
class HoneycombGeometry:
    """Torus geometry and plaquette-bond topological mapping for an L x L unit-cell cluster.

    Each unit cell (m1, m2) in [0, L)^2 contains:
    - 2 physical sites: A (sublattice 0) and B (sublattice 1)
    - 3 directed nearest-neighbor bonds (z, x, y) connecting A to B
    - 1 hexagonal plaquette centered at cell (m1, m2)
    """

    L: int

    def __post_init__(self) -> None:
        if self.L < 2:
            raise ValueError(f"Linear dimension L must be >= 2, got {self.L}")

    @property
    def num_cells(self) -> int:
        """Total unit cells N_cells = L^2."""
        return self.L * self.L

    @property
    def num_sites(self) -> int:
        """Total physical lattice sites N_sites = 2 * L^2."""
        return 2 * self.num_cells

    @property
    def num_bonds(self) -> int:
        """Total directed bonds N_bonds = 3 * L^2."""
        return 3 * self.num_cells

    def site_index(self, m1: int, m2: int, sublattice: int) -> int:
        """Map cell coordinate (m1, m2) and sublattice (0 for A, 1 for B) to site index."""
        r1 = m1 % self.L
        r2 = m2 % self.L
        return 2 * (r1 + self.L * r2) + sublattice

    def site_coordinates(self, site: int) -> tuple[int, int, int]:
        """Return (m1, m2, sublattice) given physical site index."""
        sub = site % 2
        cell = site // 2
        m1 = cell % self.L
        m2 = cell // self.L
        return m1, m2, sub

    def plaquette_bonds(self, m1: int, m2: int) -> list[tuple[int, int]]:
        """Return the 6 directed bonds (site_a, site_b) bounding hexagonal plaquette (m1, m2)."""
        L = self.L
        r1 = m1 % L
        r2 = m2 % L
        return [
            (self.site_index(r1, r2, 0), self.site_index(r1, r2, 1)),  # z at (m1, m2)
            (self.site_index(r1, (r2 + 1) % L, 0), self.site_index(r1, r2, 1)),  # x at (m1, m2+1)
            (
                self.site_index(r1, (r2 + 1) % L, 0),
                self.site_index((r1 - 1) % L, (r2 + 1) % L, 1),
            ),  # y at (m1, m2+1)
            (
                self.site_index((r1 - 1) % L, (r2 + 1) % L, 0),
                self.site_index((r1 - 1) % L, (r2 + 1) % L, 1),
            ),  # z at (m1-1, m2+1)
            (
                self.site_index((r1 - 1) % L, (r2 + 1) % L, 0),
                self.site_index((r1 - 1) % L, r2, 1),
            ),  # x at (m1-1, m2+1)
            (self.site_index(r1, r2, 0), self.site_index((r1 - 1) % L, r2, 1)),  # y at (m1, m2)
        ]

    def build_dual_graph(
        self,
    ) -> dict[tuple[int, int], list[tuple[tuple[int, int], tuple[int, int]]]]:
        """Build dual triangular graph where nodes are plaquettes and edges are shared bonds.

        Returns:
            Mapping from plaquette (m1, m2) to list of (neighbor_plaquette, shared_bond).
        """
        L = self.L
        bond_to_plaquettes: dict[tuple[int, int], list[tuple[int, int]]] = {}
        for m2 in range(L):
            for m1 in range(L):
                p = (m1, m2)
                for b in self.plaquette_bonds(m1, m2):
                    bond_to_plaquettes.setdefault(b, []).append(p)

        adj: dict[tuple[int, int], list[tuple[tuple[int, int], tuple[int, int]]]] = {
            (m1, m2): [] for m2 in range(L) for m1 in range(L)
        }
        for b, pls in bond_to_plaquettes.items():
            if len(pls) == 2:
                pA, pB = pls[0], pls[1]
                adj[pA].append((pB, b))
                adj[pB].append((pA, b))

        return adj

    def generate_branch_cut(
        self,
        p1: tuple[int, int],
        p2: tuple[int, int],
    ) -> list[tuple[int, int]]:
        """Find the minimal-length dual branch cut (string of bonds) connecting p1 and p2 using BFS."""
        if p1 == p2:
            return []

        adj = self.build_dual_graph()
        visited = {p1}
        queue: deque[tuple[tuple[int, int], list[tuple[int, int]]]] = deque([(p1, [])])

        while queue:
            curr, path_bonds = queue.popleft()
            if curr == p2:
                return path_bonds
            for nbr, bond in adj[curr]:
                if nbr not in visited:
                    visited.add(nbr)
                    queue.append((nbr, path_bonds + [bond]))

        return []

    def torus_distance(
        self,
        p1: tuple[int, int],
        p2: tuple[int, int],
    ) -> float:
        """Compute the shortest Euclidean distance between two plaquette centers on the torus."""
        # Canonicalize coordinates to fundamental torus domain [0, L)^2
        L = self.L
        p1_can = (p1[0] % L, p1[1] % L)
        p2_can = (p2[0] % L, p2[1] % L)

        # Hexagonal plaquette centers in unit-cell basis:
        # a1 = [sqrt(3), 0], a2 = [sqrt(3)/2, 3/2]
        a1 = np.array([np.sqrt(3.0), 0.0])
        a2 = np.array([np.sqrt(3.0) / 2.0, 1.5])

        min_dist_sq = float("inf")
        # Check periodic image offsets
        for d1 in (-1, 0, 1):
            for d2 in (-1, 0, 1):
                dm1 = (p2_can[0] - p1_can[0]) + d1 * L
                dm2 = (p2_can[1] - p1_can[1]) + d2 * L
                dr = dm1 * a1 + dm2 * a2
                dist_sq = float(np.dot(dr, dr))
                if dist_sq < min_dist_sq:
                    min_dist_sq = dist_sq

        return float(np.sqrt(min_dist_sq))


def create_two_vortex_gauge_configuration(
    base_gauge: GaugeConfiguration,
    p1: tuple[int, int],
    p2: tuple[int, int],
    geom: HoneycombGeometry,
) -> tuple[GaugeConfiguration, list[tuple[int, int]]]:
    """Create static two-vortex gauge configuration G(p1, p2) by inverting links along a dual branch cut.

    Returns:
        (quenched_gauge, cut_bonds)
    """
    cut_bonds = geom.generate_branch_cut(p1, p2)
    gauge = base_gauge
    for sa, sb in cut_bonds:
        gauge = flip_bond_gauge_link(gauge, site_i=sa, site_j=sb)

    return gauge, cut_bonds


def evaluate_all_plaquette_fluxes(
    gauge: GaugeConfiguration,
    geom: HoneycombGeometry,
) -> dict[tuple[int, int], int]:
    """Compute gauge-invariant Wilson loops W_p = prod_{b in partial p} u_b for all L^2 plaquettes."""
    L = geom.L
    fluxes: dict[tuple[int, int], int] = {}
    for m2 in range(L):
        for m1 in range(L):
            p = (m1, m2)
            bonds = geom.plaquette_bonds(m1, m2)
            fluxes[p] = compute_plaquette_flux(gauge, bonds)

    return fluxes


def identify_vortices(
    gauge: GaugeConfiguration,
    geom: HoneycombGeometry,
) -> list[tuple[int, int]]:
    """Return list of plaquette coordinates with non-trivial flux W_p = -1."""
    fluxes = evaluate_all_plaquette_fluxes(gauge, geom)
    return [p for p, w in fluxes.items() if w == -1]


def apply_local_spin_quench(
    gauge: GaugeConfiguration,
    site_j: int,
    bond_type: int,  # 1 (x), 2 (y), 3 (z)
    geom: HoneycombGeometry,
) -> tuple[GaugeConfiguration, tuple[int, int]]:
    """Simulate action of local spin operator sigma_j^alpha by flipping the corresponding bond gauge link.

    Returns:
        (quenched_gauge, flipped_bond)
    """
    m1, m2, sub = geom.site_coordinates(site_j)
    L = geom.L

    # Find the target bond
    if sub == 0:  # site_j is in sublattice A
        sa = site_j
        if bond_type == 3:  # z-bond
            sb = geom.site_index(m1, m2, 1)
        elif bond_type == 1:  # x-bond
            sb = geom.site_index(m1, m2 - 1, 1)
        elif bond_type == 2:  # y-bond
            sb = geom.site_index(m1 - 1, m2, 1)
        else:
            raise ValueError(f"Invalid bond_type {bond_type}, must be in {{1, 2, 3}}")
    else:  # site_j is in sublattice B
        sb = site_j
        if bond_type == 3:  # z-bond
            sa = geom.site_index(m1, m2, 0)
        elif bond_type == 1:  # x-bond
            sa = geom.site_index(m1, (m2 + 1) % L, 0)
        elif bond_type == 2:  # y-bond
            sa = geom.site_index((m1 + 1) % L, m2, 0)
        else:
            raise ValueError(f"Invalid bond_type {bond_type}, must be in {{1, 2, 3}}")

    quenched_gauge = flip_bond_gauge_link(gauge, site_i=sa, site_j=sb)
    return quenched_gauge, (sa, sb)


def categorize_spin_quench(
    gauge_initial: GaugeConfiguration,
    gauge_quenched: GaugeConfiguration,
    geom: HoneycombGeometry,
) -> dict[str, Any]:
    """Analyze the transition between initial and quenched gauge sectors."""
    vortices_init = set(identify_vortices(gauge_initial, geom))
    vortices_quenched = set(identify_vortices(gauge_quenched, geom))

    num_init = len(vortices_init)
    num_quenched = len(vortices_quenched)
    moved_from = list(vortices_init - vortices_quenched)
    moved_to = list(vortices_quenched - vortices_init)

    if num_init == num_quenched and len(moved_from) == 1 and len(moved_to) == 1:
        transition_type = "vortex_translation"
    elif num_quenched == num_init + 2:
        transition_type = "pair_creation"
        moved_from = []
        moved_to = list(vortices_quenched - vortices_init)
    elif num_quenched == num_init - 2:
        transition_type = "annihilation"
        moved_from = list(vortices_init - vortices_quenched)
        moved_to = []
    else:
        transition_type = "complex_transition"

    return {
        "transition_type": transition_type,
        "initial_vortex_count": num_init,
        "quenched_vortex_count": num_quenched,
        "initial_vortices": sorted(list(vortices_init)),
        "quenched_vortices": sorted(list(vortices_quenched)),
        "displaced_vortices": (moved_from, moved_to),
    }
