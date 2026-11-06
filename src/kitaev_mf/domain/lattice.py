"""Honeycomb lattice graph Lambda = (V, E) and multi-sublattice unit cell geometry."""

from dataclasses import dataclass
from typing import cast

import numpy as np

# Standard nearest-neighbor distance for unit Bravais lattice constant a = |a1| = |a2| = 1.0
DEFAULT_D_NN: float = 1.0 / np.sqrt(3.0)


@dataclass(frozen=True)
class Bond:
    """Directed edge <ij>_gamma in E connecting vertex i (sublattice_i) to vertex j (sublattice_j).

    Attributes:
        sublattice_i: Source vertex index i in {0, ..., num_sublattices - 1}.
        sublattice_j: Target vertex index j in {0, ..., num_sublattices - 1}.
        cell_offset: Bravais translation integer tuple (n1, n2) for R_ij = n1 * a1 + n2 * a2.
        bond_type: Directional index gamma in {1: x-bond, 2: y-bond, 3: z-bond}.
    """

    sublattice_i: int
    sublattice_j: int
    cell_offset: tuple[int, int]
    bond_type: int

    @property
    def site_i(self) -> int:
        """Vertex index i of the directed edge <ij>_gamma."""
        return self.sublattice_i

    @property
    def site_j(self) -> int:
        """Vertex index j of the directed edge <ij>_gamma."""
        return self.sublattice_j


@dataclass(frozen=True)
class Plaquette:
    """Hexagonal plaquette face on the honeycomb lattice graph.

    Attributes:
        name: Plaquette identifier (e.g. 'p0', 'p1').
        bond_indices: Indices of the 6 directed bonds forming the boundary of the plaquette.
        oriented_bonds: Optional sequence of (bond_index, orientation) steps where orientation in (+1, -1).
    """

    name: str
    bond_indices: tuple[int, ...]
    oriented_bonds: tuple[tuple[int, int], ...] = ()


@dataclass(frozen=True)
class UnitCell:
    """Unit cell geometry defining a periodic lattice graph Lambda = (V, E)."""

    name: str
    num_sublattices: int
    a1: np.ndarray
    a2: np.ndarray
    basis_positions: np.ndarray  # Shape: (num_sublattices, 2)
    bonds: tuple[Bond, ...]
    plaquettes: tuple[Plaquette, ...] = ()

    @property
    def vertices(self) -> tuple[int, ...]:
        """Ordered vertex index set V = (0, 1, ..., num_sublattices - 1)."""
        return tuple(range(self.num_sublattices))

    @property
    def edges(self) -> tuple[Bond, ...]:
        """Directed edge set E of the unit cell."""
        return self.bonds

    @property
    def reciprocal_vectors(self) -> tuple[np.ndarray, np.ndarray]:
        """Compute reciprocal lattice vectors b1, b2 satisfying a_m . b_n = 2*pi*delta_mn."""
        a_mat = np.column_stack([self.a1, self.a2])
        b_mat = 2.0 * np.pi * np.linalg.inv(a_mat).T
        return b_mat[:, 0], b_mat[:, 1]

    def bond_vector(self, bond: Bond) -> np.ndarray:
        """Cartesian vector delta_ij = r_j + n1 * a1 + n2 * a2 - r_i for a directed bond <ij>_gamma."""
        n1, n2 = bond.cell_offset
        r_i = self.basis_positions[bond.sublattice_i]
        r_j = self.basis_positions[bond.sublattice_j] + n1 * self.a1 + n2 * self.a2
        return cast(np.ndarray, r_j - r_i)

    def make_k_grid(self, n1: int, n2: int) -> np.ndarray:
        """Generate a regular (n1 * n2, 2) grid of fractional reciprocal coordinates in [0, 1)^2."""
        m1 = np.arange(n1, dtype=np.float64) / float(n1)
        m2 = np.arange(n2, dtype=np.float64) / float(n2)
        grid_1, grid_2 = np.meshgrid(m1, m2, indexing="ij")
        return np.column_stack([grid_1.ravel(), grid_2.ravel()])

    def fractional_to_cartesian_k(self, k_frac: np.ndarray) -> np.ndarray:
        """Convert fractional coordinates (..., 2) to physical Cartesian momenta k = k1*b1 + k2*b2."""
        b1, b2 = self.reciprocal_vectors
        b_mat = np.column_stack([b1, b2])
        return np.asarray(k_frac @ b_mat.T, dtype=np.float64)


def make_honeycomb_2site(d_nn: float = DEFAULT_D_NN) -> UnitCell:
    """Standard 2-sublattice primitive honeycomb unit cell with unit lattice spacing a = sqrt(3)*d_nn = 1.

    Vertex 0: A-site at (0, 0)
    Vertex 1: B-site at (0, d_nn)
    Bravais vectors: a1 = (1/2, sqrt(3)/2) * a, a2 = (-1/2, sqrt(3)/2) * a
    """
    a = np.sqrt(3.0) * d_nn
    a1 = np.array([0.5 * a, (np.sqrt(3.0) / 2.0) * a], dtype=np.float64)
    a2 = np.array([-0.5 * a, (np.sqrt(3.0) / 2.0) * a], dtype=np.float64)

    positions = np.array(
        [
            [0.0, 0.0],
            [0.0, d_nn],
        ],
        dtype=np.float64,
    )

    bonds = (
        Bond(sublattice_i=0, sublattice_j=1, cell_offset=(0, 0), bond_type=3),
        Bond(sublattice_i=0, sublattice_j=1, cell_offset=(0, -1), bond_type=1),
        Bond(sublattice_i=0, sublattice_j=1, cell_offset=(-1, 0), bond_type=2),
    )

    plaquettes = (
        Plaquette(
            name="p0",
            bond_indices=(0, 1, 2, 0, 1, 2),
            oriented_bonds=((0, 1), (1, -1), (2, 1), (0, -1), (1, 1), (2, -1)),
        ),
    )

    return UnitCell("primitive_2site", 2, a1, a2, positions, bonds, plaquettes)


def make_honeycomb_4site_rectangular(d_nn: float = DEFAULT_D_NN) -> UnitCell:
    """4-sublattice rectangular unit cell accommodating Zigzag and Stripy orders.

    Uses the identical physical length scale d_nn = 1/sqrt(3) (a = 1) as the 2-site cell:
    A1 = a1 - a2 = (1, 0)
    A2 = a1 + a2 = (0, sqrt(3))

    Vertices:
    i=0 (A1): (0, 0)
    i=1 (B1): (0, d_nn)
    i=2 (A2): (sqrt(3)/2 * d_nn, 1.5 * d_nn)
    i=3 (B2): (sqrt(3)/2 * d_nn, 2.5 * d_nn)
    """
    a1 = np.array([np.sqrt(3.0) * d_nn, 0.0], dtype=np.float64)
    a2 = np.array([0.0, 3.0 * d_nn], dtype=np.float64)

    positions = np.array(
        [
            [0.0, 0.0],
            [0.0, d_nn],
            [(np.sqrt(3.0) / 2.0) * d_nn, 1.5 * d_nn],
            [(np.sqrt(3.0) / 2.0) * d_nn, 2.5 * d_nn],
        ],
        dtype=np.float64,
    )

    bonds = (
        Bond(sublattice_i=0, sublattice_j=1, cell_offset=(0, 0), bond_type=3),
        Bond(sublattice_i=2, sublattice_j=3, cell_offset=(0, 0), bond_type=3),
        Bond(sublattice_i=0, sublattice_j=3, cell_offset=(0, -1), bond_type=1),
        Bond(sublattice_i=2, sublattice_j=1, cell_offset=(1, 0), bond_type=1),
        Bond(sublattice_i=2, sublattice_j=1, cell_offset=(0, 0), bond_type=2),
        Bond(sublattice_i=0, sublattice_j=3, cell_offset=(-1, -1), bond_type=2),
    )

    plaquettes = (
        Plaquette(
            name="p0",
            bond_indices=(0, 3, 1, 2, 4, 5),
            oriented_bonds=((0, 1), (3, -1), (1, 1), (2, -1), (4, 1), (5, -1)),
        ),
        Plaquette(
            name="p1",
            bond_indices=(0, 4, 1, 5, 3, 2),
            oriented_bonds=((0, 1), (4, -1), (1, 1), (5, -1), (3, 1), (2, -1)),
        ),
    )

    return UnitCell("rectangular_4site", 4, a1, a2, positions, bonds, plaquettes)


def make_honeycomb_supercell(n1: int, n2: int, d_nn: float = DEFAULT_D_NN) -> UnitCell:
    """Construct an arbitrary n1 x n2 honeycomb supercell with n_c = 2 * n1 * n2 sublattices.

    Parameters:
        n1: Supercell replication factor along primitive translation a1 (>= 1).
        n2: Supercell replication factor along primitive translation a2 (>= 1).
        d_nn: Nearest-neighbor bond distance (defaults to 1/sqrt(3)).

    Returns:
        UnitCell with n_c = 2 * n1 * n2 sublattices, 3 * n1 * n2 directed bonds,
        and n1 * n2 hexagonal plaquettes, guaranteeing exact Brillouin-zone folding parity.
    """
    if n1 < 1 or n2 < 1:
        raise ValueError(f"Replication factors (n1, n2) must be >= 1, got ({n1}, {n2})")

    a = np.sqrt(3.0) * d_nn
    a1_prim = np.array([0.5 * a, (np.sqrt(3.0) / 2.0) * a], dtype=np.float64)
    a2_prim = np.array([-0.5 * a, (np.sqrt(3.0) / 2.0) * a], dtype=np.float64)
    A1 = float(n1) * a1_prim
    A2 = float(n2) * a2_prim

    num_cells = n1 * n2
    num_sublattices = 2 * num_cells

    positions = np.zeros((num_sublattices, 2), dtype=np.float64)
    for m2 in range(n2):
        for m1 in range(n1):
            cell_origin = m1 * a1_prim + m2 * a2_prim
            i_A = 2 * (m1 + n1 * m2)
            i_B = i_A + 1
            positions[i_A] = cell_origin
            positions[i_B] = cell_origin + np.array([0.0, d_nn], dtype=np.float64)

    bond_map: dict[tuple[int, int, int], int] = {}
    bonds_list: list[Bond] = []
    disp = [(0, 0, 3), (0, -1, 1), (-1, 0, 2)]
    idx = 0
    for m2 in range(n2):
        for m1 in range(n1):
            i_A = 2 * (m1 + n1 * m2)
            for d1, d2, gamma in disp:
                tm1 = m1 + d1
                tm2 = m2 + d2
                r1 = tm1 % n1
                N1 = tm1 // n1
                r2 = tm2 % n2
                N2 = tm2 // n2
                i_B = 2 * (r1 + n1 * r2) + 1
                b = Bond(
                    sublattice_i=i_A,
                    sublattice_j=i_B,
                    cell_offset=(int(N1), int(N2)),
                    bond_type=gamma,
                )
                bonds_list.append(b)
                bond_map[(m1, m2, gamma)] = idx
                idx += 1

    plaquettes: list[Plaquette] = []
    for m2 in range(n2):
        for m1 in range(n1):
            b1 = bond_map[(m1, m2, 3)]
            b2 = bond_map[(m1 % n1, (m2 + 1) % n2, 1)]
            b3 = bond_map[(m1 % n1, (m2 + 1) % n2, 2)]
            b4 = bond_map[((m1 - 1) % n1, (m2 + 1) % n2, 3)]
            b5 = bond_map[((m1 - 1) % n1, (m2 + 1) % n2, 1)]
            b6 = bond_map[(m1, m2, 2)]
            ori = ((b1, 1), (b2, -1), (b3, 1), (b4, -1), (b5, 1), (b6, -1))
            plaq = Plaquette(
                name=f"p_{m1}_{m2}",
                bond_indices=(b1, b2, b3, b4, b5, b6),
                oriented_bonds=ori,
            )
            plaquettes.append(plaq)

    name = f"supercell_{num_sublattices}site_{n1}x{n2}"
    return UnitCell(
        name=name,
        num_sublattices=num_sublattices,
        a1=A1,
        a2=A2,
        basis_positions=positions,
        bonds=tuple(bonds_list),
        plaquettes=tuple(plaquettes),
    )
