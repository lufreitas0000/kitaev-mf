"""Honeycomb lattice and multi-sublattice unit cell geometry."""
from dataclasses import dataclass
from typing import Dict, List, Tuple
import numpy as np


@dataclass(frozen=True)
class Bond:
    """Directed bond on the lattice."""
    sublattice_i: int
    sublattice_j: int
    cell_offset: Tuple[int, int]
    bond_type: int  # 1: x-bond, 2: y-bond, 3: z-bond


@dataclass(frozen=True)
class UnitCell:
    """Unit cell geometry containing basis vectors and internal bonds."""
    name: str
    num_sublattices: int
    a1: np.ndarray
    a2: np.ndarray
    basis_positions: np.ndarray  # Shape: (num_sublattices, 2)
    bonds: List[Bond]

    @property
    def reciprocal_vectors(self) -> Tuple[np.ndarray, np.ndarray]:
        """Compute reciprocal lattice vectors b1, b2 satisfying ai . bj = 2*pi*delta_ij."""
        a_mat = np.column_stack([self.a1, self.a2])
        b_mat = 2.0 * np.pi * np.linalg.inv(a_mat).T
        return b_mat[:, 0], b_mat[:, 1]


def make_honeycomb_2site() -> UnitCell:
    """Standard 2-sublattice primitive honeycomb unit cell.
    
    Sublattice 0: A-site at (0, 0)
    Sublattice 1: B-site at delta_z = (0, 1/sqrt(3))
    """
    d_nn = 1.0 / np.sqrt(3.0)
    a1 = np.array([0.5, np.sqrt(3.0) / 2.0], dtype=np.float64)
    a2 = np.array([-0.5, np.sqrt(3.0) / 2.0], dtype=np.float64)
    
    positions = np.array([
        [0.0, 0.0],
        [0.0, d_nn]
    ], dtype=np.float64)
    
    # Directed bonds from A to B:
    # z-bond: cell (0, 0)
    # x-bond: cell (1, 0)
    # y-bond: cell (0, 1)
    bonds = [
        Bond(sublattice_i=0, sublattice_j=1, cell_offset=(0, 0), bond_type=3),
        Bond(sublattice_i=0, sublattice_j=1, cell_offset=(1, 0), bond_type=1),
        Bond(sublattice_i=0, sublattice_j=1, cell_offset=(0, 1), bond_type=2),
    ]
    
    return UnitCell("primitive_2site", 2, a1, a2, positions, bonds)


def make_honeycomb_4site_rectangular() -> UnitCell:
    """4-sublattice rectangular unit cell accommodating Zigzag and Stripy collinear orders.
    
    A1: (0, 0)
    B1: (0, 1)
    A2: (sqrt(3)/2, 3/2)
    B2: (sqrt(3)/2, 5/2)
    
    Vectors: A1 = (sqrt(3), 0), A2 = (0, 3)
    """
    d_nn = 1.0
    a1 = np.array([np.sqrt(3.0) * d_nn, 0.0], dtype=np.float64)
    a2 = np.array([0.0, 3.0 * d_nn], dtype=np.float64)
    
    positions = np.array([
        [0.0, 0.0],
        [0.0, d_nn],
        [np.sqrt(3.0) / 2.0 * d_nn, 1.5 * d_nn],
        [np.sqrt(3.0) / 2.0 * d_nn, 2.5 * d_nn]
    ], dtype=np.float64)
    
    bonds = [
        # z-bonds (internal): 0-1 and 2-3
        Bond(sublattice_i=0, sublattice_j=1, cell_offset=(0, 0), bond_type=3),
        Bond(sublattice_i=2, sublattice_j=3, cell_offset=(0, 0), bond_type=3),
        
        # x-bonds: 2-1 (internal) and 0-3 (offset -a2)
        Bond(sublattice_i=2, sublattice_j=1, cell_offset=(0, 0), bond_type=1),
        Bond(sublattice_i=0, sublattice_j=3, cell_offset=(0, -1), bond_type=1),
        
        # y-bonds: 2-1 (offset a1) and 0-3 (offset a1 - a2)
        Bond(sublattice_i=2, sublattice_j=1, cell_offset=(1, 0), bond_type=2),
        Bond(sublattice_i=0, sublattice_j=3, cell_offset=(1, -1), bond_type=2),
    ]
    
    return UnitCell("rectangular_4site", 4, a1, a2, positions, bonds)
