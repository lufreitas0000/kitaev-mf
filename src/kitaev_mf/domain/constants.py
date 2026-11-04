"""Fundamental mathematical matrices and algebraic constants for SO(4) Majorana representation."""

import numpy as np
from numpy.typing import NDArray

# 2x2 Pauli matrices
SIGMA_0 = np.eye(2, dtype=np.complex128)
SIGMA_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
SIGMA_Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
SIGMA_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)

# Real 2x2 building blocks
I_2 = np.eye(2, dtype=np.float64)
S_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.float64)
S_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.float64)
EPSILON = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=np.float64)  # i * sigma_y
PI_PLUS = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=np.float64)  # (I + S_Z) / 2

# Explicit real 4x4 SO(4) generators in the Majorana ordering (0, 1, 2, 3).
# Kronecker definitions: G=(I⊗ε, ε⊗Sz, ε⊗Sx), M=(Sz⊗ε, ε⊗I, Sx⊗ε).

G_MATRICES: dict[int, NDArray[np.float64]] = {
    1: np.array(
        [[0.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -1.0, 0.0]], dtype=np.float64
    ),
    2: np.array(
        [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, -1.0],
         [-1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]], dtype=np.float64
    ),
    3: np.array(
        [[0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0, 0.0],
         [0.0, -1.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0]], dtype=np.float64
    ),
}

M_MATRICES: dict[int, NDArray[np.float64]] = {
    1: np.array(
        [[0.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0, -1.0], [0.0, 0.0, 1.0, 0.0]], dtype=np.float64
    ),
    2: np.array(
        [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0],
         [-1.0, 0.0, 0.0, 0.0], [0.0, -1.0, 0.0, 0.0]], dtype=np.float64
    ),
    3: np.array(
        [[0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -1.0, 0.0],
         [0.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0]], dtype=np.float64
    ),
}

N_MATRICES: dict[int, NDArray[np.float64]] = {
    1: np.array(
        [[0.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]], dtype=np.float64
    ),
    2: np.array(
        [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 0.0],
         [-1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]], dtype=np.float64
    ),
    3: np.array(
        [[0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 0.0, 0.0],
         [0.0, 0.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0]], dtype=np.float64
    ),
}
