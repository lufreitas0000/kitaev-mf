"""Fundamental mathematical matrices and algebraic constants for SO(4) Majorana representation."""
import numpy as np

# Pauli matrices
SIGMA_0 = np.eye(2, dtype=np.complex128)
SIGMA_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
SIGMA_Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
SIGMA_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)

# 4x4 SO(4) gauge constraint matrices G^gamma (anti-self-dual 't Hooft symbols)
# Generating [G^a, G^b] = -2 * epsilon^{abc} * G^c
# T_i^gamma = (i / 4) * c_i^T * G^gamma * c_i
G_MATRICES = {
    1: np.array([
        [ 0.0,  1.0,  0.0,  0.0],
        [-1.0,  0.0,  0.0,  0.0],
        [ 0.0,  0.0,  0.0,  1.0],
        [ 0.0,  0.0, -1.0,  0.0]
    ], dtype=np.float64),
    2: np.array([
        [ 0.0,  0.0,  1.0,  0.0],
        [ 0.0,  0.0,  0.0, -1.0],
        [-1.0,  0.0,  0.0,  0.0],
        [ 0.0,  1.0,  0.0,  0.0]
    ], dtype=np.float64),
    3: np.array([
        [ 0.0,  0.0,  0.0,  1.0],
        [ 0.0,  0.0,  1.0,  0.0],
        [ 0.0, -1.0,  0.0,  0.0],
        [-1.0,  0.0,  0.0,  0.0]
    ], dtype=np.float64)
}

# 4x4 SO(4) chiral physical spin generators M^gamma (self-dual 't Hooft symbols)
# Generating [M^a, M^b] = +2 * epsilon^{abc} * M^c, and [M^a, G^b] = 0
# S_i^gamma = (i / 4) * c_i^T * M^gamma * c_i
M_MATRICES = {
    1: np.array([
        [ 0.0,  1.0,  0.0,  0.0],
        [-1.0,  0.0,  0.0,  0.0],
        [ 0.0,  0.0,  0.0, -1.0],
        [ 0.0,  0.0,  1.0,  0.0]
    ], dtype=np.float64),
    2: np.array([
        [ 0.0,  0.0,  1.0,  0.0],
        [ 0.0,  0.0,  0.0,  1.0],
        [-1.0,  0.0,  0.0,  0.0],
        [ 0.0, -1.0,  0.0,  0.0]
    ], dtype=np.float64),
    3: np.array([
        [ 0.0,  0.0,  0.0,  1.0],
        [ 0.0,  0.0, -1.0,  0.0],
        [ 0.0,  1.0,  0.0,  0.0],
        [-1.0,  0.0,  0.0,  0.0]
    ], dtype=np.float64)
}

# Kitaev representation matrices N^gamma = M^gamma - G^gamma
# S_i^gamma(Kitaev) = (i / 4) * c_i^T * N^gamma * c_i
N_MATRICES = {
    gamma: M_MATRICES[gamma] - G_MATRICES[gamma] for gamma in (1, 2, 3)
}
