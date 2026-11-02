"""Fundamental mathematical matrices and algebraic constants."""
import numpy as np

# Pauli matrices
SIGMA_0 = np.eye(2, dtype=np.complex128)
SIGMA_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
SIGMA_Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
SIGMA_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)

# 4x4 SO(4) Antisymmetric basis matrices for Majorana representation
# N^gamma matrices: sigma_i^gamma = (i/2) c_i^T N^gamma c_i
N_MATRICES = {
    1: np.array([
        [0.0,  1.0,  0.0,  0.0],
        [-1.0, 0.0,  0.0,  0.0],
        [0.0,  0.0,  0.0,  0.0],
        [0.0,  0.0,  0.0,  0.0]
    ], dtype=np.float64),
    2: np.array([
        [0.0,  0.0,  1.0,  0.0],
        [0.0,  0.0,  0.0,  0.0],
        [-1.0, 0.0,  0.0,  0.0],
        [0.0,  0.0,  0.0,  0.0]
    ], dtype=np.float64),
    3: np.array([
        [0.0,  0.0,  0.0,  1.0],
        [0.0,  0.0,  0.0,  0.0],
        [0.0,  0.0,  0.0,  0.0],
        [-1.0, 0.0,  0.0,  0.0]
    ], dtype=np.float64)
}

# G^gamma gauge constraint matrices: c_i^T G^gamma c_i = 0
G_MATRICES = {
    1: np.array([
        [0.0,  1.0,  0.0,  0.0],
        [-1.0, 0.0,  0.0,  0.0],
        [0.0,  0.0,  0.0, -1.0],
        [0.0,  0.0,  1.0,  0.0]
    ], dtype=np.float64),
    2: np.array([
        [0.0,  0.0,  1.0,  0.0],
        [0.0,  0.0,  0.0,  1.0],
        [-1.0, 0.0,  0.0,  0.0],
        [0.0, -1.0,  0.0,  0.0]
    ], dtype=np.float64),
    3: np.array([
        [0.0,  0.0,  0.0,  1.0],
        [0.0,  0.0, -1.0,  0.0],
        [0.0,  1.0,  0.0,  0.0],
        [-1.0, 0.0,  0.0,  0.0]
    ], dtype=np.float64)
}

# Physical SO(4) chiral spin generators: M^gamma = N^gamma + G^gamma
M_MATRICES = {
    gamma: N_MATRICES[gamma] + G_MATRICES[gamma] for gamma in (1, 2, 3)
}
