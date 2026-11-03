"""Fundamental mathematical matrices and algebraic constants for SO(4) Majorana representation."""
import numpy as np

# Pauli matrices
SIGMA_0 = np.eye(2, dtype=np.complex128)
SIGMA_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
SIGMA_Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
SIGMA_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)

# Canonical Kitaev matrices: sigma_i^gamma = (i / 2) c_i^T N^gamma c_i = i c_i^0 c_i^gamma
# Non-zero entries only at (0, gamma) and (gamma, 0)
N_MATRICES = {
    1: np.array([
        [ 0.0,  1.0,  0.0,  0.0],
        [-1.0,  0.0,  0.0,  0.0],
        [ 0.0,  0.0,  0.0,  0.0],
        [ 0.0,  0.0,  0.0,  0.0]
    ], dtype=np.float64),
    2: np.array([
        [ 0.0,  0.0,  1.0,  0.0],
        [ 0.0,  0.0,  0.0,  0.0],
        [-1.0,  0.0,  0.0,  0.0],
        [ 0.0,  0.0,  0.0,  0.0]
    ], dtype=np.float64),
    3: np.array([
        [ 0.0,  0.0,  0.0,  1.0],
        [ 0.0,  0.0,  0.0,  0.0],
        [ 0.0,  0.0,  0.0,  0.0],
        [-1.0,  0.0,  0.0,  0.0]
    ], dtype=np.float64)
}

# Gauge constraint matrices: c_i^T G^gamma c_i = 0
# G^gamma = N^gamma - N_tilde^gamma, with Tr(G^T G) = 4.0
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

# Physical SO(4) chiral spin generators: M^gamma = 2 * N^gamma - G^gamma
# S_i^gamma = (i / 4) c_i^T M^gamma c_i, with Tr(M^T M) = 4.0
M_MATRICES = {
    gamma: 2.0 * N_MATRICES[gamma] - G_MATRICES[gamma] for gamma in (1, 2, 3)
}
