"""Fundamental mathematical matrices and algebraic constants for SO(4) Majorana representation."""

import numpy as np

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

# 4x4 SO(4) gauge constraint matrices G^gamma (anti-self-dual 't Hooft symbols)
# Algebraic formulation:
#   G^1 = I (x) EPSILON
#   G^2 = EPSILON (x) S_Z
#   G^3 = EPSILON (x) S_X
# Algebra: [G^a, G^b] = -2 * epsilon_{abc} * G^c
G_MATRICES = {
    1: np.kron(I_2, EPSILON),
    2: np.kron(EPSILON, S_Z),
    3: np.kron(EPSILON, S_X),
}

# 4x4 SO(4) chiral physical spin generators M^gamma (self-dual 't Hooft symbols)
# Algebraic formulation:
#   M^1 = S_Z (x) EPSILON
#   M^2 = EPSILON (x) I
#   M^3 = S_X (x) EPSILON
# Algebra: [M^a, M^b] = +2 * epsilon_{abc} * M^c, with [M^a, G^b] = 0
M_MATRICES = {
    1: np.kron(S_Z, EPSILON),
    2: np.kron(EPSILON, I_2),
    3: np.kron(S_X, EPSILON),
}

# 4x4 Kitaev representation matrices N^gamma
# Algebraic formulation: N^gamma = (M^gamma + G^gamma) / 2
#   N^1 = PI_PLUS (x) EPSILON
#   N^2 = EPSILON (x) PI_PLUS
#   N^3 = (S_X (x) EPSILON + EPSILON (x) S_X) / 2
# Property: sigma_i^gamma = (i / 2) c_i^T N^gamma c_i = i c_i^0 c_i^gamma
N_MATRICES = {gamma: 0.5 * (M_MATRICES[gamma] + G_MATRICES[gamma]) for gamma in (1, 2, 3)}
