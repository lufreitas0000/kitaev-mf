"""Test exact operator relations for SO(4) Majorana and gauge matrices."""

import numpy as np

from kitaev_mf.domain.constants import (
    EPSILON,
    G_MATRICES,
    I_2,
    M_MATRICES,
    N_MATRICES,
    PI_PLUS,
    S_X,
    S_Z,
)


def test_antisymmetry() -> None:
    for gamma in (1, 2, 3):
        assert np.allclose(N_MATRICES[gamma], -N_MATRICES[gamma].T)
        assert np.allclose(G_MATRICES[gamma], -G_MATRICES[gamma].T)
        assert np.allclose(M_MATRICES[gamma], -M_MATRICES[gamma].T)


def test_pauli_kronecker_construction() -> None:
    """Verify algebraic Kronecker definitions against analytical targets."""
    # Gauge matrices G^gamma
    assert np.allclose(G_MATRICES[1], np.kron(I_2, EPSILON))
    assert np.allclose(G_MATRICES[2], np.kron(EPSILON, S_Z))
    assert np.allclose(G_MATRICES[3], np.kron(EPSILON, S_X))

    # Chiral spin matrices M^gamma
    assert np.allclose(M_MATRICES[1], np.kron(S_Z, EPSILON))
    assert np.allclose(M_MATRICES[2], np.kron(EPSILON, I_2))
    assert np.allclose(M_MATRICES[3], np.kron(S_X, EPSILON))

    # Kitaev projection matrices N^gamma
    assert np.allclose(N_MATRICES[1], np.kron(PI_PLUS, EPSILON))
    assert np.allclose(N_MATRICES[2], np.kron(EPSILON, PI_PLUS))
    n3_expected = 0.5 * (np.kron(S_X, EPSILON) + np.kron(EPSILON, S_X))
    assert np.allclose(N_MATRICES[3], n3_expected)


def test_orthogonality_and_normalization() -> None:
    for a in (1, 2, 3):
        for b in (1, 2, 3):
            # M and G are orthogonal Lie subalgebras
            trace_mg = np.trace(M_MATRICES[a].T @ G_MATRICES[b])
            assert np.isclose(trace_mg, 0.0)

            # Trace overlap: Tr(N^a^T G^b) = 2.0 * delta_ab
            trace_ng = np.trace(N_MATRICES[a].T @ G_MATRICES[b])
            expected_ng = 2.0 if a == b else 0.0
            assert np.isclose(trace_ng, expected_ng)

    for a in (1, 2, 3):
        assert np.isclose(np.trace(N_MATRICES[a].T @ N_MATRICES[a]), 2.0)
        assert np.isclose(np.trace(G_MATRICES[a].T @ G_MATRICES[a]), 4.0)
        assert np.isclose(np.trace(M_MATRICES[a].T @ M_MATRICES[a]), 4.0)


def test_so4_lie_algebra_commutations() -> None:
    eps = np.zeros((4, 4, 4))
    eps[1, 2, 3] = eps[2, 3, 1] = eps[3, 1, 2] = 1.0
    eps[1, 3, 2] = eps[3, 2, 1] = eps[2, 1, 3] = -1.0

    for a in (1, 2, 3):
        for b in (1, 2, 3):
            # [G^a, G^b] = -2 * sum_c eps_{abc} G^c
            comm_g = G_MATRICES[a] @ G_MATRICES[b] - G_MATRICES[b] @ G_MATRICES[a]
            expected_g = -2.0 * sum(eps[a, b, c] * G_MATRICES[c] for c in (1, 2, 3))
            assert np.allclose(comm_g, expected_g)

            # [M^a, M^b] = +2 * sum_c eps_{abc} M^c
            comm_m = M_MATRICES[a] @ M_MATRICES[b] - M_MATRICES[b] @ M_MATRICES[a]
            expected_m = 2.0 * sum(eps[a, b, c] * M_MATRICES[c] for c in (1, 2, 3))
            assert np.allclose(comm_m, expected_m)

            # [M^a, G^b] = 0
            comm_mg = M_MATRICES[a] @ G_MATRICES[b] - G_MATRICES[b] @ M_MATRICES[a]
            assert np.allclose(comm_mg, 0.0)
