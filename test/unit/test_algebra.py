"""Test exact operator relations for SO(4) Majorana and gauge matrices."""
import numpy as np
from kitaev_mf.domain.constants import G_MATRICES, M_MATRICES, N_MATRICES


def test_antisymmetry():
    for gamma in (1, 2, 3):
        assert np.allclose(N_MATRICES[gamma], -N_MATRICES[gamma].T)
        assert np.allclose(G_MATRICES[gamma], -G_MATRICES[gamma].T)
        assert np.allclose(M_MATRICES[gamma], -M_MATRICES[gamma].T)


def test_orthogonality_and_normalization():
    for a in (1, 2, 3):
        for b in (1, 2, 3):
            trace_mg = np.trace(M_MATRICES[a].T @ G_MATRICES[b])
            assert np.isclose(trace_mg, 0.0)
            trace_ng = np.trace(N_MATRICES[a].T @ G_MATRICES[b])
            # Tr(N^a G^b) = Tr((M^a - G^a) G^b) = - Tr(G^a G^b) = -4 delta_ab
            expected_ng = -4.0 if a == b else 0.0
            assert np.isclose(trace_ng, expected_ng)

    for a in (1, 2, 3):
        assert np.isclose(np.trace(G_MATRICES[a].T @ G_MATRICES[a]), 4.0)
        assert np.isclose(np.trace(M_MATRICES[a].T @ M_MATRICES[a]), 4.0)
        assert np.isclose(np.trace(N_MATRICES[a].T @ N_MATRICES[a]), 4.0)


def test_so4_lie_algebra_commutations():
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
