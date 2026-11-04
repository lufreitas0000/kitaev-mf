"""Test exact operator relations for SO(4) Majorana and gauge matrices."""

import numpy as np

from src.domain.constants import G_MATRICES, M_MATRICES, N_MATRICES


def test_antisymmetry():
    for gamma in (1, 2, 3):
        assert np.allclose(N_MATRICES[gamma], -N_MATRICES[gamma].T)
        assert np.allclose(G_MATRICES[gamma], -G_MATRICES[gamma].T)
        assert np.allclose(M_MATRICES[gamma], -M_MATRICES[gamma].T)


def test_orthogonality():
    for a in (1, 2, 3):
        for b in (1, 2, 3):
            trace_ng = np.trace(N_MATRICES[a].T @ G_MATRICES[b])
            assert np.isclose(trace_ng, 0.0)


def test_trace_normalization():
    for a in (1, 2, 3):
        assert np.isclose(np.trace(N_MATRICES[a].T @ N_MATRICES[a]), 2.0)
        assert np.isclose(np.trace(G_MATRICES[a].T @ G_MATRICES[a]), 4.0)
