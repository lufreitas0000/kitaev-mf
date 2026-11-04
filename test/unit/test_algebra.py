"""Test exact operator relations for SO(4) Majorana and gauge matrices."""

import numpy as np

from src.domain.constants import G_MATRICES, M_MATRICES, N_MATRICES


def test_antisymmetry() -> None:
    for gamma in (1, 2, 3):
        assert np.allclose(N_MATRICES[gamma], -N_MATRICES[gamma].T)
        assert np.allclose(G_MATRICES[gamma], -G_MATRICES[gamma].T)
        assert np.allclose(M_MATRICES[gamma], -M_MATRICES[gamma].T)


def test_cross_gram_matrix() -> None:
    indices = (1, 2, 3)
    cross_gram = np.array([
        [
            np.trace(N_MATRICES[a].T @ G_MATRICES[b])
            for b in indices
        ]
        for a in indices
    ])
    expected = 2.0 * np.eye(3)
    np.testing.assert_allclose(cross_gram, expected)


def test_trace_normalization() -> None:
    for a in (1, 2, 3):
        assert np.isclose(np.trace(N_MATRICES[a].T @ N_MATRICES[a]), 2.0)
        assert np.isclose(np.trace(G_MATRICES[a].T @ G_MATRICES[a]), 4.0)
