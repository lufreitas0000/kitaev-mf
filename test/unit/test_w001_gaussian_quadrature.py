"""Independent small Fock-space check of the proposal's single-slice kernel."""

import runpy
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy.linalg import expm


@pytest.fixture(scope="module")
def candidate() -> dict[str, Any]:
    root = Path(__file__).resolve().parents[2]
    return runpy.run_path(
        str(root / "doc/rnd/exact-perturbation-theory/code/w001_gaussian_quadrature.py")
    )


def test_single_slice_kernel_against_full_fock_exponential(candidate: dict[str, Any]) -> None:
    """Compare determinant and insertion kernels to independent 4x4 Fock matrices."""
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1.0, -1.0])
    gamma = [np.kron(x, np.eye(2)), np.kron(y, np.eye(2)), np.kron(z, x), np.kron(z, y)]
    block = np.array([[0.0, 1.0], [-1.0, 0.0]])
    a0 = np.kron(np.diag([0.7, 1.3]), block)
    random = np.random.default_rng(11).normal(size=(4, 4))
    rotation = expm(0.2 * (random - random.T))
    a1 = rotation @ np.kron(np.diag([0.8, 1.5]), block) @ rotation.T

    def hamiltonian(a: np.ndarray) -> np.ndarray:
        result = np.zeros((4, 4), dtype=complex)
        for row in range(4):
            for col in range(4):
                result += 0.25j * a[row, col] * gamma[row] @ gamma[col]
        return result

    energies, vectors = np.linalg.eigh(hamiltonian(a0))
    state = vectors[:, 0]
    w0, u0 = np.linalg.eigh(1j * a0)
    w1, u1 = np.linalg.eigh(1j * a1)
    creation = u0[:, w0 < 0] @ u0[:, w0 < 0].conj().T
    annihilation = creation.T
    vacuum_shift = -0.5 * sum(w1[w1 > 0]) - energies[0]
    for tau in [0.0, 0.3, 1.7]:
        left, right = candidate["make_LR"](w1, u1, tau)
        matrix = creation @ left + annihilation @ right
        proposed_overlap = np.exp(-tau * vacuum_shift) * np.sqrt(np.linalg.det(matrix))
        propagator = expm(-tau * (hamiltonian(a1) - energies[0] * np.eye(4)))
        exact_overlap = np.vdot(state, propagator @ state)
        proposed = 2 * left @ np.linalg.solve(matrix, annihilation)
        exact = np.array(
            [
                [np.vdot(state, a @ propagator @ b @ state) / exact_overlap for b in gamma]
                for a in gamma
            ]
        )
        # Well-conditioned dimension-four matrices; observed residuals ~1e-15.
        assert abs(proposed_overlap - exact_overlap) < 2e-13
        assert np.allclose(proposed, exact, atol=2e-13, rtol=0)


def test_candidate_rejects_ambiguous_vacuum_and_invalid_quadrature(
    candidate: dict[str, Any],
) -> None:
    with pytest.raises(ValueError, match="zero modes"):
        candidate["compute_unprojected_candidate"](3, 4, 4)
    with pytest.raises(ValueError, match="quadrature orders"):
        candidate["compute_unprojected_candidate"](4, 0, 4)
    with pytest.raises(ValueError, match="nonnegative"):
        candidate["make_LR"](np.array([-1.0, 1.0]), np.eye(2, dtype=complex), -1.0)
