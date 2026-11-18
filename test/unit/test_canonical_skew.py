"""Unit tests for canonical_skew_schur real-Schur decomposition including null blocks."""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import InhomogeneousMeanFieldState
from kitaev_mf.services.canonical_skew import ScalableSchurResult, canonical_skew_schur
from kitaev_mf.services.real_space import RealSpaceTorusCluster, compute_skew_pfaffian


def _build_sigma(energies: np.ndarray) -> np.ndarray:
    """Build block-diagonal canonical skew matrix diag(e_k * [[0, 1], [-1, 0]])."""
    size = 2 * len(energies)
    sigma = np.zeros((size, size), dtype=np.float64)
    for idx, energy in enumerate(energies):
        sigma[2 * idx, 2 * idx + 1] = float(energy)
        sigma[2 * idx + 1, 2 * idx] = -float(energy)
    return sigma


def test_canonical_skew_nonsingular_reconstruction_and_orthogonality() -> None:
    """Verify A = Q Sigma Q^T, Q^T Q = I, and ascending positive frequencies."""
    rng = np.random.default_rng(123)
    raw = rng.standard_normal((8, 8))
    q_rand, _ = np.linalg.qr(raw)
    target_energies = np.array([0.4, 1.1, 2.7, 5.2], dtype=np.float64)
    A = q_rand @ _build_sigma(target_energies) @ q_rand.T
    A = 0.5 * (A - A.T)

    res = canonical_skew_schur(A)
    assert isinstance(res, ScalableSchurResult)
    assert np.allclose(res.energies, target_energies, atol=1e-12)
    assert res.orthogonality_error < 1e-12
    assert res.reconstruction_error < 1e-12
    assert np.isclose(abs(res.det_Q), 1.0, atol=1e-12)
    assert np.allclose(res.Q.T @ res.Q, np.eye(8), atol=1e-12)
    assert np.allclose(res.Q @ _build_sigma(res.energies) @ res.Q.T, A, atol=1e-11)
    assert np.isclose(
        compute_skew_pfaffian(A),
        res.det_Q * float(np.prod(res.energies)),
        atol=1e-11,
    )

    # Swapping two basis columns flips Pfaffian sign and det_Q consistently
    swap = [1, 0] + list(range(2, 8))
    A_neg = A[np.ix_(swap, swap)]
    res_neg = canonical_skew_schur(A_neg)
    assert res_neg.det_Q == -res.det_Q
    assert np.isclose(
        compute_skew_pfaffian(A_neg),
        res_neg.det_Q * float(np.prod(res_neg.energies)),
        atol=1e-11,
    )


def test_canonical_skew_interleaved_separated_null_blocks() -> None:
    """Verify separated 1x1 Schur zero blocks are paired without distorting nonzero blocks."""
    # Place 1x1 zeros at non-adjacent indices (2 and 5) before orthogonal conjugation.
    A_base = np.zeros((6, 6), dtype=np.float64)
    A_base[0, 1], A_base[1, 0] = 3.0, -3.0
    A_base[3, 4], A_base[4, 3] = 1.25, -1.25

    rng = np.random.default_rng(77)
    q_orth, _ = np.linalg.qr(rng.standard_normal((6, 6)))
    A = q_orth @ A_base @ q_orth.T
    A = 0.5 * (A - A.T)

    res = canonical_skew_schur(A, tol=1e-12)
    assert np.allclose(res.energies, np.array([0.0, 1.25, 3.0]), atol=1e-12)
    assert res.orthogonality_error < 1e-12
    assert res.reconstruction_error < 1e-12
    assert np.isclose(abs(res.det_Q), 1.0, atol=1e-12)


def test_canonical_skew_identically_zero_matrix() -> None:
    """Verify identically zero even-dimensional matrix decomposes cleanly to zero pairs."""
    res = canonical_skew_schur(np.zeros((4, 4), dtype=np.float64))
    assert np.allclose(res.energies, np.zeros(2), atol=1e-15)
    assert res.orthogonality_error < 1e-14
    assert res.reconstruction_error == 0.0
    assert np.isclose(abs(res.det_Q), 1.0, atol=1e-14)


def test_canonical_skew_dirac_resonance_torus_null_modes() -> None:
    """Verify real-Schur null pairing handles exact Dirac zero modes on a 3 | L torus."""
    cluster = RealSpaceTorusCluster(L=3)
    params = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))
    state = InhomogeneousMeanFieldState.from_pure_kitaev(num_cells=cluster.num_cells)
    A = cluster.assemble_skew_hamiltonian(params, state)

    res = canonical_skew_schur(A, tol=1e-12)
    assert np.sum(res.energies < 1e-12) >= 2
    assert res.orthogonality_error < 1e-12
    assert res.reconstruction_error < 1e-11
    assert np.isclose(abs(res.det_Q), 1.0, atol=1e-12)


def test_canonical_skew_input_validation() -> None:
    """Verify invalid shapes, dtypes, non-skew inputs, and invalid tolerances raise ValueError."""
    with pytest.raises(ValueError, match="nonempty even-dimensional"):
        canonical_skew_schur(np.zeros((0, 0)))

    with pytest.raises(ValueError, match="nonempty even-dimensional"):
        canonical_skew_schur(np.zeros((5, 5)))

    with pytest.raises(ValueError, match="nonempty even-dimensional"):
        canonical_skew_schur(np.zeros((4, 6)))

    with pytest.raises(ValueError, match="real and finite"):
        canonical_skew_schur(np.array([[0.0, np.nan], [-np.nan, 0.0]]))

    with pytest.raises(ValueError, match="real and finite"):
        canonical_skew_schur(np.array([[0.0, 1.0j], [-1.0j, 0.0]]))

    with pytest.raises(ValueError, match="finite and positive"):
        canonical_skew_schur(np.array([[0.0, 1.0], [-1.0, 0.0]]), tol=0.0)

    with pytest.raises(ValueError, match="skew-symmetric"):
        canonical_skew_schur(np.array([[0.0, 1.0], [0.5, 0.0]]), tol=1e-12)
