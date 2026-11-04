"""Independent finite-Fock-space tests for the local Majorana algebra."""

import numpy as np
from numpy.typing import NDArray

from kitaev_mf.domain.constants import (
    EPSILON,
    G_MATRICES,
    I_2,
    M_MATRICES,
    N_MATRICES,
    S_X,
    S_Z,
)

IDENTITY_2 = np.eye(2, dtype=np.complex128)
PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PAULI_Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
PAULI_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
TOL = 1.0e-12


def _kron_all(factors: list[NDArray[np.complex128]]) -> NDArray[np.complex128]:
    result = factors[0]
    for factor in factors[1:]:
        result = np.kron(result, factor).astype(np.complex128, copy=False)
    return result


def _majorana_operators(num_modes: int) -> tuple[NDArray[np.complex128], ...]:
    """Jordan-Wigner representation with c_(2m)=Z^m X and c_(2m+1)=Z^m Y."""
    operators: list[NDArray[np.complex128]] = []
    for mode in range(num_modes):
        for local_pauli in (PAULI_X, PAULI_Y):
            factors = [PAULI_Z if site < mode else IDENTITY_2 for site in range(num_modes)]
            factors[mode] = local_pauli
            operators.append(_kron_all(factors))
    return tuple(operators)


def _bilinear(
    majoranas: tuple[NDArray[np.complex128], ...], matrix: NDArray[np.float64]
) -> NDArray[np.complex128]:
    result = np.zeros_like(majoranas[0], dtype=np.complex128)
    for mu in range(4):
        for nu in range(4):
            result += matrix[mu, nu] * (majoranas[mu] @ majoranas[nu])
    return 0.25j * result


def _skew_basis() -> tuple[NDArray[np.float64], ...]:
    matrices: list[NDArray[np.float64]] = []
    for a in range(4):
        for b in range(a + 1, 4):
            matrix = np.zeros((4, 4), dtype=np.float64)
            matrix[a, b] = 1.0
            matrix[b, a] = -1.0
            matrices.append(matrix)
    return tuple(matrices)


def _levi_civita_3() -> NDArray[np.float64]:
    epsilon = np.zeros((3, 3, 3), dtype=np.float64)
    epsilon[0, 1, 2] = epsilon[1, 2, 0] = epsilon[2, 0, 1] = 1.0
    epsilon[0, 2, 1] = epsilon[2, 1, 0] = epsilon[1, 0, 2] = -1.0
    return epsilon


def test_t00a_01_majorana_car() -> None:
    majoranas = _majorana_operators(2)
    identity = np.eye(4, dtype=np.complex128)

    assert len(majoranas) == 4
    for mu in range(4):
        assert np.allclose(majoranas[mu].conj().T, majoranas[mu], atol=TOL, rtol=0.0)
        for nu in range(4):
            anticommutator = majoranas[mu] @ majoranas[nu] + majoranas[nu] @ majoranas[mu]
            expected = 2.0 * identity if mu == nu else np.zeros_like(identity)
            assert np.allclose(anticommutator, expected, atol=TOL, rtol=0.0)


def test_t00a_02_dirac_car_and_inverse_relations() -> None:
    c0, c1, c2, c3 = _majorana_operators(2)
    d = ((c0 + 1.0j * c1) / 2.0, (c2 + 1.0j * c3) / 2.0)
    identity = np.eye(4, dtype=np.complex128)
    zero = np.zeros_like(identity)

    for a in range(2):
        for b in range(2):
            assert np.allclose(
                d[a] @ d[b].conj().T + d[b].conj().T @ d[a],
                identity if a == b else zero,
                atol=TOL,
                rtol=0.0,
            )
            assert np.allclose(d[a] @ d[b] + d[b] @ d[a], zero, atol=TOL, rtol=0.0)
            assert np.allclose(
                d[a].conj().T @ d[b].conj().T + d[b].conj().T @ d[a].conj().T,
                zero,
                atol=TOL,
                rtol=0.0,
            )

    reconstructed = (
        d[0] + d[0].conj().T,
        -1.0j * (d[0] - d[0].conj().T),
        d[1] + d[1].conj().T,
        -1.0j * (d[1] - d[1].conj().T),
    )
    for actual, expected in zip(reconstructed, (c0, c1, c2, c3), strict=True):
        assert np.allclose(actual, expected, atol=TOL, rtol=0.0)


def test_t00a_03_bilinear_hermiticity_and_commutator_identity() -> None:
    majoranas = _majorana_operators(2)
    skew_matrices = _skew_basis()
    bilinears = tuple(_bilinear(majoranas, matrix) for matrix in skew_matrices)

    for bilinear in bilinears:
        assert np.allclose(bilinear.conj().T, bilinear, atol=TOL, rtol=0.0)

    nonzero_commutator_seen = False
    for a_matrix, a_operator in zip(skew_matrices, bilinears, strict=True):
        for b_matrix, b_operator in zip(skew_matrices, bilinears, strict=True):
            matrix_commutator = a_matrix @ b_matrix - b_matrix @ a_matrix
            operator_commutator = a_operator @ b_operator - b_operator @ a_operator
            expected = 1.0j * _bilinear(majoranas, matrix_commutator)
            assert np.allclose(operator_commutator, expected, atol=TOL, rtol=0.0)
            nonzero_commutator_seen |= not np.allclose(matrix_commutator, 0.0, atol=TOL, rtol=0.0)

    assert nonzero_commutator_seen


def test_t00a_04_explicit_gauge_and_spin_contractions() -> None:
    majoranas = _majorana_operators(2)
    cyclic_pairs = {1: (2, 3), 2: (3, 1), 3: (1, 2)}

    for gamma, (gamma_plus, gamma_minus) in cyclic_pairs.items():
        gauge_contraction = 2.0 * (majoranas[0] @ majoranas[gamma] + majoranas[gamma_plus] @ majoranas[gamma_minus])
        spin_contraction = 2.0 * (majoranas[0] @ majoranas[gamma] - majoranas[gamma_plus] @ majoranas[gamma_minus])
        assert np.allclose(
            _bilinear(majoranas, G_MATRICES[gamma]),
            0.25j * gauge_contraction,
            atol=TOL,
            rtol=0.0,
        )
        assert np.allclose(
            _bilinear(majoranas, M_MATRICES[gamma]),
            0.25j * spin_contraction,
            atol=TOL,
            rtol=0.0,
        )

    e_11 = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=np.float64)
    assert np.allclose(e_11, 0.5 * (I_2 + S_Z), atol=TOL, rtol=0.0)
    assert np.allclose(N_MATRICES[1], np.kron(e_11, EPSILON), atol=TOL, rtol=0.0)
    assert np.allclose(N_MATRICES[2], np.kron(EPSILON, e_11), atol=TOL, rtol=0.0)
    expected_n3 = 0.5 * (np.kron(S_X, EPSILON) + np.kron(EPSILON, S_X))
    assert np.allclose(N_MATRICES[3], expected_n3, atol=TOL, rtol=0.0)


def test_t00a_05_matrix_lie_algebra_and_mutual_commutation() -> None:
    epsilon = _levi_civita_3()

    for alpha in (1, 2, 3):
        assert np.allclose(G_MATRICES[alpha].T, -G_MATRICES[alpha], atol=TOL, rtol=0.0)
        assert np.allclose(M_MATRICES[alpha].T, -M_MATRICES[alpha], atol=TOL, rtol=0.0)
        for beta in (1, 2, 3):
            expected_g = -2.0 * sum(
                epsilon[alpha - 1, beta - 1, gamma - 1] * G_MATRICES[gamma] for gamma in (1, 2, 3)
            )
            expected_m = 2.0 * sum(
                epsilon[alpha - 1, beta - 1, gamma - 1] * M_MATRICES[gamma] for gamma in (1, 2, 3)
            )
            commutator_g = G_MATRICES[alpha] @ G_MATRICES[beta] - G_MATRICES[beta] @ G_MATRICES[alpha]
            commutator_m = M_MATRICES[alpha] @ M_MATRICES[beta] - M_MATRICES[beta] @ M_MATRICES[alpha]
            commutator_mg = M_MATRICES[alpha] @ G_MATRICES[beta] - G_MATRICES[beta] @ M_MATRICES[alpha]
            assert np.allclose(commutator_g, expected_g, atol=TOL, rtol=0.0)
            assert np.allclose(commutator_m, expected_m, atol=TOL, rtol=0.0)
            assert np.allclose(commutator_mg, 0.0, atol=TOL, rtol=0.0)


def test_t00a_06_full_trace_gram_matrices() -> None:
    gauge_gram = np.array([[np.trace(G_MATRICES[a].T @ G_MATRICES[b]) for b in (1, 2, 3)] for a in (1, 2, 3)])
    spin_gram = np.array([[np.trace(M_MATRICES[a].T @ M_MATRICES[b]) for b in (1, 2, 3)] for a in (1, 2, 3)])
    mixed_gram = np.array([[np.trace(M_MATRICES[a].T @ G_MATRICES[b]) for b in (1, 2, 3)] for a in (1, 2, 3)])

    assert np.allclose(gauge_gram, 4.0 * np.eye(3), atol=TOL, rtol=0.0)
    assert np.allclose(spin_gram, 4.0 * np.eye(3), atol=TOL, rtol=0.0)
    assert np.allclose(mixed_gram, np.zeros((3, 3)), atol=TOL, rtol=0.0)


def test_t00a_07_even_bilinears_commute_across_sites() -> None:
    majoranas = _majorana_operators(4)
    skew_matrices = _skew_basis()
    site_zero = majoranas[:4]
    site_one = majoranas[4:]

    for matrix_a in skew_matrices:
        operator_a = _bilinear(site_zero, matrix_a)
        for matrix_b in skew_matrices:
            operator_b = _bilinear(site_one, matrix_b)
            commutator = operator_a @ operator_b - operator_b @ operator_a
            assert np.allclose(commutator, 0.0, atol=TOL, rtol=0.0)
