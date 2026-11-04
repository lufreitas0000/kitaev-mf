"""Independent finite-Fock-space tests for the local spin/gauge representation."""

import numpy as np
from numpy.typing import NDArray

from kitaev_mf.domain.constants import G_MATRICES, M_MATRICES

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


def _majorana_operators() -> tuple[NDArray[np.complex128], ...]:
    """Return four Majoranas in an independent two-mode Jordan-Wigner representation."""
    operators: list[NDArray[np.complex128]] = []
    for mode in range(2):
        for local_pauli in (PAULI_X, PAULI_Y):
            factors = [
                PAULI_Z if site < mode else IDENTITY_2
                for site in range(2)
            ]
            factors[mode] = local_pauli
            operators.append(_kron_all(factors))
    return tuple(operators)


def _bilinear(
    majoranas: tuple[NDArray[np.complex128], ...],
    matrix: NDArray[np.float64],
) -> NDArray[np.complex128]:
    result = np.zeros_like(majoranas[0], dtype=np.complex128)
    for mu in range(4):
        for nu in range(4):
            result += matrix[mu, nu] * (majoranas[mu] @ majoranas[nu])
    return 0.25j * result


def _local_operators() -> tuple[
    NDArray[np.complex128],
    dict[int, NDArray[np.complex128]],
    dict[int, NDArray[np.complex128]],
]:
    majoranas = _majorana_operators()
    parity = majoranas[0] @ majoranas[1] @ majoranas[2] @ majoranas[3]
    spin = {
        gamma: (0.5 * _bilinear(majoranas, M_MATRICES[gamma])).astype(
            np.complex128, copy=False
        )
        for gamma in (1, 2, 3)
    }
    gauge = {
        gamma: (0.5 * _bilinear(majoranas, G_MATRICES[gamma])).astype(
            np.complex128, copy=False
        )
        for gamma in (1, 2, 3)
    }
    return parity, spin, gauge


def _casimir(
    generators: dict[int, NDArray[np.complex128]],
) -> NDArray[np.complex128]:
    result = np.zeros((4, 4), dtype=np.complex128)
    for gamma in (1, 2, 3):
        result += generators[gamma] @ generators[gamma]
    return result


def test_t00b_01_local_parity_is_hermitian_involution_with_equal_sectors() -> None:
    parity, _, _ = _local_operators()
    identity = np.eye(4, dtype=np.complex128)

    assert np.allclose(parity.conj().T, parity, atol=TOL, rtol=0.0)
    assert np.allclose(parity @ parity, identity, atol=TOL, rtol=0.0)
    assert np.isclose(np.trace(parity), 0.0, atol=TOL, rtol=0.0)


def test_t00b_02_parity_commutes_with_spin_and_gauge_generators() -> None:
    parity, spin, gauge = _local_operators()

    for gamma in (1, 2, 3):
        assert np.allclose(parity @ spin[gamma] - spin[gamma] @ parity, 0.0, atol=TOL, rtol=0.0)
        assert np.allclose(parity @ gauge[gamma] - gauge[gamma] @ parity, 0.0, atol=TOL, rtol=0.0)


def test_t00b_03_spin_and_gauge_casimirs_equal_parity_expressions() -> None:
    parity, spin, gauge = _local_operators()
    identity = np.eye(4, dtype=np.complex128)

    spin_casimir = _casimir(spin)
    gauge_casimir = _casimir(gauge)
    assert np.allclose(
        spin_casimir,
        0.375 * (identity + parity),
        atol=TOL,
        rtol=0.0,
    )
    assert np.allclose(
        gauge_casimir,
        0.375 * (identity - parity),
        atol=TOL,
        rtol=0.0,
    )


def test_t00b_04_parity_sectors_are_spin_and_gauge_doublets() -> None:
    parity, spin, gauge = _local_operators()
    identity = np.eye(4, dtype=np.complex128)
    projector_plus = 0.5 * (identity + parity)
    projector_minus = 0.5 * (identity - parity)
    spin_casimir = _casimir(spin)
    gauge_casimir = _casimir(gauge)

    assert np.allclose(projector_plus @ projector_plus, projector_plus, atol=TOL, rtol=0.0)
    assert np.allclose(projector_minus @ projector_minus, projector_minus, atol=TOL, rtol=0.0)
    assert np.allclose(projector_plus @ projector_minus, 0.0, atol=TOL, rtol=0.0)
    assert np.isclose(np.trace(projector_plus), 2.0, atol=TOL, rtol=0.0)
    assert np.isclose(np.trace(projector_minus), 2.0, atol=TOL, rtol=0.0)

    assert np.allclose(spin_casimir @ projector_plus, 0.75 * projector_plus, atol=TOL, rtol=0.0)
    assert np.allclose(gauge_casimir @ projector_plus, 0.0, atol=TOL, rtol=0.0)
    assert np.allclose(spin_casimir @ projector_minus, 0.0, atol=TOL, rtol=0.0)
    assert np.allclose(gauge_casimir @ projector_minus, 0.75 * projector_minus, atol=TOL, rtol=0.0)

    for gamma in (1, 2, 3):
        assert np.allclose(gauge[gamma] @ projector_plus, 0.0, atol=TOL, rtol=0.0)
        assert np.allclose(spin[gamma] @ projector_minus, 0.0, atol=TOL, rtol=0.0)


def test_t00b_05_gauge_sign_and_dirac_occupation_parity() -> None:
    parity, spin, gauge = _local_operators()
    majoranas = _majorana_operators()
    identity = np.eye(4, dtype=np.complex128)
    epsilon = np.zeros((3, 3, 3), dtype=np.float64)
    epsilon[0, 1, 2] = epsilon[1, 2, 0] = epsilon[2, 0, 1] = 1.0
    epsilon[0, 2, 1] = epsilon[2, 1, 0] = epsilon[1, 0, 2] = -1.0

    conventional_gauge = {gamma: -gauge[gamma] for gamma in (1, 2, 3)}
    for alpha in (1, 2, 3):
        for beta in (1, 2, 3):
            spin_commutator = spin[alpha] @ spin[beta] - spin[beta] @ spin[alpha]
            gauge_commutator = (
                conventional_gauge[alpha] @ conventional_gauge[beta]
                - conventional_gauge[beta] @ conventional_gauge[alpha]
            )
            mixed_commutator = (
                spin[alpha] @ conventional_gauge[beta]
                - conventional_gauge[beta] @ spin[alpha]
            )
            expected_spin = 1.0j * sum(
                epsilon[alpha - 1, beta - 1, gamma - 1] * spin[gamma]
                for gamma in (1, 2, 3)
            )
            expected_gauge = 1.0j * sum(
                epsilon[alpha - 1, beta - 1, gamma - 1] * conventional_gauge[gamma]
                for gamma in (1, 2, 3)
            )
            assert np.allclose(spin_commutator, expected_spin, atol=TOL, rtol=0.0)
            assert np.allclose(gauge_commutator, expected_gauge, atol=TOL, rtol=0.0)
            assert np.allclose(mixed_commutator, 0.0, atol=TOL, rtol=0.0)

    d1 = 0.5 * (majoranas[0] + 1.0j * majoranas[1])
    d2 = 0.5 * (majoranas[2] + 1.0j * majoranas[3])
    number_1 = d1.conj().T @ d1
    number_2 = d2.conj().T @ d2
    dirac_parity = (identity - 2.0 * number_1) @ (identity - 2.0 * number_2)

    assert np.allclose(parity, -dirac_parity, atol=TOL, rtol=0.0)
