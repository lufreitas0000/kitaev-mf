"""Finite-Fock-space tests for physical projection and the Pauli algebra."""

import numpy as np
from numpy.typing import NDArray

from kitaev_mf.domain.constants import G_MATRICES, N_MATRICES

IDENTITY_2 = np.eye(2, dtype=np.complex128)
PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PAULI_Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
PAULI_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
TOL = 1.0e-12
COMPONENTS = (1, 2, 3)


def _kron_all(factors: list[NDArray[np.complex128]]) -> NDArray[np.complex128]:
    result = factors[0]
    for factor in factors[1:]:
        result = np.kron(result, factor).astype(np.complex128, copy=False)
    return result


def _majorana_operators(num_modes: int) -> tuple[NDArray[np.complex128], ...]:
    """Jordan-Wigner representation c_(2m)=Z^m X and c_(2m+1)=Z^m Y."""
    operators: list[NDArray[np.complex128]] = []
    for mode in range(num_modes):
        for local_pauli in (PAULI_X, PAULI_Y):
            factors = [
                PAULI_Z if site < mode else IDENTITY_2
                for site in range(num_modes)
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


def _local_operators(
    majoranas: tuple[NDArray[np.complex128], ...],
) -> tuple[
    NDArray[np.complex128],
    dict[int, NDArray[np.complex128]],
    dict[int, NDArray[np.complex128]],
]:
    parity = majoranas[0] @ majoranas[1] @ majoranas[2] @ majoranas[3]
    spin_pauli = {
        gamma: 1.0j * majoranas[0] @ majoranas[gamma]
        for gamma in COMPONENTS
    }
    gauge = {
        gamma: 0.5 * _bilinear(majoranas, G_MATRICES[gamma])
        for gamma in COMPONENTS
    }
    return parity, spin_pauli, gauge


def _levi_civita(alpha: int, beta: int, gamma: int) -> int:
    if (alpha, beta, gamma) in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
        return 1
    if (alpha, beta, gamma) in ((2, 1, 3), (3, 2, 1), (1, 3, 2)):
        return -1
    return 0


def _local_projector(
    majoranas: tuple[NDArray[np.complex128], ...],
) -> NDArray[np.complex128]:
    parity, _, _ = _local_operators(majoranas)
    identity = np.eye(parity.shape[0], dtype=np.complex128)
    return 0.5 * (identity + parity)


def _global_site_operators(
    majoranas: tuple[NDArray[np.complex128], ...],
    site: int,
) -> tuple[
    NDArray[np.complex128],
    dict[int, NDArray[np.complex128]],
    dict[int, NDArray[np.complex128]],
]:
    local_majoranas = majoranas[4 * site : 4 * site + 4]
    return _local_operators(local_majoranas)


def test_t00c_01_local_physical_projector_is_orthogonal_rank_two() -> None:
    majoranas = _majorana_operators(2)
    parity, _, _ = _local_operators(majoranas)
    projector = _local_projector(majoranas)
    identity = np.eye(4, dtype=np.complex128)

    assert np.allclose(projector, 0.5 * (identity + parity), atol=TOL, rtol=0.0)
    assert np.allclose(projector.conj().T, projector, atol=TOL, rtol=0.0)
    assert np.allclose(projector @ projector, projector, atol=TOL, rtol=0.0)
    assert np.isclose(np.trace(projector), 2.0, atol=TOL, rtol=0.0)


def test_t00c_02_global_projector_is_product_of_local_constraints() -> None:
    majoranas = _majorana_operators(4)
    identity = np.eye(16, dtype=np.complex128)
    local_projectors = []
    for site in (0, 1):
        parity, _, _ = _global_site_operators(majoranas, site)
        local_projectors.append(0.5 * (identity + parity))

    assert np.allclose(
        local_projectors[0] @ local_projectors[1],
        local_projectors[1] @ local_projectors[0],
        atol=TOL,
        rtol=0.0,
    )
    global_projector = local_projectors[0] @ local_projectors[1]
    assert np.allclose(global_projector.conj().T, global_projector, atol=TOL, rtol=0.0)
    assert np.allclose(
        global_projector @ global_projector,
        global_projector,
        atol=TOL,
        rtol=0.0,
    )
    assert np.isclose(np.trace(global_projector), 4.0, atol=TOL, rtol=0.0)

    total_parity = (
        majoranas[0] @ majoranas[1] @ majoranas[2] @ majoranas[3]
        @ majoranas[4] @ majoranas[5] @ majoranas[6] @ majoranas[7]
    )
    total_parity_projector = 0.5 * (identity + total_parity)
    assert not np.allclose(
        global_projector,
        total_parity_projector,
        atol=TOL,
        rtol=0.0,
    )


def test_t00c_03_kitaev_bilinears_obey_unprojected_clifford_relations() -> None:
    majoranas = _majorana_operators(2)
    _, spin_pauli, _ = _local_operators(majoranas)
    identity = np.eye(4, dtype=np.complex128)
    zero = np.zeros((4, 4), dtype=np.complex128)

    for gamma in COMPONENTS:
        direct = 1.0j * majoranas[0] @ majoranas[gamma]
        from_n_matrix = 2.0 * _bilinear(majoranas, N_MATRICES[gamma])
        assert np.allclose(direct, from_n_matrix, atol=TOL, rtol=0.0)
        assert np.allclose(direct.conj().T, direct, atol=TOL, rtol=0.0)
        assert np.allclose(direct @ direct, identity, atol=TOL, rtol=0.0)

    for alpha in COMPONENTS:
        for beta in COMPONENTS:
            anticommutator = (
                spin_pauli[alpha] @ spin_pauli[beta]
                + spin_pauli[beta] @ spin_pauli[alpha]
            )
            expected = 2.0 * identity if alpha == beta else zero
            assert np.allclose(anticommutator, expected, atol=TOL, rtol=0.0)


def test_t00c_04_full_space_product_identity_includes_gauge_term() -> None:
    majoranas = _majorana_operators(2)
    _, spin_pauli, gauge = _local_operators(majoranas)
    identity = np.eye(4, dtype=np.complex128)

    for alpha in COMPONENTS:
        for beta in COMPONENTS:
            expected = identity if alpha == beta else np.zeros_like(identity)
            for gamma in COMPONENTS:
                epsilon = _levi_civita(alpha, beta, gamma)
                expected += 1.0j * epsilon * spin_pauli[gamma]
                expected -= 4.0j * epsilon * gauge[gamma]
            actual = spin_pauli[alpha] @ spin_pauli[beta]
            assert np.allclose(actual, expected, atol=TOL, rtol=0.0)


def test_t00c_05_projection_recovers_local_pauli_multiplication() -> None:
    majoranas = _majorana_operators(2)
    parity, spin_pauli, _ = _local_operators(majoranas)
    identity = np.eye(4, dtype=np.complex128)
    projector = 0.5 * (identity + parity)

    for alpha in COMPONENTS:
        assert np.allclose(
            projector @ spin_pauli[alpha],
            spin_pauli[alpha] @ projector,
            atol=TOL,
            rtol=0.0,
        )
        for beta in COMPONENTS:
            expected = identity if alpha == beta else np.zeros_like(identity)
            for gamma in COMPONENTS:
                expected += (
                    1.0j
                    * _levi_civita(alpha, beta, gamma)
                    * spin_pauli[gamma]
                )
            residual = projector @ (
                spin_pauli[alpha] @ spin_pauli[beta] - expected
            ) @ projector
            assert np.allclose(residual, 0.0, atol=TOL, rtol=0.0)


def test_t00c_06_two_site_projected_operators_form_tensor_product_pauli_algebras() -> None:
    majoranas = _majorana_operators(4)
    identity = np.eye(16, dtype=np.complex128)
    site_data = [_global_site_operators(majoranas, site) for site in (0, 1)]
    local_projectors = [
        0.5 * (identity + site_data[site][0])
        for site in (0, 1)
    ]
    global_projector = local_projectors[0] @ local_projectors[1]
    projected = [
        {
            gamma: global_projector @ site_data[site][1][gamma] @ global_projector
            for gamma in COMPONENTS
        }
        for site in (0, 1)
    ]

    for site in (0, 1):
        for alpha in COMPONENTS:
            for beta in COMPONENTS:
                expected = (
                    global_projector
                    if alpha == beta
                    else np.zeros_like(global_projector)
                )
                for gamma in COMPONENTS:
                    expected += (
                        1.0j
                        * _levi_civita(alpha, beta, gamma)
                        * projected[site][gamma]
                    )
                actual = projected[site][alpha] @ projected[site][beta]
                assert np.allclose(actual, expected, atol=TOL, rtol=0.0)

    for alpha in COMPONENTS:
        for beta in COMPONENTS:
            commutator = (
                projected[0][alpha] @ projected[1][beta]
                - projected[1][beta] @ projected[0][alpha]
            )
            assert np.allclose(commutator, 0.0, atol=TOL, rtol=0.0)
