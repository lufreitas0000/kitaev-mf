"""Finite-Fock-space tests for exchange tensors and the exact quartic Hamiltonian."""

import numpy as np
from numpy.typing import NDArray

from kitaev_mf.domain.constants import G_MATRICES, N_MATRICES

IDENTITY_2 = np.eye(2, dtype=np.complex128)
PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PAULI_Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
PAULI_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
COMPONENTS = (1, 2, 3)
TOL = 1.0e-12


def _kron_all(factors: list[NDArray[np.complex128]]) -> NDArray[np.complex128]:
    result = factors[0]
    for factor in factors[1:]:
        result = np.kron(result, factor).astype(np.complex128, copy=False)
    return result


def _majorana_operators(num_modes: int) -> tuple[NDArray[np.complex128], ...]:
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


def _site_majoranas(
    majoranas: tuple[NDArray[np.complex128], ...],
    site: int,
) -> tuple[NDArray[np.complex128], ...]:
    return majoranas[4 * site : 4 * site + 4]


def _site_parity(
    majoranas: tuple[NDArray[np.complex128], ...],
) -> NDArray[np.complex128]:
    return majoranas[0] @ majoranas[1] @ majoranas[2] @ majoranas[3]


def _site_spin(
    majoranas: tuple[NDArray[np.complex128], ...],
) -> dict[int, NDArray[np.complex128]]:
    return {
        gamma: 1.0j * majoranas[0] @ majoranas[gamma]
        for gamma in COMPONENTS
    }


def _exchange_tensor(
    gamma: int,
    j_value: float,
    k_value: float,
    gamma_value: float,
    gamma_prime: float,
    dm_value: float,
) -> NDArray[np.float64]:
    gamma_plus = 1 + (gamma % 3)
    gamma_minus = 1 + ((gamma + 1) % 3)
    tensor = j_value * np.eye(3, dtype=np.float64)
    tensor[gamma - 1, gamma - 1] += k_value
    tensor[gamma_plus - 1, gamma_minus - 1] += gamma_value + dm_value
    tensor[gamma_minus - 1, gamma_plus - 1] += gamma_value - dm_value
    tensor[gamma_plus - 1, gamma - 1] += gamma_prime
    tensor[gamma - 1, gamma_plus - 1] += gamma_prime
    tensor[gamma_minus - 1, gamma - 1] += gamma_prime
    tensor[gamma - 1, gamma_minus - 1] += gamma_prime
    return tensor


def test_t00d_01_exchange_tensor_has_expected_symmetric_and_dm_parts() -> None:
    couplings = (0.7, -1.1, 0.4, -0.3, 0.9)
    j_value, k_value, gamma_value, gamma_prime, dm_value = couplings

    for gamma in COMPONENTS:
        tensor = _exchange_tensor(
            gamma, j_value, k_value, gamma_value, gamma_prime, dm_value
        )
        symmetric = 0.5 * (tensor + tensor.T)
        antisymmetric = 0.5 * (tensor - tensor.T)

        gamma_plus = 1 + (gamma % 3)
        gamma_minus = 1 + ((gamma + 1) % 3)
        cyclic_basis = (gamma_plus - 1, gamma_minus - 1, gamma - 1)
        expected_symmetric = np.zeros((3, 3), dtype=np.float64)
        expected_symmetric[np.ix_(cyclic_basis, cyclic_basis)] = np.array(
            [
                [j_value, gamma_value, gamma_prime],
                [gamma_value, j_value, gamma_prime],
                [gamma_prime, gamma_prime, j_value + k_value],
            ],
            dtype=np.float64,
        )
        expected_antisymmetric = np.zeros((3, 3), dtype=np.float64)
        expected_antisymmetric[gamma_plus - 1, gamma_minus - 1] = dm_value
        expected_antisymmetric[gamma_minus - 1, gamma_plus - 1] = -dm_value

        assert np.allclose(symmetric, expected_symmetric, atol=TOL, rtol=0.0)
        assert np.allclose(antisymmetric, expected_antisymmetric, atol=TOL, rtol=0.0)
        assert np.allclose(tensor.T, tensor - 2.0 * antisymmetric, atol=TOL, rtol=0.0)


def test_t00d_02_ordered_edge_sum_equals_single_directed_edge() -> None:
    majoranas = _majorana_operators(4)
    identity = np.eye(16, dtype=np.complex128)
    spin = [_site_spin(_site_majoranas(majoranas, site)) for site in (0, 1)]
    tensor = np.array(
        [[0.8, 0.2, -0.4], [-0.1, 0.6, 0.3], [0.5, -0.7, -0.2]],
        dtype=np.float64,
    )
    forward = np.zeros_like(identity)
    reverse = np.zeros_like(identity)

    for alpha in COMPONENTS:
        for beta in COMPONENTS:
            forward += tensor[alpha - 1, beta - 1] * spin[0][alpha] @ spin[1][beta]
            reverse += (
                tensor.T[alpha - 1, beta - 1]
                * spin[1][alpha]
                @ spin[0][beta]
            )

    ordered_half_sum = 0.5 * (forward + reverse)
    assert np.allclose(reverse, forward, atol=TOL, rtol=0.0)
    assert np.allclose(ordered_half_sum, forward, atol=TOL, rtol=0.0)


def test_t00d_03_majorana_bilinear_substitution_and_exchange_prefactor() -> None:
    majoranas = _majorana_operators(4)
    sites = [_site_majoranas(majoranas, site) for site in (0, 1)]
    spins = [_site_spin(site_majoranas) for site_majoranas in sites]
    tensor = np.array(
        [[0.4, 0.7, -0.2], [0.1, -0.6, 0.5], [-0.3, 0.8, 0.9]],
        dtype=np.float64,
    )

    majorana_bilinears: list[dict[int, NDArray[np.complex128]]] = []
    for site_majoranas, spin in zip(sites, spins, strict=True):
        operators = {
            gamma: 4.0 * _bilinear(site_majoranas, N_MATRICES[gamma])
            for gamma in COMPONENTS
        }
        for gamma in COMPONENTS:
            assert np.allclose(
                operators[gamma],
                2.0 * spin[gamma],
                atol=TOL,
                rtol=0.0,
            )
        majorana_bilinears.append(operators)

    spin_exchange = np.zeros_like(identity)
    majorana_exchange = np.zeros_like(identity)
    for alpha in COMPONENTS:
        for beta in COMPONENTS:
            spin_exchange += (
                0.5
                * tensor[alpha - 1, beta - 1]
                * spins[0][alpha]
                @ spins[1][beta]
            )
            majorana_exchange += (
                0.125
                * tensor[alpha - 1, beta - 1]
                * majorana_bilinears[0][alpha]
                @ majorana_bilinears[1][beta]
            )
    assert np.allclose(majorana_exchange, spin_exchange, atol=TOL, rtol=0.0)


def test_t00d_04_site_dependent_field_prefactor_is_exact() -> None:
    majoranas = _majorana_operators(4)
    identity = np.eye(16, dtype=np.complex128)
    fields = (
        np.array([0.2, -0.4, 0.7], dtype=np.float64),
        np.array([-0.3, 0.6, 0.1], dtype=np.float64),
    )
    spin = [
        _site_spin(_site_majoranas(majoranas, site))
        for site in (0, 1)
    ]
    direct_field = np.zeros_like(identity)
    majorana_field = np.zeros_like(identity)

    for site in (0, 1):
        site_majoranas = _site_majoranas(majoranas, site)
        for gamma in COMPONENTS:
            direct_field -= fields[site][gamma - 1] * spin[site][gamma]
            bilinear = 4.0 * _bilinear(site_majoranas, N_MATRICES[gamma])
            majorana_field -= (
                0.25 * 2.0 * fields[site][gamma - 1] * bilinear
            )

    assert np.allclose(majorana_field, direct_field, atol=TOL, rtol=0.0)


def test_t00d_05_constraint_multiplier_normalization_and_projection() -> None:
    majoranas = _majorana_operators(4)
    identity = np.eye(16, dtype=np.complex128)
    local_projectors = []
    for site in (0, 1):
        site_majoranas = _site_majoranas(majoranas, site)
        parity = _site_parity(site_majoranas)
        local_projectors.append(0.5 * (identity + parity))
    physical_projector = local_projectors[0] @ local_projectors[1]

    multiplier_term = np.zeros_like(identity)
    for site in (0, 1):
        site_majoranas = _site_majoranas(majoranas, site)
        for gamma in COMPONENTS:
            multiplier = 0.17 * (site + 1) * gamma
            gauge_generator = 0.5 * _bilinear(
                site_majoranas, G_MATRICES[gamma]
            )
            i_cGc = 4.0 * _bilinear(site_majoranas, G_MATRICES[gamma])
            assert np.allclose(
                0.25 * multiplier * i_cGc,
                2.0 * multiplier * gauge_generator,
                atol=TOL,
                rtol=0.0,
            )
            multiplier_term += 0.25 * multiplier * i_cGc

    assert np.allclose(
        physical_projector @ multiplier_term @ physical_projector,
        0.0,
        atol=TOL,
        rtol=0.0,
    )


def test_t00d_06_two_site_projected_spin_and_majorana_hamiltonians_agree() -> None:
    majoranas = _majorana_operators(4)
    identity = np.eye(16, dtype=np.complex128)
    sites = [_site_majoranas(majoranas, site) for site in (0, 1)]
    spin = [_site_spin(site_majoranas) for site_majoranas in sites]
    local_projectors = [
        0.5 * (identity + _site_parity(site_majoranas))
        for site_majoranas in sites
    ]
    physical_projector = local_projectors[0] @ local_projectors[1]

    tensor = np.array(
        [[0.9, 0.35, -0.2], [-0.15, -0.4, 0.55], [0.3, -0.1, 0.6]],
        dtype=np.float64,
    )
    fields = (
        np.array([0.13, -0.21, 0.32], dtype=np.float64),
        np.array([-0.17, 0.29, 0.08], dtype=np.float64),
    )
    multipliers = (
        np.array([0.11, -0.07, 0.19], dtype=np.float64),
        np.array([-0.23, 0.05, 0.14], dtype=np.float64),
    )

    spin_hamiltonian = np.zeros_like(identity)
    majorana_hamiltonian = np.zeros_like(identity)

    for alpha in COMPONENTS:
        for beta in COMPONENTS:
            spin_hamiltonian += 0.5 * (
                tensor[alpha - 1, beta - 1]
                * spin[0][alpha]
                @ spin[1][beta]
                + tensor.T[alpha - 1, beta - 1]
                * spin[1][alpha]
                @ spin[0][beta]
            )
            left_bilinear = 4.0 * _bilinear(sites[0], N_MATRICES[alpha])
            right_bilinear = 4.0 * _bilinear(sites[1], N_MATRICES[beta])
            reverse_left = 4.0 * _bilinear(sites[1], N_MATRICES[alpha])
            reverse_right = 4.0 * _bilinear(sites[0], N_MATRICES[beta])
            majorana_hamiltonian += 0.125 * (
                tensor[alpha - 1, beta - 1]
                * left_bilinear
                @ right_bilinear
                + tensor.T[alpha - 1, beta - 1]
                * reverse_left
                @ reverse_right
            )

    for site in (0, 1):
        for gamma in COMPONENTS:
            field_bilinear = 4.0 * _bilinear(
                sites[site], N_MATRICES[gamma]
            )
            gauge_bilinear = 4.0 * _bilinear(
                sites[site], G_MATRICES[gamma]
            )
            spin_hamiltonian -= fields[site][gamma - 1] * spin[site][gamma]
            majorana_hamiltonian -= (
                0.25 * 2.0 * fields[site][gamma - 1] * field_bilinear
            )
            majorana_hamiltonian += (
                0.25 * multipliers[site][gamma - 1] * gauge_bilinear
            )

    projected_spin = (
        physical_projector @ spin_hamiltonian @ physical_projector
    )
    projected_majorana = (
        physical_projector @ majorana_hamiltonian @ physical_projector
    )
    assert np.allclose(
        projected_majorana,
        projected_spin,
        atol=TOL,
        rtol=0.0,
    )
    assert not np.allclose(
        majorana_hamiltonian,
        spin_hamiltonian,
        atol=TOL,
        rtol=0.0,
    )
