"""Finite-Fock-space checks of the Gaussian exchange energy and double-counting term."""

import numpy as np
from numpy.typing import NDArray

from kitaev_mf.domain.constants import G_MATRICES, N_MATRICES

IDENTITY_2 = np.eye(2, dtype=np.complex128)
PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PAULI_Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
PAULI_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
COMPONENTS = (1, 2, 3)
TOL = 3.0e-12


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


def _expectation(
    state: NDArray[np.complex128],
    operator: NDArray[np.complex128],
) -> complex:
    return complex(np.vdot(state, operator @ state))


def _gaussian_state(
    majoranas: tuple[NDArray[np.complex128], ...],
) -> NDArray[np.complex128]:
    rng = np.random.default_rng(1403)
    raw = rng.standard_normal((len(majoranas), len(majoranas)))
    matrix = 0.15 * (raw - raw.T)
    hamiltonian = np.zeros_like(majoranas[0])
    for mu, c_mu in enumerate(majoranas):
        for nu, c_nu in enumerate(majoranas):
            hamiltonian += 0.25j * matrix[mu, nu] * (c_mu @ c_nu)
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    assert eigenvalues[1] - eigenvalues[0] > 1.0e-10
    return eigenvectors[:, 0]


def _covariances(
    state: NDArray[np.complex128],
    majoranas: tuple[NDArray[np.complex128], ...],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    bond = np.zeros((4, 4), dtype=np.float64)
    onsite = np.zeros((2, 4, 4), dtype=np.float64)
    for site in (0, 1):
        for mu in range(4):
            for nu in range(4):
                if mu != nu:
                    value = _expectation(
                        state,
                        1.0j * majoranas[4 * site + mu] @ majoranas[4 * site + nu],
                    )
                    assert abs(value.imag) < TOL
                    onsite[site, mu, nu] = value.real
    for mu in range(4):
        for nu in range(4):
            value = _expectation(
                state,
                1.0j * majoranas[mu] @ majoranas[4 + nu],
            )
            assert abs(value.imag) < TOL
            bond[mu, nu] = value.real
    return bond, onsite


def _spin_operators(
    majoranas: tuple[NDArray[np.complex128], ...],
) -> tuple[dict[int, NDArray[np.complex128]], dict[int, NDArray[np.complex128]]]:
    operators: list[dict[int, NDArray[np.complex128]]] = []
    for site in (0, 1):
        offset = 4 * site
        operators.append(
            {
                alpha: 1.0j * majoranas[offset] @ majoranas[offset + alpha]
                for alpha in COMPONENTS
            }
        )
    return operators[0], operators[1]


def _energy_components(
    bond: NDArray[np.float64],
    onsite: NDArray[np.float64],
    exchange: NDArray[np.float64],
    fields: NDArray[np.float64],
    multipliers: NDArray[np.float64],
) -> tuple[float, float, float, NDArray[np.float64], NDArray[np.float64]]:
    moments = np.array(
        [
            [np.trace(onsite[site].T @ N_MATRICES[a]) for a in COMPONENTS]
            for site in (0, 1)
        ],
        dtype=np.float64,
    )
    gauge_residuals = np.array(
        [
            [np.trace(G_MATRICES[a] @ onsite[site]) for a in COMPONENTS]
            for site in (0, 1)
        ],
        dtype=np.float64,
    )
    exchange_energy = 0.0
    for alpha in COMPONENTS:
        for beta in COMPONENTS:
            contraction = np.trace(
                N_MATRICES[alpha].T
                @ bond
                @ N_MATRICES[beta]
                @ bond.T
            )
            exchange_energy += exchange[alpha - 1, beta - 1] * (
                0.25 * moments[0, alpha - 1] * moments[1, beta - 1]
                - 0.5 * contraction
            )
    field_energy = -0.5 * float(np.sum(fields * moments))
    multiplier_energy = -0.25 * float(np.sum(multipliers * gauge_residuals))
    return (
        float(exchange_energy),
        field_energy,
        multiplier_energy,
        moments,
        gauge_residuals,
    )


def _quadratic_hamiltonian(
    majoranas: tuple[NDArray[np.complex128], ...],
    bond: NDArray[np.float64],
    moments: NDArray[np.float64],
    exchange: NDArray[np.float64],
    fields: NDArray[np.float64],
    multipliers: NDArray[np.float64],
) -> NDArray[np.complex128]:
    onsite_blocks = np.zeros((2, 4, 4), dtype=np.float64)
    for alpha in COMPONENTS:
        for beta in COMPONENTS:
            onsite_blocks[0] += (
                exchange[alpha - 1, beta - 1]
                * moments[1, beta - 1]
                * N_MATRICES[alpha]
            )
            onsite_blocks[1] += (
                exchange[beta - 1, alpha - 1]
                * moments[0, beta - 1]
                * N_MATRICES[alpha]
            )
    for site in (0, 1):
        for gamma in COMPONENTS:
            onsite_blocks[site] += (
                multipliers[site, gamma - 1] * G_MATRICES[gamma]
                - 2.0 * fields[site, gamma - 1] * N_MATRICES[gamma]
            )

    hopping = np.zeros((4, 4), dtype=np.float64)
    for alpha in COMPONENTS:
        for beta in COMPONENTS:
            hopping += (
                2.0
                * exchange[alpha - 1, beta - 1]
                * N_MATRICES[alpha]
                @ bond
                @ N_MATRICES[beta]
            )

    matrix = np.zeros((8, 8), dtype=np.float64)
    matrix[:4, :4] = onsite_blocks[0]
    matrix[4:, 4:] = onsite_blocks[1]
    matrix[:4, 4:] = hopping
    matrix[4:, :4] = -hopping.T

    hamiltonian = np.zeros_like(majoranas[0])
    for mu, c_mu in enumerate(majoranas):
        for nu, c_nu in enumerate(majoranas):
            hamiltonian += 0.25j * matrix[mu, nu] * (c_mu @ c_nu)
    return hamiltonian


def test_energetics_wick_exchange_energy_matches_direct_spin_expectation() -> None:
    majoranas = _majorana_operators(4)
    state = _gaussian_state(majoranas)
    bond, onsite = _covariances(state, majoranas)
    spins_i, spins_j = _spin_operators(majoranas)
    exchange = np.array(
        [[0.13, 0.21, -0.17], [-0.08, -0.23, 0.19], [0.31, -0.11, 0.07]],
        dtype=np.float64,
    )
    fields = np.zeros((2, 3), dtype=np.float64)
    multipliers = np.zeros((2, 3), dtype=np.float64)

    exchange_energy, _, _, _, _ = _energy_components(
        bond, onsite, exchange, fields, multipliers
    )
    direct_energy = 0.0
    for alpha in COMPONENTS:
        for beta in COMPONENTS:
            direct_energy += exchange[alpha - 1, beta - 1] * _expectation(
                state,
                spins_i[alpha] @ spins_j[beta],
            ).real

    assert np.isclose(exchange_energy, direct_energy, atol=TOL, rtol=0.0)


def test_energetics_double_counting_constant_reproduces_gaussian_energy() -> None:
    majoranas = _majorana_operators(4)
    state = _gaussian_state(majoranas)
    bond, onsite = _covariances(state, majoranas)
    exchange = np.array(
        [[0.13, 0.21, -0.17], [-0.08, -0.23, 0.19], [0.31, -0.11, 0.07]],
        dtype=np.float64,
    )
    fields = np.array([[0.05, -0.03, 0.02], [-0.04, 0.06, 0.01]])
    multipliers = np.array([[0.03, -0.02, 0.04], [-0.01, 0.05, -0.03]])

    exchange_energy, field_energy, multiplier_energy, moments, _ = _energy_components(
        bond, onsite, exchange, fields, multipliers
    )
    quadratic = _quadratic_hamiltonian(
        majoranas, bond, moments, exchange, fields, multipliers
    )
    quadratic_expectation = _expectation(state, quadratic)

    assert abs(quadratic_expectation.imag) < TOL
    assert np.isclose(
        quadratic_expectation.real,
        2.0 * exchange_energy + field_energy + multiplier_energy,
        atol=TOL,
        rtol=0.0,
    )
    scalar_constant = exchange_energy
    assert np.isclose(
        quadratic_expectation.real - scalar_constant,
        exchange_energy + field_energy + multiplier_energy,
        atol=TOL,
        rtol=0.0,
    )


def test_energetics_pure_kitaev_energy_normalization_per_bond_and_site() -> None:
    k_value = -1.0
    matter_correlator = -0.5248645
    gauge_correlator = 1.0
    bond_energies = []

    for gamma in COMPONENTS:
        bond = np.zeros((4, 4), dtype=np.float64)
        bond[0, 0] = matter_correlator
        bond[gamma, gamma] = gauge_correlator
        contraction = float(
            np.trace(N_MATRICES[gamma].T @ bond @ N_MATRICES[gamma] @ bond.T)
        )
        pauli_correlator = -0.5 * contraction
        bond_energies.append(k_value * pauli_correlator)

    energy_per_cell = float(sum(bond_energies))
    energy_per_site = energy_per_cell / 2.0

    assert np.allclose(bond_energies, [-0.5248645] * 3, atol=1.0e-12, rtol=0.0)
    assert np.isclose(energy_per_cell, -1.5745935, atol=1.0e-12, rtol=0.0)
    assert np.isclose(energy_per_site, -0.78729675, atol=1.0e-12, rtol=0.0)
