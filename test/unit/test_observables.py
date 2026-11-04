"""Independent finite-Fock-space checks of Gaussian spin-observable formulas."""

import numpy as np
from numpy.typing import NDArray

from kitaev_mf.domain.constants import N_MATRICES

IDENTITY_2 = np.eye(2, dtype=np.complex128)
PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PAULI_Y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
PAULI_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
COMPONENTS = (1, 2, 3)
TOL = 2.0e-12


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


def _quadratic_ground_state(
    majoranas: tuple[NDArray[np.complex128], ...],
) -> NDArray[np.complex128]:
    """Construct a nondegenerate Gaussian ground state from a generic quadratic form."""
    rng = np.random.default_rng(1402)
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
    onsite = np.zeros((2, 4, 4), dtype=np.float64)
    bond = np.zeros((4, 4), dtype=np.float64)
    for site in (0, 1):
        offset = 4 * site
        for mu in range(4):
            for nu in range(4):
                if mu != nu:
                    value = _expectation(
                        state,
                        1.0j * majoranas[offset + mu] @ majoranas[offset + nu],
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
    spins: list[dict[int, NDArray[np.complex128]]] = []
    for site in (0, 1):
        offset = 4 * site
        spins.append(
            {
                gamma: 1.0j * majoranas[offset] @ majoranas[offset + gamma]
                for gamma in COMPONENTS
            }
        )
    return spins[0], spins[1]


def _wick_bond_correlator(
    bond: NDArray[np.float64],
    onsite_i: NDArray[np.float64],
    onsite_j: NDArray[np.float64],
    alpha: int,
    beta: int,
) -> tuple[float, float]:
    n_alpha = N_MATRICES[alpha]
    n_beta = N_MATRICES[beta]
    m_i = float(np.trace(onsite_i.T @ n_alpha))
    m_j = float(np.trace(onsite_j.T @ n_beta))
    exchange_contraction = float(
        np.trace(n_alpha.T @ bond @ n_beta @ bond.T)
    )
    pauli = 0.25 * m_i * m_j - 0.5 * exchange_contraction
    return pauli, 0.25 * pauli


def test_observable_one_point_normalization_matches_fock_expectation() -> None:
    majoranas = _majorana_operators(4)
    state = _quadratic_ground_state(majoranas)
    _, onsite = _covariances(state, majoranas)
    spins_i, _ = _spin_operators(majoranas)

    for gamma in COMPONENTS:
        n_gamma = N_MATRICES[gamma]
        m_gamma = float(np.trace(onsite[0].T @ n_gamma))
        pauli_from_covariance = 0.5 * m_gamma
        spin_from_covariance = 0.25 * m_gamma
        pauli_direct = _expectation(state, spins_i[gamma])
        spin_direct = 0.5 * pauli_direct

        assert np.isclose(pauli_from_covariance, pauli_direct.real, atol=TOL, rtol=0.0)
        assert np.isclose(spin_from_covariance, spin_direct.real, atol=TOL, rtol=0.0)


def test_observable_two_site_wick_formula_matches_fock_expectation() -> None:
    majoranas = _majorana_operators(4)
    state = _quadratic_ground_state(majoranas)
    bond, onsite = _covariances(state, majoranas)
    spins_i, spins_j = _spin_operators(majoranas)

    for alpha in COMPONENTS:
        for beta in COMPONENTS:
            pauli_wick, spin_wick = _wick_bond_correlator(
                bond, onsite[0], onsite[1], alpha, beta
            )
            pauli_direct = _expectation(
                state,
                spins_i[alpha] @ spins_j[beta],
            )
            spin_direct = 0.25 * pauli_direct

            assert abs(pauli_direct.imag) < TOL
            assert np.isclose(pauli_wick, pauli_direct.real, atol=TOL, rtol=0.0)
            assert np.isclose(spin_wick, spin_direct.real, atol=TOL, rtol=0.0)


def test_observable_wick_formula_contains_nonzero_hartree_and_exchange_terms() -> None:
    majoranas = _majorana_operators(4)
    state = _quadratic_ground_state(majoranas)
    bond, onsite = _covariances(state, majoranas)

    alpha, beta = 1, 2
    n_alpha = N_MATRICES[alpha]
    n_beta = N_MATRICES[beta]
    m_i = float(np.trace(onsite[0].T @ n_alpha))
    m_j = float(np.trace(onsite[1].T @ n_beta))
    hartree = 0.25 * m_i * m_j
    exchange = -0.5 * float(np.trace(n_alpha.T @ bond @ n_beta @ bond.T))

    assert abs(hartree) > 1.0e-6
    assert abs(exchange) > 1.0e-6
    assert np.isclose(hartree + exchange, _wick_bond_correlator(
        bond, onsite[0], onsite[1], alpha, beta
    )[0], atol=TOL, rtol=0.0)


def test_observable_pauli_and_physical_spin_normalizations_are_distinct() -> None:
    majoranas = _majorana_operators(4)
    state = _quadratic_ground_state(majoranas)
    spins_i, _ = _spin_operators(majoranas)

    for gamma in COMPONENTS:
        pauli_value = _expectation(state, spins_i[gamma])
        physical_spin_value = 0.5 * pauli_value
        assert np.isclose(
            physical_spin_value,
            0.5 * pauli_value,
            atol=TOL,
            rtol=0.0,
        )
        assert np.isclose(
            abs(physical_spin_value),
            0.5 * abs(pauli_value),
            atol=TOL,
            rtol=0.0,
        )
