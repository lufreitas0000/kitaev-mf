"""Independent numerical check of primitive/rectangular Bloch folding."""

import numpy as np
from numpy.typing import NDArray

ComplexMatrix = NDArray[np.complex128]


def _primitive_bloch(k: NDArray[np.float64], coupling: float = 1.0) -> ComplexMatrix:
    """Hermitian matter Majorana generator iA(k), in the (A, B) basis."""
    k1, k2 = k
    structure_factor = (
        1.0
        + np.exp(-2.0j * np.pi * k1)
        + np.exp(-2.0j * np.pi * k2)
    )
    hopping = 2.0 * coupling * structure_factor
    return np.array(
        [[0.0, 1.0j * hopping], [-1.0j * np.conjugate(hopping), 0.0]],
        dtype=np.complex128,
    )


def _rectangular_bloch(
    kappa: NDArray[np.float64], coupling: float = 1.0
) -> ComplexMatrix:
    """Hermitian matter Majorana generator iA(kappa), ordered (A0, B0, A1, B1)."""
    kappa1, kappa2 = kappa
    hopping = np.zeros((2, 2), dtype=np.complex128)

    # (A-coset, B-coset, rectangular-cell offset).
    bonds = (
        (0, 0, (0, 0)),   # z: A0 -> B0
        (1, 1, (0, 0)),   # z: A1 -> B1
        (0, 1, (0, -1)),  # x: A0 -> B1
        (1, 0, (1, 0)),   # x: A1 -> B0
        (1, 0, (0, 0)),   # y: A1 -> B0
        (0, 1, (-1, -1)), # y: A0 -> B1
    )

    for row, column, (m1, m2) in bonds:
        phase = np.exp(2.0j * np.pi * (kappa1 * m1 + kappa2 * m2))
        hopping[row, column] += 2.0 * coupling * phase

    matrix = np.zeros((4, 4), dtype=np.complex128)
    for a_coset in range(2):
        for b_coset in range(2):
            a_index = 2 * a_coset
            b_index = 2 * b_coset + 1
            matrix[a_index, b_index] = 1.0j * hopping[a_coset, b_coset]
            matrix[b_index, a_index] = -1.0j * np.conjugate(
                hopping[a_coset, b_coset]
            )
    return matrix


def _folding_unitary(k1_a: float) -> ComplexMatrix:
    """Transform folded primitive momenta to the rectangular-cell basis."""
    phase = np.exp(2.0j * np.pi * k1_a)
    identity = np.eye(2, dtype=np.complex128)
    return np.block(
        [[identity, identity], [phase * identity, -phase * identity]]
    ) / np.sqrt(2.0)


def _folded_momenta(kappa: NDArray[np.float64]) -> tuple[NDArray[np.float64], ...]:
    """Return the two primitive momenta folded onto rectangular momentum kappa."""
    kappa1, kappa2 = kappa
    k_a = np.array(
        [(kappa1 + kappa2) / 2.0, (kappa2 - kappa1) / 2.0],
        dtype=np.float64,
    )
    k_b = np.mod(k_a + 0.5, 1.0)
    return k_a, k_b


def test_folding_unitary_is_unitary() -> None:
    for kappa in (
        np.array([0.13, 0.27]),
        np.array([0.0, 0.0]),
        np.array([0.5, 0.0]),
        np.array([0.5, 0.5]),
    ):
        k_a, _ = _folded_momenta(kappa)
        unitary = _folding_unitary(float(k_a[0]))
        np.testing.assert_allclose(
            unitary.conj().T @ unitary,
            np.eye(4, dtype=np.complex128),
            atol=1e-13,
            rtol=0.0,
        )


def test_rectangular_bloch_matrix_is_unitarily_folded_primitive_matrix() -> None:
    rng = np.random.default_rng(20261104)
    generic_points = rng.random((20, 2))
    boundary_points = np.array(
        [
            [0.0, 0.0],
            [0.5, 0.0],
            [0.0, 0.5],
            [0.5, 0.5],
            [0.0, 0.23],
            [0.5, 0.23],
            [0.37, 0.0],
            [0.37, 0.5],
        ],
        dtype=np.float64,
    )

    for kappa in np.vstack((generic_points, boundary_points)):
        k_a, k_b = _folded_momenta(kappa)
        unitary = _folding_unitary(float(k_a[0]))
        folded = np.zeros((4, 4), dtype=np.complex128)
        folded[:2, :2] = _primitive_bloch(k_a)
        folded[2:, 2:] = _primitive_bloch(k_b)

        reconstructed = unitary @ folded @ unitary.conj().T
        np.testing.assert_allclose(
            _rectangular_bloch(kappa),
            reconstructed,
            atol=2e-12,
            rtol=2e-12,
        )
