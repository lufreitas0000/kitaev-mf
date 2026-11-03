"""Test ModelParameters and bond-directional exchange matrices J_{<ij>_gamma}."""

import numpy as np

from kitaev_mf.domain.parameters import ModelParameters


def test_exchange_matrix_z_bond_matches_manuscript() -> None:
    """Verify J_{<ij>_z} matches Freitas et al. (2024) Eq. (1) & Eq. (13) with Gamma'."""
    params = ModelParameters(J=0.1, K=-1.0, Gamma=0.2, Gamma_prime=0.05, D=0.03)
    j_z = params.exchange_matrix(3)
    expected_z = np.array(
        [
            [0.1, 0.2 + 0.03, 0.05],
            [0.2 - 0.03, 0.1, 0.05],
            [0.05, 0.05, 0.1 - 1.0],
        ],
        dtype=np.float64,
    )
    assert np.allclose(j_z, expected_z, atol=1e-15)


def test_exchange_matrix_cyclic_permutation() -> None:
    """Verify x and y bond exchange matrices follow by cyclic permutation of (1, 2, 3)."""
    params = ModelParameters(J=0.15, K=-1.0, Gamma=0.25, Gamma_prime=-0.04, D=0.02)
    j_x = params.exchange_matrix(1)
    j_y = params.exchange_matrix(2)
    j_z = params.exchange_matrix(3)

    p_cyc = np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    assert np.allclose(j_x, p_cyc.T @ j_z @ p_cyc, atol=1e-15)
    assert np.allclose(j_y, p_cyc @ j_z @ p_cyc.T, atol=1e-15)


def test_anisotropic_kitaev_couplings() -> None:
    """Verify bond-dependent tuple couplings K = (Kx, Ky, Kz) as in Eq. (9) of Freitas et al. (2024)."""
    params = ModelParameters(K=(-1.0, -1.0, -1.2))
    assert np.isclose(params.exchange_matrix(1)[0, 0], -1.0, atol=1e-15)
    assert np.isclose(params.exchange_matrix(2)[1, 1], -1.0, atol=1e-15)
    assert np.isclose(params.exchange_matrix(3)[2, 2], -1.2, atol=1e-15)
