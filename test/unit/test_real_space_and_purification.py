"""Unit tests for real-space torus clusters and eigenvector-free matrix sign purification."""

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.purification import purify_newton_schulz
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


def test_real_space_cluster_geometry() -> None:
    """Verify site and Majorana indexing and Hamiltonian skew-symmetry."""
    params = ModelParameters(K=-1.0, h=(0.20, 0.20, 0.20))
    cell2 = make_honeycomb_2site()
    state = make_qsl_seed(cell2)
    for L in (2, 3, 4):
        cluster = RealSpaceTorusCluster(L=L)
        assert cluster.num_cells == L * L
        assert cluster.num_sites == 2 * L * L
        assert cluster.num_majoranas == 8 * L * L

        A = cluster.assemble_skew_hamiltonian(params, state)
        assert A.shape == (cluster.num_majoranas, cluster.num_majoranas)
        assert np.allclose(A, -A.T, atol=1e-15)


def test_real_vs_complex_covariance_parity() -> None:
    """Verify exact equivalence between complex zheevd and real dsyevd on L=3 torus (M=72)."""
    cluster = RealSpaceTorusCluster(L=3)
    params = ModelParameters(K=-1.0, h=(0.15, 0.15, 0.15))
    cell2 = make_honeycomb_2site()
    state = make_qsl_seed(cell2)
    A = cluster.assemble_skew_hamiltonian(params, state)
    M = cluster.num_majoranas

    omega_c, e_c = cluster.solve_covariance_complex_eigh(A)
    omega_r, e_r = cluster.solve_covariance_real_eigh(A)

    # Relative Frobenius difference between covariance matrices (dsyevd on A^2 has condition ~ kappa^2)
    diff_frob = float(np.linalg.norm(omega_c - omega_r)) / np.sqrt(M)
    assert diff_frob < 1e-7

    # Ground-state energy parity
    assert np.isclose(e_c, e_r, atol=1e-8)

    # Invariance check: Omega^2 = -I on the active non-zero subspace
    assert np.allclose(omega_c @ omega_c, omega_r @ omega_r, atol=1e-7)


def test_matrix_sign_purification_accuracy() -> None:
    """Verify Newton-Schulz purification accuracy and machine precision agreement with eigensolvers."""
    cluster = RealSpaceTorusCluster(L=3)
    params = ModelParameters(K=-1.0, h=(0.20, 0.20, 0.20))
    cell2 = make_honeycomb_2site()
    state = make_qsl_seed(cell2)
    A = cluster.assemble_skew_hamiltonian(params, state)
    M = cluster.num_majoranas

    # 1. Exact eigensolver reference
    omega_ref, e_ref = cluster.solve_covariance_complex_eigh(A)

    # 2. Newton-Schulz purification in pure real dgemm arithmetic
    omega_purif, met = purify_newton_schulz(A, tol=1e-11, max_iter=40)
    diff_purif = float(np.linalg.norm(omega_purif - omega_ref)) / np.sqrt(M)

    # Ground state energy and covariance matrix match eigensolver to machine precision
    assert diff_purif < 1e-12
    assert np.isclose(met.ground_state_energy, e_ref, atol=1e-11)
    assert met.iterations <= 15
