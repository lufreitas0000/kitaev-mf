"""Unit test suite for Phase 3 Task T3.1.2: Bogoliubov Overlaps & Thouless Matrix Elements.

Audits:
- Hypothesis H_BOG_01: Canonical Bogoliubov unitarity and anti-symmetry relations.
- Hypothesis H_BOG_02: Thouless vacuum overlap vs. full Fock-space exact diagonalization.
- Hypothesis H_BOG_03: Majorana transition matrix elements and unitary sum rule.
- Hypothesis H_BOG_04: Particle-hole resolution of singular transformations (det(X) = 0).
- Physical gauge-quench Bogoliubov spectrum on honeycomb clusters.
"""

import numpy as np
from scipy.linalg import schur

from kitaev_mf.services.phase_3.bogoliubov import (
    compute_bogoliubov_transformation,
    compute_majorana_matrix_elements,
    extract_mode_wavefunctions,
    verify_bogoliubov_unitarity,
)
from kitaev_mf.services.phase_3.gauge_theory import (
    assemble_matter_skew_matrix,
    build_honeycomb_torus_cluster,
    build_uniform_gauge_configuration,
    diagonalize_matter_hamiltonian,
    flip_bond_gauge_link,
)


def _get_mode_phi(A: np.ndarray) -> np.ndarray:
    """Helper to extract standardized mode wavefunctions phi from skew-symmetric A."""
    N = A.shape[0] // 2
    T, Q = schur(A, output="real")
    for k in range(N):
        if T[2 * k, 2 * k + 1] < 0:
            Q[:, [2 * k, 2 * k + 1]] = Q[:, [2 * k + 1, 2 * k]]
    return extract_mode_wavefunctions(Q)


def test_h_bog_01_canonical_unitarity_relations() -> None:
    """Verify Hypothesis H_BOG_01: Bogoliubov blocks X, Y satisfy canonical identities to < 1e-14."""
    np.random.seed(42)
    N = 4
    # Random skew-symmetric matrices
    M1 = np.random.randn(2 * N, 2 * N)
    A1 = M1 - M1.T
    M2 = np.random.randn(2 * N, 2 * N)
    A2 = M2 - M2.T

    phi_a = _get_mode_phi(A1)
    phi_b = _get_mode_phi(A2)

    bog = compute_bogoliubov_transformation(phi_a, phi_b)
    is_valid, errors = verify_bogoliubov_unitarity(bog.X, bog.Y, tol=1e-13)

    assert is_valid is True
    assert errors["max_err_XX_YY"] < 1e-14
    assert errors["max_err_XdagX_YtYstar"] < 1e-14
    assert errors["max_err_XYt_YXt"] < 1e-14


def test_h_bog_02_thouless_vacuum_overlap_vs_fock_exact() -> None:
    """Verify Hypothesis H_BOG_02: Thouless vacuum overlap sqrt(|det(X)|) matches Fock-space exact."""
    # Seed 42 gives det(Qa) = det(Qb) = -1 => det(W) = +1 (regular transformation)
    np.random.seed(42)
    a0 = np.array([[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0]], dtype=np.complex128)
    a1 = np.array([[0, 0, 1, 0], [0, 0, 0, -1], [0, 0, 0, 0], [0, 0, 0, 0]], dtype=np.complex128)
    a0_dag = a0.conj().T
    a1_dag = a1.conj().T

    M1 = np.random.randn(4, 4)
    A1 = M1 - M1.T
    M2 = np.random.randn(4, 4)
    A2 = M2 - M2.T

    phi_a = _get_mode_phi(A1)
    phi_b = _get_mode_phi(A2)
    bog = compute_bogoliubov_transformation(phi_a, phi_b)
    assert bog.is_singular is False

    # Construct b-mode operators in Fock space
    X, Y = bog.X, bog.Y
    b0 = X[0, 0] * a0 + X[0, 1] * a1 + Y[0, 0] * a0_dag + Y[0, 1] * a1_dag
    b1 = X[1, 0] * a0 + X[1, 1] * a1 + Y[1, 0] * a0_dag + Y[1, 1] * a1_dag
    b0_dag = b0.conj().T
    b1_dag = b1.conj().T

    # Find Fock vacuum |0_b> via SVD nullspace
    M_null = np.vstack([b0, b1])
    _, _, vh = np.linalg.svd(M_null)
    psi_0b = vh[-1].conj()
    psi_0b /= np.linalg.norm(psi_0b)

    # Fock vacuum |0_a> is basis vector [1, 0, 0, 0]
    psi_0a = np.array([1, 0, 0, 0], dtype=np.complex128)
    exact_overlap_mag = float(np.abs(np.dot(psi_0b.conj(), psi_0a)))

    # Thouless predicted overlap
    thouless_overlap_mag = bog.vacuum_overlap_magnitude
    assert np.isclose(exact_overlap_mag, thouless_overlap_mag, atol=1e-12)
    assert exact_overlap_mag > 0.1

    # Verify Thouless state reconstruction of |0_a> in the b-basis
    assert bog.Z_b is not None
    overlap_complex = np.dot(psi_0b.conj(), psi_0a)
    psi_thouless = overlap_complex * (psi_0b + bog.Z_b[0, 1] * (b0_dag @ b1_dag @ psi_0b))
    assert np.allclose(psi_thouless, psi_0a, atol=1e-13)


def test_h_bog_03_majorana_matrix_elements_and_sum_rule() -> None:
    """Verify Hypothesis H_BOG_03: Transition matrix elements satisfy completeness sum rule."""
    np.random.seed(42)
    a0 = np.array([[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0]], dtype=np.complex128)
    a1 = np.array([[0, 0, 1, 0], [0, 0, 0, -1], [0, 0, 0, 0], [0, 0, 0, 0]], dtype=np.complex128)
    a0_dag = a0.conj().T
    a1_dag = a1.conj().T

    M1 = np.random.randn(4, 4)
    A1 = M1 - M1.T
    M2 = np.random.randn(4, 4)
    A2 = M2 - M2.T

    phi_a = _get_mode_phi(A1)
    phi_b = _get_mode_phi(A2)
    bog = compute_bogoliubov_transformation(phi_a, phi_b)

    X, Y = bog.X, bog.Y
    b0 = X[0, 0] * a0 + X[0, 1] * a1 + Y[0, 0] * a0_dag + Y[0, 1] * a1_dag
    b1 = X[1, 0] * a0 + X[1, 1] * a1 + Y[1, 0] * a0_dag + Y[1, 1] * a1_dag
    b0_dag = b0.conj().T
    b1_dag = b1.conj().T

    M_null = np.vstack([b0, b1])
    _, _, vh = np.linalg.svd(M_null)
    psi_0b = vh[-1].conj()
    psi_0b /= np.linalg.norm(psi_0b)

    # Initial state: one-fermion state |psi_init> = a0_dag |0_a>
    psi_init = a0_dag @ np.array([1, 0, 0, 0], dtype=np.complex128)

    # Majorana operator c_0
    c0 = 2.0 * (
        phi_a[0, 0].conj() * a0
        + phi_a[0, 0] * a0_dag
        + phi_a[1, 0].conj() * a1
        + phi_a[1, 0] * a1_dag
    )

    c_psi = c0 @ psi_init
    # Exact completeness norm: || c_0 |psi> ||^2 must be 1.0 since {c, c} = 2 => c^2 = 1
    total_norm = float(np.linalg.norm(c_psi) ** 2)
    assert np.isclose(total_norm, 1.0, atol=1e-12)

    # Compare analytical matrix elements
    overlap_phase = np.dot(psi_0b.conj(), np.array([1, 0, 0, 0], dtype=np.complex128))
    phase_normalized = overlap_phase / abs(overlap_phase)

    res = compute_majorana_matrix_elements(
        bog=bog,
        phi_a=phi_a,
        phi_b=phi_b,
        site_j=0,
        initial_parity=-1,
        occupied_mode=0,
        overlap_phase=phase_normalized,
    )

    # Check vacuum component
    mat_elem_fock_0 = np.dot(psi_0b.conj(), c_psi)
    mat_elem_analytic_0 = res["amplitudes"][()]
    assert np.isclose(mat_elem_fock_0, mat_elem_analytic_0, atol=1e-12)

    # Check 2-particle component: <2_b| = <0_b| b1 b0
    psi_2b = b0_dag @ b1_dag @ psi_0b
    mat_elem_fock_2 = np.dot(psi_2b.conj(), c_psi)
    mat_elem_analytic_2 = res["amplitudes"][(0, 1)]
    assert np.isclose(mat_elem_fock_2, mat_elem_analytic_2, atol=1e-12)

    # Completeness sum rule: sum of squares over all even channels must be 1.0
    assert np.isclose(res["sum_of_squares"], 1.0, atol=1e-12)

    # Also verify initial_parity = +1 (vacuum ground state |0_a>)
    psi_init_even = np.array([1, 0, 0, 0], dtype=np.complex128)
    c_psi_even = c0 @ psi_init_even
    res_even = compute_majorana_matrix_elements(
        bog=bog,
        phi_a=phi_a,
        phi_b=phi_b,
        site_j=0,
        initial_parity=+1,
        overlap_phase=phase_normalized,
    )
    psi_1b_0 = b0_dag @ psi_0b
    psi_1b_1 = b1_dag @ psi_0b
    fock_even_0 = np.dot(psi_1b_0.conj(), c_psi_even)
    fock_even_1 = np.dot(psi_1b_1.conj(), c_psi_even)
    assert np.isclose(fock_even_0, res_even["amplitudes"][(0,)], atol=1e-12)
    assert np.isclose(fock_even_1, res_even["amplitudes"][(1,)], atol=1e-12)
    assert np.isclose(res_even["sum_of_squares"], 1.0, atol=1e-12)


def test_h_bog_04_particle_hole_resolution_of_singular_transformation() -> None:
    """Verify Hypothesis H_BOG_04: Particle-hole interchange resolves det(X) = 0 and matches Fock space."""
    # Seed 101 gives det(Qa) = -1, det(Qb) = +1 => det(W) = -1 (singular transformation)
    np.random.seed(101)
    a0 = np.array([[0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0]], dtype=np.complex128)
    a1 = np.array([[0, 0, 1, 0], [0, 0, 0, -1], [0, 0, 0, 0], [0, 0, 0, 0]], dtype=np.complex128)
    a0_dag = a0.conj().T
    a1_dag = a1.conj().T

    M1 = np.random.randn(4, 4)
    A1 = M1 - M1.T
    M2 = np.random.randn(4, 4)
    A2 = M2 - M2.T

    phi_a = _get_mode_phi(A1)
    phi_b = _get_mode_phi(A2)
    bog = compute_bogoliubov_transformation(phi_a, phi_b)

    # Must be detected as singular
    assert bog.is_singular is True
    assert bog.vacuum_overlap_magnitude == 0.0

    # Build b-operators and Fock vacuum |0_b>
    X, Y = bog.X, bog.Y
    b0 = X[0, 0] * a0 + X[0, 1] * a1 + Y[0, 0] * a0_dag + Y[0, 1] * a1_dag
    b1 = X[1, 0] * a0 + X[1, 1] * a1 + Y[1, 0] * a0_dag + Y[1, 1] * a1_dag
    b0_dag = b0.conj().T
    b1_dag = b1.conj().T

    M_null = np.vstack([b0, b1])
    _, _, vh = np.linalg.svd(M_null)
    psi_0b = vh[-1].conj()
    psi_0b /= np.linalg.norm(psi_0b)

    # Fock vacuum |0_a>
    psi_0a = np.array([1, 0, 0, 0], dtype=np.complex128)
    # Vacuum overlap vanishes
    assert np.isclose(float(np.abs(np.dot(psi_0b.conj(), psi_0a))), 0.0, atol=1e-12)

    # Single-quasiparticle state for the singular mode mu: |1_b> = b_mu_dag |0_b>
    mu = bog.singular_mode
    assert mu in (0, 1)
    b_mu_dag = b0_dag if mu == 0 else b1_dag
    psi_1b = b_mu_dag @ psi_0b
    psi_1b /= np.linalg.norm(psi_1b)
    exact_1b_overlap = float(np.abs(np.dot(psi_1b.conj(), psi_0a)))
    assert exact_1b_overlap > 0.1

    # Predicted overlap from particle-hole transformed det(X')
    assert bog.singular_overlap_magnitude is not None
    assert np.isclose(exact_1b_overlap, bog.singular_overlap_magnitude, atol=1e-12)

    # Verify matrix elements in singular transformation against exact Fock space
    # 1. Even initial parity (|0_a>): connects to odd quasiparticles of tilde{b}
    c0 = 2.0 * (
        phi_a[0, 0].conj() * a0
        + phi_a[0, 0] * a0_dag
        + phi_a[1, 0].conj() * a1
        + phi_a[1, 0] * a1_dag
    )
    c_psi_0a = c0 @ psi_0a
    overlap_complex_sing = np.dot(psi_1b.conj(), psi_0a)
    phase_sing = overlap_complex_sing / abs(overlap_complex_sing)

    res_sing_even = compute_majorana_matrix_elements(
        bog=bog,
        phi_a=phi_a,
        phi_b=phi_b,
        site_j=0,
        initial_parity=+1,
        overlap_phase=phase_sing,
    )
    psi_tilde_1_0 = (b0_dag if mu != 0 else np.eye(4, dtype=np.complex128)) @ (
        psi_1b if mu != 0 else psi_0b
    )
    psi_tilde_1_1 = (b1_dag if mu != 1 else np.eye(4, dtype=np.complex128)) @ (
        psi_1b if mu != 1 else psi_0b
    )
    fock_sing_0 = np.dot(psi_tilde_1_0.conj(), c_psi_0a)
    fock_sing_1 = np.dot(psi_tilde_1_1.conj(), c_psi_0a)
    assert np.isclose(fock_sing_0, res_sing_even["amplitudes"][(0,)], atol=1e-12)
    assert np.isclose(fock_sing_1, res_sing_even["amplitudes"][(1,)], atol=1e-12)
    assert np.isclose(res_sing_even["sum_of_squares"], 1.0, atol=1e-12)

    # 2. Odd initial parity (|1_a> = a0_dag |0_a>): connects to even quasiparticles of tilde{b}
    psi_1a = a0_dag @ psi_0a
    c_psi_1a = c0 @ psi_1a
    res_sing_odd = compute_majorana_matrix_elements(
        bog=bog,
        phi_a=phi_a,
        phi_b=phi_b,
        site_j=0,
        initial_parity=-1,
        occupied_mode=0,
        overlap_phase=phase_sing,
    )
    fock_sing_odd_0 = np.dot(psi_1b.conj(), c_psi_1a)
    assert np.isclose(fock_sing_odd_0, res_sing_odd["amplitudes"][()], atol=1e-12)
    psi_tilde_2 = (b0_dag if mu == 1 else b1_dag) @ psi_0b
    fock_sing_odd_2 = np.dot(psi_tilde_2.conj(), c_psi_1a)
    assert np.isclose(fock_sing_odd_2, res_sing_odd["amplitudes"][(0, 1)], atol=1e-12)
    assert np.isclose(res_sing_odd["sum_of_squares"], 1.0, atol=1e-12)


def test_h_bog_05_physical_quench_on_honeycomb_cluster() -> None:
    """Verify Bogoliubov transformation across physical sudden spin quench on honeycomb torus."""
    L = 3
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    gauge_0 = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    kappa = 0.15
    K = 1.0
    A_0 = assemble_matter_skew_matrix(gauge_0, K=K, kappa=kappa, triplets=triplets)
    diag_0 = diagonalize_matter_hamiltonian(A_0, gauge=gauge_0, boundary="torus", L=L)

    # Perform sudden spin quench: flip bond 0
    sa = bonds[0].site_a
    sb = bonds[0].site_b
    gauge_quenched = flip_bond_gauge_link(gauge_0, site_i=sa, site_j=sb)
    A_quenched = assemble_matter_skew_matrix(gauge_quenched, K=K, kappa=kappa, triplets=triplets)
    diag_quenched = diagonalize_matter_hamiltonian(
        A_quenched, gauge=gauge_quenched, boundary="torus", L=L
    )

    phi_0 = extract_mode_wavefunctions(diag_0.Q)
    phi_quenched = extract_mode_wavefunctions(diag_quenched.Q)

    # Bogoliubov transformation between clean vacuum and quenched 2-vortex sector
    bog = compute_bogoliubov_transformation(phi_0, phi_quenched)
    is_valid, errors = verify_bogoliubov_unitarity(bog.X, bog.Y, tol=1e-12)

    assert is_valid is True
    assert errors["max_err_XX_YY"] < 1e-13
    assert errors["max_err_XdagX_YtYstar"] < 1e-13
    assert errors["max_err_XYt_YXt"] < 1e-13

    # Overlap magnitude is bounded in [0, 1]
    assert 0.0 <= bog.vacuum_overlap_magnitude <= 1.0
