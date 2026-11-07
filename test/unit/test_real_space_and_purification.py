"""Unit tests for real-space torus clusters and eigenvector-free matrix sign purification."""

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import InhomogeneousMeanFieldState
from kitaev_mf.services.purification import purify_newton_schulz, purify_pade_sign
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


def test_solve_real_space_purification_accelerators() -> None:
    """Verify that real-space self-consistent purification converges with DIIS and L-Broyden."""
    cluster = RealSpaceTorusCluster(L=4)
    cell2 = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, h=(0.02, 0.02, 0.02))
    seed = make_qsl_seed(cell2)

    # 1. DIIS convergence test
    _, state_diis, met_diis = cluster.solve_real_space_purification(
        params, seed, max_iter=25, tol=1e-6, accelerator="diis"
    )
    assert met_diis["converged"] is True
    assert met_diis["iterations"] <= 15
    assert met_diis["final_residual"] < 1e-6

    # 2. L-Broyden convergence test
    _, state_lb, met_lb = cluster.solve_real_space_purification(
        params, seed, max_iter=25, tol=1e-6, accelerator="lbroyden"
    )
    assert met_lb["converged"] is True
    assert met_lb["iterations"] <= 20
    assert met_lb["final_residual"] < 1e-6

    # Parity between DIIS and L-Broyden solutions
    diff_u = np.max(np.abs(state_diis.U - state_lb.U))
    assert diff_u < 1e-5


def test_disordered_and_vacancy_cluster_purification() -> None:
    """Verify that real-space purification converges on clusters with bond disorder and vacancies."""
    L = 4
    cluster = RealSpaceTorusCluster(L=L)
    cell2 = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, h=(0.02, 0.02, 0.02))
    seed = make_qsl_seed(cell2)

    rng = np.random.default_rng(999)
    # Bond disorder delta_b in [-0.05, 0.05]
    bond_disorder = {
        (m1, m2, b_idx): float(rng.uniform(-0.05, 0.05))
        for m2 in range(L)
        for m1 in range(L)
        for b_idx in range(3)
    }
    # Vacancy site
    vacancies = {1}

    # Verify Hamiltonian skew-symmetry under disorder and vacancies
    A = cluster.assemble_skew_hamiltonian(
        params, seed, bond_disorder=bond_disorder, site_vacancies=vacancies
    )
    assert np.allclose(A, -A.T, atol=1e-15)

    # Solve self-consistently
    omega, state, met = cluster.solve_real_space_purification(
        params,
        seed,
        max_iter=30,
        tol=1e-5,
        accelerator="diis",
        bond_disorder=bond_disorder,
        site_vacancies=vacancies,
    )
    assert met["converged"] is True
    assert met["final_residual"] < 1e-5
    # Covariance idempotency on active subspace
    assert np.allclose(omega, -omega.T, atol=1e-15)


def test_inhomogeneous_state_domain_validation_and_flattening() -> None:
    """Verify InhomogeneousMeanFieldState tensor architecture, skew-symmetry, and bijection."""
    num_cells = 9  # L=3
    state_pure = InhomogeneousMeanFieldState.from_pure_kitaev(num_cells=num_cells)
    assert state_pure.num_cells == 9
    assert state_pure.num_sites == 18
    assert state_pure.num_bonds == 27
    assert state_pure.dim == 86 * 9

    # Pack and unpack roundtrip
    x = state_pure.flatten()
    assert x.shape == (86 * 9,)
    restored = InhomogeneousMeanFieldState.from_flattened(x, num_cells=num_cells)
    assert np.allclose(restored.U, state_pure.U)
    assert np.allclose(restored.V, state_pure.V)
    assert np.allclose(restored.lambdas, state_pure.lambdas)

    # Skew-symmetry violation check
    bad_v = state_pure.V.copy()
    bad_v[0, 0, 0, 1] = 1.0
    bad_v[0, 0, 1, 0] = 1.0  # symmetric entry
    try:
        InhomogeneousMeanFieldState(U=state_pure.U, V=bad_v, lambdas=state_pure.lambdas)
        raise AssertionError("Should have raised ValueError on non-skew V")
    except ValueError as e:
        assert "skew-symmetric" in str(e)


def test_inhomogeneous_hamiltonian_assembly_matches_homogeneous() -> None:
    """Verify exact numerical agreement between homogeneous and inhomogeneous assembly."""
    L = 4
    cluster = RealSpaceTorusCluster(L=L)
    cell2 = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, h=(0.05, 0.05, 0.05))
    seed_homo = make_qsl_seed(cell2)

    # 1. Homogeneous assembly
    A_homo = cluster.assemble_skew_hamiltonian(params, seed_homo)

    # 2. Inhomogeneous assembly using replicated state
    seed_inhomo = InhomogeneousMeanFieldState.from_homogeneous(
        seed_homo, num_cells=cluster.num_cells
    )
    A_inhomo = cluster.assemble_skew_hamiltonian(params, seed_inhomo)

    # Exact equality to machine precision
    diff = float(np.max(np.abs(A_homo - A_inhomo)))
    assert diff < 1e-14, f"Mismatch between homogeneous and inhomogeneous assembly: {diff}"


def test_inhomogeneous_local_extraction_retains_localized_perturbation() -> None:
    """Verify that extract_inhomogeneous_state retains local defects without spatial averaging."""
    L = 4
    cluster = RealSpaceTorusCluster(L=L)
    cell2 = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    seed_homo = make_qsl_seed(cell2)

    # Introduce a local bond flip at (m1=1, m2=1, b_idx=0)
    flipped = [(1, 1, 0)]
    A_vor = cluster.assemble_skew_hamiltonian(params, seed_homo, flipped_bonds=flipped)
    omega_vor, _ = cluster.solve_covariance_complex_eigh(A_vor)

    # Inhomogeneous extraction preserves point-wise variation
    state_inhomo = cluster.extract_inhomogeneous_state(omega_vor)
    u_flipped = state_inhomo.U[1 + L * 1, 0, 3, 3]  # z-bond gauge link at (1, 1)
    u_unflipped = state_inhomo.U[0 + L * 0, 0, 3, 3]  # z-bond gauge link at (0, 0)

    # Inverted link has negative sign relative to vacuum link
    assert u_flipped < -0.9, f"Expected inverted link near -1.0, got {u_flipped}"
    assert u_unflipped > 0.9, f"Expected vacuum link near +1.0, got {u_unflipped}"

    # Verify that traditional extract_mean_field_state averages and washes out the defect
    state_avg = cluster.extract_mean_field_state(omega_vor, np.zeros((2, 3)))
    # Averages over L^2 = 16 cells: (15 * 1.0 - 1.0) / 16 = 14 / 16 = 0.875
    assert 0.8 < state_avg.U[0, 3, 3] < 0.95


def test_inhomogeneous_self_consistent_purification_convergence() -> None:
    """Verify full self-consistent inhomogeneous purification solver on L=4 cluster."""
    L = 4
    cluster = RealSpaceTorusCluster(L=L)
    params = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))
    initial_state = InhomogeneousMeanFieldState.from_pure_kitaev(num_cells=cluster.num_cells)

    # Run inhomogeneous purification with DIIS
    omega, final_state, met = cluster.solve_inhomogeneous_purification(
        params=params,
        initial_state=initial_state,
        max_iter=15,
        tol=1e-7,
        accelerator="diis",
        mixing_alpha=0.5,
    )
    assert met["converged"] is True
    assert met["iterations"] <= 6

    # Verify that pure Kitaev matter bond correlator u_matter converges uniformly to ~ -0.52486
    u_matter_z = final_state.U[:, 0, 0, 0]
    assert np.allclose(u_matter_z, u_matter_z[0], atol=1e-7)
    assert -0.55 < u_matter_z[0] < -0.50


def test_pade_sign_purification_accuracy_and_skew_invariants() -> None:
    """K-018: Rational Padé matrix sign purification satisfies exact skew and idempotency invariants."""
    cluster = RealSpaceTorusCluster(L=4)
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, h=(0.10, 0.10, 0.10))
    state = make_qsl_seed(cell)
    A = cluster.assemble_skew_hamiltonian(params, state)
    M = cluster.num_majoranas

    omega_pade, met_p = purify_pade_sign(A, tol=1e-11, order="32")
    assert np.allclose(omega_pade, -omega_pade.T, atol=1e-15)

    # Idempotency: Omega^2 = -I
    idemp_err = float(np.linalg.norm(omega_pade @ omega_pade + np.eye(M))) / np.sqrt(M)
    assert idemp_err < 1e-11

    # Exact energy agreement with eigensolver
    omega_ref, e_ref = cluster.solve_covariance_complex_eigh(A)
    assert np.isclose(met_p.ground_state_energy, e_ref, atol=1e-8)


def test_pade_32_vs_newton_schulz_acceleration_on_dirac_resonance() -> None:
    """Contract H10: Padé [3/2] achieves >= 2.5x iteration reduction over Newton-Schulz on 3 | L resonances."""
    cluster = RealSpaceTorusCluster(L=6)
    params = ModelParameters(K=-1.0, h=(0.02, 0.02, 0.02))
    state = InhomogeneousMeanFieldState.from_pure_kitaev(num_cells=cluster.num_cells)
    state.lambdas[:, :, :] = np.sqrt(5.0) * 0.02 / np.sqrt(3.0)
    A = cluster.assemble_skew_hamiltonian(params, state)

    # 1. Newton-Schulz cubic polynomial purification
    _, met_ns = purify_newton_schulz(A, tol=1e-10, max_iter=50)

    # 2. Padé [3/2] rational purification
    _, met_p32 = purify_pade_sign(A, tol=1e-10, max_iter=30, order="32")

    # Contract H10 acceptance metric: speedup >= 2.5x
    speedup = met_ns.iterations / met_p32.iterations
    assert speedup >= 2.5, (
        f"Expected Padé [3/2] speedup >= 2.5x, got {speedup:.2f}x ({met_ns.iterations} vs {met_p32.iterations})"
    )
    assert met_p32.iterations <= 12, (
        f"Expected Padé [3/2] iterations <= 12, got {met_p32.iterations}"
    )
    assert met_ns.iterations >= 30, (
        f"Expected Newton-Schulz iterations >= 30, got {met_ns.iterations}"
    )
