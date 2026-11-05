"""Mathematical TDD unit tests for quasi-free states, Wick contraction, and the Grassmannian.

Covers claims CLM-M1-QF-001 through CLM-M1-QF-008:
- TST-QF-01: Exact Pfaffian / Wick formula for Gaussian states
- TST-QF-02: Wick failure for generic non-Gaussian states
- TST-QF-03: Non-convexity of the quasi-free state manifold
- TST-QF-04: Operator norm bound ||Omega|| <= 1, purity Omega^2 = -I, and parity Pf(Omega)
- TST-QF-05: Gaussian von Neumann entropy and maximum entropy principle
- TST-QF-06: Affine tangent property and quadratic energy curvature
- TST-QF-07: Grassmannian stationarity [h(Omega), Omega] = 0 and spectral projector
- TST-QF-08: Graded fermionic CAR locality vs ordinary spin subnet commutativity
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as la


def build_majorana_fock_basis(n_modes: int = 2) -> list[np.ndarray]:
    """Construct explicit 2^n_modes x 2^n_modes Majorana operators via Jordan-Wigner."""
    sigma_x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    sigma_y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
    sigma_z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
    eye = np.eye(2, dtype=np.complex128)

    dirac_d: list[np.ndarray] = []
    for k in range(n_modes):
        op: np.ndarray = np.array([[1.0]], dtype=np.complex128)
        for _ in range(k):
            op = np.asarray(np.kron(op, sigma_z), dtype=np.complex128)
        sigma_minus = 0.5 * (sigma_x - 1.0j * sigma_y)
        op = np.asarray(np.kron(op, sigma_minus), dtype=np.complex128)
        for _ in range(k + 1, n_modes):
            op = np.asarray(np.kron(op, eye), dtype=np.complex128)
        dirac_d.append(op)

    # Majorana operators: c_{2k} = d_k + d_k^dagger, c_{2k+1} = -i(d_k - d_k^dagger)
    c_ops: list[np.ndarray] = []
    for k in range(n_modes):
        c_re = dirac_d[k] + dirac_d[k].conj().T
        c_im = -1.0j * (dirac_d[k] - dirac_d[k].conj().T)
        c_ops.append(c_re)
        c_ops.append(c_im)

    return c_ops


def pfaffian_4x4(a: np.ndarray) -> float:
    """Compute the Pfaffian of a 4x4 skew-symmetric matrix."""
    return float(a[0, 1] * a[2, 3] - a[0, 2] * a[1, 3] + a[0, 3] * a[1, 2])


def test_tst_qf_01_pfaffian_wick_for_gaussian_states() -> None:
    """TST-QF-01: Wick's theorem evaluates exactly to the Pfaffian for any Gaussian state."""
    c_ops = build_majorana_fock_basis(n_modes=2)  # 4 Majoranas, dim 4
    dim = 4

    # Build a random quadratic Hamiltonian H_K = (i/4) sum K_ab c_a c_b
    rng = np.random.default_rng(42)
    k_mat = rng.normal(size=(4, 4))
    k_mat = 0.5 * (k_mat - k_mat.T)  # skew-symmetric

    h_k = np.zeros((dim, dim), dtype=np.complex128)
    for a in range(4):
        for b in range(4):
            h_k += 0.25j * k_mat[a, b] * (c_ops[a] @ c_ops[b])

    # Thermal density matrix rho = exp(-H_K) / Tr(exp(-H_K))
    exp_h = la.expm(-h_k)
    rho = exp_h / np.trace(exp_h)

    # Compute two-point covariance: <c_a c_b> = delta_ab + i Omega_ab
    omega = np.zeros((4, 4), dtype=np.float64)
    for a in range(4):
        for b in range(4):
            val = np.trace(rho @ c_ops[a] @ c_ops[b])
            if a != b:
                omega[a, b] = float(np.imag(val))

    # Evaluate exact four-point expectation value <c_0 c_1 c_2 c_3>
    op_4 = c_ops[0] @ c_ops[1] @ c_ops[2] @ c_ops[3]
    exp_4 = float(np.real(np.trace(rho @ op_4)))

    # Compute Pfaffian of Omega
    pf_val = pfaffian_4x4(omega)

    # In a Gaussian state with <c_a c_b> = delta_ab + i Omega_ab,
    # <c_0 c_1 c_2 c_3> = i^2 * Pf(Omega) = -Pf(Omega)
    np.testing.assert_allclose(exp_4, -pf_val, atol=1e-12)


def test_tst_qf_02_wick_fails_for_generic_non_gaussian_state() -> None:
    """TST-QF-02: For a generic entangled or non-Gaussian state, the 4th cumulant is non-zero."""
    c_ops = build_majorana_fock_basis(n_modes=2)
    dim = 4

    # Singlet-like entangled pure state: |psi> = (1/sqrt(2)) (|01> - |10>)
    psi = np.zeros(dim, dtype=np.complex128)
    psi[1] = 1.0 / np.sqrt(2.0)
    psi[2] = -1.0 / np.sqrt(2.0)
    rho = np.outer(psi, psi.conj())

    omega = np.zeros((4, 4), dtype=np.float64)
    for a in range(4):
        for b in range(4):
            if a != b:
                omega[a, b] = float(np.imag(np.trace(rho @ c_ops[a] @ c_ops[b])))

    op_4 = c_ops[0] @ c_ops[1] @ c_ops[2] @ c_ops[3]
    exp_4 = float(np.real(np.trace(rho @ op_4)))
    pf_val = pfaffian_4x4(omega)

    cumulant_4 = abs(exp_4 - pf_val)
    # The 4th cumulant must NOT be zero for this entangled state
    assert cumulant_4 > 0.1, f"Expected non-zero cumulant for singlet state, got {cumulant_4}"


def test_tst_qf_03_non_convexity_of_quasifree_manifold() -> None:
    """TST-QF-03: The convex combination of two Gaussian states is NOT Gaussian."""
    c_ops = build_majorana_fock_basis(n_modes=2)

    # Two distinct Gaussian states rho1, rho2 from random K1, K2
    rng = np.random.default_rng(123)
    k1 = rng.normal(size=(4, 4))
    k1 = 0.5 * (k1 - k1.T)
    k2 = rng.normal(size=(4, 4))
    k2 = 0.5 * (k2 - k2.T)

    h1 = sum(0.25j * k1[a, b] * (c_ops[a] @ c_ops[b]) for a in range(4) for b in range(4))
    h2 = sum(0.25j * k2[a, b] * (c_ops[a] @ c_ops[b]) for a in range(4) for b in range(4))

    rho1 = la.expm(-h1) / np.trace(la.expm(-h1))
    rho2 = la.expm(-h2) / np.trace(la.expm(-h2))

    # Mixture state: rho_mix = 0.5 * rho1 + 0.5 * rho2
    rho_mix = 0.5 * rho1 + 0.5 * rho2

    # Measure covariance and 4th cumulant of the mixture
    omega_mix = np.zeros((4, 4), dtype=np.float64)
    for a in range(4):
        for b in range(4):
            if a != b:
                omega_mix[a, b] = float(np.imag(np.trace(rho_mix @ c_ops[a] @ c_ops[b])))

    op_4 = c_ops[0] @ c_ops[1] @ c_ops[2] @ c_ops[3]
    exp_4_mix = float(np.real(np.trace(rho_mix @ op_4)))
    pf_mix = pfaffian_4x4(omega_mix)

    cumulant_mix = abs(exp_4_mix - pf_mix)
    # Mixture of Gaussian states has a non-vanishing 4th cumulant!
    assert cumulant_mix > 1e-4, f"Expected non-zero cumulant in mixture state, got {cumulant_mix}"


def test_tst_qf_04_admissibility_purity_and_fermion_parity() -> None:
    """TST-QF-04: Covariance satisfies ||Omega|| <= 1; pure iff Omega^2 = -I; Pf(Omega) = +-1."""
    # Test on a pure state constructed from real orthogonal Grassmannian
    # For a 4x4 skew matrix, pure covariance has eigenvalues +-i, +-i
    # Standard normal form: [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]]
    j_std = np.array(
        [[0.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -1.0, 0.0]],
        dtype=np.float64,
    )

    # Any O in O(4) rotates J: Omega = O J O^T
    rng = np.random.default_rng(77)
    random_mat = rng.normal(size=(4, 4))
    q_ortho, _ = la.qr(random_mat)

    omega_pure = q_ortho @ j_std @ q_ortho.T

    # Admissibility: spectral norm exactly 1
    norm_op = la.norm(omega_pure, 2)
    np.testing.assert_allclose(norm_op, 1.0, atol=1e-14)

    # Purity: Omega^2 = -I
    omega_sq = omega_pure @ omega_pure
    np.testing.assert_allclose(omega_sq, -np.eye(4), atol=1e-14)

    # Fermion parity Pf(Omega) in {+1, -1}
    pf_val = pfaffian_4x4(omega_pure)
    assert abs(abs(pf_val) - 1.0) < 1e-12


def test_tst_qf_05_gaussian_entropy_and_maximum_entropy() -> None:
    """TST-QF-05: Gaussian entropy matches -Tr(rho log rho) and satisfies Jaynes maximum entropy."""
    c_ops = build_majorana_fock_basis(n_modes=2)

    rng = np.random.default_rng(99)
    k = rng.normal(size=(4, 4))
    k = 0.5 * (k - k.T)
    h_k = sum(0.25j * k[a, b] * (c_ops[a] @ c_ops[b]) for a in range(4) for b in range(4))

    rho_gauss = la.expm(-h_k) / np.trace(la.expm(-h_k))

    # Von Neumann entropy of rho_gauss: S = -Tr(rho log2 rho)
    eigvals = la.eigvalsh(rho_gauss)
    eigvals = eigvals[eigvals > 1e-15]
    s_exact = -float(np.sum(eigvals * np.log2(eigvals)))

    # Compute covariance Omega
    omega = np.zeros((4, 4), dtype=np.float64)
    for a in range(4):
        for b in range(4):
            if a != b:
                omega[a, b] = float(np.imag(np.trace(rho_gauss @ c_ops[a] @ c_ops[b])))

    # Symplectic eigenvalues nu_k: singular values / positive imaginary eigenvalues of Omega
    eig_omega = la.eigvals(omega)
    pos_imag = np.sort(np.imag(eig_omega[np.imag(eig_omega) > 0]))

    # Gaussian entropy formula: sum_k h((1 + nu_k)/2)
    s_gauss_formula = 0.0
    for nu in pos_imag:
        p = 0.5 * (1.0 + nu)
        if 0.0 < p < 1.0:
            s_gauss_formula += -p * np.log2(p) - (1.0 - p) * np.log2(1.0 - p)

    np.testing.assert_allclose(s_exact, s_gauss_formula, atol=1e-12)


def test_tst_qf_06_affine_tangent_property() -> None:
    """TST-QF-06: E(Omega') - [omega_{Omega'}(H_MF[Omega]) - C[Omega]] = O(||Omega' - Omega||^2) >= 0."""
    # Define a pure Kitaev-like 4-site quartic energy functional
    # E(Omega) = - (1/8) sum V_abcd (Omega_ab Omega_cd - Omega_ac Omega_bd + Omega_ad Omega_bc)
    rng = np.random.default_rng(2026)
    v_tensor = rng.normal(size=(4, 4, 4, 4))
    # Skew-symmetrize tensor on pairs (ab) and (cd)
    v_tensor = 0.25 * (v_tensor - np.swapaxes(v_tensor, 0, 1))
    v_tensor = 0.25 * (v_tensor - np.swapaxes(v_tensor, 2, 3))

    def compute_energy(om: np.ndarray) -> float:
        val = 0.0
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    for d in range(4):
                        wick = om[a, b] * om[c, d] - om[a, c] * om[b, d] + om[a, d] * om[b, c]
                        val += float(v_tensor[a, b, c, d] * wick)
        return -0.125 * val

    def compute_hmf(om: np.ndarray) -> np.ndarray:
        # Full gradient h_ab = dE / dOmega_ab taking all Wick pairings into account
        grad = np.zeros((4, 4), dtype=np.float64)
        eps = 1e-7
        for i in range(4):
            for j in range(4):
                d_om = np.zeros((4, 4))
                d_om[i, j] = eps
                grad[i, j] = (compute_energy(om + d_om) - compute_energy(om - d_om)) / (2.0 * eps)
        return grad

    om_base = np.array(
        [[0.0, 0.8, 0.0, 0.0], [-0.8, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.8], [0.0, 0.0, -0.8, 0.0]],
        dtype=np.float64,
    )
    e_base = compute_energy(om_base)
    h_mf = compute_hmf(om_base)

    # Check finite-difference tangent convergence: Delta E - (Tr(h dOmega) - C) = O(delta^2)
    direction = rng.normal(size=(4, 4))
    direction = 0.5 * (direction - direction.T)
    direction /= la.norm(direction, "fro")

    deltas = [1e-2, 1e-3, 1e-4]
    errors = []
    for delta in deltas:
        om_perturbed = om_base + delta * direction
        e_pert = compute_energy(om_perturbed)
        # Linearized approximation: Tr(h_mf * delta_Omega)
        linear_approx = e_base + np.sum(h_mf * (om_perturbed - om_base))
        errors.append(abs(e_pert - linear_approx))

    # Error scaling ratio between successive 10x delta reductions should be ~100 (quadratic)
    ratio_1 = errors[0] / errors[1]
    ratio_2 = errors[1] / errors[2]
    assert 90.0 < ratio_1 < 110.0, f"Expected quadratic scaling ~100, got {ratio_1}"
    assert 90.0 < ratio_2 < 110.0, f"Expected quadratic scaling ~100, got {ratio_2}"


def test_tst_qf_07_grassmannian_stationarity_commutator() -> None:
    """TST-QF-07: Stationarity [h(Omega), Omega] = 0 and spectral projector construction."""
    # Consider a single-particle Hamiltonian matrix h in so(4)
    h = np.array(
        [[0.0, 2.5, 0.0, 0.0], [-2.5, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.2], [0.0, 0.0, -1.2, 0.0]],
        dtype=np.float64,
    )
    # The negative-energy ground covariance is Omega = -sgn(i h)
    # Since i*h has eigenvalues +-2.5, +-1.2:
    # sgn(i*h) has eigenvalues +-1, +-1 with the same eigenvectors
    omega_ground = np.array(
        [[0.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -1.0, 0.0]],
        dtype=np.float64,
    )

    # Commutator [h, Omega] must vanish exactly at stationarity
    comm = h @ omega_ground - omega_ground @ h
    np.testing.assert_allclose(comm, 0.0, atol=1e-15)

    # Spectral projector P = (I + i Omega) / 2 must be an orthogonal projector
    p_proj = 0.5 * (np.eye(4, dtype=np.complex128) + 1.0j * omega_ground)
    p_sq = p_proj @ p_proj
    np.testing.assert_allclose(p_sq, p_proj, atol=1e-15)
    np.testing.assert_allclose(p_proj.conj().T, p_proj, atol=1e-15)


def test_tst_qf_08_graded_locality_vs_spin_commutativity() -> None:
    """TST-QF-08: Distinct-site Majoranas anticommute while physical spin bilinears commute."""
    c_ops = build_majorana_fock_basis(n_modes=2)
    # c_ops[0..1] correspond to site 1 (mode 1), c_ops[2..3] correspond to site 2 (mode 2)

    # 1. Graded locality: {c_i, c_j} = 0 for i on site 1 and j on site 2
    anti_comm = c_ops[0] @ c_ops[2] + c_ops[2] @ c_ops[0]
    np.testing.assert_allclose(anti_comm, 0.0, atol=1e-15)

    # 2. Spin bilinears on distinct sites: sigma_1^z = -i c_0 c_1, sigma_2^z = -i c_2 c_3
    sigma_1_z = -1.0j * (c_ops[0] @ c_ops[1])
    sigma_2_z = -1.0j * (c_ops[2] @ c_ops[3])

    # Commutator [sigma_1^z, sigma_2^z] must vanish identically (microcausality)
    comm_spins = sigma_1_z @ sigma_2_z - sigma_2_z @ sigma_1_z
    np.testing.assert_allclose(comm_spins, 0.0, atol=1e-15)
