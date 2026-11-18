"""Unit tests for advanced numerical solver paths, purification accelerators, and domain helpers (K-018, K-021, K-023).

Verifies:
1. Halley Padé [2/1] rational sign purification (order="21").
2. Adiabatic warm-start commutator acceleration in Newton-Schulz and Padé purification.
3. SC-MF fixed-point solver with dual Lehmann preconditioning (SolverConfig.precondition_dual=True).
4. Damped Picard accelerator in real-space cluster purification.
5. InhomogeneousMeanFieldState spatial coarse-graining via to_homogeneous_mean().
6. Skew-Pfaffian singular rank deficiency and odd-dimension branches returning 0.0.
7. Single-sublattice wave function support in compute_vortex_mode_core_weight.
8. Lattice Bond property helpers and fractional_to_cartesian_k conversion.
"""

import numpy as np
import pytest

import kitaev_mf.services.solver as solver_mod
from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import InhomogeneousMeanFieldState, MeanFieldState
from kitaev_mf.services.purification import (
    purify_newton_schulz,
    purify_pade_sign,
)
from kitaev_mf.services.real_space import (
    RealSpaceTorusCluster,
    compute_skew_pfaffian,
)
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def test_halley_pade_21_rational_purification() -> None:
    """Verify Padé [2/1] Halley rational iteration (order='21') for matrix sign purification."""
    # Build 8x8 random skew-symmetric matrix
    rng = np.random.default_rng(42)
    B = rng.standard_normal((8, 8))
    A = B - B.T

    omega_21, metrics_21 = purify_pade_sign(A, order="21", tol=1e-10)
    omega_32, _ = purify_pade_sign(A, order="32", tol=1e-10)

    # Both Padé approximations must converge to the exact same spectral sign
    assert np.allclose(omega_21, omega_32, atol=1e-8)
    # Exact skew-symmetry and projector idempotency Omega^2 = -I
    assert np.allclose(omega_21, -omega_21.T, atol=1e-12)
    assert np.allclose(omega_21 @ omega_21, -np.eye(8), atol=1e-9)
    assert metrics_21.final_idempotency_error < 1e-8


def test_adiabatic_warm_start_commutator_acceleration() -> None:
    """Verify adiabatic warm-start commutator gradient moves closer to A2 and works in outer SCF."""
    cluster = RealSpaceTorusCluster(L=4)
    cell2 = make_honeycomb_2site(cluster.d_nn)
    state = MeanFieldState.from_pure_kitaev(cell2, u_gauge=1.0, u_matter=-0.5)

    params1 = ModelParameters(K=-1.0, h=(0.04, 0.04, 0.04))
    params2 = ModelParameters(K=-1.0, h=(0.06, 0.04, 0.03))
    A1 = cluster.assemble_skew_hamiltonian(params1, state)
    A2 = cluster.assemble_skew_hamiltonian(params2, state)

    omega_1, _ = purify_newton_schulz(A1, tol=1e-10)
    # Verify non-trivial commutator gradient A2 + omega_1 A2 omega_1 != 0
    grad_comm = A2 + omega_1 @ A2 @ omega_1
    assert float(np.linalg.norm(grad_comm)) > 0.05

    omega_ref, _ = cluster.solve_covariance_complex_eigh(A2)
    omega_warm, _ = purify_newton_schulz(A2, warm_start=omega_1, tol=1e-10)
    omega_p_warm, m_p_warm = purify_pade_sign(A2, warm_start=omega_1, order="32", tol=1e-10)

    # Single inner warm-start yields exact skew idempotent projector closer to A2's ground state
    # and with lower Rayleigh-Ritz energy on A2 than the un-rotated omega_1.
    for om in (omega_warm, omega_p_warm):
        assert np.allclose(om, -om.T, atol=1e-12)
        assert np.allclose(om @ om, -np.eye(cluster.num_majoranas), atol=1e-8)
        assert np.linalg.norm(om - omega_ref) < np.linalg.norm(omega_1 - omega_ref)
        e_warm = 0.25 * float(np.sum(A2 * om))
        e_prev = 0.25 * float(np.sum(A2 * omega_1))
        assert e_warm < e_prev
    assert m_p_warm.final_idempotency_error < 1e-8

    # A second warm-started tangent step from omega_warm contracts error toward omega_ref further.
    omega_warm2, _ = purify_newton_schulz(A2, warm_start=omega_warm, tol=1e-10)
    assert np.linalg.norm(omega_warm2 - omega_ref) < np.linalg.norm(omega_warm - omega_ref)

    # Verify SCF on A1 warm-started vs cold-start converges cleanly under outer DIIS (warm_start=False).
    _, scf_state, met_scf = cluster.solve_real_space_purification(
        params=params2,
        initial_state=state,
        max_iter=25,
        tol=1e-6,
        accelerator="diis",
        warm_start=False,
    )
    assert met_scf["converged"] is True
    assert met_scf["final_residual"] < 1e-6
    assert np.allclose(scf_state.V, -np.swapaxes(scf_state.V, -1, -2), atol=1e-12)


def test_dual_lehmann_preconditioning_in_sc_mf() -> None:
    """Verify converged accuracy and residual drop under precondition_dual=True."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, h=(0.04, 0.04, 0.04))
    config = SolverConfig(
        max_iterations=30,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        k_grid_shape=(4, 4),
        method="diis",
        precondition_dual=True,
    )
    init_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5)
    state, metrics = solve_self_consistent_mean_field(cell, params, init_state, config)

    assert state.U.shape == (3, 4, 4)
    assert metrics.converged is True
    assert metrics.final_gauge_residual < 1e-8
    assert metrics.final_correlator_residual < 1e-8
    assert metrics.iterations > 0


def test_dual_lehmann_indefinite_hessian_fallback(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify indefinite finite-difference Hessian safely falls back to damped gradient ascent."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, h=(0.04, 0.04, 0.04))
    config = SolverConfig(
        max_iterations=5,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        k_grid_shape=(4, 4),
        method="fp",
        lambda_step=0.2,
        precondition_dual=True,
    )
    init_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5)

    def _indefinite_hessian(*args: object, **kwargs: object) -> np.ndarray:
        return -np.eye(6, dtype=np.float64)

    monkeypatch.setattr(solver_mod, "compute_lehmann_dual_hessian", _indefinite_hessian)
    st_fb, metrics_fb = solve_self_consistent_mean_field(cell, params, init_state, config)

    config_plain = SolverConfig(
        max_iterations=5,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        k_grid_shape=(4, 4),
        method="fp",
        lambda_step=0.2,
        precondition_dual=False,
    )
    st_plain, metrics_plain = solve_self_consistent_mean_field(
        cell, params, init_state, config_plain
    )
    assert metrics_fb.iterations == 5
    assert np.allclose(metrics_fb.residual_history, metrics_plain.residual_history, atol=1e-14)
    assert np.allclose(st_fb.lambdas, st_plain.lambdas, atol=1e-14)


def test_purification_input_validation() -> None:
    """Verify purification and solver loops reject invalid max_iter, order, and accelerator."""
    A = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=np.float64)
    with pytest.raises(ValueError, match="max_iter must be >= 1"):
        purify_newton_schulz(A, max_iter=0)
    with pytest.raises(ValueError, match="max_iter must be >= 1"):
        purify_pade_sign(A, max_iter=0)
    with pytest.raises(ValueError, match="order must be '32' or '21'"):
        purify_pade_sign(A, order="42")  # type: ignore[arg-type]

    cluster = RealSpaceTorusCluster(L=2)
    cell2 = make_honeycomb_2site(cluster.d_nn)
    init_state = MeanFieldState.from_pure_kitaev(cell2)
    params = ModelParameters(K=-1.0)
    with pytest.raises(ValueError, match="max_iter must be >= 1"):
        cluster.solve_real_space_purification(params, init_state, max_iter=0)
    with pytest.raises(ValueError, match="Unknown accelerator"):
        cluster.solve_real_space_purification(params, init_state, accelerator="cg")  # type: ignore[arg-type]


def test_damped_picard_real_space_purification() -> None:
    """Verify damped Picard accelerator branch in solve_real_space_purification."""
    cluster = RealSpaceTorusCluster(L=4)
    cell2 = make_honeycomb_2site(cluster.d_nn)
    params = ModelParameters(K=-1.0, h=(0.02, 0.02, 0.02))
    init_state = MeanFieldState.from_pure_kitaev(cell2, u_gauge=1.0, u_matter=-0.5)

    omega, state, metrics = cluster.solve_real_space_purification(
        params=params,
        initial_state=init_state,
        max_iter=5,
        tol=1e-4,
        accelerator="picard",
        mixing_alpha=0.3,
    )
    assert omega.shape == (cluster.num_majoranas, cluster.num_majoranas)
    assert metrics["iterations"] > 0


def test_inhomogeneous_state_to_homogeneous_mean() -> None:
    """Verify InhomogeneousMeanFieldState.to_homogeneous_mean() spatial averaging."""
    cluster = RealSpaceTorusCluster(L=4)
    cell2 = make_honeycomb_2site(cluster.d_nn)
    homo = MeanFieldState.from_pure_kitaev(cell2, u_gauge=1.0, u_matter=-0.5)
    inhomo = InhomogeneousMeanFieldState.from_homogeneous(homo, num_cells=cluster.num_cells)

    averaged = inhomo.to_homogeneous_mean()
    assert np.allclose(averaged.U, homo.U, atol=1e-12)
    assert np.allclose(averaged.V, homo.V, atol=1e-12)
    assert np.allclose(averaged.lambdas, homo.lambdas, atol=1e-12)


def test_skew_pfaffian_singular_and_odd_branches() -> None:
    """Verify compute_skew_pfaffian returns 0.0 for odd dimensions or singular skew matrices."""
    # Odd dimension
    A_odd = np.zeros((5, 5))
    assert compute_skew_pfaffian(A_odd) == 0.0

    # Singular even-dimensional skew matrix (zero row/col)
    A_sing = np.zeros((6, 6))
    A_sing[0, 1] = 1.0
    A_sing[1, 0] = -1.0
    # Remaining 4 rows/cols are identically zero -> rank deficient
    assert compute_skew_pfaffian(A_sing) == 0.0


def test_compute_skew_pfaffian_pivoting_and_oracle() -> None:
    """Verify complete-pivoting swaps, negative Pfaffian sign, and Pf^2 == det oracle."""
    # Explicit 4x4 skew matrix with zero at (0, 1) to force p != k, q != k+1 swaps.
    A4 = np.array(
        [
            [0.0, 0.0, 1.3, -2.1],
            [0.0, 0.0, 0.7, 3.5],
            [-1.3, -0.7, 0.0, 0.4],
            [2.1, -3.5, -0.4, 0.0],
        ],
        dtype=np.float64,
    )
    pf4_exact = A4[0, 1] * A4[2, 3] - A4[0, 2] * A4[1, 3] + A4[0, 3] * A4[1, 2]
    assert np.isclose(compute_skew_pfaffian(A4), pf4_exact, atol=1e-13)

    # Row-col transposition flips Pfaffian sign cleanly
    perm = [1, 0, 2, 3]
    A4_swapped = A4[np.ix_(perm, perm)]
    assert np.isclose(compute_skew_pfaffian(A4_swapped), -pf4_exact, atol=1e-13)

    # Random dense 12x12 skew-symmetric matrices: verify Pf(A)^2 == det(A) and sign flip under swap
    rng = np.random.default_rng(2026)
    for _ in range(5):
        B = rng.standard_normal((12, 12))
        M = B - B.T
        pf = compute_skew_pfaffian(M)
        det = float(np.linalg.det(M))
        assert np.isclose(pf * pf, det, rtol=1e-10, atol=1e-10)
        swap = [1, 0] + list(range(2, 12))
        assert np.isclose(compute_skew_pfaffian(M[np.ix_(swap, swap)]), -pf, rtol=1e-10)


def test_vortex_mode_core_weight_single_sublattice_support() -> None:
    """Verify compute_vortex_mode_core_weight with len(psi) == num_total_sites (matter sector)."""
    cluster = RealSpaceTorusCluster(L=4)
    psi_matter = np.zeros(cluster.num_sites)
    # Localize wave function at site 0
    psi_matter[0] = 1.0

    w_core, w_uniform, ratio = cluster.compute_vortex_mode_core_weight(
        psi=psi_matter,
        p1=(1, 1),
        p2=(2, 2),
    )
    assert 0.0 <= w_core <= 1.0
    assert w_uniform > 0.0
    assert ratio >= 0.0


def test_lattice_bond_helpers_and_k_conversion() -> None:
    """Verify Bond.site_i, Bond.site_j property helpers and fractional_to_cartesian_k."""
    cell = make_honeycomb_2site()
    bond0 = cell.bonds[0]
    assert bond0.site_i == bond0.sublattice_i
    assert bond0.site_j == bond0.sublattice_j

    # Reciprocal momentum conversion
    k_cart = cell.fractional_to_cartesian_k(np.array([0.0, 0.0]))
    assert np.allclose(k_cart, np.array([0.0, 0.0]))
