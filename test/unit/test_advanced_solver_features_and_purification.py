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
    """Verify adiabatic warm-start commutator gradient delta Y = -tau (A + Omega_0 A Omega_0)."""
    cluster = RealSpaceTorusCluster(L=4)
    cell2 = make_honeycomb_2site(cluster.d_nn)
    state = MeanFieldState.from_pure_kitaev(cell2, u_gauge=1.0, u_matter=-0.5)

    params1 = ModelParameters(K=-1.0, h=(0.04, 0.04, 0.04))
    A1 = cluster.assemble_skew_hamiltonian(params1, state)

    # Baseline cold start
    omega_1, _ = purify_newton_schulz(A1, tol=1e-10)

    # Warm-started purification executes the commutator gradient branch
    omega_1_warm, metrics_warm = purify_newton_schulz(A1, warm_start=omega_1, tol=1e-10)
    assert np.allclose(omega_1_warm, -omega_1_warm.T, atol=1e-12)
    assert np.allclose(omega_1_warm @ omega_1_warm, -np.eye(cluster.num_majoranas), atol=1e-8)
    assert metrics_warm.final_idempotency_error < 1e-8

    # Test warm start with Padé rational solver as well
    omega_p_warm, m_p_warm = purify_pade_sign(A1, warm_start=omega_1, order="32", tol=1e-10)
    assert np.allclose(omega_p_warm, -omega_p_warm.T, atol=1e-12)
    assert np.allclose(omega_p_warm @ omega_p_warm, -np.eye(cluster.num_majoranas), atol=1e-8)
    assert m_p_warm.final_idempotency_error < 1e-8


def test_dual_lehmann_preconditioning_in_sc_mf() -> None:
    """Verify execution of dual Lehmann preconditioned Newton step in solve_self_consistent_mean_field."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, h=(0.04, 0.04, 0.04))
    config = SolverConfig(
        max_iterations=15,
        tol_correlator=1e-6,
        tol_gauge=1e-6,
        k_grid_shape=(4, 4),
        method="fp",
        precondition_dual=True,
    )
    init_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5)
    state, metrics = solve_self_consistent_mean_field(cell, params, init_state, config)

    assert state.U.shape == (3, 4, 4)
    assert metrics.iterations > 0


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
