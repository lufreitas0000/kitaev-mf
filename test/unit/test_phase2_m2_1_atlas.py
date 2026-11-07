"""Test suite for Phase 2 Milestone M2.1 Task T2.1.1: Real-Space Observable Atlas."""

from pathlib import Path

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


def test_t2_1_1_zero_field_qsl_invariants() -> None:
    """Verify zero-field pure Kitaev invariants on real torus (L=6, N=72 sites)."""
    L = 6
    cluster = RealSpaceTorusCluster(L=L)
    cell2 = make_honeycomb_2site()
    seed_state = make_qsl_seed(cell2)
    params = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))

    omega, final_state, met = cluster.solve_real_space_purification(
        params=params,
        initial_state=seed_state,
        max_iter=15,
        tol=1e-8,
        accelerator="lbroyden",
        mixing_alpha=0.5,
        lbroyden_history_size=10,
    )

    assert met["converged"]
    assert met["iterations"] <= 5

    # 1. Zero magnetization at h=0
    m_111 = float(np.sum(final_state.V[:, 0, 1:]) / (2 * np.sqrt(3.0)))
    assert abs(m_111) < 1e-12

    # 2. Vanishing Lagrange multipliers at h=0
    assert np.linalg.norm(final_state.lambdas) < 1e-12

    # 3. Energy per site within finite-size correction of thermodynamic U_bar = 0.524865
    breakdown = compute_gaussian_energy(final_state, cell2, params)
    # E_site = -1.5 * U_bar -> ~ -0.787
    assert -0.80 < breakdown.per_site < -0.77

    # 4. Plaquette flux W_p == +1.0 in the vortex-free ground state
    u_z = final_state.U[0, 3, 3]
    u_x = final_state.U[1, 1, 1]
    u_y = final_state.U[2, 2, 2]
    w_p = (u_x * u_y * u_z) ** 2
    assert np.isclose(w_p, 1.0, atol=1e-10)
    assert np.isclose(u_x, 1.0, atol=1e-10)
    assert np.isclose(u_y, 1.0, atol=1e-10)
    assert np.isclose(u_z, 1.0, atol=1e-10)


def test_t2_1_1_symmetry_and_c3_invariance() -> None:
    """Verify exact C3 symmetry (m_x == m_y == m_z) and sublattice equality in field h || [111]."""
    L = 6
    cluster = RealSpaceTorusCluster(L=L)
    cell2 = make_honeycomb_2site()
    seed_state = make_qsl_seed(cell2)
    h_mag = 0.15
    h_vec = (h_mag / np.sqrt(3.0), h_mag / np.sqrt(3.0), h_mag / np.sqrt(3.0))
    params = ModelParameters(K=-1.0, h=h_vec)

    omega, final_state, met = cluster.solve_real_space_purification(
        params=params,
        initial_state=seed_state,
        max_iter=30,
        tol=1e-6,
        accelerator="lbroyden",
        mixing_alpha=0.5,
        lbroyden_history_size=10,
    )

    assert met["converged"]

    # Extract sublattice magnetizations
    m_A = np.array([final_state.V[0, 0, 1], final_state.V[0, 0, 2], final_state.V[0, 0, 3]])
    m_B = np.array([final_state.V[1, 0, 1], final_state.V[1, 0, 2], final_state.V[1, 0, 3]])

    # Sublattice equality
    assert np.allclose(m_A, m_B, atol=1e-8)

    # C3 rotational symmetry on spin (m_x == m_y == m_z)
    assert np.isclose(m_A[0], m_A[1], atol=1e-8)
    assert np.isclose(m_A[1], m_A[2], atol=1e-8)

    # C3 rotational symmetry on bond correlators (U[0, 3, 3] == U[1, 1, 1] == U[2, 2, 2])
    u_z = final_state.U[0, 3, 3]
    u_x = final_state.U[1, 1, 1]
    u_y = final_state.U[2, 2, 2]
    assert np.isclose(u_x, u_y, atol=1e-8)
    assert np.isclose(u_y, u_z, atol=1e-8)

    # C3 on matter bond correlators (U[0, 0, 0] == U[1, 0, 0] == U[2, 0, 0])
    u0_z = final_state.U[0, 0, 0]
    u0_x = final_state.U[1, 0, 0]
    u0_y = final_state.U[2, 0, 0]
    assert np.isclose(u0_x, u0_y, atol=1e-8)
    assert np.isclose(u0_y, u0_z, atol=1e-8)

    # Positive magnetization along [111]
    m_111 = float(np.sum(m_A) / np.sqrt(3.0))
    assert 0.0 < m_111 < 1.0


def test_t2_1_1_time_reversal_odd_parity() -> None:
    """Verify time-reversal odd parity m(-h) == -m(h)."""
    L = 6
    cluster = RealSpaceTorusCluster(L=L)
    cell2 = make_honeycomb_2site()
    seed_state = make_qsl_seed(cell2)

    h_mag = 0.10
    h_pos = (h_mag / np.sqrt(3.0), h_mag / np.sqrt(3.0), h_mag / np.sqrt(3.0))
    h_neg = (-h_mag / np.sqrt(3.0), -h_mag / np.sqrt(3.0), -h_mag / np.sqrt(3.0))

    _, state_pos, met_pos = cluster.solve_real_space_purification(
        params=ModelParameters(K=-1.0, h=h_pos),
        initial_state=seed_state,
        max_iter=30,
        tol=1e-6,
        accelerator="lbroyden",
        mixing_alpha=0.5,
        lbroyden_history_size=10,
    )
    _, state_neg, met_neg = cluster.solve_real_space_purification(
        params=ModelParameters(K=-1.0, h=h_neg),
        initial_state=seed_state,
        max_iter=30,
        tol=1e-6,
        accelerator="lbroyden",
        mixing_alpha=0.5,
        lbroyden_history_size=10,
    )

    assert met_pos["converged"] and met_neg["converged"]

    m_pos = state_pos.V[0, 0, 1:]
    m_neg = state_neg.V[0, 0, 1:]
    assert np.allclose(m_neg, -m_pos, atol=1e-8)


def test_t2_1_1_dataset_integrity() -> None:
    """Verify exported benchmark dataset formatting, convergence, and physical bounds."""
    dat_path = Path(
        "/home/lucas/Projects/kitaev-mf/data/processed/benchmarks/m2_1/real_space_observable_atlas.dat"
    )
    assert dat_path.exists()

    data = np.loadtxt(dat_path)
    assert data.shape[0] == 26
    assert data.shape[1] == 13

    # Columns: h, e_site, m_111, chi_111, m_x, m_y, m_z, u_00, u_gg, w_p, lambda_norm, conv, iters
    h = data[:, 0]
    e_site = data[:, 1]
    m_111 = data[:, 2]
    w_p = data[:, 9]
    conv = data[:, 11]

    # Every single point in reference benchmark must be converged
    assert np.all(conv == 1)

    # Monotonicity of field
    assert np.all(np.diff(h) > 0)

    # Zero-field plaquette flux is +1
    assert np.isclose(w_p[0], 1.0, atol=1e-10)

    # Bounds on magnetization
    assert np.all((m_111 >= 0.0) & (m_111 <= 1.0 + 1e-12))

    # Magnetization is monotonically non-decreasing
    assert np.all(np.diff(m_111) >= -1e-10)

    # Energy is monotonically non-increasing with field
    assert np.all(np.diff(e_site) <= 1e-10)
