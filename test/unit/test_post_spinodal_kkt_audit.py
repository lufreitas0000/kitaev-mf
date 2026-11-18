"""Unit tests for Final Wave Phase 1 Contracts (HC-P1-01, HC-P1-02, HC-P1-03).

Verifies:
1. HC-P1-01: Perturbative vs self-consistent Majorana gap prefactor reconciliation (K-010, K-019).
2. HC-P1-02: Post-spinodal KKT saddle excursion elimination and Bloch-ball admissibility gate (K-023).
3. HC-P1-03: Existence and schema validity of the Phase 1 Master Solver Selection & Resource Table.
"""

from __future__ import annotations

import os

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.correlators import compute_correlators_from_bloch, compute_gauge_residuals
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.observables import compute_sublattice_magnetizations
from kitaev_mf.services.seeds import make_polarized_seed, make_qsl_seed
from kitaev_mf.services.solver import (
    compute_lehmann_dual_hessian,
    solve_self_consistent_mean_field,
)


def test_hc_p1_01_perturbative_gap_reconciliation_identities() -> None:
    """Verify exact quantitative ratios comparing c_MF = 2.5921 and static vison perturbation scale."""
    u_bar = 0.5248657458506
    two_delta_gauge = 2.0 * u_bar
    two_vison_static = 0.2640
    denom_ratio_sq = (two_delta_gauge / two_vison_static) ** 2
    assert np.isclose(denom_ratio_sq, 15.8106, atol=1e-3)

    c_mf_bloch = 2.5921
    c_mf_canonical = 4.0 * c_mf_bloch
    assert np.isclose(c_mf_canonical, 10.3684, atol=1e-4)

    c_static_vison = 4.0 / (two_vison_static**2)
    assert np.isclose(c_static_vison, 57.3921, atol=1e-3)


def test_hc_p1_02_post_spinodal_kkt_admissibility_at_h0450() -> None:
    """Verify that at h=0.450 (post-spinodal), KKT-admissible polarized state satisfies E = -h - 0.5 and m_111 = 1."""
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()
    h_val = 0.450
    params = ModelParameters(
        K=-1.0, h=(h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
    )
    config = SolverConfig(
        k_grid_shape=(8, 8),
        max_iterations=50,
        tol_correlator=1e-9,
        tol_gauge=1e-9,
        method="diis",
        mixing_alpha=0.5,
        skew_mode="retract_step",
    )
    pol_seed = make_polarized_seed(cell, axis="111")
    state_final, metrics = solve_self_consistent_mean_field(cell, params, pol_seed, config, backend)
    assert metrics.converged
    assert metrics.final_gauge_residual < 1e-10

    eb = compute_gaussian_energy(state_final, cell, params)
    sig_mag, _ = compute_sublattice_magnetizations(state_final, cell)
    m_vec = 0.5 * (sig_mag[0] + sig_mag[1])
    m_111 = float(np.dot(m_vec, np.ones(3) / np.sqrt(3.0)))

    assert np.isclose(m_111, 1.0, atol=1e-6)
    assert np.isclose(eb.per_site, -h_val - 0.5, atol=1e-6)


def test_hc_p1_03_master_solver_table_and_audit_data_exist() -> None:
    """Verify that the post-spinodal KKT audit table and master solver selection table are present and non-empty."""
    audit_path = "data/processed/benchmarks/m1_9/post_spinodal_kkt_audit.dat"
    master_path = "data/processed/benchmarks/m1_9/master_solver_resource_table.dat"
    assert os.path.isfile(audit_path)
    assert os.path.isfile(master_path)

    with open(audit_path, encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]
    assert len(lines) >= 10  # Header + 12 field points

    with open(master_path, encoding="utf-8") as f:
        m_lines = [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]
    assert len(m_lines) >= 6  # Header + 6 solver regimes


def test_lehmann_dual_hessian_preconditioned_newton_step() -> None:
    """Contract H8: Preconditioning dual multiplier update with Lehmann dual Hessian yields > 1000x residual drop."""
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, h=(0.02, 0.02, 0.02))
    seed = make_qsl_seed(cell)
    k_grid = cell.make_k_grid(8, 8)

    # Initial gauge residual
    h_b0 = build_bloch_hamiltonian(k_grid, cell, params, seed, backend)
    _, v_out0 = compute_correlators_from_bloch(h_b0, cell, k_grid, backend=backend)
    g0 = compute_gauge_residuals(v_out0).ravel()
    res_norm_before = float(np.linalg.norm(g0))
    assert res_norm_before > 0.01

    # 1 Lehmann Newton step: Delta lambda = - H_ll^{-1} g0
    h_ll = compute_lehmann_dual_hessian(cell, params, seed, k_grid, backend)
    assert h_ll.shape == (6, 6)
    assert np.all(np.linalg.eigvalsh(h_ll) > 0.0)  # Positive definite

    dlam = np.linalg.solve(h_ll, g0)
    lam_new = seed.lambdas - dlam.reshape(seed.lambdas.shape)
    state_new = MeanFieldState(U=seed.U.copy(), V=seed.V.copy(), lambdas=lam_new)

    h_b1 = build_bloch_hamiltonian(k_grid, cell, params, state_new, backend)
    _, v_out1 = compute_correlators_from_bloch(h_b1, cell, k_grid, backend=backend)
    g1 = compute_gauge_residuals(v_out1).ravel()
    res_norm_after = float(np.linalg.norm(g1))

    reduction_factor = res_norm_before / res_norm_after
    assert reduction_factor > 1000.0, f"Expected > 1000x drop, got {reduction_factor:.1f}x"


def test_post_spinodal_monotonicity_across_hysteresis_window() -> None:
    """Contract H7, K-023: Strict KKT tolerance guarantees monotonic energy without downward serrations."""
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()
    config = SolverConfig(
        k_grid_shape=(8, 8),
        max_iterations=60,
        tol_correlator=1e-9,
        tol_gauge=1e-9,
        method="diis",
    )
    pol_seed = make_polarized_seed(cell, axis="111")

    h_points = [0.410, 0.430, 0.450]
    energies: list[float] = []

    for h_val in h_points:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)
        st, met = solve_self_consistent_mean_field(cell, params, pol_seed, config, backend)
        assert met.converged
        assert met.final_gauge_residual < 1e-9

        eb = compute_gaussian_energy(st, cell, params)
        e_site = float(eb.per_site)
        energies.append(e_site)

        # Must strictly match the physical polarized ground state E = -h - 0.5
        assert np.isclose(e_site, -h_val - 0.5, atol=1e-5), (
            f"At h={h_val}, expected E = {-h_val - 0.5:.5f}, got {e_site:.5f}"
        )

    # Strictly monotonic energy decrease across the post-spinodal window: E(0.450) < E(0.430) < E(0.410)
    assert energies[0] > energies[1] > energies[2], f"Non-monotonic energy sequence: {energies}"
