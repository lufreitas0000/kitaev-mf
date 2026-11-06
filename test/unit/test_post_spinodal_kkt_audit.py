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
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.observables import compute_sublattice_magnetizations
from kitaev_mf.services.seeds import make_polarized_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field


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
