"""Unit tests for Wave 2B (UOW-03 / HC-02): Gutzwiller D_j=+1 VMC stability across SU(2)_R & flux sectors."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.sparse import csr_matrix

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.qsl_bridge import SO4_L
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.gutzwiller_projection import (
    bca_confidence_interval,
    torus_bonds,
    torus_plaquettes,
    torus_site,
)
from kitaev_mf.services.gutzwiller_stability import (
    FluxSectorName,
    build_flux_sector_bonds,
    build_site_su2_gauge_orthogonal,
    compute_chain_diagnostics,
    evaluate_flux_sector_gutzwiller,
    evaluate_n_decoupling_su2_breaking,
    evaluate_su2_gauge_stability,
    make_staggered_su2_patterns,
    verify_thouless_su2_block_diagonality,
)
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field

REPO_ROOT = Path(__file__).resolve().parents[2]


def _solve_reference_states() -> tuple[MeanFieldState, MeanFieldState]:
    cell = make_honeycomb_2site()
    cfg = SolverConfig(k_grid_shape=(24, 24), tol_correlator=1e-10, tol_gauge=1e-10)
    st0, met0 = solve_self_consistent_mean_field(
        cell, ModelParameters(K=-1.0), make_qsl_seed(cell), cfg
    )
    assert met0.converged
    c25 = 0.25 / float(np.sqrt(3.0))
    st25, met25 = solve_self_consistent_mean_field(
        cell, ModelParameters(K=-1.0, h=(c25, c25, c25)), st0, cfg
    )
    assert met25.converged
    return st0, st25


def _exact_spin_ed_sector_ground_energy(nx: int, ny: int, target_wp: list[int]) -> float:
    """Independent exact 2^N spin-Hamiltonian sector ground energy per site at h=0."""
    n_sites = 2 * nx * ny
    dim = 1 << n_sites
    idx = np.arange(dim, dtype=np.int32)
    spins_x = [np.where((idx >> s) & 1, 1.0, -1.0) for s in range(n_sites)]
    phases_y = [np.where((idx >> s) & 1, 1.0j, -1.0j) for s in range(n_sites)]

    def _op_csr(factors: list[tuple[int, str]]) -> csr_matrix:
        flip = 0
        ph = np.ones(dim, dtype=np.complex128)
        for s, ax in factors:
            if ax == "x":
                ph *= spins_x[s]
            elif ax == "z":
                flip |= 1 << s
            else:
                flip |= 1 << s
                ph *= phases_y[s]
        return csr_matrix((ph, (idx, idx ^ flip)), shape=(dim, dim))

    h_mat = csr_matrix((dim, dim), dtype=np.complex128)
    for u, v, ax in torus_bonds(nx, ny):
        h_mat = h_mat - _op_csr([(u, ax), (v, ax)])

    eye_sp = csr_matrix((np.ones(dim, dtype=np.complex128), (idx, idx)), shape=(dim, dim))
    proj = eye_sp
    for fac, w_val in zip(torus_plaquettes(nx, ny), target_wp, strict=True):
        proj = proj @ (0.5 * (eye_sp + float(w_val) * _op_csr(list(fac))))

    sector_dim = dim >> (nx * ny - 1)
    rng = np.random.default_rng(17)
    basis_mat = proj @ rng.standard_normal((dim, sector_dim + 16))
    u_svd, s_vals, _ = np.linalg.svd(basis_mat, full_matrices=False)
    rank = int(np.sum(s_vals > 1e-10))
    assert rank == sector_dim
    q_sub = u_svd[:, :rank]
    h_sub = q_sub.conj().T @ (h_mat @ q_sub)
    return float(np.linalg.eigvalsh(h_sub)[0]) / float(n_sites)


def test_su2_gauge_rodrigues_and_thouless_block_diagonality_theorem() -> None:
    """Verify exact Rodrigues exponential vs scipy.linalg.expm and U_cf^dag O U_cf = diag(U_O, U_O*)."""
    rng = np.random.default_rng(2026)
    thetas = rng.uniform(-np.pi, np.pi, size=(6, 3))
    thetas[0] = np.zeros(3)
    thetas[1] = np.array([np.pi, 0.0, 0.0])
    thetas[2] = np.array([2.0 * np.pi, 0.0, 0.0])

    o_mat = build_site_su2_gauge_orthogonal(thetas)
    for s in range(6):
        gen = thetas[s, 0] * SO4_L["x"] + thetas[s, 1] * SO4_L["y"] + thetas[s, 2] * SO4_L["z"]
        np.testing.assert_allclose(
            o_mat[4 * s : 4 * s + 4, 4 * s : 4 * s + 4], expm(gen), atol=1e-14
        )
    np.testing.assert_allclose(o_mat.T @ o_mat, np.eye(24), atol=1e-14)
    assert np.linalg.det(o_mat) == pytest.approx(1.0, abs=1e-13)
    np.testing.assert_allclose(o_mat[8:12, 8:12], -np.eye(4), atol=1e-14)

    for ref in (np.zeros(6, dtype=np.int8), np.array([0, 1, 1, 0, 1, 0], dtype=np.int8)):
        chk = verify_thouless_su2_block_diagonality(thetas, ref)
        assert chk.max_off_block < 1e-14
        assert chk.max_conj_err < 1e-14
        assert chk.max_unitarity_err < 1e-14
        assert chk.max_det_err < 1e-14


def test_staggered_su2_r_gutzwiller_exact_and_vmc_invariance_n8_and_n18_h025() -> None:
    """Verify exact cond, bitstring E_loc, fidelity, and 50-chain VMC invariance across staggered SU(2)_R."""
    _, st25 = _solve_reference_states()
    results = evaluate_su2_gauge_stability(
        2,
        2,
        0.25,
        st25,
        n_chains=50,
        burn_in_sweeps=15,
        samples_per_chain=40,
        seed_base=3100,
        independent_vmc_seeds=True,
    )
    assert len(results) == 5

    for r in results:
        if r.pattern_name != "identity":
            assert r.matrix_rotation_norm > 1.0
            assert r.z_eff_diff_norm > 1.0
        assert r.cond_rel_diff < 1e-12
        assert r.max_probe_amp_diff < 1e-12
        assert r.max_probe_eloc_diff < 1e-12
        assert r.exact_fidelity is not None
        assert abs(1.0 - r.exact_fidelity) < 1e-13
        assert r.delta_e_exact is not None and r.delta_e_exact < 1e-13
        assert r.delta_m111_exact is not None and r.delta_m111_exact < 1e-13
        assert r.delta_wp_exact is not None and r.delta_wp_exact < 1e-13
        assert r.diagnostics.split_rhat < 1.05
        assert r.diagnostics.acc_rate_1flip > 0.1
        assert r.diagnostics.acc_rate_2flip > 0.1
        assert r.diagnostics.tau_int_wolff == pytest.approx(0.5 * r.diagnostics.g_sokal, abs=1e-14)
        assert r.diagnostics.n_eff_sokal == pytest.approx(
            (50.0 * 40.0) / (2.0 * r.diagnostics.tau_int_wolff), rel=1e-12
        )
        ci_width = r.vmc_e_ci_high - r.vmc_e_ci_low
        assert r.vmc_delta_e_vs_base < 1.5 * max(0.01, ci_width)

    # N=18 (3x3) probe-level exact invariance + VMC check under worst-case pi flips
    pats18 = make_staggered_su2_patterns(18, seed=91)
    res18 = evaluate_su2_gauge_stability(
        3,
        3,
        0.25,
        st25,
        patterns={"identity": pats18["identity"], "pi_flips_xyz": pats18["pi_flips_xyz"]},
        n_chains=12,
        burn_in_sweeps=8,
        samples_per_chain=16,
        seed_base=3900,
        independent_vmc_seeds=True,
    )
    assert res18[1].matrix_rotation_norm > 10.0
    assert res18[1].cond_rel_diff < 1e-12
    assert res18[1].max_probe_amp_diff < 1e-11
    assert res18[1].max_probe_eloc_diff < 1e-11


def test_n_decoupling_breaks_continuous_su2_r_while_preserving_z2_and_pure_spin_r() -> None:
    """Explicit F1 Lakatosian contrast: N^a breaks continuous SU(2)_R -> Z_2, whereas O A O^T & R^a preserve SU(2)_R."""
    _, st25 = _solve_reference_states()
    res = evaluate_n_decoupling_su2_breaking(2, 2, 0.25, st25)
    assert abs(1.0 - res.direct_conjugation_fidelity) < 1e-13
    assert abs(1.0 - res.pure_spin_r_decoupling_fidelity) < 1e-13
    assert abs(1.0 - res.n_decoupling_z2_fidelity) < 1e-13
    assert res.n_decoupling_staggered_fidelity < 0.85


def test_flux_sector_zero_variance_ed_recovery_and_2x3_pi_flux_parity_repair() -> None:
    """Verify exact flux-sector ED match, zero variance, N=8 pi-flux recovery, and 2x3 pi-flux 1-QP parity repair."""
    st0, _ = _solve_reference_states()

    res8_pi = evaluate_flux_sector_gutzwiller(
        2, 2, 0.0, st0, "pi_flux", n_chains=16, burn_in_sweeps=8, samples_per_chain=20
    )
    assert not res8_pi.parity_repaired
    assert res8_pi.e_proj_exact == pytest.approx(-float(np.sqrt(3.0)) / 2.0, abs=1e-11)
    assert res8_pi.vmc_e_mean == pytest.approx(-float(np.sqrt(3.0)) / 2.0, abs=1e-10)
    assert res8_pi.plaquette_flux_exact == pytest.approx(-1.0, abs=1e-12)
    assert res8_pi.exact_energy_variance is not None and res8_pi.exact_energy_variance < 1e-13
    assert res8_pi.vmc_e_loc_std < 1e-12
    assert (
        res8_pi.vmc_e_ci_low
        == res8_pi.vmc_e_ci_high
        == pytest.approx(res8_pi.vmc_e_mean, abs=1e-12)
    )

    res8_2v = evaluate_flux_sector_gutzwiller(
        2, 2, 0.0, st0, "two_vison", n_chains=12, burn_in_sweeps=6, samples_per_chain=16
    )
    assert res8_2v.plaquette_flux_exact == pytest.approx(0.0, abs=1e-12)
    assert res8_2v.exact_energy_variance is not None and res8_2v.exact_energy_variance < 1e-13
    assert res8_2v.vmc_e_loc_std < 1e-12

    e12_0_ed = _exact_spin_ed_sector_ground_energy(2, 3, [+1] * 6)
    e12_pi_ed = _exact_spin_ed_sector_ground_energy(2, 3, [-1] * 6)
    touch_idx = [
        pi
        for pi, fac in enumerate(torus_plaquettes(2, 3))
        if {torus_site(2, 3, 0, 0, 0), torus_site(2, 3, 0, 0, 1)} <= {s for s, _ in fac}
    ]
    e12_2v_ed = _exact_spin_ed_sector_ground_energy(
        2, 3, [-1 if i in touch_idx else +1 for i in range(6)]
    )

    sector_cases: tuple[tuple[FluxSectorName, float, float, bool], ...] = (
        ("zero_flux", e12_0_ed, +1.0, False),
        ("two_vison", e12_2v_ed, +1.0 / 3.0, False),
        ("pi_flux", e12_pi_ed, -1.0, True),
    )
    for sector, ed_val, wp_exp, expect_repair in sector_cases:
        res12 = evaluate_flux_sector_gutzwiller(
            2, 3, 0.0, st0, sector, n_chains=12, burn_in_sweeps=6, samples_per_chain=16
        )
        assert res12.parity_repaired is expect_repair
        assert res12.e_proj_exact == pytest.approx(ed_val, abs=1e-10)
        assert res12.vmc_e_mean == pytest.approx(ed_val, abs=1e-10)
        assert res12.plaquette_flux_exact == pytest.approx(wp_exp, abs=1e-12)
        assert res12.exact_energy_variance is not None and res12.exact_energy_variance < 1e-13
        assert res12.vmc_e_loc_std < 1e-12

    with pytest.raises(ValueError, match="topologically forbidden"):
        build_flux_sector_bonds(3, 3, "pi_flux")

    res18_0 = evaluate_flux_sector_gutzwiller(
        3, 3, 0.0, st0, "zero_flux", n_chains=10, burn_in_sweeps=6, samples_per_chain=16
    )
    assert res18_0.vmc_e_mean == pytest.approx(-0.793972367896, abs=1e-9)
    assert res18_0.vmc_e_loc_std < 1e-12

    res18_2v = evaluate_flux_sector_gutzwiller(
        3, 3, 0.0, st0, "two_vison", n_chains=10, burn_in_sweeps=6, samples_per_chain=16
    )
    assert res18_2v.parity_repaired
    assert res18_2v.vmc_e_loc_std < 1e-12
    assert res18_2v.vmc_e_mean > res18_0.vmc_e_mean


def test_autocorrelation_geyer_ips_split_rhat_and_bca_calibration() -> None:
    """Verify AR(1) exact autocorrelation recovery, split-Rhat trap detection, and BCa CI calibration."""
    rng = np.random.default_rng(4242)
    phi = 0.6
    m_chains, n_steps = 60, 200
    ar1 = np.empty((m_chains, n_steps), dtype=np.float64)
    sigma_eps = float(np.sqrt(1.0 - phi * phi))
    x = rng.normal(0.0, 1.0, size=m_chains)
    for t in range(n_steps):
        x = phi * x + sigma_eps * rng.normal(0.0, 1.0, size=m_chains)
        ar1[:, t] = x
    diag_ar1 = compute_chain_diagnostics(ar1, acc_rate_1flip=0.5, acc_rate_2flip=0.4)
    assert 3.0 < diag_ar1.g_ips_geyer < 5.0
    assert 3.0 < diag_ar1.g_sokal < 5.0
    assert diag_ar1.tau_int_wolff == pytest.approx(0.5 * diag_ar1.g_sokal, abs=1e-14)
    assert diag_ar1.split_rhat < 1.05
    assert abs(diag_ar1.mcse_chains - diag_ar1.mcse_autocorr) / diag_ar1.mcse_autocorr < 0.35

    trapped = rng.normal(0.0, 0.1, size=(20, 80)) + np.linspace(-2.0, 2.0, 20)[:, None]
    diag_trap = compute_chain_diagnostics(trapped)
    assert diag_trap.split_rhat > 1.2

    with pytest.raises(ValueError, match="finite"):
        bca_confidence_interval(np.array([1.0, np.nan, 2.0]))
    covered = 0
    n_trials = 120
    for trial in range(n_trials):
        sample = rng.exponential(scale=1.0, size=50)
        lo, hi = bca_confidence_interval(sample, n_boot=1000, seed=5000 + trial)
        if lo <= 1.0 <= hi:
            covered += 1
    assert 0.88 <= covered / float(n_trials) <= 1.0


def test_wave2b_exported_benchmark_provenance_and_invariants() -> None:
    """Verify ACID header, SHA-256 sidecar, zero-variance rows, and SU(2)_R invariance in exported .dat."""
    dat_path = REPO_ROOT / "data/processed/benchmarks/wave2/hc02_gutzwiller_su2_flux_stability.dat"
    meta_path = dat_path.with_suffix(".meta.json")
    assert dat_path.is_file(), f"Missing benchmark table: {dat_path}"
    assert meta_path.is_file(), f"Missing benchmark sidecar: {meta_path}"

    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    assert meta["sha256"] == hashlib.sha256(dat_path.read_bytes()).hexdigest()
    assert meta["unit_convention"] == "pauli_k_minus_one_field_111"

    lines = dat_path.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("# Benchmark: hc02_gutzwiller_su2_flux_stability.dat")
    assert lines[1].startswith("# Generated:")
    assert lines[2].startswith("# Git Commit:")
    data_rows = [ln.split() for ln in lines if ln.strip() and not ln.startswith("#")]
    assert len(data_rows) == 20

    for row in data_rows:
        kind = row[0]
        probe_eloc_diff = float(row[12])
        e_loc_std = float(row[17])
        split_rhat = float(row[23])
        assert split_rhat < 1.05
        assert probe_eloc_diff < 1e-11
        if kind == "flux":
            assert e_loc_std < 1e-12
        elif row[11] != "nan":
            assert abs(1.0 - float(row[11])) < 1e-6
