"""Unit tests for Wave 2C (UOW-05 / Q-L2-1): Two-vison binding V(R,L), Gutzwiller projection, and disorder stability."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import InhomogeneousMeanFieldState, MeanFieldState
from kitaev_mf.services.disorder_stability_wave2c import (
    evaluate_bond_disordered_two_vison_ensemble,
    evaluate_vacancy_pair_flux_binding,
    evaluate_vacancy_stability_and_p0_parity,
    probe_single_vacancy_kkt_level_crossing,
)
from kitaev_mf.services.gutzwiller_projection import (
    _exact_spin_observables,
    _integrated_autocorr,
    assemble_torus_majorana_matrix,
    bca_confidence_interval,
    build_thouless_data,
    exact_wavefunction_from_thouless,
)
from kitaev_mf.services.gutzwiller_stability import (
    evaluate_flux_sector_gutzwiller,
    holonomy_matched_ref,
)
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field
from kitaev_mf.services.vison_binding_wave2c import (
    compare_oscillatory_mzm_scaling,
    evaluate_deconfinement_vs_confinement,
    evaluate_exact_projected_two_vison_h0,
)
from kitaev_mf.services.vortex_pair_scaling import TwoVisonPointResult

REPO_ROOT = Path(__file__).resolve().parents[2]


def _solve_h0_state() -> MeanFieldState:
    cell = make_honeycomb_2site()
    cfg = SolverConfig(k_grid_shape=(24, 24), tol_correlator=1e-10, tol_gauge=1e-10)
    st0, met = solve_self_consistent_mean_field(
        cell, ModelParameters(K=-1.0), make_qsl_seed(cell), cfg
    )
    assert met.converged
    return st0


def test_uow04_f1_holonomy_ref_zero_mode_and_energy_twist_ranking() -> None:
    """UOW-04 F1: Deterministic 1D yz-loop holonomy ref, 3|L MZM ground state, & energy ranking at L>=4."""
    st0 = _solve_h0_state()

    # 3x3 two_vison with allow_zero_mode=True recovers exact 2^18 Spin-Lanczos ED -0.762134236785 at bc=(1,1)
    res_3x3_zm = evaluate_flux_sector_gutzwiller(
        3,
        3,
        0.0,
        st0,
        "two_vison",
        n_chains=8,
        burn_in_sweeps=4,
        samples_per_chain=12,
        allow_zero_mode=True,
    )
    assert res_3x3_zm.bc_twist == (1, 1)
    assert not res_3x3_zm.parity_repaired
    assert res_3x3_zm.vmc_e_mean == pytest.approx(-0.762134236785, abs=1e-9)
    assert res_3x3_zm.vmc_e_loc_std < 1e-12

    # L=4 (N=32) two_vison: deterministic holonomy ref succeeds in 1 step & ranks by energy -> bc=(1,1)
    res_4x4 = evaluate_flux_sector_gutzwiller(
        4,
        4,
        0.0,
        st0,
        "two_vison",
        n_chains=2,
        burn_in_sweeps=2,
        samples_per_chain=4,
        allow_zero_mode=True,
    )
    assert res_4x4.bc_twist == (1, 1)
    assert res_4x4.vmc_e_mean == pytest.approx(-0.7836309302, abs=1e-8)
    assert res_4x4.vmc_e_loc_std < 1e-12

    # Direct holonomy check: ref has parity flips exactly on columns where -bc_y * prod u_b == -1
    ref_4 = holonomy_matched_ref(4, 4, (1, 1), [(0, 0, "z")])
    assert int(np.sum(ref_4)) == 3  # x=1,2,3 have prod=+1 => eta=-1; x=0 has prod=-1 => eta=+1


def test_uow04_f2_matter_mode_filtering_preserves_zero_flux_plaquettes() -> None:
    """UOW-04 F2: 1-QP parity repair filters for c^0 matter modes, preventing b-gauge W_p contamination."""
    st0 = _solve_h0_state()
    params0 = ModelParameters(K=-1.0)
    a_mat = assemble_torus_majorana_matrix(2, 2, params0, st0, bc_twist=(1, 1), flipped_bonds=[])
    det_ref = holonomy_matched_ref(2, 2, (1, 1), [])
    tdata = build_thouless_data(
        a_mat, 8, excite_lowest_qp=True, candidate_refs=(det_ref,), allow_zero_mode=False
    )
    assert tdata is not None
    psi = exact_wavefunction_from_thouless(tdata[0], tdata[1])
    assert psi is not None
    _, _, wp_mean, _ = _exact_spin_observables(2, 2, 0.0, psi)
    assert wp_mean == pytest.approx(1.0, abs=1e-12)


def test_uow04_f3_bca_confidence_midp_and_uncapped_sokal() -> None:
    """UOW-04 F3: Validate BCa confidence bounds, Mid-P tie symmetry, and uncapped Sokal autocorr."""
    sample = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float64)
    with pytest.raises(ValueError, match="confidence"):
        bca_confidence_interval(sample, confidence=-0.5)
    with pytest.raises(ValueError, match="confidence"):
        bca_confidence_interval(sample, confidence=1.0)

    # Reflection symmetry under Mid-P on discrete tied data
    a0 = np.zeros(50, dtype=np.float64)
    a0[0] = 1.0
    a1 = 1.0 - a0
    lo0, hi0 = bca_confidence_interval(a0, n_boot=2000, seed=77, mid_p=True)
    lo1, hi1 = bca_confidence_interval(a1, n_boot=2000, seed=77, mid_p=True)
    assert abs((1.0 - hi0) - lo1) < 1e-12
    assert abs((1.0 - lo0) - hi1) < 1e-12

    # Uncapped Sokal autocorr on strongly correlated AR(1) phi=0.95 (exact g = 39.0)
    rng = np.random.default_rng(9876)
    phi = 0.95
    n_steps = 4000
    x = np.empty(n_steps, dtype=np.float64)
    cur = 0.0
    sig = float(np.sqrt(1.0 - phi * phi))
    for t in range(n_steps):
        cur = phi * cur + sig * float(rng.standard_normal())
        x[t] = cur
    g_uncapped = _integrated_autocorr(x)
    assert g_uncapped > 28.0


def test_exact_projected_two_vison_restores_lieb_positivity_and_deconfinement() -> None:
    """Exact Zero-Variance V_proj(R,L) restores Lieb positivity on 3|L (L=6) and favors deconfinement."""
    st0 = _solve_h0_state()
    e0_p: dict[int, tuple[float, tuple[int, int]]] = {}
    e0_mf: dict[int, tuple[float, InhomogeneousMeanFieldState]] = {}

    # Resonant L=6 (N=72): unprojected scMFT turns negative at R=2,3, whereas V_proj stays > 0!
    pts_l6 = [
        evaluate_exact_projected_two_vison_h0(
            6, r, st0, n_zero_var_probes=3, e0_proj_cache=e0_p, e0_scmft_cache=e0_mf
        )
        for r in (1, 2, 3)
    ]
    assert all(p.zero_var_probe_std < 1e-12 for p in pts_l6)
    assert all(p.v_proj > 0.20 for p in pts_l6)
    assert pts_l6[1].v_scmft < 0.0 and pts_l6[2].v_scmft < 0.0

    # Non-resonant L=4, 5: Zero Variance and deconfinement model selection check
    pts_nr = [
        evaluate_exact_projected_two_vison_h0(
            L, r, st0, n_zero_var_probes=3, e0_proj_cache=e0_p, e0_scmft_cache=e0_mf
        )
        for L in (4, 5)
        for r in range(1, L // 2 + 1)
    ]
    assert all(p.zero_var_probe_std < 1e-12 for p in pts_nr)

    # Synthetic deconfined plateau test (V ~ 0.3072 across N_cut = 1..6) vs sigma_0 = 0.20 string
    deconf = evaluate_deconfinement_vs_confinement(
        [1, 2, 3, 4, 5, 6],
        [1, 2, 3, 4, 5, 6],
        [0.263, 0.295, 0.290, 0.306, 0.305, 0.307],
    )
    assert deconf.delta_bic_deconf_over_conf > 15.0
    assert abs(deconf.sigma_eff_unconstrained) < 0.02
    assert deconf.delta_v_rel_err_vs_kitaev < 0.05


def test_oscillatory_mzm_envelope_resolves_wave1_intervalley_ambiguity() -> None:
    """Lahtinen oscillatory envelope |cos(2*pi*R/3 + phi)| exp(-R/xi) beats monotone exp by Delta_BIC >> 10."""
    w1_path = REPO_ROOT / "data/processed/benchmarks/wave1/hc03_multi_l_vison_scaling.dat"
    raw_rows = [
        ln.split()
        for ln in w1_path.read_text(encoding="utf-8").splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    xis: list[float] = []
    for kap in (0.04, 0.06, 0.10):
        pts = [
            TwoVisonPointResult(
                L=int(r[0]),
                resonance_class=r[1],
                separation=int(r[2]),
                kappa_static=float(r[3]),
                h_field=float(r[4]),
                converged=bool(int(r[5])),
                iterations=int(r[6]),
                max_gauge_residual=float(r[7]),
                defect_count=int(r[8]),
                min_abs_pfaffian_wp=float(r[9]),
                e_total=float(r[10]),
                v_pair=float(r[11]),
                omega_mzm_can=float(r[12]),
                delta_bulk_can=float(r[13]),
                e0_over_delta=float(r[14]),
                omega_over_delta=float(r[15]),
                mzm_core_weight=float(r[16]),
                bulk_core_weight=float(r[17]),
            )
            for r in raw_rows
            if abs(float(r[3]) - kap) < 1e-6 and float(r[4]) == 0.0
        ]
        comp = compare_oscillatory_mzm_scaling(pts)
        assert comp.delta_bic_osc_exp_vs_mono_exp > 15.0
        assert comp.leave_lmax_rel_err_osc < comp.leave_lmax_rel_err_mono
        assert comp.leave_lmin_rel_err_osc < comp.leave_lmin_rel_err_mono
        xis.append(comp.xi_osc_linear)
    # Strictly contracting correlation length with increasing chiral gap kappa
    assert xis[0] > xis[1] > xis[2] > 0.5


def test_inhomogeneous_bond_disordered_two_vison_trapping() -> None:
    """Quenched bond disorder delta <= 0.20 preserves 100% two-vison trapping and KKT gauge neutrality."""
    rows = evaluate_bond_disordered_two_vison_ensemble(
        6, 3, (0.0, 0.10, 0.20), kappa_static=0.06, h_field=0.0, n_seeds=3, seed_base=202611
    )
    assert len(rows) == 3
    for r in rows:
        assert r.n_converged == r.n_seeds
        assert r.n_trapped_2v == r.n_seeds
        assert r.max_gauge_residual < 1e-10
        assert r.min_abs_pfaffian_wp > 0.999
        assert r.max_pade_vs_eigh_diff < 1e-11
        assert r.pade32_iters_mean < r.ns_iters_mean
        assert r.v_pair_ci_low <= r.v_pair_mean <= r.v_pair_ci_high
        assert r.core_enrichment_mean > 2.5


def test_vacancy_p0_pfaffian_parity_positivity_diamond_and_wcm_flux_binding() -> None:
    """Verify odd/even N_act O(4) P_0 Pfaffian parity, KKT crossing obstruction, Pade advantage & WCM binding."""
    # 1. Sublattice-balanced clusters (n_vac in {0, 2}) with bond disorder delta=0.05 on L=4 converge cleanly
    vrows_bal = evaluate_vacancy_stability_and_p0_parity(
        4, (0.0, 0.10), (0, 2), delta_disorder=0.05, seed=202611
    )
    assert len(vrows_bal) == 4
    for v in vrows_bal:
        assert v.converged_pade
        assert v.max_gauge_residual < 1e-10
        assert v.e_double_counting_err < 1e-11
        assert v.pade_vs_eigh_max_diff < 1e-11
        assert v.p0_pfaffian_parity_ok
        assert v.max_positivity_diamond_err < 1e-13
        assert v.max_site_spectral_norm <= 1.0 + 1e-12
        if v.h_field > 0.0:
            assert v.pade32_purify_iters < v.ns_purify_iters
            assert abs(v.pfaffian_kmf) > 0.5
            assert v.pfaffian_v2 == pytest.approx(v.expected_p0_sign * v.pfaffian_kmf, rel=1e-10)

    # 2. Clean resonant L=6 odd-active-site single vacancy (N_act=71) converges cleanly & flips Pf across P_0!
    vrows_odd = evaluate_vacancy_stability_and_p0_parity(6, (0.10,), (1,), delta_disorder=0.0)
    v_odd = vrows_odd[0]
    assert v_odd.n_active_sites == 71 and v_odd.expected_p0_sign == -1
    assert v_odd.converged_pade
    assert v_odd.max_gauge_residual < 1e-10
    assert v_odd.e_double_counting_err < 1e-11
    assert v_odd.pade32_purify_iters < v_odd.ns_purify_iters
    assert abs(v_odd.pfaffian_kmf) > 0.5
    assert v_odd.pfaffian_v2 == pytest.approx(-v_odd.pfaffian_kmf, rel=1e-10)
    assert v_odd.p0_pfaffian_parity_ok
    assert v_odd.max_positivity_diamond_err < 1e-13

    # 3. Lakatosian counterexample: deterministic bulk-state C_3 scan shows T=0 KKT crossing at lam ~ 2h_gamma
    cross = probe_single_vacancy_kkt_level_crossing(L=4, h_field=0.10, delta_lambda=0.01)
    assert cross.g_z_below < -0.35 and cross.g_z_above > 0.35
    assert cross.gauge_jump > 0.85
    assert 1e-6 < cross.min_active_sv < 1e-4
    cross_rep = probe_single_vacancy_kkt_level_crossing(L=4, h_field=0.10, delta_lambda=0.01)
    assert cross_rep.g_z_below == pytest.approx(cross.g_z_below, abs=1e-14)
    assert cross_rep.g_z_above == pytest.approx(cross.g_z_above, abs=1e-14)
    assert cross_rep.gauge_jump == pytest.approx(cross.gauge_jump, abs=1e-14)

    # 4. Willans-Chalker-Moessner / Kao-Knolle Eq. (8) A-B vacancy-pair bound-flux energy < 0
    bind_l8 = evaluate_vacancy_pair_flux_binding(8)
    assert bind_l8.converged_zero and bind_l8.converged_bound
    assert bind_l8.max_gauge_residual < 1e-11
    assert -0.06 < bind_l8.e_bind_per_vacancy < -0.015


def test_wave2c_exported_benchmarks_acid_provenance_and_invariants() -> None:
    """Verify ACID headers, Git Dirty tags (UOW-04 F4), SHA-256 sidecars, and physical invariants."""
    bdir = REPO_ROOT / "data/processed/benchmarks/wave2"
    vison_dat = bdir / "q_l2_1_two_vison_binding_projection.dat"
    dis_dat = bdir / "q_l2_1_disorder_vacancy_stability.dat"
    for dat in (vison_dat, dis_dat):
        meta_path = dat.with_suffix(".meta.json")
        assert dat.is_file(), f"Missing benchmark file: {dat}"
        assert meta_path.is_file(), f"Missing benchmark sidecar: {meta_path}"
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        assert meta["sha256"] == hashlib.sha256(dat.read_bytes()).hexdigest()
        assert "git_dirty" in meta and isinstance(meta["git_dirty"], bool)
        assert "base_git_sha" in meta and len(meta["base_git_sha"]) == 40
        lines = dat.read_text(encoding="utf-8").splitlines()
        assert lines[0].startswith(f"# Benchmark: {dat.name}")
        assert any(ln.startswith("# Git Dirty:") for ln in lines[:8])

    # Physical invariants inside vison binding benchmark
    v_meta = json.loads(vison_dat.with_suffix(".meta.json").read_text(encoding="utf-8"))
    v_rows = [
        ln.split()
        for ln in vison_dat.read_text(encoding="utf-8").splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    proj_rows = [r for r in v_rows if r[0] == "proj_h0"]
    assert len(proj_rows) == 14
    for r in proj_rows:
        assert float(r[9]) > 0.0  # Strict Lieb positivity of V_proj across all L=4..8 incl. L=6
        assert float(r[12]) < 1e-12  # Certified Zero Variance
    for kap_key in ("0.04", "0.06", "0.10"):
        osc = v_meta["parameters"]["oscillatory_mzm_by_kappa"][kap_key]
        assert osc["delta_bic_osc_exp_vs_mono_exp"] > 15.0
        assert osc["v_bind_r1_lmax"] < 0.0

    # Physical invariants inside disorder/vacancy benchmark
    d_rows = [
        ln.split()
        for ln in dis_dat.read_text(encoding="utf-8").splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    bond_rows = [r for r in d_rows if r[0] == "bond_2v"]
    vac_rows = [r for r in d_rows if r[0] == "vacancy"]
    bind_rows = [r for r in d_rows if r[0] == "vac_bind"]
    assert len(bond_rows) == 10
    assert len(vac_rows) == 15
    assert len(bind_rows) == 3
    for r in bond_rows:
        assert int(r[5]) == int(r[4]) and int(r[6]) == int(r[4])  # 100% converged & trapped
        assert float(r[7]) < 1e-10
    for r in vac_rows:
        assert int(r[6]) == 1  # Exact O(4) P_0 Pfaffian parity holds on all 15 rows
        L_val, h_val, delta_val, n_vac = int(r[1]), float(r[2]), float(r[3]), int(r[4])
        if L_val == 4 and h_val > 0.0 and delta_val > 0.0 and n_vac == 1:
            # Documented Lakatosian pure-state T=0 KKT level-crossing obstruction on unpaired gapped vacancy
            assert int(r[5]) == 0 and float(r[7]) > 0.4
        else:
            assert int(r[5]) == 1 and float(r[7]) < 1e-10
    for r in bind_rows:
        assert int(r[5]) == 1 and float(r[11]) < 0.0  # negative vacancy flux binding
