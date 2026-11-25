"""ACID benchmark exporters for Wave-2C (UOW-05 / Q-L2-1) two-vison binding and disorder stability."""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import subprocess
from pathlib import Path

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import InhomogeneousMeanFieldState
from kitaev_mf.services.disorder_stability_wave2c import (
    evaluate_bond_disordered_two_vison_ensemble,
    evaluate_vacancy_pair_flux_binding,
    evaluate_vacancy_stability_and_p0_parity,
)
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field
from kitaev_mf.services.vison_binding_wave2c import (
    ProjectedTwoVisonPoint,
    compare_oscillatory_mzm_scaling,
    evaluate_deconfinement_vs_confinement,
    evaluate_exact_projected_two_vison_h0,
)
from kitaev_mf.services.vortex_pair_scaling import TwoVisonPointResult, measure_two_vison_point


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def _git_provenance(repo_root: Path) -> tuple[str, bool]:
    sha = subprocess.check_output(
        ["git", "-C", str(repo_root), "rev-parse", "HEAD"], text=True
    ).strip()
    dirty_out = subprocess.check_output(
        ["git", "-C", str(repo_root), "status", "--porcelain"], text=True
    ).strip()
    return sha, bool(dirty_out)


def export_wave2c_vison_binding_benchmark(repo_root: Path) -> tuple[Path, Path]:
    """Export Wave-2C Q-L2-1 isolated two-vison binding V(R,L) & oscillatory MZM benchmark."""
    sha, git_dirty = _git_provenance(repo_root)
    out_dir = repo_root / "data" / "processed" / "benchmarks" / "wave2"
    dat_path = out_dir / "q_l2_1_two_vison_binding_projection.dat"
    meta_path = dat_path.with_suffix(".meta.json")

    cell = make_honeycomb_2site()
    cfg = SolverConfig(k_grid_shape=(24, 24), tol_correlator=1e-10, tol_gauge=1e-10)
    st0, met0 = solve_self_consistent_mean_field(
        cell, ModelParameters(K=-1.0), make_qsl_seed(cell), cfg
    )
    if not met0.converged:
        raise RuntimeError("h=0 scMFT solve failed to converge.")

    # Section A: Exact Zero-Variance Gutzwiller V_proj(R, L) vs. unprojected V_scMFT(R, L) at h=0
    e0_p_cache: dict[int, tuple[float, tuple[int, int]]] = {}
    e0_mf_cache: dict[int, tuple[float, InhomogeneousMeanFieldState]] = {}
    proj_rows: list[ProjectedTwoVisonPoint] = []
    for L in (4, 5, 6, 7, 8):
        for r_sep in range(1, L // 2 + 1):
            proj_rows.append(
                evaluate_exact_projected_two_vison_h0(
                    L,
                    r_sep,
                    st0,
                    n_zero_var_probes=4,
                    e0_proj_cache=e0_p_cache,
                    e0_scmft_cache=e0_mf_cache,
                    seed=202605 + 100 * L,
                )
            )

    # Deconfinement vs. confinement on non-resonant L=8 (both projected and scMFT)
    l8_rows = [r for r in proj_rows if r.L == 8]
    deconf_proj_l8 = evaluate_deconfinement_vs_confinement(
        [r.separation for r in l8_rows],
        [r.n_cut for r in l8_rows],
        [r.v_proj for r in l8_rows],
    )
    deconf_proj_nonres = evaluate_deconfinement_vs_confinement(
        [r.separation for r in proj_rows if r.resonance_class == "3!|L" and r.L >= 5],
        [r.n_cut for r in proj_rows if r.resonance_class == "3!|L" and r.L >= 5],
        [r.v_proj for r in proj_rows if r.resonance_class == "3!|L" and r.L >= 5],
    )

    # Section B: Chiral kappa_static in {0.04, 0.06, 0.10} oscillatory MZM & binding sweep (L=8..16)
    chiral_cache: dict[
        tuple[int, float, float], tuple[float, float, InhomogeneousMeanFieldState]
    ] = {}
    chiral_rows: list[TwoVisonPointResult] = []
    osc_summaries: dict[str, dict[str, float]] = {}
    for kap in (0.04, 0.06, 0.10):
        k_rows: list[TwoVisonPointResult] = []
        for L in (8, 10, 12, 14, 15):
            seps = list(range(1, L // 2 + 1)) if L <= 12 else [1, 2, 3, 4, 6, L // 2]
            tol_l = 1e-12 if L <= 12 else 1e-10
            for r_sep in seps:
                pt = measure_two_vison_point(
                    L, r_sep, kappa_static=kap, h_field=0.0, e0_cache=chiral_cache, tol=tol_l
                )
                k_rows.append(pt)
                chiral_rows.append(pt)
        osc = compare_oscillatory_mzm_scaling(k_rows)
        l_max_k = max(row.L for row in k_rows)
        v_asymp = max(
            (row.v_pair for row in k_rows if row.L == l_max_k and row.separation == l_max_k // 2),
            default=0.0,
        )
        v_r1 = next((row.v_pair for row in k_rows if row.L == l_max_k and row.separation == 1), 0.0)
        osc_summaries[f"{kap:.2f}"] = {
            "n_points": float(osc.n_points),
            "xi_osc_linear": osc.xi_osc_linear,
            "phi_osc_linear": osc.phi_osc_linear,
            "xi_osc_log": osc.xi_osc_log,
            "alpha_osc_log": osc.alpha_osc_log,
            "delta_bic_osc_exp_vs_mono_exp": osc.delta_bic_osc_exp_vs_mono_exp,
            "delta_aicc_osc_exp_vs_mono_exp": osc.delta_aicc_osc_exp_vs_mono_exp,
            "delta_bic_osc_exp_vs_osc_pow_log": osc.delta_bic_osc_exp_vs_osc_pow_log,
            "delta_aicc_osc_exp_vs_osc_pow_log": osc.delta_aicc_osc_exp_vs_osc_pow_log,
            "leave_lmin_rel_err_osc": osc.leave_lmin_rel_err_osc,
            "leave_lmax_rel_err_osc": osc.leave_lmax_rel_err_osc,
            "leave_lmin_rel_err_mono": osc.leave_lmin_rel_err_mono,
            "leave_lmax_rel_err_mono": osc.leave_lmax_rel_err_mono,
            "v_r1_lmax": v_r1,
            "v_plateau_lmax": v_asymp,
            "v_bind_r1_lmax": v_r1 - v_asymp,
        }

    lines = [
        f"# Benchmark: {dat_path.name}",
        f"# Generated: {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
        f"# Git Commit: {sha[:7]} ({sha})",
        f"# Git Dirty: {str(git_dirty).lower()}",
        "# Generator: kitaev_mf.services.wave2c_benchmark.export_wave2c_vison_binding_benchmark",
        "# Description: Wave-2C (UOW-05 / Q-L2-1) exact Gutzwiller V_proj(R,L) vs scMFT V_scMFT(R,L) and chiral oscillatory MZM scaling",
        "# Columns: kind  L  res_class  R_sep  kappa_static  bc_x  bc_y  qp_exc  cond_or_iters  v_proj_or_vpair  v_scmft_or_e0  parity_shift_or_om_over_d  zv_std_or_max_g  min_pf_or_core_w",
        "#" + "-" * 130,
    ]
    for r in proj_rows:
        lines.append(
            f"{'proj_h0':8s} {r.L:3d} {r.resonance_class:5s} {r.separation:3d} {0.0:8.4f} "
            f"{r.bc_2v_proj[0]:3d} {r.bc_2v_proj[1]:3d} {r.qp_exc_2v:2d} {r.cond_2v:10.3f} "
            f"{r.v_proj:13.6f} {r.v_scmft:13.6f} {r.parity_shift_v:+13.6f} "
            f"{r.zero_var_probe_std:11.3e} {float(r.n_zero_var_probes):8.1f}"
        )
    for c in chiral_rows:
        lines.append(
            f"{'chiral':8s} {c.L:3d} {c.resonance_class:5s} {c.separation:3d} {c.kappa_static:8.4f} "
            f"{1:3d} {1:3d} {0:2d} {float(c.iterations):10.3f} "
            f"{c.v_pair:13.6f} {c.e_total:13.6f} {c.omega_over_delta:+13.6e} "
            f"{c.max_gauge_residual:11.3e} {c.mzm_core_weight:8.5f}"
        )

    _atomic_write(dat_path, "\n".join(lines) + "\n")
    meta = {
        "schema_version": 1,
        "dataset": "q_l2_1_two_vison_binding_projection",
        "source_repo": "kitaev-mf",
        "base_git_sha": sha,
        "git_sha": sha,
        "git_dirty": git_dirty,
        "sha256": hashlib.sha256(dat_path.read_bytes()).hexdigest(),
        "unit_convention": "pauli_k_minus_one_canonical_4HB_Delta_12sqrt3_kappa_static",
        "geometry": "Periodic LxL honeycomb tori L in {4..8} (exact Gutzwiller) and L in {8..15} (chiral scMFT)",
        "parameters": {
            "uow_id": "UOW-05",
            "n_proj_rows": len(proj_rows),
            "n_chiral_rows": len(chiral_rows),
            "deconfinement_proj_L8": {
                "v_inf_plateau": deconf_proj_l8.v_inf_plateau,
                "delta_v_extracted": deconf_proj_l8.delta_v_extracted,
                "delta_v_rel_err_vs_kitaev": deconf_proj_l8.delta_v_rel_err_vs_kitaev,
                "sigma_eff_unconstrained": deconf_proj_l8.sigma_eff_unconstrained,
                "delta_bic_deconf_over_conf": deconf_proj_l8.delta_bic_deconf_over_conf,
                "bayes_factor_10": deconf_proj_l8.bayes_factor_10,
            },
            "deconfinement_proj_nonresonant_L5_L7_L8": {
                "n_points": deconf_proj_nonres.n_points,
                "v_inf_plateau": deconf_proj_nonres.v_inf_plateau,
                "delta_v_extracted": deconf_proj_nonres.delta_v_extracted,
                "sigma_eff_unconstrained": deconf_proj_nonres.sigma_eff_unconstrained,
                "delta_bic_deconf_over_conf": deconf_proj_nonres.delta_bic_deconf_over_conf,
            },
            "oscillatory_mzm_by_kappa": osc_summaries,
        },
    }
    _atomic_write(meta_path, json.dumps(meta, indent=2, sort_keys=True) + "\n")
    return dat_path, meta_path


def export_wave2c_disorder_vacancy_benchmark(repo_root: Path) -> tuple[Path, Path]:
    """Export Wave-2C Q-L2-1 inhomogeneous bond-disorder, vacancy stability & flux-binding benchmark."""
    sha, git_dirty = _git_provenance(repo_root)
    out_dir = repo_root / "data" / "processed" / "benchmarks" / "wave2"
    dat_path = out_dir / "q_l2_1_disorder_vacancy_stability.dat"
    meta_path = dat_path.with_suffix(".meta.json")

    deltas = (0.00, 0.05, 0.10, 0.15, 0.20)
    ens_l6 = evaluate_bond_disordered_two_vison_ensemble(
        6, 3, deltas, kappa_static=0.06, h_field=0.0, n_seeds=8, seed_base=202611
    )
    ens_l8 = evaluate_bond_disordered_two_vison_ensemble(
        8, 4, deltas, kappa_static=0.06, h_field=0.0, n_seeds=6, seed_base=202651
    )
    vac_l4 = evaluate_vacancy_stability_and_p0_parity(
        4, (0.0, 0.05, 0.10), (0, 1, 2), delta_disorder=0.05, seed=202611
    )
    vac_l6 = evaluate_vacancy_stability_and_p0_parity(
        6, (0.0, 0.10), (0, 1, 2), delta_disorder=0.0, seed=202621
    )
    bind_rows = [evaluate_vacancy_pair_flux_binding(L) for L in (6, 8, 10)]

    lines = [
        f"# Benchmark: {dat_path.name}",
        f"# Generated: {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
        f"# Git Commit: {sha[:7]} ({sha})",
        f"# Git Dirty: {str(git_dirty).lower()}",
        "# Generator: kitaev_mf.services.wave2c_benchmark.export_wave2c_disorder_vacancy_benchmark",
        "# Description: Wave-2C (UOW-05 / Q-L2-1) inhomogeneous bond-disordered two-vison ensembles, dilute vacancy Pade vs NS + O(4) P_0 parity, and WCM flux binding",
        "# Columns: kind  L  h_or_kap  delta  n_vac_or_seeds  conv_ok  aux_ok  max_g  kappa_act  pade_it  ns_it  val_main  ci_low_or_pf_kmf  ci_high_or_pf_v2  aux_metric1  aux_metric2",
        "#" + "-" * 135,
    ]
    for e in (*ens_l6, *ens_l8):
        lines.append(
            f"{'bond_2v':9s} {e.L:3d} {e.kappa_static:8.4f} {e.delta_disorder:6.3f} {e.n_seeds:3d} "
            f"{e.n_converged:3d} {e.n_trapped_2v:3d} {e.max_gauge_residual:10.3e} {e.kappa_active_mean:10.3e} "
            f"{e.pade32_iters_mean:6.2f} {e.ns_iters_mean:6.2f} {e.v_pair_mean:12.6f} "
            f"{e.v_pair_ci_low:12.6f} {e.v_pair_ci_high:12.6f} {e.omega_mzm_mean:12.5e} {e.min_abs_pfaffian_wp:9.6f}"
        )
    for v in (*vac_l4, *vac_l6):
        lines.append(
            f"{'vacancy':9s} {v.L:3d} {v.h_field:8.4f} {v.delta_disorder:6.3f} {v.n_vacancies:3d} "
            f"{int(v.converged_pade):3d} {int(v.p0_pfaffian_parity_ok):3d} {v.max_gauge_residual:10.3e} "
            f"{v.kappa_active:10.3e} {float(v.pade32_purify_iters):6.2f} {float(v.ns_purify_iters):6.2f} "
            f"{v.e_total:12.6f} {v.pfaffian_kmf:+12.5f} {v.pfaffian_v2:+12.5f} "
            f"{v.pade_vs_eigh_max_diff:12.5e} {v.max_positivity_diamond_err:9.2e}"
        )
    for b in bind_rows:
        lines.append(
            f"{'vac_bind':9s} {b.L:3d} {0.0:8.4f} {0.0:6.3f} {2:3d} "
            f"{int(b.converged_zero and b.converged_bound):3d} {b.cut_length:3d} {b.max_gauge_residual:10.3e} "
            f"{0.0:10.3e} {0.0:6.2f} {0.0:6.2f} {b.e_bind_per_vacancy:12.6f} "
            f"{b.e_zero_flux:12.6f} {b.e_bound_flux:12.6f} {float(b.vac_a_site):12.1f} {float(b.vac_b_site):9.1f}"
        )

    _atomic_write(dat_path, "\n".join(lines) + "\n")
    meta = {
        "schema_version": 1,
        "dataset": "q_l2_1_disorder_vacancy_stability",
        "source_repo": "kitaev-mf",
        "base_git_sha": sha,
        "git_sha": sha,
        "git_dirty": git_dirty,
        "sha256": hashlib.sha256(dat_path.read_bytes()).hexdigest(),
        "unit_convention": "pauli_k_minus_one_inhomogeneous_D_80L2",
        "geometry": "Periodic LxL honeycomb tori L in {4, 6, 8, 10} with quenched bond disorder and non-magnetic vacancies",
        "parameters": {
            "uow_id": "UOW-05",
            "n_bond_2v_rows": len(ens_l6) + len(ens_l8),
            "n_vacancy_rows": len(vac_l4) + len(vac_l6),
            "n_vac_bind_rows": len(bind_rows),
            "vacancy_flux_binding_by_L": {
                str(b.L): {
                    "resonance_class": b.resonance_class,
                    "e_bind_per_vacancy": b.e_bind_per_vacancy,
                    "cut_length": b.cut_length,
                    "max_gauge_residual": b.max_gauge_residual,
                }
                for b in bind_rows
            },
            "kao_knolle_thermo_reference_e_bind": -0.0268,
        },
    }
    _atomic_write(meta_path, json.dumps(meta, indent=2, sort_keys=True) + "\n")
    return dat_path, meta_path


__all__ = [
    "export_wave2c_disorder_vacancy_benchmark",
    "export_wave2c_vison_binding_benchmark",
]
