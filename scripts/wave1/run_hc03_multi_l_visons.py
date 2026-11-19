#!/usr/bin/env python3
"""Wave 1 HC-03 runner: Stratified multi-L two-vison energetics, chiral NNN bridge, and MZM scaling."""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import InhomogeneousMeanFieldState
from kitaev_mf.services.purification import purify_newton_schulz
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field
from kitaev_mf.services.vortex_pair_scaling import (
    TwoVisonPointResult,
    build_covariant_nnn_matrix,
    compare_exp_vs_powerlaw_scaling,
    measure_two_vison_point,
    solve_inhomogeneous_torus_sc,
)


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def _git_sha(root: Path) -> str:
    return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()


def _verify_kwg_bridge(root: Path) -> dict[str, dict[str, float]]:
    """Verify exact single-particle static chiral torus match against kitaev-walk-geometry."""
    ext_path = (
        root / "data" / "external" / "kitaev_walk_geometry" / "chiral_vison_tmatrix_scaling.dat"
    )
    kwg_ref: dict[tuple[float, int], float] = {}
    if ext_path.exists():
        for raw in ext_path.read_text(encoding="utf-8").splitlines():
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("model "):
                continue
            parts = line.split()
            if parts[0] == "covariant_nnn":
                kwg_ref[(float(parts[1]), int(parts[3]))] = float(parts[5])

    cell = make_honeycomb_2site()
    cfg = SolverConfig(k_grid_shape=(24, 24), tol_correlator=1e-10, tol_gauge=1e-10)
    st0, _ = solve_self_consistent_mean_field(
        cell, ModelParameters(K=-1.0), make_qsl_seed(cell), cfg
    )

    out: dict[str, dict[str, float]] = {}
    for kappa_kwg, L in ((0.08, 12), (0.08, 18), (0.12, 12), (0.12, 18)):
        kappa_static = 0.5 * kappa_kwg
        cl = RealSpaceTorusCluster(L)
        cut = [(0, 0, 0)]
        signs = np.ones((L, L, 3), dtype=np.float64)
        signs[0, 0, 0] = -1.0
        a_tot = cl.assemble_skew_hamiltonian(
            ModelParameters(K=-1.0), st0, flipped_bonds=cut
        ) + build_covariant_nnn_matrix(cl, kappa_static, signs)
        w2 = np.linalg.eigvalsh(0.5 * ((-a_tot @ a_tot) + (-a_tot @ a_tot).T))
        eps_min = float(np.sqrt(max(float(np.min(w2)), 0.0)))
        delta_can = 12.0 * float(np.sqrt(3.0)) * kappa_static
        e0_over_delta = eps_min / delta_can
        ref = kwg_ref.get((round(kappa_kwg, 6), L), float("nan"))
        out[f"kKWG_{kappa_kwg:.2f}_L{L}"] = {
            "kappa_static": kappa_static,
            "mf_e0_over_delta": e0_over_delta,
            "kwg_e0_torus_over_delta": ref,
            "abs_diff": abs(e0_over_delta - ref) if not np.isnan(ref) else float("nan"),
        }
    return out


def _verify_purify_vs_eigh() -> float:
    cl = RealSpaceTorusCluster(12)
    params = ModelParameters(K=-1.0, h=(0.20 / float(np.sqrt(3.0)),) * 3)
    pure = InhomogeneousMeanFieldState.from_pure_kitaev(cl.num_cells)
    seed_2v, _ = cl.create_unpinned_vortex_pair_seed(pure, (2, 2), (8, 2))
    om_e, st_e, _, _, _ = solve_inhomogeneous_torus_sc(
        cl, params, seed_2v, reference_clean_state=pure, kappa_static=0.0, use_eigh=True, tol=1e-11
    )
    a_conv = cl.assemble_skew_hamiltonian(params, st_e)
    om_p, _ = purify_newton_schulz(a_conv, tol=1e-13, max_iter=35)
    return float(np.max(np.abs(om_e - om_p)))


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    out_dir = root / "data" / "processed" / "benchmarks" / "wave1"
    dat_path = out_dir / "hc03_multi_l_vison_scaling.dat"
    sha = _git_sha(root)

    cache: dict[tuple[int, float, float], tuple[float, float, InhomogeneousMeanFieldState]] = {}
    all_rows: list[TwoVisonPointResult] = []
    summaries_by_kappa: dict[str, dict[str, float]] = {}

    for kappa in (0.04, 0.06, 0.10):
        k_rows: list[TwoVisonPointResult] = []
        for L in (8, 10, 12, 15, 16, 18):
            tol_L = 1e-12 if L <= 12 else 1e-10
            seps = list(range(1, L // 2 + 1)) if L <= 12 else [1, 2, 4, 6, L // 2]
            for r_sep in seps:
                row = measure_two_vison_point(
                    L,
                    r_sep,
                    kappa_static=kappa,
                    h_field=0.0,
                    e0_cache=cache,
                    use_eigh=True,
                    tol=tol_L,
                )
                k_rows.append(row)
                all_rows.append(row)

        if abs(kappa - 0.06) < 1e-12:
            for r_sep in (2, 5, 8, 10):
                row = measure_two_vison_point(
                    20,
                    r_sep,
                    kappa_static=kappa,
                    h_field=0.0,
                    e0_cache=cache,
                    use_eigh=True,
                    tol=1e-10,
                )
                k_rows.append(row)
                all_rows.append(row)

        comp = compare_exp_vs_powerlaw_scaling(k_rows)
        comp_res = compare_exp_vs_powerlaw_scaling(
            [row for row in k_rows if row.resonance_class == "3|L"]
        )
        comp_nonres = compare_exp_vs_powerlaw_scaling(
            [row for row in k_rows if row.resonance_class == "3!|L"]
        )
        summaries_by_kappa[f"{kappa:.2f}"] = {
            "n_points": float(comp.n_points),
            "xi_fit": comp.xi_fit,
            "alpha_fit": comp.alpha_fit,
            "rss_exp": comp.rss_exp,
            "rss_pow": comp.rss_pow,
            "delta_bic_exp_over_pow": comp.delta_bic_exp_over_pow,
            "delta_aicc_exp_over_pow": comp.delta_aicc_exp_over_pow,
            "leave_lmin_rel_err": comp.leave_lmin_rel_err,
            "leave_lmax_rel_err": comp.leave_lmax_rel_err,
            "xi_fit_3divL": comp_res.xi_fit,
            "xi_fit_3ndivL": comp_nonres.xi_fit,
        }

    # Regime A: Pure [111] Zeeman self-consistent visons without artificial kappa_static
    for h_val, l_vals in ((0.0, (8, 12)), (0.20, (8, 10, 12, 16))):
        for L in l_vals:
            row = measure_two_vison_point(
                L,
                L // 2,
                kappa_static=0.0,
                h_field=h_val,
                e0_cache=cache,
                use_eigh=True,
                tol=1e-11,
            )
            all_rows.append(row)

    kwg_bridge = _verify_kwg_bridge(root)
    purify_eigh_diff = _verify_purify_vs_eigh()

    lines = [
        f"# Benchmark: {dat_path.name}",
        f"# Generated: {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
        f"# Git Commit: {sha[:7]} ({sha})",
        "# Description: HC-03 stratified multi-L (3|L vs 3!|L) two-vison energetics, chiral NNN bridge, and MZM scaling",
        "# Summary: "
        + ", ".join(
            f"k={k}:xi={v['xi_fit']:.4f},dBIC={v['delta_bic_exp_over_pow']:.2f},lmin={v['leave_lmin_rel_err']:.3f},lmax={v['leave_lmax_rel_err']:.3f}"
            for k, v in summaries_by_kappa.items()
        ),
        "# Columns: L  res_class  R_sep  kappa_static  h_field  conv  iters  max_g  defs  min_abs_pf_wp  e_total  V_pair  omega_mzm_can  delta_bulk_can  e0_over_delta  omega_over_delta  mzm_core_w  bulk_core_w",
        "#" + "-" * 78,
    ]
    for r in all_rows:
        lines.append(
            f"{r.L:3d} {r.resonance_class:5s} {r.separation:3d} {r.kappa_static:8.4f} {r.h_field:6.2f} "
            f"{int(r.converged):1d} {r.iterations:3d} {r.max_gauge_residual:10.3e} {r.defect_count:2d} "
            f"{r.min_abs_pfaffian_wp:10.6f} {r.e_total:14.8f} {r.v_pair:12.6f} "
            f"{r.omega_mzm_can:12.6e} {r.delta_bulk_can:10.6f} {r.e0_over_delta:12.6e} "
            f"{r.omega_over_delta:12.6e} {r.mzm_core_weight:8.5f} {r.bulk_core_weight:8.5f}"
        )

    _atomic_write(dat_path, "\n".join(lines) + "\n")
    meta = {
        "schema_version": 1,
        "dataset": "hc03_multi_l_vison_scaling",
        "source_repo": "kitaev-mf",
        "git_sha": sha,
        "sha256": hashlib.sha256(dat_path.read_bytes()).hexdigest(),
        "unit_convention": "pauli_k_minus_one_canonical_4HB_Delta_12sqrt3_kappa_static",
        "geometry": "Periodic LxL honeycomb torus (stratified 3|L vs 3!|L)",
        "parameters": {
            "scaling_by_kappa": summaries_by_kappa,
            "kwg_static_bridge": kwg_bridge,
            "purify_vs_eigh_omega_max_diff_L12_h020": purify_eigh_diff,
        },
    }
    _atomic_write(
        dat_path.with_suffix(".meta.json"), json.dumps(meta, indent=2, sort_keys=True) + "\n"
    )
    print(
        f"Exported HC-03 benchmark -> {dat_path} ({len(all_rows)} rows, purify diff={purify_eigh_diff:.2e})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
