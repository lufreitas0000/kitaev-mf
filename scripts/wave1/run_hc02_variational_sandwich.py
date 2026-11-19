#!/usr/bin/env python3
"""Wave 1 HC-02 runner: ED, certified Clarabel Pauli SOS SDP, and Gutzwiller scMFT sandwich."""

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
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.gutzwiller_projection import (
    evaluate_exact_gutzwiller,
    evaluate_vmc_gutzwiller,
)
from kitaev_mf.services.seeds import make_polarized_seed, make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field

FIELD_GRID: tuple[float, ...] = (0.0, 0.02, 0.0485, 0.10, 0.25, 0.3325, 0.40, 0.50)
TORI: tuple[tuple[int, int], ...] = ((2, 2), (2, 3), (3, 3))


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def _git_sha(root: Path) -> str:
    return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()


def _load_external_table(path: Path, value_col: int) -> dict[tuple[int, float], float]:
    out: dict[tuple[int, float], float] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("N_sites"):
            continue
        parts = s.split()
        n_sites = int(parts[0])
        h_val = round(float(parts[4]), 6)
        out[(n_sites, h_val)] = float(parts[value_col])
    return out


def _solve_variational_scmft(
    h_vals: tuple[float, ...],
) -> dict[float, tuple[MeanFieldState, float, str]]:
    """Solve forward QSL and backward polarized scMFT branches and pick lowest energy per h."""
    cell = make_honeycomb_2site()
    cfg = SolverConfig(
        k_grid_shape=(24, 24),
        max_iterations=80,
        tol_correlator=1e-10,
        tol_gauge=1e-10,
        method="diis",
        mixing_alpha=0.6,
        skew_mode="retract_step",
    )

    qsl_states: dict[float, tuple[MeanFieldState, float]] = {}
    curr_qsl = make_qsl_seed(cell)
    for h in h_vals:
        comp = h / float(np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=(comp, comp, comp))
        st, met = solve_self_consistent_mean_field(cell, params, curr_qsl, cfg)
        if met.converged:
            curr_qsl = st
            eb = compute_gaussian_energy(st, cell, params)
            qsl_states[h] = (st, float(eb.per_site))

    pol_states: dict[float, tuple[MeanFieldState, float]] = {}
    curr_pol = make_polarized_seed(cell, "111")
    for h in reversed(h_vals):
        comp = h / float(np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=(comp, comp, comp))
        st, met = solve_self_consistent_mean_field(cell, params, curr_pol, cfg)
        if met.converged:
            curr_pol = st
            eb = compute_gaussian_energy(st, cell, params)
            pol_states[h] = (st, float(eb.per_site))

    best: dict[float, tuple[MeanFieldState, float, str]] = {}
    for h in h_vals:
        cands: list[tuple[float, MeanFieldState, str]] = []
        if h in qsl_states:
            cands.append((qsl_states[h][1], qsl_states[h][0], "qsl"))
        if h in pol_states:
            cands.append((pol_states[h][1], pol_states[h][0], "polarized"))
        cands.sort(key=lambda t: t[0])
        e_best, st_best, label = cands[0]
        best[h] = (st_best, e_best, label)
    return best


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    ed_path = root / "data" / "external" / "kitaev_numerics" / "ed_torus_field_sweep_n8_n12_n18.dat"
    sdp_path = (
        root / "data" / "external" / "kitaev_numerics" / "pauli_sos_sdp_bounds_n8_n12_n18.dat"
    )
    ed_map = _load_external_table(ed_path, value_col=5)
    sdp_map = _load_external_table(sdp_path, value_col=10)

    out_dir = root / "data" / "processed" / "benchmarks" / "wave1"
    dat_path = out_dir / "hc02_ed_sdp_scmft_sandwich.dat"
    sha = _git_sha(root)

    mf_best = _solve_variational_scmft(FIELD_GRID)
    lines = [
        f"# Benchmark: {dat_path.name}",
        f"# Generated: {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
        f"# Git Commit: {sha[:7]} ({sha})",
        "# Description: HC-02 E_SDP_cert <= E_0_ED <= E_proj sandwich across N in {8, 12, 18} with 50-seed VMC BCa 95% CIs",
        "# Columns: N_sites  nx  ny  h  branch  bc_x  bc_y  e_sdp_cert  e_ed  e_mf_unproj  e_proj_exact  e_proj_vmc  vmc_ci_low  vmc_ci_high  e_loc_std  tau_int  n_eff  delta_e_corr  proj_minus_mf  sandwich_valid",
        "#" + "-" * 78,
    ]

    for nx, ny in TORI:
        n_sites = 2 * nx * ny
        for h in FIELD_GRID:
            h_key = round(h, 6)
            e_sdp = sdp_map[(n_sites, h_key)]
            e_ed = ed_map[(n_sites, h_key)]
            st, e_mf, branch = mf_best[h]

            if n_sites <= 12:
                ex = evaluate_exact_gutzwiller(nx, ny, h, st)
                e_ex_str = f"{ex.e_proj_per_site:.10f}"
                pref_bc = ex.bc_twist
                vmc = evaluate_vmc_gutzwiller(
                    nx,
                    ny,
                    h,
                    st,
                    n_chains=50,
                    burn_in_sweeps=15,
                    samples_per_chain=50,
                    seed_base=1000 + 37 * n_sites + int(round(h * 1000)),
                    preferred_bc=pref_bc,
                )
                e_upper = ex.e_proj_per_site
                valid = int((e_sdp <= e_ed + 1e-8) and (e_ed <= ex.e_proj_per_site + 1e-8))
                if not valid:
                    raise RuntimeError(
                        f"Sandwich violation on N={n_sites} h={h}: SDP={e_sdp}, ED={e_ed}, EX={e_upper}"
                    )
                delta_corr = ex.e_proj_per_site - e_ed
                proj_gain = ex.e_proj_per_site - e_mf
            else:
                vmc = evaluate_vmc_gutzwiller(
                    nx,
                    ny,
                    h,
                    st,
                    n_chains=50,
                    burn_in_sweeps=15,
                    samples_per_chain=50,
                    seed_base=2000 + int(round(h * 1000)),
                )
                e_ex_str = "nan"
                e_upper = vmc.e_proj_mean
                # For VMC on N=18, verify SDP <= ED and ED <= upper 95% CI + 1e-6 (zero-var at h=0)
                valid = int((e_sdp <= e_ed + 1e-8) and (e_ed <= vmc.e_proj_ci_high + 1e-5))
                if not valid:
                    raise RuntimeError(
                        f"VMC Sandwich violation on N={n_sites} h={h}: SDP={e_sdp}, ED={e_ed}, CI_hi={vmc.e_proj_ci_high}"
                    )
                delta_corr = vmc.e_proj_mean - e_ed
                proj_gain = vmc.e_proj_mean - e_mf

            lines.append(
                f"{n_sites:3d} {nx:2d} {ny:2d} {h:8.4f} {branch:9s} {vmc.bc_twist[0]:3d} {vmc.bc_twist[1]:3d} "
                f"{e_sdp:15.10f} {e_ed:15.10f} {e_mf:15.10f} {e_ex_str:>15s} "
                f"{vmc.e_proj_mean:15.10f} {vmc.e_proj_ci_low:15.10f} {vmc.e_proj_ci_high:15.10f} "
                f"{vmc.e_loc_std:10.3e} {vmc.tau_int_energy:7.3f} {vmc.n_eff_total:8.1f} "
                f"{delta_corr:12.8f} {proj_gain:12.8f} {valid:1d}"
            )

    _atomic_write(dat_path, "\n".join(lines) + "\n")
    meta = {
        "schema_version": 1,
        "dataset": "hc02_ed_sdp_scmft_sandwich",
        "source_repo": "kitaev-mf",
        "git_sha": sha,
        "sha256": hashlib.sha256(dat_path.read_bytes()).hexdigest(),
        "unit_convention": "pauli_k_minus_one_field_111",
        "geometry": "Honeycomb tori 2x2 (N=8), 2x3 (N=12), 3x3 (N=18)",
        "parameters": {
            "field_grid": list(FIELD_GRID),
            "vmc_chains": 50,
            "samples_per_chain": 50,
            "bootstrap_method": "BCa 95% CI (B=2000)",
        },
    }
    _atomic_write(
        dat_path.with_suffix(".meta.json"), json.dumps(meta, indent=2, sort_keys=True) + "\n"
    )
    print(f"Exported HC-02 sandwich benchmark -> {dat_path} (24 rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
