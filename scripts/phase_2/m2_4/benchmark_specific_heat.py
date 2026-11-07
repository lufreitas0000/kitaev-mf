#!/usr/bin/env python3
"""Benchmark script for Task T2.4.1: Two-Peak Specific Heat & Fractionalization Entropy Stages."""

import datetime
import subprocess
from pathlib import Path

import numpy as np

from kitaev_mf.services.phase_2.thermodynamics import (
    SpecificHeatCalculator,
    ThermodynamicBayesianFitter,
)


def get_git_revision() -> str:
    """Retrieve the current Git commit hash."""
    try:
        rev = (
            subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL)
            .decode("ascii")
            .strip()
        )
        return rev
    except Exception:
        return "unknown"


def main() -> None:
    print("=" * 78)
    print("Task T2.4.1: Two-Peak Specific Heat C_v(T) & Fractionalization Entropy Stages")
    print("=" * 78)

    calc = SpecificHeatCalculator(delta_v=0.137, k_grid=(32, 32))
    temps = np.logspace(-2.5, 0.8, 40)

    print(f"\n[1/2] Computing thermodynamic curve across {len(temps)} temperatures...")
    res = calc.compute_thermodynamic_curve(temps)

    print(f"  Low-T Vison Peak:      T_L = {res.t_low:.4f} |K|, C_v(T_L) = {res.c_v_low:.4f}")
    print(f"  High-T Majorana Peak:  T_H = {res.t_high:.4f} |K|, C_v(T_H) = {res.c_v_high:.4f}")
    print(f"  Peak Ratio:            T_H / T_L = {res.peak_ratio:.2f} (Target >= 5.0)")
    print(
        f"  Intermediate Dip:      T_mid = {res.t_mid:.4f} |K|, C_v(T_mid) = {res.c_v_mid:.4f}, "
        f"S(T_mid) = {res.s_mid:.5f} (Target 1/2 ln 2 = {0.5 * np.log(2):.5f})"
    )
    print(
        f"  High-T Entropy:        T_max = {res.t_max:.4f} |K|, S(T_max) = {res.s_max:.5f} "
        f"(Target ln 2 = {np.log(2):.5f})"
    )

    print("\n[2/2] Performing Bayesian model selection (Two-Peak vs Single-Peak)...")
    cv_values = np.array([p.c_v_total for p in res.points])
    fitter = ThermodynamicBayesianFitter()
    bayes_res = fitter.fit_and_compare(temps, cv_values)

    print(
        f"  n={bayes_res.n_samples}, RSS0={bayes_res.rss_null:.4e}, RSS1={bayes_res.rss_cand:.4e}, "
        f"Delta_BIC={bayes_res.delta_bic:.2f}, BF={bayes_res.bayes_factor:.2e}, "
        f"FavorsTwoPeak={bayes_res.favors_two_peak}"
    )

    # Save ACID dataset
    out_dir = Path("data/processed/benchmarks/m2_4")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "specific_heat_benchmark.dat"

    git_hash = get_git_revision()
    now_utc = datetime.datetime.now(datetime.timezone.utc).isoformat()

    header = (
        f"# Kitaev Model Phase 2 Benchmark: Task T2.4.1 Specific Heat & Entropy Stages\n"
        f"# Generated: {now_utc}\n"
        f"# Git Commit: {git_hash}\n"
        f"# Model: Pure Ferromagnetic Kitaev K=-1.0, zero field h=0.0\n"
        f"# Vison Gap Delta_v: {calc.delta_v:.4f} |K|\n"
        f"# Peak Analysis: T_L={res.t_low:.5f}, T_H={res.t_high:.5f}, Peak_Ratio={res.peak_ratio:.2f}\n"
        f"# Entropy Stages: S(T_mid={res.t_mid:.4f})={res.s_mid:.5f} (0.5 ln 2={0.5*np.log(2):.5f}), "
        f"S(T_max={res.t_max:.4f})={res.s_max:.5f} (ln 2={np.log(2):.5f})\n"
        f"# Bayesian Model Selection Summary: Delta_BIC={bayes_res.delta_bic:.2f}, BF={bayes_res.bayes_factor:.2e}, "
        f"FavorsTwoPeak={bayes_res.favors_two_peak}\n"
        f"# Columns: Temperature  C_v_vison  C_v_maj  C_v_total  S_vison  S_maj  S_total\n"
    )

    data_rows = []
    for p in res.points:
        data_rows.append(
            f"{p.temperature:12.6e} {p.c_v_vison:12.6e} {p.c_v_maj:12.6e} {p.c_v_total:12.6e} "
            f"{p.s_vison:12.6e} {p.s_maj:12.6e} {p.s_total:12.6e}\n"
        )

    with open(out_file, "w", encoding="utf-8") as f:
        f.write(header)
        f.writelines(data_rows)

    print(f"\n[OK] Benchmark dataset successfully saved to {out_file}\n")


if __name__ == "__main__":
    main()
