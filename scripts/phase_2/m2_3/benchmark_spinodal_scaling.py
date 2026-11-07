#!/usr/bin/env python3
"""Benchmark script for Task T2.3.2: Critical Exponents & Scaling at the Polarized Breakdown Transition.

Traces high-resolution forward continuation approaching the spinodal breakdown h_c2+ = 0.400,
computes differential susceptibility chi_111(h) = dm/dh, and performs Bayesian model selection
comparing universal mean-field spinodal divergence (gamma_sp = 0.50) against regular crossover.
"""

import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from kitaev_mf.services.phase_2.spinodal_scaling import (
    SpinodalBayesianFitter,
    SpinodalScalingCalculator,
)


def get_git_commit() -> str:
    """Retrieve the current short git commit hash."""
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
        )
        return out.decode("ascii").strip()
    except Exception:
        return "unknown"


def main() -> None:
    output_dir = Path("data/processed/benchmarks/m2_3")
    output_dir.mkdir(parents=True, exist_ok=True)
    out_file = output_dir / "spinodal_scaling_benchmark.dat"

    print("=" * 78)
    print("Task T2.3.2: Critical Exponents & Scaling at Polarized Breakdown Transition")
    print("=" * 78)

    calc = SpinodalScalingCalculator(k_grid=(24, 24), h_spinodal=0.400)
    fitter = SpinodalBayesianFitter()

    h_fine = np.array([
        0.350, 0.360, 0.370, 0.375, 0.380,
        0.385, 0.390, 0.392, 0.394, 0.396, 0.398,
    ])

    print(f"\n[1/3] Tracing forward continuation across {len(h_fine)} field points...")
    spinodal_points = calc.run_forward_continuation(h_fine)

    for pt in spinodal_points:
        print(f"  h={pt.h_mag:.3f}: m_111={pt.m_111:.5f}, iters={pt.iterations:2d}, conv={pt.converged}")

    print("\n[2/3] Evaluating differential susceptibility series chi_111 = dm/dh...")
    susc_series = calc.compute_susceptibility_series(spinodal_points)

    all_rows = []
    for sp in susc_series:
        print(f"  h_mid={sp.h_mid:.4f} (t={sp.distance_to_spinodal:.4f}): chi_111={sp.chi_111:.4f}")
        all_rows.append(f"{sp.h_mid:8.4f}  {sp.distance_to_spinodal:8.4f}  {sp.chi_111:10.4f}")

    print("\n[3/3] Performing Bayesian model selection (gamma=0.50 vs constant crossover)...")
    bayes_res = fitter.fit_and_compare(susc_series)

    print(
        f"  Fitted gamma: {bayes_res.fitted_gamma:.4f} "
        f"(rel error vs 0.50: {bayes_res.gamma_relative_error*100:.2f}%)"
    )
    print(
        f"  n={bayes_res.n_samples}, RSS0={bayes_res.rss_null:.4e}, RSS1={bayes_res.rss_candidate:.4e}, "
        f"Delta_BIC={bayes_res.delta_bic:.2f}, BF={bayes_res.bayes_factor:.2e}, "
        f"FavorsSpinodal={bayes_res.favors_spinodal_divergence}"
    )

    git_hash = get_git_commit()
    timestamp = datetime.now(timezone.utc).isoformat()
    tmp_file = out_file.with_suffix(".tmp")

    with open(tmp_file, "w", encoding="utf-8") as f:
        f.write("# Task T2.3.2: Critical Exponents & Scaling at Polarized Breakdown Transition\n")
        f.write(f"# Git Commit: {git_hash}\n")
        f.write(f"# Generated: {timestamp}\n")
        f.write("# Model: Pure Ferromagnetic Kitaev K=-1.0 under [111] Magnetic Field\n")
        f.write(
            f"# Bayesian Model Selection Summary: fitted_gamma={bayes_res.fitted_gamma:.4f}, "
            f"RSS_null={bayes_res.rss_null:.4e}, RSS_cand={bayes_res.rss_candidate:.4e}, "
            f"Delta_BIC={bayes_res.delta_bic:.2f}, BayesFactor={bayes_res.bayes_factor:.2e}, "
            f"FavorsSpinodal={bayes_res.favors_spinodal_divergence}\n"
        )
        f.write("#\n")
        f.write("# h_mid     t=(h_sp-h)    chi_111\n")
        for row in all_rows:
            f.write(f"{row}\n")

    tmp_file.replace(out_file)
    print(f"\n[OK] Benchmark dataset successfully saved to {out_file}")


if __name__ == "__main__":
    main()
