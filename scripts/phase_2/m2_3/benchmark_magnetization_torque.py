#!/usr/bin/env python3
"""Benchmark script for Task T2.3.1: Global Field-Directional Magnetization Surface & Torque Discontinuities.

Evaluates field-directional magnetization surfaces, magnetic torque profiles,
and computes Bayesian evidence comparing anisotropic Kitaev torque against isotropic null.
"""

import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from kitaev_mf.services.phase_2.magnetization_torque import (
    MagnetizationTorqueCalculator,
    TorqueBayesianFitter,
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
    out_file = output_dir / "magnetization_torque_benchmark.dat"

    print("=" * 78)
    print("Task T2.3.1: Global Field-Directional Magnetization Surface & Torque Discontinuities")
    print("=" * 78)

    calc = MagnetizationTorqueCalculator(k_grid=(16, 16))
    fitter = TorqueBayesianFitter()

    thetas_deg = np.linspace(0.0, 90.0, 13)
    # Fields spanning weak, intermediate, and near-critical regimes
    fields = [0.10, 0.20, 0.35]

    all_rows = []
    all_thetas_rad = []
    all_taus_proj = []

    for h_mag in fields:
        print(f"\n--- Field Magnitude h={h_mag:.2f} | Angular Sweep [001] -> [111] -> [110] ---")
        points = calc.sweep_angles(h_mag, thetas_deg)

        for pt in points:
            all_thetas_rad.append(pt.theta_rad)
            all_taus_proj.append(pt.tau_projection)

            print(
                f"  theta={pt.theta_deg:5.1f} deg: ||m||={pt.m_magnitude:6.4f}, "
                f"tau_proj={pt.tau_projection:+10.6f}, tau_norm={pt.tau_norm:9.6f}, "
                f"E_site={pt.energy_per_site:9.5f}"
            )

            all_rows.append(
                f"{pt.h_mag:6.3f}  {pt.theta_deg:6.2f}  {pt.theta_rad:8.5f}  "
                f"{pt.m_magnitude:8.5f}  {pt.tau_projection:+12.6e}  {pt.tau_norm:12.6e}  "
                f"{pt.energy_per_site:10.5f}"
            )

    thetas_arr = np.array(all_thetas_rad)
    taus_arr = np.array(all_taus_proj)
    h_arr = np.array([float(r.split()[0]) for r in all_rows])
    bayes_res = fitter.fit_and_compare(thetas_arr, taus_arr, h_magnitudes=h_arr)

    print("\n--- Bayesian Model Selection Summary (Anisotropic vs Isotropic) ---")
    print(
        f"n={bayes_res.n_samples}, RSS0={bayes_res.rss_null:.4e}, RSS1={bayes_res.rss_candidate:.4e}, "
        f"Delta_BIC={bayes_res.delta_bic:.2f}, BF={bayes_res.bayes_factor:.2e}, "
        f"FavorsAnisotropic={bayes_res.favors_anisotropic_torque}"
    )

    git_hash = get_git_commit()
    timestamp = datetime.now(timezone.utc).isoformat()
    tmp_file = out_file.with_suffix(".tmp")

    with open(tmp_file, "w", encoding="utf-8") as f:
        f.write(
            "# Task T2.3.1: Global Field-Directional Magnetization Surface & Torque Discontinuities\n"
        )
        f.write(f"# Git Commit: {git_hash}\n")
        f.write(f"# Generated: {timestamp}\n")
        f.write("# Model: Pure Ferromagnetic Kitaev K=-1.0 under Directional Magnetic Field\n")
        f.write(
            f"# Bayesian Model Selection Summary: RSS_null={bayes_res.rss_null:.4e}, "
            f"RSS_cand={bayes_res.rss_candidate:.4e}, Delta_BIC={bayes_res.delta_bic:.2f}, "
            f"BayesFactor={bayes_res.bayes_factor:.2e}, FavorsAnisotropic={bayes_res.favors_anisotropic_torque}\n"
        )
        f.write("#\n")
        f.write(
            "# h_mag   theta_deg  theta_rad     ||m||       tau_proj       tau_norm       E_site\n"
        )
        for row in all_rows:
            f.write(f"{row}\n")

    tmp_file.replace(out_file)
    print(f"\n[OK] Benchmark dataset successfully saved to {out_file}")


if __name__ == "__main__":
    main()
