#!/usr/bin/env python3
"""Benchmark Script for Task T3.W1.2: MZM Thermodynamic Scaling and Sharpness.

Exports ACID-compliant benchmark table with metadata provenance:
- Table 1: data/processed/benchmarks/phase_3/mzm_hybridization_scaling.dat
- Table 2: data/processed/benchmarks/phase_3/mzm_thermodynamic_sharpness.dat
"""

import datetime
import os
import subprocess
from pathlib import Path

from kitaev_mf.services.phase_3.mzm_thermodynamic_scaling import (
    compute_thermodynamic_zero_bias_sharpness,
    measure_mzm_hybridization_scaling,
)


def get_git_revision() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()
    except Exception:
        return "unknown"


def main() -> None:
    git_rev = get_git_revision()
    now_iso = datetime.datetime.now(datetime.timezone.utc).isoformat()
    out_dir = Path(__file__).resolve().parents[2] / "data/processed/benchmarks/phase_3"
    os.makedirs(out_dir, exist_ok=True)

    # 1. Hybridization scaling across commensurate sizes L in {6, 12, 18, 24}
    print("Measuring MZM hybridization scaling...")
    data_hyb, fit_hyb = measure_mzm_hybridization_scaling(L_values=[6, 12, 18, 24])

    path_hyb = os.path.join(out_dir, "mzm_hybridization_scaling.dat")
    with open(path_hyb, "w") as f:
        f.write("# ACID DATA PROVENANCE HEADER\n")
        f.write("# Task: T3.W1.2 MZM Hybridization Thermodynamic Scaling\n")
        f.write(f"# Date: {now_iso}\n")
        f.write(f"# Git-Commit: {git_rev}\n")
        f.write(
            f"# Fit: xi_M = {fit_hyb.xi_M:.4f} +/- {fit_hyb.xi_M_stderr:.4f} a, R^2 = {fit_hyb.r_squared:.4f}\n"
        )
        f.write(
            f"# Intercept: alpha = {fit_hyb.intercept_alpha:.4f}, p_value = {fit_hyb.p_value:.4e}\n"
        )
        f.write("# Schema: L dy R_units_a delta_eps_0\n")
        for d in data_hyb:
            f.write(f"{d.L:4d} {d.dy:4d} {d.R_units_a:10.4f} {d.delta_eps_0:14.8e}\n")
    print(f"Exported {path_hyb}")

    # 2. Thermodynamic sharpness and subgap weight across L in [6, 8, 12, 16, 20]
    print("Computing zero-bias sharpness across L...")
    results_sharp = compute_thermodynamic_zero_bias_sharpness(L_values=[6, 8, 12, 16, 20])

    path_sharp = os.path.join(out_dir, "mzm_thermodynamic_sharpness.dat")
    with open(path_sharp, "w") as f:
        f.write("# ACID DATA PROVENANCE HEADER\n")
        f.write("# Task: T3.W1.2 Zero-Bias Sharpness and Thermodynamic Convergence\n")
        f.write(f"# Date: {now_iso}\n")
        f.write(f"# Git-Commit: {git_rev}\n")
        f.write("# Parameters: K=1.0, kappa=0.15, boundary=obc\n")
        f.write(
            "# Schema: L boundary omega_MZM tot_weight subgap_weight sharpness_ratio bulk_min_E bulk_subgap_weight\n"
        )
        for r in results_sharp:
            f.write(
                f"{r.L:4d} {r.boundary:>8s} {r.omega_MZM:10.6f} {r.total_spectral_weight:10.6f} "
                f"{r.subgap_weight:10.6f} {r.sharpness_ratio:10.6f} {r.bulk_min_energy:10.6f} "
                f"{r.bulk_subgap_weight:14.8e}\n"
            )
    print(f"Exported {path_sharp}")


if __name__ == "__main__":
    main()
