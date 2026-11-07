#!/usr/bin/env python3
"""Benchmark script for Task T2.3.3: Global Phase Diagram Mapping & First-Order Boundary Anisotropy.

Sweeps the first-order thermodynamic boundary h_c*(theta) across 13 orientations in the (1, 1, 0) plane
connecting [001] -> [111] -> [110], and computes Bayesian evidence comparing anisotropic modulation
against isotropic phase boundary hypotheses.
"""

import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from kitaev_mf.services.phase_2.phase_diagram import (
    PhaseBoundaryBayesianFitter,
    PhaseBoundaryCalculator,
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
    out_file = output_dir / "phase_boundary_benchmark.dat"

    print("=" * 78)
    print("Task T2.3.3: Global Phase Diagram Mapping & First-Order Boundary Anisotropy")
    print("=" * 78)

    calc = PhaseBoundaryCalculator(k_grid=(16, 16))
    fitter = PhaseBoundaryBayesianFitter()

    thetas_deg = np.linspace(0.0, 90.0, 13)
    print(f"\n[1/2] Sweeping first-order critical field h_c*(theta) across {len(thetas_deg)} angles...")
    points = calc.sweep_boundary(thetas_deg)

    all_rows = []
    h_c_list = []
    for pt in points:
        h_c_list.append(pt.h_c_star)
        print(
            f"  theta={pt.theta_deg:5.1f} deg (rad={pt.theta_rad:6.4f}): "
            f"h_c^*={pt.h_c_star:.5f}, E_cross={pt.energy_crossing:.5f}"
        )
        all_rows.append(
            f"{pt.theta_deg:6.2f}  {pt.theta_rad:8.5f}  {pt.h_c_star:10.5f}  {pt.energy_crossing:10.5f}  "
            f"{pt.qsl_energy_330:10.5f}  {pt.pol_energy_330:10.5f}  "
            f"{pt.qsl_energy_340:10.5f}  {pt.pol_energy_340:10.5f}"
        )

    h_c_arr = np.array(h_c_list)
    print("\n[2/2] Performing Bayesian model selection (Anisotropic vs Isotropic boundary)...")
    bayes_res = fitter.fit_and_compare(thetas_deg, h_c_arr)

    print(f"  h_base = {bayes_res.h_base:.5f}, Delta_h = {bayes_res.delta_h:+.5f}")
    print(f"  theta_max = {bayes_res.theta_max_deg:.2f} deg (expected 54.74 deg)")
    print(
        f"  n={bayes_res.n_samples}, RSS0={bayes_res.rss_null:.4e}, RSS1={bayes_res.rss_candidate:.4e}, "
        f"Delta_BIC={bayes_res.delta_bic:.2f}, BF={bayes_res.bayes_factor:.2e}, "
        f"FavorsAnisotropic={bayes_res.favors_anisotropic_boundary}"
    )

    git_hash = get_git_commit()
    timestamp = datetime.now(timezone.utc).isoformat()
    tmp_file = out_file.with_suffix(".tmp")

    with open(tmp_file, "w", encoding="utf-8") as f:
        f.write("# Task T2.3.3: Global Phase Diagram Mapping & First-Order Boundary Anisotropy\n")
        f.write(f"# Git Commit: {git_hash}\n")
        f.write(f"# Generated: {timestamp}\n")
        f.write("# Model: Pure Ferromagnetic Kitaev K=-1.0 in (1, 1, 0) Field Tilt Plane\n")
        f.write(
            f"# Bayesian Model Selection Summary: h_base={bayes_res.h_base:.5f}, "
            f"Delta_h={bayes_res.delta_h:+.5f}, theta_max={bayes_res.theta_max_deg:.2f}, "
            f"RSS_null={bayes_res.rss_null:.4e}, RSS_cand={bayes_res.rss_candidate:.4e}, "
            f"Delta_BIC={bayes_res.delta_bic:.2f}, BayesFactor={bayes_res.bayes_factor:.2e}, "
            f"FavorsAnisotropic={bayes_res.favors_anisotropic_boundary}\n"
        )
        f.write("#\n")
        f.write("# theta_deg  theta_rad     h_c_star     E_cross     E_qsl_330   E_pol_330   E_qsl_340   E_pol_340\n")
        for row in all_rows:
            f.write(f"{row}\n")

    tmp_file.replace(out_file)
    print(f"\n[OK] Benchmark dataset successfully saved to {out_file}")


if __name__ == "__main__":
    main()
