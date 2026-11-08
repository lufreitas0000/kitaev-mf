#!/usr/bin/env python3
"""Benchmark script for Task T2.2.1: Vortex Separation Potential V(R) & Deconfinement.

Evaluates inter-vortex potential V(R) across Euclidean separations on honeycomb tori,
verifies gauge flux pinning (W_p = -1), and conducts Bayesian model comparison
falsifying linear confinement in favor of topological vison deconfinement.
"""

import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.phase_2.vortex_potential import (
    DeconfinementBayesianFitter,
    VortexSeparationPotentialCalculator,
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
    output_dir = Path("data/processed/benchmarks/m2_2")
    output_dir.mkdir(parents=True, exist_ok=True)
    out_file = output_dir / "vortex_potential_benchmark.dat"

    print("=" * 78)
    print("Task T2.2.1: Vortex Separation Potential V(R) & Topological Deconfinement")
    print("=" * 78)

    clusters = [5, 7, 8]
    fields = [0.0, 0.05]
    all_rows = []

    fitter = DeconfinementBayesianFitter()
    summary_results = []

    for L in clusters:
        print(f"\n--- Cluster L={L} (N={2 * L**2} sites, {L**2} plaquettes) ---")
        calc = VortexSeparationPotentialCalculator(L=L)

        for h_mag in fields:
            hx = hy = hz = h_mag / np.sqrt(3.0)
            params = ModelParameters(K=-1.0, h=(hx, hy, hz))

            pts = calc.scan_separation_profile(params)
            r_arr = np.array([pt.R for pt in pts])
            v_arr = np.array([pt.V_R for pt in pts])
            cuts = np.array([pt.cut_length for pt in pts])

            res = fitter.fit_and_compare(r_arr, v_arr, cut_lengths=cuts)
            summary_results.append((L, h_mag, res))

            print(
                f"  h={h_mag:.2f}: n={res.n_samples}, V_inf={res.v_inf:.4f}, "
                f"Delta_v={res.vison_gap:.4f} (err: {res.vison_gap_relative_error * 100:.1f}%), "
                f"Delta_BIC={res.delta_bic:.2f}, BF={res.bayes_factor:.2e}, "
                f"Deconfined={res.favors_deconfinement}"
            )

            for pt in pts:
                all_rows.append(
                    f"{L:4d}  {h_mag:6.3f}  {pt.p1[0]:3d}  {pt.p1[1]:3d}  "
                    f"{pt.p2[0]:3d}  {pt.p2[1]:3d}  {pt.cut_length:4d}  "
                    f"{pt.R:8.4f}  {pt.E_flux_free:12.6f}  {pt.E_2vortex:12.6f}  "
                    f"{pt.V_R:10.6f}  {pt.w_p1:6.2f}  {pt.w_p2:6.2f}"
                )

    git_hash = get_git_commit()
    timestamp = datetime.now(timezone.utc).isoformat()

    with open(out_file, "w", encoding="utf-8") as f:
        f.write("# Task T2.2.1: Vortex Separation Potential V(R) & Topological Deconfinement\n")
        f.write(f"# Git Commit: {git_hash}\n")
        f.write(f"# Generated: {timestamp}\n")
        f.write("# Model: Pure Ferromagnetic Kitaev K=-1.0 under [111] Magnetic Field\n")
        f.write("# Theoretical Exact Vison Gap: Delta_v = 0.1536 |K|, 2*Delta_v = 0.3072\n")
        f.write("# Bayesian Model Selection Summary:\n")
        for L, h_mag, res in summary_results:
            f.write(
                f"#   L={L}, h={h_mag:.2f}: RSS_null={res.rss_null_confinement:.6f}, "
                f"RSS_cand={res.rss_candidate_deconfinement:.6f}, Delta_BIC={res.delta_bic:.2f}, "
                f"BayesFactor={res.bayes_factor:.2e}, VisonGap={res.vison_gap:.4f}, "
                f"FavorsDeconfinement={res.favors_deconfinement}\n"
            )
        f.write("#\n")
        f.write(
            "# L      h    p1_1  p1_2  p2_1  p2_2  ncut         R           E_0          E_2v         V_R    w_p1    w_p2\n"
        )
        for row in all_rows:
            f.write(row + "\n")

    print(f"\n[OK] Benchmark completed successfully. Saved to {out_file}")


if __name__ == "__main__":
    main()
