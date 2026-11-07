#!/usr/bin/env python3
"""Benchmark script for Task T2.2.2: Majorana Zero Modes Bound to Isolated Vortices.

Evaluates in-gap spectral isolation, core localization enrichment, and Bayesian model
selection of exponential finite-size MZM decoupling on honeycomb tori.
"""

from datetime import datetime, timezone
from pathlib import Path
import subprocess
import numpy as np

from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.phase_2.mzm_localization import (
    MZMLocalizationCalculator,
    MZMSplittingBayesianFitter,
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
    out_file = output_dir / "mzm_localization_benchmark.dat"

    print("=" * 78)
    print("Task T2.2.2: Majorana Zero Modes Localization & Hybridization Decoupling")
    print("=" * 78)

    clusters = [4, 5, 7, 8]
    fields = [0.02, 0.05, 0.10]
    all_rows = []

    fitter = MZMSplittingBayesianFitter()
    bayesian_summaries = []

    for h_mag in fields:
        hx = hy = hz = h_mag / np.sqrt(3.0)
        params = ModelParameters(K=-1.0, h=(hx, hy, hz))
        print(f"\n--- Field h={h_mag:.2f} parallel to [111] ---")

        l_list = []
        eps_list = []

        for L in clusters:
            calc = MZMLocalizationCalculator(L=L)
            metrics = calc.evaluate_mzm_pair(params)

            l_list.append(float(L))
            eps_list.append(metrics.eps_mzm)

            p_d0 = metrics.hop_profile.get(0, 0.0)
            p_d1 = metrics.hop_profile.get(1, 0.0)
            p_d2 = metrics.hop_profile.get(2, 0.0)

            print(
                f"  L={L:2d} (N={2*L**2:3d}): eps_MZM={metrics.eps_mzm:.6f}, "
                f"eps_bulk={metrics.eps_bulk:.6f}, in-gap ratio={metrics.in_gap_ratio:.4f}, "
                f"W_core={metrics.w_core:.4f}, ratio={metrics.enrichment_ratio:.2f}x"
            )

            all_rows.append(
                f"{L:4d}  {h_mag:6.3f}  {metrics.R:8.4f}  {metrics.eps_mzm:10.6f}  "
                f"{metrics.eps_bulk:10.6f}  {metrics.in_gap_ratio:10.4f}  "
                f"{metrics.w_core:8.4f}  {metrics.w_uniform:8.4f}  {metrics.enrichment_ratio:8.2f}  "
                f"{p_d0:10.6f}  {p_d1:10.6f}  {p_d2:10.6f}"
            )

        # Bayesian fit for this field
        res = fitter.fit_and_compare(np.array(l_list), np.array(eps_list))
        bayesian_summaries.append((h_mag, res))
        print(
            f"  >> Bayesian Fit: Delta_BIC={res.delta_bic:.2f}, BF={res.bayes_factor:.2e}, "
            f"xi={res.decay_length_xi:.2f} unit cells, FavorsDecay={res.favors_exponential_decay}"
        )

    git_hash = get_git_commit()
    timestamp = datetime.now(timezone.utc).isoformat()

    with open(out_file, "w", encoding="utf-8") as f:
        f.write("# Task T2.2.2: Majorana Zero Modes Localization & Hybridization Decoupling\n")
        f.write(f"# Git Commit: {git_hash}\n")
        f.write(f"# Generated: {timestamp}\n")
        f.write("# Model: Pure Ferromagnetic Kitaev K=-1.0 under [111] Magnetic Field\n")
        f.write("# Bayesian Model Selection Summary (Exponential Decoupling vs Constant):\n")
        for h_mag, res in bayesian_summaries:
            f.write(
                f"#   h={h_mag:.2f}: RSS_null={res.rss_null_constant:.6f}, "
                f"RSS_cand={res.rss_candidate_exponential:.6f}, Delta_BIC={res.delta_bic:.2f}, "
                f"BayesFactor={res.bayes_factor:.2e}, xi={res.decay_length_xi:.2f}, "
                f"FavorsDecay={res.favors_exponential_decay}\n"
            )
        f.write("#\n")
        f.write(
            "# L      h         R    eps_MZM   eps_bulk   in_gap_rat   W_core   W_unif   ratio     P_d0        P_d1        P_d2\n"
        )
        for row in all_rows:
            f.write(row + "\n")

    print(f"\n[OK] Benchmark completed successfully. Saved to {out_file}")


if __name__ == "__main__":
    main()
