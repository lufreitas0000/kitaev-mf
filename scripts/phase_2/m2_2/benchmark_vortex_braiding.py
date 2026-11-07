#!/usr/bin/env python3
"""Benchmark script for Task T2.2.3: Vortex Braiding Phases & Non-Abelian Holonomy.

Evaluates geometric Berry phases along closed dual vortex loops, verifies topological
pi statistical holonomy, and performs Bayesian model selection against trivial holonomy.
"""

from datetime import datetime, timezone
from pathlib import Path
import subprocess
import numpy as np

from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.phase_2.vortex_braiding import (
    BraidingBayesianFitter,
    VortexBraidingCalculator,
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
    out_file = output_dir / "vortex_braiding_benchmark.dat"

    print("=" * 78)
    print("Task T2.2.3: Vortex Braiding Phases & Non-Abelian Statistical Holonomy")
    print("=" * 78)

    # Focus benchmarks on resonant clusters L=0 mod 3 (L=6) where the momentum grid
    # samples the Dirac points and isolated in-gap Majorana zero modes exist,
    # and scan across the low-to-intermediate topological field window.
    clusters = [6]
    fields = [0.02, 0.03, 0.05, 0.10]
    all_rows = []
    delta_gammas = []

    fitter = BraidingBayesianFitter()

    for h_mag in fields:
        hx = hy = hz = h_mag / np.sqrt(3.0)
        params = ModelParameters(K=-1.0, h=(hx, hy, hz))
        print(f"\n--- Field h={h_mag:.2f} parallel to [111] ---")

        for L in clusters:
            calc = VortexBraidingCalculator(L=L)
            res = calc.evaluate_braiding_holonomy(params, N_sub=4)
            delta_gammas.append(res.delta_gamma)

            print(
                f"  L={L:2d} (N={2*L**2:3d}): gamma_enc={res.gamma_enclosing:.4f} rad, "
                f"gamma_non={res.gamma_non_enclosing:.4f} rad, Delta_gamma={res.delta_gamma:.4f} rad "
                f"({res.delta_gamma_over_pi:.4f} pi, err: {res.relative_error_from_pi*100:.1f}%)"
            )

            all_rows.append(
                f"{L:4d}  {h_mag:6.3f}  {res.gamma_enclosing:10.4f}  "
                f"{res.gamma_non_enclosing:10.4f}  {res.delta_gamma:10.4f}  "
                f"{res.delta_gamma_over_pi:10.4f}  {res.relative_error_from_pi:10.4f}"
            )

    bayes_res = fitter.fit_and_compare(np.array(delta_gammas))
    print("\n--- Bayesian Model Selection Summary (pi Holonomy vs Trivial 0) ---")
    print(
        f"n={bayes_res.n_samples}, RSS0={bayes_res.rss_null:.4f}, RSS1={bayes_res.rss_candidate:.4f}, "
        f"Delta_BIC={bayes_res.delta_bic:.2f}, BF={bayes_res.bayes_factor:.2e}, "
        f"FavorsPi={bayes_res.favors_pi_holonomy}"
    )

    git_hash = get_git_commit()
    timestamp = datetime.now(timezone.utc).isoformat()

    with open(out_file, "w", encoding="utf-8") as f:
        f.write("# Task T2.2.3: Vortex Braiding Phases & Non-Abelian Statistical Holonomy\n")
        f.write(f"# Git Commit: {git_hash}\n")
        f.write(f"# Generated: {timestamp}\n")
        f.write("# Model: Pure Ferromagnetic Kitaev K=-1.0 under [111] Magnetic Field\n")
        f.write(
            f"# Bayesian Model Selection Summary: RSS_null={bayes_res.rss_null:.4f}, "
            f"RSS_cand={bayes_res.rss_candidate:.4f}, Delta_BIC={bayes_res.delta_bic:.2f}, "
            f"BayesFactor={bayes_res.bayes_factor:.2e}, FavorsPi={bayes_res.favors_pi_holonomy}\n"
        )
        f.write("#\n")
        f.write(
            "# L      h       gamma_enc   gamma_non   Delta_gamma  Delta_g/pi   rel_error\n"
        )
        for row in all_rows:
            f.write(row + "\n")

    print(f"\n[OK] Benchmark completed successfully. Saved to {out_file}")


if __name__ == "__main__":
    main()
