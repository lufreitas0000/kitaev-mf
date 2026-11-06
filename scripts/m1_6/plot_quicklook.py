#!/usr/bin/env python3
"""Generate quick-look exploratory diagnostic figures for M1.6 weak-field results.

Usage:
    python3 scripts/m1_6/plot_quicklook.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BUILD_DIR = Path("build/rnd/m1_6")
FIGURES_DIR = BUILD_DIR / "figures"
LOG_PATH = BUILD_DIR / "runs_weak_field.jsonl"


def plot_quicklook_figures() -> None:
    """Generate diagnostic PNG figures from runs_weak_field.jsonl."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    if not LOG_PATH.is_file():
        print(f"Log file not found: {LOG_PATH}")
        return

    records = []
    with LOG_PATH.open("r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                records.append(json.loads(line))

    if not records:
        print("No records found in log.")
        return

    # Style configuration
    plt.style.use(
        "seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default"
    )
    plt.rcParams.update({"font.size": 11, "lines.linewidth": 2.0, "lines.markersize": 6})

    # Figure 1: Convergence trajectory comparison on QSL seed at h = 0.05
    plt.figure(figsize=(8, 5))
    qsl_runs = [r for r in records if r["h"] == 0.05 and r["seed"] == "qsl" and r["delta"] == 0.01]
    for r in qsl_runs:
        solver_name = r["solver"].upper()
        hist = r.get("residual_history", [])
        if hist:
            plt.semilogy(
                range(1, len(hist) + 1),
                hist,
                marker="o",
                label=f"{solver_name} (n_phi={r['n_phi']})",
            )

    plt.xlabel("Outer Iteration")
    plt.ylabel(r"Residual Norm $\|R_k\|_\infty$")
    plt.title(r"Solver Convergence on Perturbed QSL Seed ($\delta=0.01, h=0.05$)")
    plt.axhline(1e-8, color="black", linestyle="--", alpha=0.6, label="Tolerance (1e-8)")
    plt.legend()
    plt.tight_layout()
    fig1_path = FIGURES_DIR / "fig1_solver_convergence_efficiency.png"
    plt.savefig(fig1_path, dpi=200)
    plt.close()
    print(f"Generated -> {fig1_path}")

    # Figure 2: Basin energy and flux separation (QSL vs Polarized)
    plt.figure(figsize=(9, 4.5))
    plt.subplot(1, 2, 1)
    for seed_name, color, marker in [("qsl", "teal", "o"), ("polarized", "crimson", "s")]:
        runs = [
            r
            for r in records
            if r["seed"] == seed_name and r["solver"] == "diis" and r["delta"] == 0.0
        ]
        h_vals = [r["h"] for r in runs]
        e_vals = [r["e_site"] for r in runs]
        plt.plot(h_vals, e_vals, color=color, marker=marker, label=f"{seed_name.upper()} basin")
    plt.xlabel("Magnetic Field h / |K|")
    plt.ylabel(r"Ground State Energy $E_{\mathrm{site}}$")
    plt.title("Energy Separation")
    plt.legend()

    plt.subplot(1, 2, 2)
    for seed_name, color, marker in [("qsl", "teal", "o"), ("polarized", "crimson", "s")]:
        runs = [
            r
            for r in records
            if r["seed"] == seed_name and r["solver"] == "diis" and r["delta"] == 0.0
        ]
        h_vals = [r["h"] for r in runs]
        wp_vals = [r["w_p"] for r in runs]
        plt.plot(h_vals, wp_vals, color=color, marker=marker, label=f"{seed_name.upper()} basin")
    plt.xlabel("Magnetic Field h / |K|")
    plt.ylabel(r"Plaquette Flux $W_p$")
    plt.title(r"Plaquette Wilson Loop $W_p$")
    plt.legend()
    plt.tight_layout()

    fig2_path = FIGURES_DIR / "fig2_weak_field_basin_separation.png"
    plt.savefig(fig2_path, dpi=200)
    plt.close()
    print(f"Generated -> {fig2_path}")

    # Figure 3: Multiplier scaling norm ||lambda|| vs h (Testing PH-01)
    plt.figure(figsize=(7, 4.5))
    qsl_diis = [
        r for r in records if r["seed"] == "qsl" and r["solver"] == "diis" and r["delta"] == 0.0
    ]
    h_arr = np.array([r["h"] for r in qsl_diis], dtype=np.float64)
    lam_arr = np.array([r["norm_lam"] for r in qsl_diis], dtype=np.float64)

    if len(h_arr) >= 2:
        plt.loglog(
            h_arr, lam_arr, "o-", color="darkblue", label=r"Measured $\|\boldsymbol{\lambda}^*\|$"
        )
        # Ideal quadratic scaling reference
        quad_ref = lam_arr[0] * (h_arr / h_arr[0]) ** 2
        plt.loglog(h_arr, quad_ref, "--", color="gray", label=r"Quadratic Reference $\sim h^2$")
        plt.xlabel("Magnetic Field h")
        plt.ylabel(r"Gauge Multiplier Norm $\|\boldsymbol{\lambda}^*\|$")
        plt.title(r"Weak-Field Gauge Multiplier Scaling (PH-01)")
        plt.legend()
        plt.tight_layout()
        fig3_path = FIGURES_DIR / "fig3_multiplier_scaling.png"
        plt.savefig(fig3_path, dpi=200)
        plt.close()
        print(f"Generated -> {fig3_path}")


if __name__ == "__main__":
    plot_quicklook_figures()
