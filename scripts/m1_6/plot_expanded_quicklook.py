#!/usr/bin/env python3
"""Generate expanded Phase 2 diagnostic figures for M1.6.

Usage:
    python3 scripts/m1_6/plot_expanded_quicklook.py
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

BUILD_DIR = Path("build/rnd/m1_6")
FIGURES_DIR = BUILD_DIR / "figures"
LOG_SWEEP = BUILD_DIR / "runs_basin_sweep.jsonl"
LOG_MESH = BUILD_DIR / "runs_mesh_sensitivity.jsonl"


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    records = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                records.append(json.loads(line))
    return records


def plot_phase2_figures() -> None:
    """Generate diagnostic figures fig4, fig5, fig6, fig7."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    sweep_records = load_jsonl(LOG_SWEEP)
    mesh_records = load_jsonl(LOG_MESH)

    plt.style.use(
        "seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default"
    )
    plt.rcParams.update({"font.size": 11, "lines.linewidth": 2.0, "lines.markersize": 5})

    # --- Figure 4: Hysteresis and Phase Diagram ---
    if sweep_records:
        fwd_runs = [r for r in sweep_records if r.get("sweep") == "fwd_contin"]
        bwd_runs = [r for r in sweep_records if r.get("sweep") == "bwd_contin"]

        fwd_runs.sort(key=lambda r: float(r["h"]))
        bwd_runs.sort(key=lambda r: float(r["h"]))

        plt.figure(figsize=(12, 4))

        # Panel 1: Energy E_site
        plt.subplot(1, 3, 1)
        plt.plot(
            [r["h"] for r in fwd_runs],
            [r["e_site"] for r in fwd_runs],
            "o-",
            color="teal",
            label="Forward (0 -> 5)",
        )
        plt.plot(
            [r["h"] for r in bwd_runs],
            [r["e_site"] for r in bwd_runs],
            "s--",
            color="crimson",
            label="Backward (5 -> 0)",
        )
        plt.xlabel("Field h / |K|")
        plt.ylabel(r"Energy $E_{\mathrm{site}}$")
        plt.title("Ground-State Energy")
        plt.legend()

        # Panel 2: Magnetization m_111
        plt.subplot(1, 3, 2)
        plt.plot(
            [r["h"] for r in fwd_runs],
            [r["m_111"] for r in fwd_runs],
            "o-",
            color="teal",
            label="Forward",
        )
        plt.plot(
            [r["h"] for r in bwd_runs],
            [r["m_111"] for r in bwd_runs],
            "s--",
            color="crimson",
            label="Backward",
        )
        plt.xlabel("Field h / |K|")
        plt.ylabel(r"Magnetization $m_{[111]}$")
        plt.title(r"Magnetization along $[111]$")
        plt.legend()

        # Panel 3: Wilson loop W_p
        plt.subplot(1, 3, 3)
        plt.plot(
            [r["h"] for r in fwd_runs],
            [r["w_p"] for r in fwd_runs],
            "o-",
            color="teal",
            label="Forward",
        )
        plt.plot(
            [r["h"] for r in bwd_runs],
            [r["w_p"] for r in bwd_runs],
            "s--",
            color="crimson",
            label="Backward",
        )
        plt.xlabel("Field h / |K|")
        plt.ylabel(r"Plaquette Flux $W_p$")
        plt.title(r"Topological Flux $W_p$ (PH-02)")
        plt.legend()

        plt.tight_layout()
        fig4_path = FIGURES_DIR / "fig4_hysteresis_phase_diagram.png"
        plt.savefig(fig4_path, dpi=200)
        plt.close()
        print(f"Generated -> {fig4_path}")

        # --- Figure 5: Cross-Phase Seed Entrapment ---
        plt.figure(figsize=(10, 4.5))
        plt.subplot(1, 2, 1)
        for seed_tag, color, marker in [
            ("qsl_exact", "teal", "o"),
            ("pol_exact", "crimson", "s"),
            ("pol_pert_large", "darkorange", "^"),
        ]:
            runs = [
                r for r in sweep_records if r.get("sweep") == "static" and r.get("seed") == seed_tag
            ]
            runs.sort(key=lambda r: float(r["h"]))
            plt.plot(
                [r["h"] for r in runs if r["h"] <= 1.0],
                [r["e_site"] for r in runs if r["h"] <= 1.0],
                color=color,
                marker=marker,
                label=seed_tag,
            )
        plt.xlabel("Field h / |K|")
        plt.ylabel(r"$E_{\mathrm{site}}$")
        plt.title("Seed Energy Comparison (h <= 1.0)")
        plt.legend()

        plt.subplot(1, 2, 2)
        for seed_tag, color, marker in [
            ("qsl_exact", "teal", "o"),
            ("pol_exact", "crimson", "s"),
            ("pol_pert_large", "darkorange", "^"),
        ]:
            runs = [
                r for r in sweep_records if r.get("sweep") == "static" and r.get("seed") == seed_tag
            ]
            runs.sort(key=lambda r: float(r["h"]))
            plt.plot(
                [r["h"] for r in runs if r["h"] <= 1.0],
                [r["w_p"] for r in runs if r["h"] <= 1.0],
                color=color,
                marker=marker,
                label=seed_tag,
            )
        plt.xlabel("Field h / |K|")
        plt.ylabel(r"Plaquette Flux $W_p$")
        plt.title("Wilson Loop Entrapment (PH-04)")
        plt.legend()

        plt.tight_layout()
        fig5_path = FIGURES_DIR / "fig5_cross_phase_entrapment.png"
        plt.savefig(fig5_path, dpi=200)
        plt.close()
        print(f"Generated -> {fig5_path}")

        # --- Figure 6: High-Field Classical Asymptotic Limit (PH-05) ---
        plt.figure(figsize=(8, 5))
        bwd_high = [r for r in bwd_runs if r["h"] >= 1.0]
        if bwd_high:
            h_arr = np.array([r["h"] for r in bwd_high])
            e_meas = np.array([r["e_site"] for r in bwd_high])
            # Exact classical limit for unit-vector spins: E_classical = -h - 0.5
            e_class = -h_arr - 0.5
            diff_class = np.maximum(np.abs(e_meas - e_class), 1e-16)

            plt.semilogy(
                h_arr, diff_class, "o-", color="purple", label=r"$|E_{\mathrm{site}} - (-h - 0.5)|$"
            )
            # Quadratic reference ~ 1/h^2
            ref_quad = diff_class[0] * (h_arr[0] / h_arr) ** 2
            plt.loglog(
                h_arr,
                ref_quad,
                "--",
                color="gray",
                label=r"Asymptotic $\mathcal{O}(h^{-2})$ Reference",
            )

            plt.xlabel("Magnetic Field h")
            plt.ylabel(r"Departure from Classical Limit $|E - E_{\mathrm{class}}|$")
            plt.title(r"High-Field Asymptotic Saturation (PH-05)")
            plt.legend()
            plt.tight_layout()
            fig6_path = FIGURES_DIR / "fig6_high_field_classical_limit.png"
            plt.savefig(fig6_path, dpi=200)
            plt.close()
            print(f"Generated -> {fig6_path}")

    # --- Figure 7: Mesh Sensitivity & Timing (AH-03, AH-04) ---
    if mesh_records:
        plt.figure(figsize=(10, 4.5))

        # Subplot 1: ms per map evaluation vs Nk^2
        plt.subplot(1, 2, 1)
        diis_mesh = [r for r in mesh_records if r["solver"] == "diis" and r["h"] == 0.05]
        diis_mesh.sort(key=lambda r: int(r["nk"]))
        nk_vals = [r["nk"] for r in diis_mesh]
        ms_vals = [r["ms_per_phi"] for r in diis_mesh]

        plt.plot(nk_vals, ms_vals, "o-", color="darkgreen")
        plt.xlabel(r"Mesh Dimension $N_k$")
        plt.ylabel("Wall Time per Map Evaluation (ms)")
        plt.title(r"Cost Scaling $t_\Phi(N_k)$ on CPU")

        # Subplot 2: Energy vs Nk
        plt.subplot(1, 2, 2)
        plt.plot(nk_vals, [r["e_site"] for r in diis_mesh], "s-", color="navy")
        plt.xlabel(r"Mesh Dimension $N_k$")
        plt.ylabel(r"Ground State Energy $E_{\mathrm{site}}$")
        plt.title(r"Finite-Size Energy Convergence (AH-04)")

        plt.tight_layout()
        fig7_path = FIGURES_DIR / "fig7_mesh_scaling_and_timing.png"
        plt.savefig(fig7_path, dpi=200)
        plt.close()
        print(f"Generated -> {fig7_path}")


if __name__ == "__main__":
    plot_phase2_figures()
