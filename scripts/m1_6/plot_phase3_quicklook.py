"""Generate Phase 3 quicklook diagnostic plots (Figures 8, 9, 10, 11).

Inputs:
- data/processed/benchmarks/m1_6/critical_zoom_111.dat
- data/processed/benchmarks/m1_6/majorana_gap_scaling.dat
- data/processed/benchmarks/m1_6/supercell_comparison.dat
Outputs (build/rnd/m1_6/figures/):
- fig8_critical_zoom_crossing_and_spinodal.png
- fig9_majorana_gap_power_law.png
- fig10_supercell_energy_parity.png
- fig11_contraction_factor_critical_slowing.png
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

BUILD_DIR = Path("build/rnd/m1_6")
FIGURES_DIR = BUILD_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def plot_phase3_figures() -> None:
    # -------------------------------------------------------------------------
    # Figure 8: Precision Critical Crossing and Forward Spinodal
    # -------------------------------------------------------------------------
    zoom_path = Path("data/processed/benchmarks/m1_6/critical_zoom_111.dat")
    if zoom_path.is_file():
        data_zoom: list[dict[str, Any]] = []
        with zoom_path.open("r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                data_zoom.append(
                    {
                        "h": float(parts[0]),
                        "branch": parts[1],
                        "conv": parts[2] == "True",
                        "n_phi": int(parts[3]),
                        "e_site": float(parts[5]),
                        "m_111": float(parts[6]),
                        "w_p": float(parts[7]),
                        "mu": float(parts[9]),
                    }
                )

        fwd_runs = [d for d in data_zoom if d["branch"] == "fwd_qsl"]
        bwd_runs = [d for d in data_zoom if d["branch"] == "bwd_pol"]

        fwd_runs.sort(key=lambda d: float(d["h"]))
        bwd_runs.sort(key=lambda d: float(d["h"]))

        plt.figure(figsize=(11, 4.5))

        # Subplot 1: Energy & Delta E
        plt.subplot(1, 2, 1)
        h_fwd = [d["h"] for d in fwd_runs]
        e_fwd = [d["e_site"] for d in fwd_runs]
        h_bwd = [d["h"] for d in bwd_runs]
        e_bwd = [d["e_site"] for d in bwd_runs]

        plt.plot(h_fwd, e_fwd, "-", color="teal", lw=2, label="QSL Branch (Forward)")
        plt.plot(h_bwd, e_bwd, "--", color="crimson", lw=2, label="Polarized Branch (Backward)")
        plt.axvline(
            0.3325, color="black", linestyle=":", lw=1.5, label=r"$h_c^* = 0.3325$ (Crossing)"
        )
        plt.axvline(
            0.4000, color="orange", linestyle="--", lw=1.5, label=r"$h_{c2}^+ = 0.400$ (Spinodal)"
        )

        plt.xlabel("Magnetic Field $h$")
        plt.ylabel(r"Energy per Site $E_{\mathrm{site}}$")
        plt.title("Precision Critical Energy Landscape (M1.6 Phase 3)")
        plt.legend(fontsize=9)

        # Subplot 2: Wilson loop flux W_p
        plt.subplot(1, 2, 2)
        w_fwd = [d["w_p"] for d in fwd_runs]
        w_bwd = [d["w_p"] for d in bwd_runs]
        plt.plot(h_fwd, w_fwd, "-", color="teal", lw=2, label="QSL Branch $W_p$")
        plt.plot(h_bwd, w_bwd, "--", color="crimson", lw=2, label="Polarized Branch $W_p$")
        plt.axvline(0.3325, color="black", linestyle=":", lw=1.5)
        plt.axvline(0.4000, color="orange", linestyle="--", lw=1.5)

        plt.xlabel("Magnetic Field $h$")
        plt.ylabel(r"Plaquette Flux $W_p$")
        plt.title(r"Topological Flux Collapse at Spinodal $h_{c2}^+$")
        plt.legend(fontsize=9)

        plt.tight_layout()
        fig8_path = FIGURES_DIR / "fig8_critical_zoom_crossing_and_spinodal.png"
        plt.savefig(fig8_path, dpi=200)
        plt.close()
        print(f"Generated -> {fig8_path}")

        # ---------------------------------------------------------------------
        # Figure 11: Contraction Factor & Critical Slowing Down
        # ---------------------------------------------------------------------
        plt.figure(figsize=(8, 4.5))
        # Filter valid mu values for fwd branch
        valid_mu = [d for d in fwd_runs if not np.isnan(float(d["mu"])) and float(d["h"]) < 0.400]
        h_mu = [d["h"] for d in valid_mu]
        mu_vals = [d["mu"] for d in valid_mu]
        n_phi_vals = [d["n_phi"] for d in valid_mu]

        ax1 = plt.gca()
        color1 = "tab:blue"
        ax1.set_xlabel("Magnetic Field $h$")
        ax1.set_ylabel(r"DIIS Contraction Factor $\mu_{\mathrm{final}}$", color=color1)
        ax1.plot(h_mu, mu_vals, "o-", color=color1, lw=1.5, ms=4, label=r"Contraction $\mu$")
        ax1.tick_params(axis="y", labelcolor=color1)

        ax2 = ax1.twinx()
        color2 = "tab:red"
        ax2.set_ylabel(r"Map Evaluations $N_\Phi$", color=color2)
        ax2.plot(h_mu, n_phi_vals, "s--", color=color2, lw=1.5, ms=4, label=r"Iterations $N_\Phi$")
        ax2.tick_params(axis="y", labelcolor=color2)

        plt.axvline(
            0.400, color="orange", linestyle="--", lw=1.5, label=r"Spinodal $h_{c2}^+ = 0.400$"
        )
        plt.title(r"Critical Slowing Down & Contraction Factor $\mu \to 1^-$")
        plt.tight_layout()
        fig11_path = FIGURES_DIR / "fig11_contraction_factor_critical_slowing.png"
        plt.savefig(fig11_path, dpi=200)
        plt.close()
        print(f"Generated -> {fig11_path}")

    # -------------------------------------------------------------------------
    # Figure 9: Majorana Gap Cubic Scaling
    # -------------------------------------------------------------------------
    gap_path = Path("data/processed/benchmarks/m1_6/majorana_gap_scaling.dat")
    if gap_path.is_file():
        h_gap, gap_k, gap_glob = [], [], []
        with gap_path.open("r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                h_val = float(parts[0])
                if h_val > 0.0:
                    h_gap.append(h_val)
                    gap_k.append(float(parts[1]))
                    gap_glob.append(float(parts[2]))

        h_arr = np.array(h_gap)
        gap_arr = np.array(gap_k)

        plt.figure(figsize=(8, 5))
        plt.loglog(
            h_arr, gap_arr, "o", color="darkgreen", ms=6, label=r"Measured Dirac Gap $\Delta_K(h)$"
        )

        # Fitted line
        fit_idx = (h_arr >= 0.005) & (h_arr <= 0.150)
        log_h = np.log(h_arr[fit_idx])
        log_g = np.log(gap_arr[fit_idx])
        p, log_c = np.polyfit(log_h, log_g, deg=1)
        c = np.exp(log_c)

        h_dense = np.logspace(np.log10(h_arr[0]), np.log10(h_arr[-1]), 100)
        gap_fit = c * (h_dense**p)
        plt.loglog(
            h_dense,
            gap_fit,
            "-",
            color="forestgreen",
            lw=2,
            label=f"Fit: $\\Delta(h) = {c:.2f} h^{{{p:.3f}}}$",
        )

        # Pure cubic reference
        gap_cubic = gap_arr[0] * (h_dense / h_arr[0]) ** 3
        plt.loglog(
            h_dense,
            gap_cubic,
            "--",
            color="gray",
            lw=1.5,
            label=r"Kitaev Chiral $\mathcal{O}(h^3)$ Scaling",
        )

        plt.xlabel("Magnetic Field $h$")
        plt.ylabel(r"Majorana Excitation Gap $\Delta(h)$")
        plt.title(r"Topological Majorana Mass Gap Scaling $\Delta_K \propto h^3$ (PH-08)")
        plt.legend(fontsize=10)
        plt.grid(True, which="both", ls=":", alpha=0.5)

        plt.tight_layout()
        fig9_path = FIGURES_DIR / "fig9_majorana_gap_power_law.png"
        plt.savefig(fig9_path, dpi=200)
        plt.close()
        print(f"Generated -> {fig9_path}")

    # -------------------------------------------------------------------------
    # Figure 10: Supercell Energy Parity & Nematicity
    # -------------------------------------------------------------------------
    super_path = Path("data/processed/benchmarks/m1_6/supercell_comparison.dat")
    if super_path.is_file():
        super_records: list[dict[str, Any]] = []
        with super_path.open("r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                super_records.append(
                    {
                        "h": float(parts[0]),
                        "cell": parts[1],
                        "seed": parts[2],
                        "conv": parts[3] == "True",
                        "e_site": float(parts[7]),
                        "nem": float(parts[9]),
                        "w_p": float(parts[10]),
                    }
                )

        rec_2_qsl = [r for r in super_records if r["cell"] == "2site" and r["seed"] == "qsl"]
        rec_4_qsl = [r for r in super_records if r["cell"] == "4site" and r["seed"] == "qsl"]
        rec_2_pol = [r for r in super_records if r["cell"] == "2site" and r["seed"] == "pol"]
        rec_4_pol = [r for r in super_records if r["cell"] == "4site" and r["seed"] == "pol"]

        plt.figure(figsize=(10, 4.5))

        # Subplot 1: Energy per site parity
        plt.subplot(1, 2, 1)
        plt.plot(
            [r["h"] for r in rec_2_qsl],
            [r["e_site"] for r in rec_2_qsl],
            "o-",
            color="teal",
            label="2-site QSL",
        )
        plt.plot(
            [r["h"] for r in rec_4_qsl],
            [r["e_site"] for r in rec_4_qsl],
            "s--",
            color="cyan",
            label="4-site QSL",
        )
        plt.plot(
            [r["h"] for r in rec_2_pol],
            [r["e_site"] for r in rec_2_pol],
            "^-",
            color="crimson",
            label="2-site Pol",
        )
        plt.plot(
            [r["h"] for r in rec_4_pol],
            [r["e_site"] for r in rec_4_pol],
            "v--",
            color="orange",
            label="4-site Pol",
        )

        plt.xlabel("Magnetic Field $h$")
        plt.ylabel(r"Energy per Site $E_{\mathrm{site}}$")
        plt.title("2-Site vs 4-Site Scale Parity (PH-09)")
        plt.legend(fontsize=9)

        # Subplot 2: Energy difference and nematicity
        plt.subplot(1, 2, 2)
        h_common = [r["h"] for r in rec_2_qsl]
        e_diff_qsl = [
            abs(r2["e_site"] - r4["e_site"]) for r2, r4 in zip(rec_2_qsl, rec_4_qsl, strict=True)
        ]
        nem_4 = [r["nem"] for r in rec_4_qsl]

        plt.semilogy(
            h_common,
            e_diff_qsl,
            "o-",
            color="purple",
            label=r"$|E_{\mathrm{site}}^{(2)} - E_{\mathrm{site}}^{(4)}|$",
        )
        plt.semilogy(
            h_common, nem_4, "s--", color="olive", label=r"4-site Nematicity $N_{\mathrm{nem}}$"
        )

        plt.xlabel("Magnetic Field $h$")
        plt.ylabel("Absolute Discrepancy")
        plt.title(r"Symmetry Parity & Invariance ($< 10^{-4}$)")
        plt.legend(fontsize=9)

        plt.tight_layout()
        fig10_path = FIGURES_DIR / "fig10_supercell_energy_parity.png"
        plt.savefig(fig10_path, dpi=200)
        plt.close()
        print(f"Generated -> {fig10_path}")


if __name__ == "__main__":
    plot_phase3_figures()
