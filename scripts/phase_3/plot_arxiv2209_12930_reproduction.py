#!/usr/bin/env python3
"""Quantitative Reproduction and Visualization of arXiv:2209.12930 Figures 3 & 4.

Generates publication-quality figures and structured numerical datasets:
- Figure 3: Three-tier local STS dynamical spin correlation spectra (Core vs Ring vs Bulk).
- Figure 4a: Radial exponential decay of MZM wavefunction density and low-energy spectral weight.
- Figure 4b: 2D real-space honeycomb lattice conductance map with localized vortex core peaks.
- Structured tabular data files (.dat) for PGFPlots / LaTeX integration.

Author: Lu Freitas
Date: 2026-11-07
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Add project root to path if needed
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from kitaev_mf.services.phase_3.mzm_spectroscopy import (  # noqa: E402
    Figure3Spectrogram,
    Figure4SpatialDecay,
    TwoVortexMZMSystem,
    build_two_vortex_mzm_system,
    reproduce_figure_3_spectrogram,
    reproduce_figure_4_spatial_decay,
)

ARTIFACT_DIR = Path("/home/lucas/.gemini/antigravity/brain/c1344f54-6bcf-45cd-b35a-ebf1dd03b975")
FIGURES_DIR = PROJECT_ROOT / "doc" / "labbook" / "figures" / "phase_3"


def generate_figure_3_assets(
    system: TwoVortexMZMSystem,
    out_dir: Path = FIGURES_DIR,
) -> tuple[Path, Path, Path]:
    """Generate Figure 3 spectrogram plots (.png, .pdf) and structured dataset (.dat)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    omega_grid = np.linspace(0.0, 4.0, 400)
    fig3: Figure3Spectrogram = reproduce_figure_3_spectrogram(
        system, omega_grid=omega_grid, sigma=0.06
    )

    # 1. Export structured tabular data (.dat)
    dat_path = out_dir / "fig3_data.dat"
    with dat_path.open("w", encoding="utf-8") as f:
        f.write("# arXiv:2209.12930 Figure 3 Reproduction Data\n")
        f.write(f"# Parameters: L={system.L}, K=1.0, kappa=0.15, sigma=0.06\n")
        f.write(
            f"# Core Site: {fig3.core_site}, Ring Site: {fig3.ring_site}, Bulk Site: {fig3.bulk_site}\n"
        )
        f.write(
            f"# Lowest Core Energy: {fig3.core_lowest_energy:.6f}, Lowest Bulk Energy: {fig3.bulk_lowest_energy:.6f}\n"
        )
        f.write(f"# {'omega/K':<12} {'S_core':<16} {'S_ring':<16} {'S_bulk':<16}\n")
        for w, c, r, b in zip(
            fig3.omega,
            fig3.core_spectrum.spectral_density,
            fig3.ring_spectrum.spectral_density,
            fig3.bulk_spectrum.spectral_density,
            strict=False,
        ):
            f.write(f"{w:12.6f} {c:16.8e} {r:16.8e} {b:16.8e}\n")

    # 2. Render Figure 3 Plot
    fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=300)

    ax.plot(
        fig3.omega,
        fig3.core_spectrum.spectral_density,
        color="#D9381E",
        linewidth=2.4,
        label="Core: Tip on vortex boundary ($r = 0$)",
    )
    ax.plot(
        fig3.omega,
        fig3.ring_spectrum.spectral_density,
        color="#1E78B4",
        linewidth=2.0,
        linestyle="--",
        label="Ring: First neighbor ($r \\approx 1.73 a$)",
    )
    ax.plot(
        fig3.omega,
        fig3.bulk_spectrum.spectral_density,
        color="#2E2E2E",
        linewidth=2.0,
        linestyle="-.",
        label="Bulk: Far from vortices ($r \\gg \\xi_M$)",
    )

    # Key physics annotations
    # Core peak
    core_peak_idx = int(np.argmax(fig3.core_spectrum.spectral_density))
    w_peak = fig3.omega[core_peak_idx]
    s_peak = fig3.core_spectrum.spectral_density[core_peak_idx]
    ax.annotate(
        f"MZM Zero-Bias Peak\n$\\omega = {w_peak:.3f}|K|$",
        xy=(w_peak, s_peak),
        xytext=(w_peak + 0.45, s_peak * 0.85),
        arrowprops={"facecolor": "#D9381E", "edgecolor": "#D9381E", "width": 1.5, "headwidth": 7},
        fontsize=10,
        fontweight="bold",
        color="#9E200C",
    )

    # Bulk threshold onset
    ax.axvline(
        x=fig3.bulk_lowest_energy,
        color="#707070",
        linestyle=":",
        linewidth=1.3,
        alpha=0.85,
    )
    ax.annotate(
        f"Bulk Vison Threshold\n$\\Delta E_{{2v}} = {fig3.bulk_lowest_energy:.3f}|K|$",
        xy=(fig3.bulk_lowest_energy, 0.5),
        xytext=(fig3.bulk_lowest_energy + 0.35, 2.5),
        arrowprops={"facecolor": "#505050", "edgecolor": "#505050", "width": 1.2, "headwidth": 6},
        fontsize=9.5,
        color="#333333",
    )

    ax.set_xlim(0.0, 3.5)
    ax.set_ylim(-0.2, max(s_peak * 1.15, 14.0))
    ax.set_xlabel(r"Excitation Frequency $\omega / |K|$", fontsize=12)
    ax.set_ylabel(
        r"Local Dynamical Spin Correlation $S_{jj}^{\gamma\gamma}(\omega)\; [1/|K|]$", fontsize=12
    )
    ax.set_title(
        r"$\bf{arXiv:2209.12930\ Fig.\ 3\ Reproduction:}$ Local STS Spin Correlations"
        f"\n($L={system.L}$, $K=1.0$, $\\kappa={system.bulk_gap_estimate / (6 * np.sqrt(3)):.2f}$, $\\sigma=0.06$)",
        fontsize=12,
        pad=10,
    )
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(
        frameon=True, facecolor="#FCFAF5", edgecolor="#CCCCCC", fontsize=10.5, loc="upper right"
    )

    fig.tight_layout()

    png_path = out_dir / "fig3_sts_spectrogram.png"
    pdf_path = out_dir / "fig3_sts_spectrogram.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)

    return png_path, pdf_path, dat_path


def generate_figure_4_assets(
    system: TwoVortexMZMSystem,
    out_dir: Path = FIGURES_DIR,
) -> tuple[Path, Path, Path]:
    """Generate Figure 4 spatial decay & 2D lattice map (.png, .pdf) and dataset (.dat)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fig4: Figure4SpatialDecay = reproduce_figure_4_spatial_decay(
        system, omega_cutoff=0.35, max_fit_distance=3.5
    )

    # 1. Export structured tabular data (.dat)
    dat_path = out_dir / "fig4_data.dat"
    with dat_path.open("w", encoding="utf-8") as f:
        f.write("# arXiv:2209.12930 Figure 4 Reproduction Data\n")
        f.write(f"# Parameters: L={system.L}, K=1.0, kappa=0.15\n")
        f.write(f"# Fitted Localization Length xi_M: {fig4.fitted_xi_M:.6f} a\n")
        f.write(f"# Fit Slope: {fig4.fit_slope:.6f}, R^2: {fig4.fit_r_squared:.6f}\n")
        f.write(f"# {'r/a':<12} {'rho_MZM':<16} {'I_low':<16} {'rho_fit':<16}\n")
        for r, rho, w in zip(
            fig4.radial_distances, fig4.mzm_densities, fig4.low_energy_weights, strict=False
        ):
            fit_val = float(np.exp(fig4.fit_slope * r + np.log(fig4.mzm_densities[0])))
            f.write(f"{r:12.6f} {rho:16.8e} {w:16.8e} {fit_val:16.8e}\n")

    # 2. Render Figure 4 Combined Plot (Panel A: Radial Decay, Panel B: 2D Real-Space Map)
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 5.8), dpi=300)

    # --- Panel A: Radial Exponential Decay ---
    r_vals = fig4.radial_distances
    rho_vals = fig4.mzm_densities
    w_vals = fig4.low_energy_weights

    ax_a.semilogy(
        r_vals,
        rho_vals,
        "o",
        color="#D9381E",
        markersize=7.5,
        markeredgecolor="#7A1808",
        label=r"MZM Density $\rho_{\mathrm{MZM}}(r) = |\phi_0|^2 + |\phi_1|^2$",
    )
    ax_a.semilogy(
        r_vals,
        w_vals,
        "s",
        color="#1E78B4",
        markersize=6.5,
        markeredgecolor="#0E446B",
        alpha=0.8,
        label=r"Low-Energy Conductance $I_{\mathrm{low}}(r)$ ($\omega \leq 0.35$)",
    )

    # Regression fit line
    r_fit = np.linspace(0.0, 3.5, 100)
    # intercept from first point
    intercept = np.log(rho_vals[0])
    rho_fit_curve = np.exp(fig4.fit_slope * r_fit + intercept)
    ax_a.plot(
        r_fit,
        rho_fit_curve,
        "--",
        color="#B02510",
        linewidth=2.0,
        label=f"Fit: $\\xi_M = {fig4.fitted_xi_M:.2f}a$ ($R^2 = {fig4.fit_r_squared:.3f}$)",
    )

    # Theoretical bulk decay length estimate
    xi_theory = 1.0 / (2.0 * np.sqrt(3.0) * (system.bulk_gap_estimate / (6 * np.sqrt(3))))
    ax_a.text(
        0.05,
        0.12,
        f"Fitted $\\xi_M = {fig4.fitted_xi_M:.2f} a$\n"
        f"Theory $\\xi_{{\\mathrm{{bulk}}}} \\approx \\frac{{\\sqrt{{3}}K}}{{6\\kappa}} = {xi_theory:.2f} a$\n"
        f"Fit $R^2 = {fig4.fit_r_squared:.3f}$",
        transform=ax_a.transAxes,
        fontsize=10.5,
        bbox={"boxstyle": "round,pad=0.5", "facecolor": "#FCFAF5", "edgecolor": "#DDDDDD"},
    )

    ax_a.set_xlim(-0.2, 3.8)
    ax_a.set_ylim(min(rho_vals) * 0.4, max(rho_vals) * 2.5)
    ax_a.set_xlabel(r"Euclidean Distance $r / a$ from Vortex Core", fontsize=11.5)
    ax_a.set_ylabel(r"Local Spectral Weight & Density (log scale)", fontsize=11.5)
    ax_a.set_title(r"$\bf{(a)}$ Spatial Radial Decay of MZM Wavefunction", fontsize=12)
    ax_a.grid(True, linestyle=":", alpha=0.6)
    ax_a.legend(
        frameon=True, facecolor="#FCFAF5", edgecolor="#CCCCCC", fontsize=9.5, loc="upper right"
    )

    # --- Panel B: 2D Real-Space Honeycomb Conductance Map ---
    geom = system.geom

    # Honeycomb lattice vectors for regular hexagons
    a1 = np.array([np.sqrt(3.0) / 2.0, 1.5])
    a2 = np.array([-np.sqrt(3.0) / 2.0, 1.5])
    dz = np.array([0.0, 1.0])

    site_positions = np.zeros((geom.num_sites, 2))
    for s in range(geom.num_sites):
        m1, m2, sub = geom.site_coordinates(s)
        site_positions[s] = m1 * a1 + m2 * a2 + (dz if sub == 1 else np.zeros(2))

    # Draw nearest-neighbor bonds
    for bond in system.bonds:
        pos_a = site_positions[bond.site_a]
        pos_b = site_positions[bond.site_b]
        # Only draw intra-cell or short bonds (skip periodic wrap-around links)
        if np.linalg.norm(pos_a - pos_b) < 1.5:
            ax_b.plot(
                [pos_a[0], pos_b[0]],
                [pos_a[1], pos_b[1]],
                color="#CCCCCC",
                linewidth=1.2,
                zorder=1,
            )

    # Plot sites colored by MZM density
    sc = ax_b.scatter(
        site_positions[:, 0],
        site_positions[:, 1],
        c=system.mzm_density,
        cmap="YlOrRd",
        s=90,
        edgecolors="#444444",
        linewidths=0.8,
        zorder=2,
    )
    cbar = fig.colorbar(sc, ax=ax_b, fraction=0.046, pad=0.04)
    cbar.set_label(r"Majorana Zero Mode Density $\rho_{\mathrm{MZM}}(\mathbf{r}_j)$", fontsize=10.5)

    # Calculate and mark vortex plaquette centers
    def plaquette_center(p: tuple[int, int]) -> np.ndarray:
        m1, m2 = p
        # Average coordinates of bounding sites
        bonds = geom.plaquette_bonds(m1, m2)
        p_sites = list({s for b in bonds for s in b})
        # Use positions relative to first site to avoid torus wrap jumps
        site_positions[p_sites[0]]
        coords = []
        for s in p_sites:
            pos = site_positions[s]
            coords.append(pos)
        return np.mean(coords, axis=0)

    center_p1 = plaquette_center(system.p1)
    center_p2 = plaquette_center(system.p2)

    ax_b.scatter(
        [center_p1[0]],
        [center_p1[1]],
        marker="*",
        s=300,
        color="#228B22",
        edgecolors="black",
        linewidths=1.2,
        label=f"Vortex $p_1={system.p1}$",
        zorder=3,
    )
    ax_b.scatter(
        [center_p2[0]],
        [center_p2[1]],
        marker="*",
        s=300,
        color="#8B008B",
        edgecolors="black",
        linewidths=1.2,
        label=f"Vortex $p_2={system.p2}$",
        zorder=3,
    )

    ax_b.set_aspect("equal", adjustable="box")
    ax_b.set_xlabel(r"$x / a$", fontsize=11.5)
    ax_b.set_ylabel(r"$y / a$", fontsize=11.5)
    ax_b.set_title(r"$\bf{(b)}$ 2D Real-Space MZM Conductance Map", fontsize=12)
    ax_b.legend(
        frameon=True, facecolor="#FCFAF5", edgecolor="#CCCCCC", fontsize=9.5, loc="lower left"
    )

    fig.suptitle(
        r"$\bf{arXiv:2209.12930\ Fig.\ 4\ Reproduction:}$ MZM Localization & Real-Space Topology",
        fontsize=13.5,
        y=0.98,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))

    png_path = out_dir / "fig4_mzm_spatial_decay.png"
    pdf_path = out_dir / "fig4_mzm_spatial_decay.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)

    return png_path, pdf_path, dat_path


def main() -> int:
    """Generate all Phase 3 publication figures, benchmarks, and data files."""
    print("=" * 70)
    print("Executing Phase 3 Reproduction: arXiv:2209.12930 Figures 3 & 4")
    print("=" * 70)

    print("Building two-vortex MZM system (L=6, p1=(1,1), p2=(4,4), K=1.0, kappa=0.15)...")
    system = build_two_vortex_mzm_system(L=6, p1=(1, 1), p2=(4, 4), K=1.0, kappa=0.15)
    print(
        f"Mid-gap MZM energies: eps_0 = {system.mzm_energies[0]:.4f}, eps_1 = {system.mzm_energies[1]:.4f}"
    )
    print(f"Bulk gap estimate:    Delta_bulk = {system.bulk_gap_estimate:.4f}")
    print(f"In-gap ratio:         eps_1 / eps_2 = {system.in_gap_ratio:.4f}")

    print("\nGenerating Figure 3 STS Spectrogram assets...")
    f3_png, f3_pdf, f3_dat = generate_figure_3_assets(system)
    print(f"  Saved: {f3_png} ({f3_png.stat().st_size / 1024:.1f} kB)")
    print(f"  Saved: {f3_pdf} ({f3_pdf.stat().st_size / 1024:.1f} kB)")
    print(f"  Saved: {f3_dat} ({f3_dat.stat().st_size / 1024:.1f} kB)")

    print("\nGenerating Figure 4 Spatial Decay & 2D Lattice Map assets...")
    f4_png, f4_pdf, f4_dat = generate_figure_4_assets(system)
    print(f"  Saved: {f4_png} ({f4_png.stat().st_size / 1024:.1f} kB)")
    print(f"  Saved: {f4_pdf} ({f4_pdf.stat().st_size / 1024:.1f} kB)")
    print(f"  Saved: {f4_dat} ({f4_dat.stat().st_size / 1024:.1f} kB)")

    # Copy to artifact directory if available
    if ARTIFACT_DIR.is_dir():
        for p in [f3_png, f4_png]:
            target = ARTIFACT_DIR / p.name
            shutil.copy(p, target)
            print(f"  Mirrored to artifact: {target}")

    print("\nPhase 3 visualization assets generated successfully.")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
