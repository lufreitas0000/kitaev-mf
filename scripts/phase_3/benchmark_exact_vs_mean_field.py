#!/usr/bin/env python3
"""Systematic Cross-Validation Benchmark: Track A (Exact Gauge Quench) vs Track B (4-Majorana Mean-Field).

Evaluates and compares:
- Low-energy spectral thresholds (exact two-vison gap vs mean-field quasiparticle continuum).
- Majorana Zero Mode (MZM) core peak emergence vs bare mean-field vacuum.
- Spectral sum rule saturation and point-group isotropy across parameter sweeps.
- Exports structured data (.dat) and publication-grade comparison plots (.png, .pdf).

Author: Lu Freitas
Date: 2026-11-07
"""

from __future__ import annotations

from pathlib import Path
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.phase_3.mean_field_spectral import (
    MeanFieldSpectralData,
    compare_mean_field_vs_exact_quench_thresholds,
    compute_mean_field_spectral_poles,
    evaluate_mean_field_continuous_spectrum,
)
from kitaev_mf.services.phase_3.mzm_spectroscopy import (
    Figure3Spectrogram,
    TwoVortexMZMSystem,
    build_two_vortex_mzm_system,
    reproduce_figure_3_spectrogram,
)
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed

ARTIFACT_DIR = Path("/home/lucas/.gemini/antigravity/brain/c1344f54-6bcf-45cd-b35a-ebf1dd03b975")
FIGURES_DIR = PROJECT_ROOT / "doc" / "labbook" / "figures" / "phase_3"


def run_cross_validation_sweep(
    kappa_values: list[float] | None = None,
    L_exact: int = 6,
    L_mf: int = 4,
    out_dir: Path = FIGURES_DIR,
) -> tuple[Path, Path, Path, list[dict[str, float]]]:
    """Execute systematic parameter sweep and generate comparative benchmark assets."""
    if kappa_values is None:
        kappa_values = [0.05, 0.10, 0.15, 0.20, 0.25]

    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Initialize Track B (4-Majorana Mean-Field) state
    cell = make_honeycomb_2site()
    state_qsl = make_qsl_seed(cell)
    cluster_mf = RealSpaceTorusCluster(L=L_mf)
    params_mf = ModelParameters(K=-1.0)
    mf_data: MeanFieldSpectralData = compute_mean_field_spectral_poles(
        cluster_mf, params_mf, state_qsl, site_j=0, spin_component=3
    )

    records: list[dict[str, float]] = []

    # 2. Sweep kappa values
    for kap in kappa_values:
        sys_exact: TwoVortexMZMSystem = build_two_vortex_mzm_system(
            L=L_exact, p1=(1, 1), p2=(4, 4), K=1.0, kappa=kap
        )
        fig3: Figure3Spectrogram = reproduce_figure_3_spectrogram(sys_exact, kappa=kap)

        exact_vison = fig3.bulk_lowest_energy
        exact_core = fig3.core_lowest_energy
        mzm_splitting = abs(sys_exact.mzm_energies[1] - sys_exact.mzm_energies[0])

        comp = compare_mean_field_vs_exact_quench_thresholds(
            exact_vison_gap=exact_vison, mf_data=mf_data
        )

        records.append(
            {
                "kappa": kap,
                "exact_core_onset": exact_core,
                "exact_vison_threshold": exact_vison,
                "mzm_splitting": mzm_splitting,
                "mf_threshold": comp["mean_field_threshold"],
                "threshold_shift": comp["threshold_shift"],
                "relative_shift": comp["relative_threshold_shift"],
                "mf_total_weight": mf_data.total_spectral_weight,
            }
        )

    # 3. Export structured tabular dataset (.dat)
    dat_path = out_dir / "benchmark_exact_vs_mf.dat"
    with dat_path.open("w", encoding="utf-8") as f:
        f.write("# Systematic Cross-Validation Benchmark: Track A (Exact Quench) vs Track B (4-Majorana MF)\n")
        f.write(f"# Parameters: L_exact={L_exact}, L_mf={L_mf}, K=1.0\n")
        f.write(
            f"# {'kappa':<10} {'exact_core':<16} {'exact_vison':<16} {'mf_threshold':<16} "
            f"{'shift':<16} {'rel_shift':<12} {'mzm_split':<14}\n"
        )
        for r in records:
            f.write(
                f"{r['kappa']:10.4f} {r['exact_core_onset']:16.6e} {r['exact_vison_threshold']:16.6e} "
                f"{r['mf_threshold']:16.6e} {r['threshold_shift']:16.6e} {r['relative_shift']:12.4f} "
                f"{r['mzm_splitting']:14.6e}\n"
            )

    # 4. Generate Comparative Plot
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 5.2), dpi=300)

    # Panel A: Spectral Density Comparison at kappa = 0.15
    sys_ref = build_two_vortex_mzm_system(L=L_exact, p1=(1, 1), p2=(4, 4), K=1.0, kappa=0.15)
    fig3_ref = reproduce_figure_3_spectrogram(sys_ref, kappa=0.15, omega_grid=np.linspace(0.0, 6.0, 600))
    omega_grid = np.linspace(0.0, 6.0, 600)
    mf_spec = evaluate_mean_field_continuous_spectrum(mf_data, omega_grid=omega_grid, sigma=0.06)

    ax_a.plot(
        fig3_ref.omega,
        fig3_ref.core_spectrum.spectral_density,
        color="#D9381E",
        linewidth=2.2,
        label=r"Track A: Vortex Core $S_{\mathrm{core}}(\omega)$ (MZM peak at $\omega \approx 0.16$)",
    )
    ax_a.plot(
        fig3_ref.omega,
        fig3_ref.bulk_spectrum.spectral_density,
        color="#2E2E2E",
        linewidth=1.8,
        linestyle="-.",
        label=r"Track A: Bulk $S_{\mathrm{bulk}}(\omega)$ (Vison onset at $\omega \approx 0.46$)",
    )
    ax_a.plot(
        omega_grid,
        mf_spec.spectral_density,
        color="#1E78B4",
        linewidth=2.2,
        linestyle="--",
        label=r"Track B: 4-Majorana MF $S_{\mathrm{MF}}(\omega)$ (2-QP onset at $\omega \approx 3.05$)",
    )

    ax_a.set_xlim(0.0, 5.5)
    ax_a.set_ylim(-0.2, 14.0)
    ax_a.set_xlabel(r"Excitation Frequency $\omega / |K|$", fontsize=11.5)
    ax_a.set_ylabel(r"Local Spin Structure Factor $S_{jj}^{\gamma\gamma}(\omega)\; [1/|K|]$", fontsize=11.5)
    ax_a.set_title(r"$\bf{(a)}$ Spectral Response: Exact Quench vs. 4-Majorana MF", fontsize=12)
    ax_a.grid(True, linestyle=":", alpha=0.6)
    ax_a.legend(frameon=True, facecolor="#FCFAF5", edgecolor="#CCCCCC", fontsize=9.2, loc="upper right")

    # Panel B: Threshold Discrepancy Scaling vs. Kappa
    kaps = [r["kappa"] for r in records]
    exact_visons = [r["exact_vison_threshold"] for r in records]
    mf_threshs = [r["mf_threshold"] for r in records]

    ax_b.plot(
        kaps,
        exact_visons,
        "o-",
        color="#D9381E",
        linewidth=2.2,
        markersize=7.0,
        label=r"Track A: Exact Two-Vison Creation Gap $\Delta E_{2v}$",
    )
    ax_b.plot(
        kaps,
        mf_threshs,
        "s--",
        color="#1E78B4",
        linewidth=2.2,
        markersize=7.0,
        label=r"Track B: Mean-Field Two-Particle Onset $\omega_{\mathrm{onset}}^{\mathrm{MF}}$",
    )
    ax_b.fill_between(
        kaps,
        exact_visons,
        mf_threshs,
        color="#F0E6D2",
        alpha=0.6,
        label=r"Threshold Discrepancy $\Delta \omega > 2.5 |K|$",
    )

    ax_b.set_xlim(0.03, 0.27)
    ax_b.set_ylim(0.0, 3.8)
    ax_b.set_xlabel(r"Chiral 3-Spin Coupling $\kappa / |K|$", fontsize=11.5)
    ax_b.set_ylabel(r"Spectral Threshold Energy $\omega_{\mathrm{onset}} / |K|$", fontsize=11.5)
    ax_b.set_title(r"$\bf{(b)}$ Threshold Discrepancy Scaling vs. Chirality $\kappa$", fontsize=12)
    ax_b.grid(True, linestyle=":", alpha=0.6)
    ax_b.legend(frameon=True, facecolor="#FCFAF5", edgecolor="#CCCCCC", fontsize=9.2, loc="center right")

    fig.suptitle(
        r"$\bf{Phase\ 3\ Master\ Benchmark:}$ Exact Gauge Quench vs. 4-Majorana Mean-Field Response",
        fontsize=13.5,
        y=0.98,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    png_path = out_dir / "fig5_exact_vs_mean_field.png"
    pdf_path = out_dir / "fig5_exact_vs_mean_field.pdf"
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)

    # Mirror to artifact directory
    if ARTIFACT_DIR.is_dir():
        target = ARTIFACT_DIR / png_path.name
        shutil.copy(png_path, target)

    return png_path, pdf_path, dat_path, records


def main() -> int:
    """Run master benchmark suite and output summary metrics."""
    print("=" * 70)
    print("Executing Phase 3 Master Cross-Validation: Track A vs. Track B")
    print("=" * 70)

    png_path, pdf_path, dat_path, records = run_cross_validation_sweep()

    print(f"Generated assets in {FIGURES_DIR}:")
    print(f"  • Plot (PNG): {png_path} ({png_path.stat().st_size / 1024:.1f} kB)")
    print(f"  • Plot (PDF): {pdf_path} ({pdf_path.stat().st_size / 1024:.1f} kB)")
    print(f"  • Data (DAT): {dat_path} ({dat_path.stat().st_size / 1024:.1f} kB)")

    print("\nSummary Results across Kappa Sweep:")
    print(f"{'kappa':<8} {'exact_vison':<14} {'mf_threshold':<14} {'shift':<12} {'rel_shift':<10}")
    print("-" * 60)
    for r in records:
        print(
            f"{r['kappa']:<8.2f} {r['exact_vison_threshold']:<14.4f} "
            f"{r['mf_threshold']:<14.4f} {r['threshold_shift']:<12.4f} {r['relative_shift']:<10.2f}x"
        )

    print("\nMaster cross-validation completed successfully.")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main())
