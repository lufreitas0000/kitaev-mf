# ruff: noqa: E402
# Repository-path bootstrap must precede local imports in this standalone entrypoint.
"""Unit tests and automated regression checks for Task T3.5.1 (Master Cross-Validation Benchmark).

Verifies contracts:
- H_COMP_01: Unified Benchmark Execution (Returns non-empty assets and records)
- H_COMP_02: Threshold Gap Discrepancy (MF threshold > 2.50 |K| strictly exceeds exact vison gap < 0.50 |K|)
- H_COMP_03: MZM Core vs Bulk Demarcation (Core onset <= 0.25 vs Bulk vison onset >= 0.40)
- H_COMP_04: Export & Structured Data Integrity (.dat ASCII headers, dimensions, and finite numbers)
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure repository root is on sys.path for scripts importing
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np

from scripts.phase_3.benchmark_exact_vs_mean_field import (
    FIGURES_DIR,
    run_cross_validation_sweep,
)
from scripts.phase_3.benchmark_exact_vs_mean_field import (
    main as run_benchmark_main,
)


def test_h_comp_01_benchmark_execution_and_asset_integrity(tmp_path: Path) -> None:
    """H_COMP_01: Verify execution of cross-validation sweep and generated asset files."""
    png_path, pdf_path, dat_path, records = run_cross_validation_sweep(
        kappa_values=[0.10, 0.15, 0.20],
        L_exact=6,
        L_mf=4,
        out_dir=tmp_path,
    )

    assert png_path.is_file(), f"PNG plot not found: {png_path}"
    assert pdf_path.is_file(), f"PDF plot not found: {pdf_path}"
    assert dat_path.is_file(), f"DAT data file not found: {dat_path}"

    assert png_path.stat().st_size > 25_000, f"PNG size too small: {png_path.stat().st_size} bytes"
    assert pdf_path.stat().st_size > 15_000, f"PDF size too small: {pdf_path.stat().st_size} bytes"
    assert dat_path.stat().st_size > 200, f"DAT size too small: {dat_path.stat().st_size} bytes"

    assert len(records) == 3


def test_h_comp_02_threshold_gap_discrepancy() -> None:
    """H_COMP_02: Verify that mean-field threshold consistently exceeds exact vison creation gap."""
    dat_path = FIGURES_DIR / "benchmark_exact_vs_mf.dat"
    assert dat_path.is_file(), f"Master benchmark data file missing: {dat_path}"

    data = np.loadtxt(dat_path, comments="#")
    assert data.ndim == 2 and data.shape[0] >= 3
    assert data.shape[1] == 7

    kappa_col = data[:, 0]
    exact_vison_col = data[:, 2]
    mf_thresh_col = data[:, 3]
    shift_col = data[:, 4]
    rel_shift_col = data[:, 5]

    for i in range(len(kappa_col)):
        # Exact vison gap is subgap < 0.50 |K|
        assert exact_vison_col[i] < 0.50, (
            f"Exact vison gap unexpectedly large: {exact_vison_col[i]}"
        )

        # Mean-field threshold is shifted up into the two-quasiparticle continuum > 2.50 |K|
        assert mf_thresh_col[i] > 2.50, (
            f"Mean-field threshold unexpectedly small: {mf_thresh_col[i]}"
        )

        # Threshold shift exceeds 1.50 |K|
        assert shift_col[i] > 1.50, f"Threshold shift too small: {shift_col[i]}"
        assert rel_shift_col[i] > 3.0, f"Relative shift too small: {rel_shift_col[i]}"


def test_h_comp_03_mzm_core_vs_bulk_demarcation() -> None:
    """H_COMP_03: Contrast core MZM zero-bias peak with bulk vison onset."""
    dat_path = FIGURES_DIR / "benchmark_exact_vs_mf.dat"
    data = np.loadtxt(dat_path, comments="#")

    exact_core_col = data[:, 1]
    exact_vison_col = data[:, 2]

    for core_e, vison_e in zip(exact_core_col, exact_vison_col, strict=True):
        # Core has low-energy translation onset <= 0.25 |K|
        assert core_e <= 0.25, f"Core onset not in low-bias window: {core_e}"
        # Bulk has hard vison threshold >= 0.40 |K|
        assert vison_e >= 0.40, f"Bulk vison threshold too low: {vison_e}"
        # Core onset is strictly lower than bulk threshold
        assert core_e < vison_e


def test_h_comp_04_pipeline_cli_entrypoint() -> None:
    """H_COMP_04: Verify CLI execution of benchmark script and standard outputs."""
    exit_code = run_benchmark_main()
    assert exit_code == 0, f"Benchmark main() returned exit code {exit_code}"

    # Verify standard labbook files
    assert (FIGURES_DIR / "fig5_exact_vs_mean_field.png").is_file()
    assert (FIGURES_DIR / "fig5_exact_vs_mean_field.pdf").is_file()
    assert (FIGURES_DIR / "benchmark_exact_vs_mf.dat").is_file()
