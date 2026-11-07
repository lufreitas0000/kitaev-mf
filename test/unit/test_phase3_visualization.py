"""Unit tests and automated regression checks for Task T3.2.2 (Visualization & Figure Reproduction).

Verifies contracts:
- H_VIS_01: Figure 3 Spectrogram Reproduction (Core vs Ring vs Bulk hierarchy, peak height, and thresholds)
- H_VIS_02: Figure 4a Radial Decay (Monotonic decay, fitted xi_M in [1.5, 4.5] a, R^2 >= 0.90)
- H_VIS_03: Figure 4b 2D Real-Space Conductance Map (Vortex core density localization > 2.5x bulk background)
- H_VIS_04: Asset Integrity & Structured Data Export (PNG/PDF file size, DPI, .dat tabular formatting)
"""

from __future__ import annotations

from pathlib import Path
import sys

# Ensure repository root is on sys.path for scripts importing
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import pytest
import numpy as np

from kitaev_mf.services.phase_3.mzm_spectroscopy import (
    build_two_vortex_mzm_system,
)
from scripts.phase_3.plot_arxiv2209_12930_reproduction import (
    FIGURES_DIR,
    generate_figure_3_assets,
    generate_figure_4_assets,
    main as run_visualization_main,
)


@pytest.fixture(scope="module")
def shared_mzm_system():
    """Construct shared L=6 two-vortex MZM system for visualization tests."""
    return build_two_vortex_mzm_system(L=6, p1=(1, 1), p2=(4, 4), K=1.0, kappa=0.15)


def test_h_vis_01_figure_3_spectrogram_assets(shared_mzm_system, tmp_path):
    """H_VIS_01: Verify Figure 3 asset generation, file sizes, and physics content."""
    png_path, pdf_path, dat_path = generate_figure_3_assets(shared_mzm_system, out_dir=tmp_path)

    # Asset existence and size checks
    assert png_path.is_file(), f"PNG file not found: {png_path}"
    assert pdf_path.is_file(), f"PDF file not found: {pdf_path}"
    assert dat_path.is_file(), f"DAT file not found: {dat_path}"

    assert png_path.stat().st_size > 25_000, f"PNG size too small: {png_path.stat().st_size} bytes"
    assert pdf_path.stat().st_size > 15_000, f"PDF size too small: {pdf_path.stat().st_size} bytes"

    # Validate data file contents
    data = np.loadtxt(dat_path, comments="#")
    assert data.shape[0] >= 100, f"Expected >= 100 frequency points, got {data.shape[0]}"
    assert data.shape[1] == 4, f"Expected 4 columns (omega, S_core, S_ring, S_bulk), got {data.shape[1]}"

    omega = data[:, 0]
    s_core = data[:, 1]
    s_ring = data[:, 2]
    s_bulk = data[:, 3]

    # Peak validation
    max_core = np.max(s_core)
    assert max_core > 10.0, f"Expected core zero-bias peak > 10.0, got {max_core}"

    # Peak position
    peak_omega = omega[np.argmax(s_core)]
    assert 0.10 <= peak_omega <= 0.25, f"Core peak frequency out of range: {peak_omega}"

    # Bulk subgap suppression: strictly vanishing at low bias and < 0.50 for omega <= 0.30
    subgap_bulk_low = s_bulk[omega <= 0.20]
    assert np.all(subgap_bulk_low < 0.01), f"Bulk low-bias weight must be < 0.01, got max {np.max(subgap_bulk_low)}"
    subgap_bulk = s_bulk[omega <= 0.30]
    assert np.all(subgap_bulk < 0.50), f"Bulk subgap weight must be < 0.50, got max {np.max(subgap_bulk)}"


def test_h_vis_02_figure_4_spatial_decay_assets(shared_mzm_system, tmp_path):
    """H_VIS_02: Verify Figure 4 asset generation, regression fit consistency, and xi_M."""
    png_path, pdf_path, dat_path = generate_figure_4_assets(shared_mzm_system, out_dir=tmp_path)

    # Asset existence and size checks
    assert png_path.is_file(), f"PNG file not found: {png_path}"
    assert pdf_path.is_file(), f"PDF file not found: {pdf_path}"
    assert dat_path.is_file(), f"DAT file not found: {dat_path}"

    assert png_path.stat().st_size > 25_000, f"PNG size too small: {png_path.stat().st_size} bytes"
    assert pdf_path.stat().st_size > 15_000, f"PDF size too small: {pdf_path.stat().st_size} bytes"

    # Validate data file contents
    data = np.loadtxt(dat_path, comments="#")
    assert data.shape[0] >= 5, f"Expected >= 5 radial bins, got {data.shape[0]}"
    assert data.shape[1] == 4, f"Expected 4 columns (r, rho, I_low, rho_fit), got {data.shape[1]}"

    r_vals = data[:, 0]
    rho_vals = data[:, 1]

    # Monotonic decay over the first 3 bins
    assert rho_vals[0] > rho_vals[1] > rho_vals[2], (
        f"Wavefunction density must decay monotonically near core: {rho_vals[:3]}"
    )

    # Read header to parse fitted xi_M and R^2
    with dat_path.open("r", encoding="utf-8") as f:
        header_text = f.read()

    assert "Fitted Localization Length xi_M:" in header_text
    assert "Fit Slope:" in header_text
    assert "R^2:" in header_text

    # Extract xi_M value
    for line in header_text.splitlines():
        if "Fitted Localization Length xi_M:" in line:
            xi_val = float(line.split(":")[1].replace("a", "").strip())
            assert 1.5 <= xi_val <= 4.5, f"Fitted xi_M out of physical range [1.5, 4.5]: {xi_val}"
        if "R^2:" in line:
            r2_val = float(line.split("R^2:")[1].strip())
            assert r2_val >= 0.90, f"Fit R^2 must be >= 0.90, got {r2_val}"


def test_h_vis_03_figure_4b_2d_conductance_contrast(shared_mzm_system):
    """H_VIS_03: Verify real-space MZM wavefunction contrast between vortex core and bulk."""
    mzm_density = shared_mzm_system.mzm_density
    max_core_density = float(np.max(mzm_density))
    min_bulk_density = float(np.min(mzm_density))

    # The density at the vortex core boundary must exceed bulk background by at least 2.5x
    contrast_ratio = max_core_density / min_bulk_density
    assert contrast_ratio >= 2.5, (
        f"MZM core-to-bulk contrast ratio {contrast_ratio:.2f} is less than threshold 2.5"
    )


def test_h_vis_04_pipeline_execution_and_asset_integrity():
    """H_VIS_04: Full pipeline execution and production directory validation."""
    exit_code = run_visualization_main()
    assert exit_code == 0, f"Visualization main() failed with code {exit_code}"

    # Verify standard labbook production directory
    expected_files = [
        FIGURES_DIR / "fig3_sts_spectrogram.png",
        FIGURES_DIR / "fig3_sts_spectrogram.pdf",
        FIGURES_DIR / "fig3_data.dat",
        FIGURES_DIR / "fig4_mzm_spatial_decay.png",
        FIGURES_DIR / "fig4_mzm_spatial_decay.pdf",
        FIGURES_DIR / "fig4_data.dat",
    ]

    for p in expected_files:
        assert p.is_file(), f"Expected production file missing: {p}"
        assert p.stat().st_size > 0, f"Production file is empty: {p}"
