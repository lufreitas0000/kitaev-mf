"""Unit tests verifying integrity and formatting of exported benchmark data files."""

from pathlib import Path
import numpy as np

BENCHMARK_DIR = Path("data/processed/benchmarks")

EXPECTED_BENCHMARKS = [
    "k006_quadrature.dat",
    "hellmann_feynman_scaling.dat",
    "susceptibility_curve.dat",
    "klein_duality_splitting.dat",
    "solver_convergence.dat",
    "dispersion_path.dat",
    "parseval_sum_rule.dat",
    "wick_cumulant_check.dat",
]


def test_benchmark_files_exist_and_have_provenance() -> None:
    """All 8 expected benchmark files must exist with valid header metadata."""
    for b_name in EXPECTED_BENCHMARKS:
        path = BENCHMARK_DIR / b_name
        assert path.is_file(), f"Missing benchmark file: {path}"
        lines = path.read_text(encoding="utf-8").strip().splitlines()
        assert len(lines) > 5, f"Benchmark {b_name} too short"

        # Check provenance headers
        assert lines[0].startswith(f"# Benchmark: {b_name}")
        assert lines[1].startswith("# Generated:")
        assert lines[2].startswith("# Git Commit:")
        assert lines[3].startswith("# Description:")
        assert lines[4].startswith("# Columns:")


def test_parseval_sum_rule_data_integrity() -> None:
    """Parseval-Frobenius residuals must be within 1e-14 across all meshes."""
    path = BENCHMARK_DIR / "parseval_sum_rule.dat"
    data = np.loadtxt(path, comments="#")
    assert data.shape[0] >= 6
    residuals = data[:, 3]
    np.testing.assert_allclose(residuals, 0.0, atol=1e-14)


def test_dispersion_path_data_integrity() -> None:
    """Band dispersion must match bandwidth 3*|K| at Gamma."""
    path = BENCHMARK_DIR / "dispersion_path.dat"
    data = np.loadtxt(path, comments="#")
    assert data.shape[0] > 100
    # First row is Gamma (0,0): e1 should be -1.5, e8 should be +1.5
    np.testing.assert_allclose(data[0, 3], -1.5, atol=1e-4)
    np.testing.assert_allclose(data[0, 10], 1.5, atol=1e-4)
