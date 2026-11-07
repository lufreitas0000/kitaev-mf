"""Unit tests for Task T2.3.2: Critical Exponents & Scaling at Polarized Breakdown Transition."""

from pathlib import Path

import numpy as np
import pytest

from kitaev_mf.services.phase_2.spinodal_scaling import (
    SpinodalBayesianFitter,
    SpinodalPoint,
    SpinodalScalingCalculator,
    SusceptibilityPoint,
)


@pytest.fixture
def spinodal_calc() -> SpinodalScalingCalculator:
    return SpinodalScalingCalculator(k_grid=(12, 12), h_spinodal=0.400)


def test_t2_3_2_forward_continuation_convergence(spinodal_calc: SpinodalScalingCalculator) -> None:
    """Verify that forward continuation converges smoothly on approach to the spinodal limit."""
    h_test = np.array([0.360, 0.380, 0.395])
    points = spinodal_calc.run_forward_continuation(h_test)

    assert len(points) == 3
    for pt in points:
        assert pt.converged is True
        assert pt.m_111 > 0.15

    # Magnetization must be strictly monotonically increasing with field
    assert points[0].m_111 < points[1].m_111 < points[2].m_111


def test_t2_3_2_compute_susceptibility_series(spinodal_calc: SpinodalScalingCalculator) -> None:
    """Verify compute_susceptibility_series differential calculation."""
    pts = [
        SpinodalPoint(h_mag=0.350, m_111=0.200, converged=True, iterations=10),
        SpinodalPoint(h_mag=0.360, m_111=0.212, converged=True, iterations=10),
        SpinodalPoint(h_mag=0.370, m_111=0.228, converged=True, iterations=10),
    ]
    susc = spinodal_calc.compute_susceptibility_series(pts)
    assert len(susc) == 2
    assert np.isclose(susc[0].chi_111, 1.20)
    assert np.isclose(susc[1].chi_111, 1.60)
    assert np.isclose(susc[0].distance_to_spinodal, 0.400 - 0.355)


def test_t2_3_2_spinodal_divergence_exponent() -> None:
    """Verify that the fitted spinodal exponent gamma_sp matches 0.50 within 15% error."""
    benchmark_path = Path("data/processed/benchmarks/m2_3/spinodal_scaling_benchmark.dat")
    assert benchmark_path.exists()

    series = []
    with open(benchmark_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip().startswith("#") and line.strip():
                parts = line.strip().split()
                series.append(
                    SusceptibilityPoint(
                        h_mid=float(parts[0]),
                        distance_to_spinodal=float(parts[1]),
                        chi_111=float(parts[2]),
                    )
                )

    fitter = SpinodalBayesianFitter()
    res = fitter.fit_and_compare(series)

    assert res.n_samples == 10
    # Pre-registered threshold: error <= 0.15 (15%)
    assert res.gamma_relative_error <= 0.15, (
        f"Fitted exponent {res.fitted_gamma:.4f} deviated from 0.50 by {res.gamma_relative_error * 100:.2f}%"
    )
    assert 0.45 <= res.fitted_gamma <= 0.60


def test_t2_3_2_bayesian_model_selection() -> None:
    """Verify that Bayesian model comparison decisively favors mean-field spinodal divergence (Delta_BIC > 10)."""
    benchmark_path = Path("data/processed/benchmarks/m2_3/spinodal_scaling_benchmark.dat")
    assert benchmark_path.exists()

    series = []
    with open(benchmark_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip().startswith("#") and line.strip():
                parts = line.strip().split()
                series.append(
                    SusceptibilityPoint(
                        h_mid=float(parts[0]),
                        distance_to_spinodal=float(parts[1]),
                        chi_111=float(parts[2]),
                    )
                )

    fitter = SpinodalBayesianFitter()
    res = fitter.fit_and_compare(series)

    assert res.delta_bic > 10.0, f"Expected Delta_BIC > 10, got {res.delta_bic:.2f}"
    assert res.bayes_factor > 148.0, f"Expected Bayes factor > 148, got {res.bayes_factor:.2e}"
    assert res.favors_spinodal_divergence is True


def test_t2_3_2_synthetic_bayesian_fitter() -> None:
    """Verify SpinodalBayesianFitter on an exact synthetic 1/2 power-law series."""
    t_vals = np.linspace(0.005, 0.050, 10)
    synthetic_chi = 0.25 * (t_vals ** (-0.50))

    series = [
        SusceptibilityPoint(h_mid=float(0.400 - t), distance_to_spinodal=float(t), chi_111=float(c))
        for t, c in zip(t_vals, synthetic_chi, strict=True)
    ]

    fitter = SpinodalBayesianFitter()
    res = fitter.fit_and_compare(series)

    assert res.favors_spinodal_divergence is True
    assert np.isclose(res.fitted_gamma, 0.50, atol=1e-3)
    assert res.delta_bic > 20.0


def test_t2_3_2_benchmark_dataset_acid_provenance() -> None:
    """Verify ACID provenance header and format of spinodal scaling benchmark table."""
    benchmark_path = Path("data/processed/benchmarks/m2_3/spinodal_scaling_benchmark.dat")
    assert benchmark_path.exists()

    with open(benchmark_path, encoding="utf-8") as f:
        lines = f.readlines()

    header_lines = [line for line in lines if line.strip().startswith("#")]
    data_lines = [line for line in lines if not line.strip().startswith("#") and line.strip()]

    assert len(header_lines) >= 6
    assert any("Git Commit:" in line for line in header_lines)
    assert any("Generated:" in line for line in header_lines)
    assert any("Bayesian Model Selection Summary:" in line for line in header_lines)
    assert len(data_lines) == 10

    for dl in data_lines:
        tokens = dl.strip().split()
        assert len(tokens) == 3
        h_mid = float(tokens[0])
        t_val = float(tokens[1])
        chi_val = float(tokens[2])
        assert 0.35 <= h_mid <= 0.40
        assert 0.0 < t_val <= 0.05
        assert chi_val > 0.0
