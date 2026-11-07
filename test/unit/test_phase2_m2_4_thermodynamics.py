"""Unit tests for Task T2.4.1: Two-Peak Specific Heat C_v(T) & Fractionalization Entropy Stages."""

from pathlib import Path

import numpy as np
import pytest

from kitaev_mf.services.phase_2.thermodynamics import (
    SpecificHeatCalculator,
    ThermodynamicBayesianFitter,
)


@pytest.fixture
def heat_calc() -> SpecificHeatCalculator:
    return SpecificHeatCalculator(delta_v=0.137, k_grid=(24, 24))


def test_t2_4_1_two_peak_structure(heat_calc: SpecificHeatCalculator) -> None:
    """Verify that specific heat C_v(T) possesses two well-separated peaks with T_H / T_L >= 5.0."""
    temps = np.logspace(-2.5, 0.8, 35)
    res = heat_calc.compute_thermodynamic_curve(temps)

    assert res.t_low < res.t_high
    assert res.peak_ratio >= 5.0, (
        f"Expected peak ratio >= 5.0, got {res.peak_ratio:.2f} (T_L={res.t_low:.4f}, T_H={res.t_high:.4f})"
    )
    assert 0.03 <= res.t_low <= 0.10, f"Expected T_L in [0.03, 0.10], got {res.t_low:.4f}"
    assert 0.40 <= res.t_high <= 0.80, f"Expected T_H in [0.40, 0.80], got {res.t_high:.4f}"
    assert res.c_v_low > 0.15
    assert res.c_v_high > 0.10


def test_t2_4_1_intermediate_entropy_plateau(heat_calc: SpecificHeatCalculator) -> None:
    """Verify that entropy plateaus at S(T_mid) ~ 1/2 ln 2 within 5.0% tolerance."""
    temps = np.logspace(-2.5, 0.8, 35)
    res = heat_calc.compute_thermodynamic_curve(temps)

    s_target = 0.5 * np.log(2.0)
    rel_err = abs(res.s_mid - s_target) / s_target

    assert rel_err <= 0.05, (
        f"Expected intermediate plateau entropy within 5% of 1/2 ln 2 ({s_target:.5f}), "
        f"got {res.s_mid:.5f} (rel error {rel_err * 100:.2f}%)"
    )


def test_t2_4_1_high_temperature_entropy_saturation(heat_calc: SpecificHeatCalculator) -> None:
    """Verify that asymptotic entropy saturates to ln 2 within 2.0% tolerance."""
    temps = np.logspace(-2.5, 0.8, 35)
    res = heat_calc.compute_thermodynamic_curve(temps)

    s_target = np.log(2.0)
    rel_err = abs(res.s_max - s_target) / s_target

    assert rel_err <= 0.02, (
        f"Expected high-T entropy within 2% of ln 2 ({s_target:.5f}), "
        f"got {res.s_max:.5f} (rel error {rel_err * 100:.2f}%)"
    )


def test_t2_4_1_bayesian_model_selection() -> None:
    """Verify that Bayesian model comparison decisively favors two-peak over single-peak model."""
    benchmark_path = Path("data/processed/benchmarks/m2_4/specific_heat_benchmark.dat")
    assert benchmark_path.exists()

    data = np.loadtxt(benchmark_path)
    temps = data[:, 0]
    cv_total = data[:, 3]

    fitter = ThermodynamicBayesianFitter()
    bayes_res = fitter.fit_and_compare(temps, cv_total)

    assert bayes_res.n_samples == 40
    assert bayes_res.delta_bic > 10.0, f"Expected Delta_BIC > 10.0, got {bayes_res.delta_bic:.2f}"
    assert bayes_res.bayes_factor > 148.0, f"Expected BF > 148.0, got {bayes_res.bayes_factor:.2e}"
    assert bayes_res.favors_two_peak is True


def test_t2_4_1_synthetic_bayesian_fitter() -> None:
    """Verify ThermodynamicBayesianFitter on an exact synthetic double-Schottky model."""
    temps = np.logspace(-2.5, 0.8, 30)
    # Synthetic double Schottky with well-separated gaps D1=0.08, D2=1.2
    x1 = 0.08 / temps
    cv1 = 0.25 * (x1**2) * np.exp(-x1) / ((1.0 + np.exp(-x1)) ** 2)
    x2 = 1.20 / temps
    cv2 = 0.20 * (x2**2) * np.exp(-x2) / ((1.0 + np.exp(-x2)) ** 2)
    synthetic_cv = cv1 + cv2

    fitter = ThermodynamicBayesianFitter()
    bayes_res = fitter.fit_and_compare(temps, synthetic_cv)

    assert bayes_res.favors_two_peak is True
    assert bayes_res.delta_bic > 50.0


def test_t2_4_1_benchmark_dataset_acid_provenance() -> None:
    """Verify ACID provenance header and formatting of the specific heat dataset."""
    benchmark_path = Path("data/processed/benchmarks/m2_4/specific_heat_benchmark.dat")
    assert benchmark_path.exists()

    with open(benchmark_path, encoding="utf-8") as f:
        lines = f.readlines()

    header_lines = [line for line in lines if line.strip().startswith("#")]
    data_lines = [line for line in lines if not line.strip().startswith("#") and line.strip()]

    assert len(header_lines) >= 6
    assert any("Git Commit:" in line for line in header_lines)
    assert any("Generated:" in line for line in header_lines)
    assert any("Bayesian Model Selection Summary:" in line for line in header_lines)
    assert any("Columns:" in line for line in header_lines)
    assert len(data_lines) == 40

    for dl in data_lines:
        tokens = dl.strip().split()
        assert len(tokens) == 7, f"Expected 7 columns, got {len(tokens)}: {dl}"
