"""Unit tests for Task T2.3.3: Global Phase Diagram Mapping & First-Order Boundary Anisotropy."""

from pathlib import Path

import numpy as np
import pytest

from kitaev_mf.services.phase_2.phase_diagram import (
    PhaseBoundaryBayesianFitter,
    PhaseBoundaryCalculator,
)


@pytest.fixture
def boundary_calc() -> PhaseBoundaryCalculator:
    return PhaseBoundaryCalculator(k_grid=(12, 12))


def test_t2_3_3_critical_field_bracket(boundary_calc: PhaseBoundaryCalculator) -> None:
    """Verify that the first-order thermodynamic crossing h_c* lies inside [0.330, 0.335]."""
    for th_deg in (0.0, 54.7356, 90.0):
        pt = boundary_calc.find_critical_field(th_deg)
        assert 0.330 <= pt.h_c_star <= 0.335, (
            f"Expected h_c* in [0.330, 0.335] at {th_deg} deg, got {pt.h_c_star:.5f}"
        )
        assert pt.qsl_energy_330 < pt.pol_energy_330  # QSL is ground state at h=0.330
        assert pt.qsl_energy_340 > pt.pol_energy_340  # Polarized is ground state at h=0.340


def test_t2_3_3_maximum_stability_at_111() -> None:
    """Verify that the critical field is maximized near [111] (theta ~ 54.74 deg)."""
    benchmark_path = Path("data/processed/benchmarks/m2_3/phase_boundary_benchmark.dat")
    assert benchmark_path.exists()

    data = np.loadtxt(benchmark_path)
    thetas_deg = data[:, 0]
    h_c_values = data[:, 2]

    fitter = PhaseBoundaryBayesianFitter()
    res = fitter.fit_and_compare(thetas_deg, h_c_values)

    # Maximum must be within 3.0 degrees of 54.74 deg (grid resolution is 7.5 deg)
    assert abs(res.theta_max_deg - 54.74) <= 3.0, (
        f"Expected max at 54.74 deg, got {res.theta_max_deg:.2f} deg"
    )
    assert res.delta_h > 0.0, (
        "Expected positive Delta_h indicating enhanced QSL stability along [111]"
    )


def test_t2_3_3_bayesian_model_selection() -> None:
    """Verify that Bayesian model comparison decisively favors anisotropic boundary modulation."""
    benchmark_path = Path("data/processed/benchmarks/m2_3/phase_boundary_benchmark.dat")
    assert benchmark_path.exists()

    data = np.loadtxt(benchmark_path)
    thetas_deg = data[:, 0]
    h_c_values = data[:, 2]

    fitter = PhaseBoundaryBayesianFitter()
    res = fitter.fit_and_compare(thetas_deg, h_c_values)

    assert res.n_samples == 13
    assert res.delta_bic > 10.0, f"Expected Delta_BIC > 10, got {res.delta_bic:.2f}"
    assert res.bayes_factor > 148.0, f"Expected Bayes factor > 148, got {res.bayes_factor:.2e}"
    assert res.favors_anisotropic_boundary is True


def test_t2_3_3_synthetic_bayesian_fitter() -> None:
    """Verify PhaseBoundaryBayesianFitter on an exact synthetic modulated boundary."""
    thetas = np.linspace(0.0, 90.0, 13)
    th_rad = np.radians(thetas)
    basis = (np.sin(th_rad) ** 4 * np.cos(th_rad) ** 2) / (4.0 / 27.0)
    synthetic_h_c = 0.33150 + 0.00100 * basis

    fitter = PhaseBoundaryBayesianFitter()
    res = fitter.fit_and_compare(thetas, synthetic_h_c)

    assert res.favors_anisotropic_boundary is True
    assert np.isclose(res.h_base, 0.33150, atol=1e-5)
    assert np.isclose(res.delta_h, 0.00100, atol=1e-5)
    assert res.delta_bic > 20.0


def test_t2_3_3_benchmark_dataset_acid_provenance() -> None:
    """Verify ACID provenance header and formatting of the phase boundary dataset."""
    benchmark_path = Path("data/processed/benchmarks/m2_3/phase_boundary_benchmark.dat")
    assert benchmark_path.exists()

    with open(benchmark_path, encoding="utf-8") as f:
        lines = f.readlines()

    header_lines = [line for line in lines if line.strip().startswith("#")]
    data_lines = [line for line in lines if not line.strip().startswith("#") and line.strip()]

    assert len(header_lines) >= 6
    assert any("Git Commit:" in line for line in header_lines)
    assert any("Generated:" in line for line in header_lines)
    assert any("Bayesian Model Selection Summary:" in line for line in header_lines)
    assert len(data_lines) == 13

    for dl in data_lines:
        tokens = dl.strip().split()
        assert len(tokens) == 8
        th_deg = float(tokens[0])
        h_c = float(tokens[2])
        assert 0.0 <= th_deg <= 90.0
        assert 0.330 <= h_c <= 0.335
