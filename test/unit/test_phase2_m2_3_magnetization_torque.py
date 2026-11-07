"""Unit tests for Task T2.3.1: Global Field-Directional Magnetization Surface & Torque Discontinuities."""

from pathlib import Path
import numpy as np
import pytest

from kitaev_mf.services.phase_2.magnetization_torque import (
    MagnetizationTorqueCalculator,
    TorqueBayesianFitter,
)


@pytest.fixture
def torque_calc() -> MagnetizationTorqueCalculator:
    return MagnetizationTorqueCalculator(k_grid=(12, 12))


def test_t2_3_1_high_symmetry_nodes(torque_calc: MagnetizationTorqueCalculator) -> None:
    """Verify that magnetic torque vanishes at high-symmetry axes [001], [111], and [110]."""
    # High-symmetry angles in degrees
    # [001]: 0 deg
    # [111]: arccos(1/sqrt(3)) ~ 54.7356 deg
    # [110]: 90 deg
    symm_angles = [0.0, float(np.degrees(np.arccos(1.0 / np.sqrt(3.0)))), 90.0]
    h_mag = 0.20

    for th_deg in symm_angles:
        pt, _ = torque_calc.evaluate_point(h_mag, th_deg)
        assert pt.tau_norm < 1e-5, f"Torque should vanish at symmetry angle {th_deg} deg, got {pt.tau_norm}"
        assert abs(pt.tau_projection) < 1e-5


def test_t2_3_1_torque_sign_inversion(torque_calc: MagnetizationTorqueCalculator) -> None:
    """Verify that projected torque inverts sign across the [111] nodal axis."""
    th_111 = float(np.degrees(np.arccos(1.0 / np.sqrt(3.0))))
    h_mag = 0.35

    pt_before, _ = torque_calc.evaluate_point(h_mag, th_111 - 15.0)
    pt_after, _ = torque_calc.evaluate_point(h_mag, th_111 + 15.0)

    # In [001] to [111] quadrant, tau_proj > 0; in [111] to [110] quadrant, tau_proj < 0
    assert pt_before.tau_projection > 5e-4, f"Expected positive torque before [111], got {pt_before.tau_projection}"
    assert pt_after.tau_projection < -5e-4, f"Expected negative torque after [111], got {pt_after.tau_projection}"


def test_t2_3_1_bayesian_model_selection() -> None:
    """Verify that Bayesian model comparison decisively favors anisotropic torque (Delta_BIC > 10)."""
    benchmark_path = Path("data/processed/benchmarks/m2_3/magnetization_torque_benchmark.dat")
    assert benchmark_path.exists()

    data = np.loadtxt(benchmark_path)
    h_arr = data[:, 0]
    th_arr = data[:, 2]
    tau_arr = data[:, 4]

    fitter = TorqueBayesianFitter()
    bayes_res = fitter.fit_and_compare(th_arr, tau_arr, h_magnitudes=h_arr)

    assert bayes_res.n_samples == 39
    assert bayes_res.delta_bic > 10.0, f"Expected Delta_BIC > 10, got {bayes_res.delta_bic}"
    assert bayes_res.bayes_factor > 148.0, f"Expected Bayes factor > 148, got {bayes_res.bayes_factor}"
    assert bayes_res.favors_anisotropic_torque is True


def test_t2_3_1_synthetic_bayesian_fitter() -> None:
    """Verify TorqueBayesianFitter on synthetic data with known node structure."""
    thetas = np.linspace(0.0, np.pi / 2.0, 13)
    basis = np.sin(2.0 * thetas) * (3.0 * np.cos(thetas)**2 - 1.0)
    synthetic_tau = 0.002 * basis

    fitter = TorqueBayesianFitter()
    res = fitter.fit_and_compare(thetas, synthetic_tau)

    assert res.favors_anisotropic_torque is True
    assert res.delta_bic > 20.0
    assert np.isclose(res.amplitude_c, 0.002, rtol=1e-3)


def test_t2_3_1_benchmark_dataset_acid_provenance() -> None:
    """Verify ACID provenance header and formatting of the benchmark dataset."""
    benchmark_path = Path("data/processed/benchmarks/m2_3/magnetization_torque_benchmark.dat")
    assert benchmark_path.exists()

    with open(benchmark_path, encoding="utf-8") as f:
        lines = f.readlines()

    header_lines = [l for l in lines if l.strip().startswith("#")]
    data_lines = [l for l in lines if not l.strip().startswith("#") and l.strip()]

    assert len(header_lines) >= 6
    assert any("Git Commit:" in l for l in header_lines)
    assert any("Generated:" in l for l in header_lines)
    assert any("Bayesian Model Selection Summary:" in l for l in header_lines)
    assert len(data_lines) == 39  # 3 fields x 13 angles

    for dl in data_lines:
        tokens = dl.strip().split()
        assert len(tokens) == 7
        h_val = float(tokens[0])
        th_deg = float(tokens[1])
        m_norm = float(tokens[3])
        assert h_val in (0.10, 0.20, 0.35)
        assert 0.0 <= th_deg <= 90.0
        assert m_norm > 0.0
