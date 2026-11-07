"""Unit tests for Task T2.2.3: Vortex Braiding Phases & Non-Abelian Holonomy."""

from pathlib import Path

import numpy as np
import pytest

from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.phase_2.vortex_braiding import (
    BraidingBayesianFitter,
    VortexBraidingCalculator,
)


@pytest.fixture
def topological_params() -> ModelParameters:
    hx = hy = hz = 0.05 / np.sqrt(3.0)
    return ModelParameters(K=-1.0, h=(hx, hy, hz))


def test_t2_2_3_non_enclosing_loop_trivial_berry_phase(topological_params: ModelParameters) -> None:
    """Verify that transporting a vortex around an empty region yields trivial Berry phase (0 mod 2pi)."""
    calc = VortexBraidingCalculator(L=6)
    res = calc.evaluate_braiding_holonomy(topological_params, N_sub=4)

    # In radians, 0 mod 2pi means either near 0 or near 2*pi
    dist_to_zero = min(res.gamma_non_enclosing, 2.0 * np.pi - res.gamma_non_enclosing)
    assert dist_to_zero < 0.05, (
        f"Expected trivial Berry phase near 0, got {res.gamma_non_enclosing:.4f}"
    )


def test_t2_2_3_enclosing_loop_pi_holonomy(topological_params: ModelParameters) -> None:
    """Verify that winding around a vortex imparts topological pi statistical holonomy."""
    calc = VortexBraidingCalculator(L=6)
    res = calc.evaluate_braiding_holonomy(topological_params, N_sub=4)

    # Both independent phases must be physically meaningful:
    dist_non_to_zero = min(res.gamma_non_enclosing, 2.0 * np.pi - res.gamma_non_enclosing)
    assert dist_non_to_zero < 0.05, f"gamma_non should be near 0, got {res.gamma_non_enclosing:.4f}"

    dist_enc_to_pi = min(abs(res.gamma_enclosing - np.pi), abs(res.gamma_enclosing - 3.0 * np.pi))
    assert dist_enc_to_pi < 0.50, f"gamma_enc should be near pi, got {res.gamma_enclosing:.4f}"

    assert res.relative_error_from_pi <= 0.18, (
        f"L=6: Holonomy error {res.relative_error_from_pi * 100:.1f}% exceeded 18% threshold "
        f"(Delta_gamma = {res.delta_gamma:.4f} rad)"
    )


def test_t2_2_3_braiding_bayesian_model_selection() -> None:
    """Verify that Bayesian model selection decisively favors pi holonomy over trivial 0 (Delta_BIC > 10)."""
    # Load all benchmark points from dataset
    benchmark_path = Path("data/processed/benchmarks/m2_2/vortex_braiding_benchmark.dat")
    assert benchmark_path.exists()

    diffs = []
    with open(benchmark_path, encoding="utf-8") as f:
        for line in f:
            if not line.strip().startswith("#") and line.strip():
                tokens = line.strip().split()
                diffs.append(float(tokens[4]))

    fitter = BraidingBayesianFitter()
    bayes_res = fitter.fit_and_compare(np.array(diffs))

    assert bayes_res.n_samples == 4
    assert bayes_res.delta_bic > 10.0, f"Expected Delta_BIC > 10, got {bayes_res.delta_bic}"
    assert bayes_res.bayes_factor > 148.0, (
        f"Expected Bayes factor > 148, got {bayes_res.bayes_factor}"
    )
    assert bayes_res.favors_pi_holonomy is True


def test_t2_2_3_synthetic_bayesian_fitter() -> None:
    """Verify BraidingBayesianFitter on an exact synthetic distribution near pi."""
    synthetic_diffs = np.array([3.14, 3.15, 3.13, 3.16])
    fitter = BraidingBayesianFitter()
    res = fitter.fit_and_compare(synthetic_diffs)

    assert res.favors_pi_holonomy is True
    assert res.delta_bic > 20.0
    assert np.isclose(res.mean_delta_gamma, np.pi, atol=0.02)


def test_t2_2_3_benchmark_dataset_acid_provenance() -> None:
    """Verify existence, header metadata provenance, and numerical integrity of benchmark table."""
    benchmark_path = Path("data/processed/benchmarks/m2_2/vortex_braiding_benchmark.dat")
    assert benchmark_path.exists(), f"Benchmark file {benchmark_path} does not exist"

    with open(benchmark_path, encoding="utf-8") as f:
        lines = f.readlines()

    header_lines = [line for line in lines if line.strip().startswith("#")]
    data_lines = [line for line in lines if not line.strip().startswith("#") and line.strip()]

    assert len(header_lines) >= 6
    assert any("Git Commit:" in line for line in header_lines)
    assert any("Generated:" in line for line in header_lines)
    assert any("Bayesian Model Selection Summary:" in line for line in header_lines)
    assert len(data_lines) >= 4  # Resonant cluster L=6 across fields h in [0.02, 0.10]

    for dl in data_lines:
        tokens = dl.strip().split()
        assert len(tokens) == 7
        # Tokens: L, h, gamma_enc, gamma_non, Delta_gamma, Delta_g/pi, rel_error
        L = int(tokens[0])
        h = float(tokens[1])
        delta_gamma = float(tokens[4])
        rel_error = float(tokens[6])

        assert L in (5, 6, 7)
        assert h in (0.02, 0.03, 0.05, 0.10)
        if L >= 6:
            assert 2.5 <= delta_gamma <= 3.8
            assert rel_error <= 0.20
