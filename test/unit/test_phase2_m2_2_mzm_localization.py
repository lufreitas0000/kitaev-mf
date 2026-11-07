"""Unit tests for Task T2.2.2: Majorana Zero Modes Bound to Isolated Vortices."""

from pathlib import Path
import numpy as np
import pytest

from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.phase_2.mzm_localization import (
    MZMLocalizationCalculator,
    MZMSplittingBayesianFitter,
)


@pytest.fixture
def topological_params() -> ModelParameters:
    hx = hy = hz = 0.05 / np.sqrt(3.0)
    return ModelParameters(K=-1.0, h=(hx, hy, hz))


def test_t2_2_2_mzm_in_gap_spectral_isolation(topological_params: ModelParameters) -> None:
    """Verify that the vortex-bound MZM mode is isolated deep inside the bulk spectral gap."""
    calc = MZMLocalizationCalculator(L=8)
    metrics = calc.evaluate_mzm_pair(topological_params)

    assert metrics.eps_bulk > 0.5, f"Expected non-zero bulk gap, got {metrics.eps_bulk}"
    assert metrics.eps_mzm > 0.0
    assert metrics.in_gap_ratio <= 0.25, (
        f"In-gap ratio {metrics.in_gap_ratio:.4f} exceeded 0.25 threshold"
    )


def test_t2_2_2_mzm_core_localization_enrichment(topological_params: ModelParameters) -> None:
    """Verify that MZM wavefunction density is enriched by at least 2.0x on vortex core boundary sites."""
    for L in (5, 7, 8):
        calc = MZMLocalizationCalculator(L=L)
        metrics = calc.evaluate_mzm_pair(topological_params)

        assert metrics.enrichment_ratio >= 2.0, (
            f"Cluster L={L} core enrichment ratio {metrics.enrichment_ratio:.2f}x is below 2.0x"
        )
        assert metrics.w_core > metrics.w_uniform


def test_t2_2_2_mzm_hop_distance_decay(topological_params: ModelParameters) -> None:
    """Verify that MZM probability density is strongly localized at d_hop=0 and decays into the bulk."""
    calc = MZMLocalizationCalculator(L=8)
    metrics = calc.evaluate_mzm_pair(topological_params)

    p0 = metrics.hop_profile.get(0, 0.0)
    p_far = metrics.hop_profile.get(max(metrics.hop_profile.keys()), 0.0)

    assert p0 > 0.02, f"Expected core density p0 > 0.02, got {p0:.4f}"
    assert p0 > 3.0 * p_far, f"Core density {p0:.4f} should exceed 3x far bulk density {p_far:.4f}"


def test_t2_2_2_mzm_bayesian_exponential_decoupling(topological_params: ModelParameters) -> None:
    """Verify that Bayesian model selection decisively favors exponential MZM decoupling (Delta_BIC > 10)."""
    clusters = [4, 5, 7, 8]
    eps_list = []

    for L in clusters:
        calc = MZMLocalizationCalculator(L=L)
        metrics = calc.evaluate_mzm_pair(topological_params)
        eps_list.append(metrics.eps_mzm)

    fitter = MZMSplittingBayesianFitter()
    res = fitter.fit_and_compare(np.array(clusters, dtype=float), np.array(eps_list))

    assert res.n_samples == 4
    assert res.delta_bic > 10.0, f"Expected Delta_BIC > 10, got {res.delta_bic}"
    assert res.bayes_factor > 148.0, f"Expected Bayes factor > 148, got {res.bayes_factor}"
    assert res.favors_exponential_decay is True
    assert 2.0 < res.decay_length_xi < 10.0, f"Unphysical decay length: {res.decay_length_xi}"


def test_t2_2_2_mzm_bayesian_synthetic_decay() -> None:
    """Verify MZMSplittingBayesianFitter with an exact analytical exponential sequence."""
    l_vals = np.array([4.0, 6.0, 8.0, 10.0])
    xi_true = 3.5
    eps_synthetic = 0.5 * np.exp(-l_vals / xi_true)

    fitter = MZMSplittingBayesianFitter()
    res = fitter.fit_and_compare(l_vals, eps_synthetic)

    assert res.favors_exponential_decay is True
    assert res.delta_bic > 10.0
    assert np.isclose(res.decay_length_xi, xi_true, rtol=1e-3)


def test_t2_2_2_benchmark_dataset_acid_provenance() -> None:
    """Verify existence, header metadata provenance, and numerical integrity of benchmark table."""
    benchmark_path = Path("data/processed/benchmarks/m2_2/mzm_localization_benchmark.dat")
    assert benchmark_path.exists(), f"Benchmark file {benchmark_path} does not exist"

    with open(benchmark_path, encoding="utf-8") as f:
        lines = f.readlines()

    header_lines = [line for line in lines if line.strip().startswith("#")]
    data_lines = [line for line in lines if not line.strip().startswith("#") and line.strip()]

    assert len(header_lines) >= 8
    assert any("Git Commit:" in line for line in header_lines)
    assert any("Generated:" in line for line in header_lines)
    assert any("Bayesian Model Selection Summary" in line for line in header_lines)
    assert len(data_lines) == 12  # 4 clusters x 3 fields

    for dl in data_lines:
        tokens = dl.strip().split()
        assert len(tokens) == 12
        # Tokens: L, h, R, eps_MZM, eps_bulk, in_gap_rat, W_core, W_unif, ratio, P_d0, P_d1, P_d2
        L = int(tokens[0])
        h = float(tokens[1])
        eps_mzm = float(tokens[3])
        eps_bulk = float(tokens[4])
        in_gap_ratio = float(tokens[5])
        ratio = float(tokens[8])

        assert L in (4, 5, 7, 8)
        assert h in (0.02, 0.05, 0.10)
        assert eps_mzm > 0.0
        assert eps_bulk > 0.5
        assert in_gap_ratio <= 0.28
        assert ratio >= 2.0
