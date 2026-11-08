"""Unit tests for Task T2.2.1: Vortex Separation Potential V(R) & Topological Deconfinement."""

from pathlib import Path

import numpy as np
import pytest

from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.phase_2.vortex_potential import (
    DeconfinementBayesianFitter,
    VortexSeparationPotentialCalculator,
)


@pytest.fixture
def kitaev_params() -> ModelParameters:
    return ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))


def test_t2_2_1_vortex_gauge_flux_pinning(kitaev_params: ModelParameters) -> None:
    """Verify that inverted branch cut strictly pins W_p = -1 at p1, p2 and W_p = +1 elsewhere."""
    calc = VortexSeparationPotentialCalculator(L=5)
    p1 = (0, 0)
    p2 = (2, 0)

    point = calc.evaluate_vortex_pair(p1, p2, kitaev_params)

    # Core flux verification
    assert np.isclose(point.w_p1, -1.0, atol=1e-6)
    assert np.isclose(point.w_p2, -1.0, atol=1e-6)
    assert point.cut_length == 2
    assert point.R > 0.0

    # Anti-vacuity check: verify zero flux leakage onto remaining L^2 - 2 plaquettes
    cut = calc.geom.generate_vortex_pair_cut(p1, p2)
    Av = calc.geom.assemble_skew_hamiltonian(kitaev_params, calc.qsl_seed, flipped_bonds=cut)
    omegav, _ = calc.geom.solve_covariance_real_eigh(Av)
    fluxes = calc.geom.compute_plaquette_fluxes(omegav)

    vortex_plaquettes = {k for k, v in fluxes.items() if v < 0.0}
    flux_free_plaquettes = {k for k, v in fluxes.items() if v > 0.9}

    assert vortex_plaquettes == {p1, p2}
    assert len(flux_free_plaquettes) == calc.geom.num_cells - 2


def test_t2_2_1_vortex_pair_energy_positivity(kitaev_params: ModelParameters) -> None:
    """Verify that isolated vortex pairs strictly cost positive physical energy (V(R) > 0)."""
    calc = VortexSeparationPotentialCalculator(L=5)
    points = calc.scan_separation_profile(kitaev_params)

    assert len(points) >= 4
    for pt in points:
        assert pt.V_R > 0.0, f"Expected V(R) > 0, got {pt.V_R} at R={pt.R}"
        assert pt.E_2vortex > pt.E_flux_free


def test_t2_2_1_vison_deconfinement_bayesian_model_selection(
    kitaev_params: ModelParameters,
) -> None:
    """Verify that Bayesian model selection decisively rejects linear confinement (Delta_BIC > 10, BF > 148)."""
    calc = VortexSeparationPotentialCalculator(L=7)
    points = calc.scan_separation_profile(kitaev_params)

    r_arr = np.array([pt.R for pt in points])
    v_arr = np.array([pt.V_R for pt in points])
    cuts = np.array([pt.cut_length for pt in points])

    fitter = DeconfinementBayesianFitter()
    res = fitter.fit_and_compare(r_arr, v_arr, cut_lengths=cuts)

    assert res.n_samples >= 7
    assert res.delta_bic > 10.0, f"Expected Delta_BIC > 10, got {res.delta_bic}"
    assert res.bayes_factor > 148.0, f"Expected Bayes factor > 148, got {res.bayes_factor}"
    assert res.favors_deconfinement is True


def test_t2_2_1_kitaev_analytical_vison_gap_agreement(kitaev_params: ModelParameters) -> None:
    """Verify that the numerical asymptotic vison gap matches Kitaev's analytical value within 10%."""
    calc = VortexSeparationPotentialCalculator(L=7)
    points = calc.scan_separation_profile(kitaev_params)

    r_arr = np.array([pt.R for pt in points])
    v_arr = np.array([pt.V_R for pt in points])
    cuts = np.array([pt.cut_length for pt in points])

    fitter = DeconfinementBayesianFitter()
    res = fitter.fit_and_compare(r_arr, v_arr, cut_lengths=cuts)

    # Kitaev exact analytical value: Delta_v = 0.1536
    assert res.vison_gap_relative_error < 0.10, (
        f"Vison gap {res.vison_gap:.4f} deviated by {res.vison_gap_relative_error * 100:.1f}% "
        f"from analytical 0.1536"
    )


def test_t2_2_1_unconstrained_effective_string_tension_vanishing(
    kitaev_params: ModelParameters,
) -> None:
    """Verify that unconstrained linear string tension is vanishingly small (sigma <= 0.05)."""
    calc = VortexSeparationPotentialCalculator(L=7)
    points = calc.scan_separation_profile(kitaev_params)

    r_arr = np.array([pt.R for pt in points])
    v_arr = np.array([pt.V_R for pt in points])
    cuts = np.array([pt.cut_length for pt in points])

    fitter = DeconfinementBayesianFitter()
    res = fitter.fit_and_compare(r_arr, v_arr, cut_lengths=cuts)

    # Bare exchange energy per bond is ~1.05; screened string tension should be <= 0.05
    assert abs(res.unconstrained_sigma) <= 0.05, (
        f"Unconstrained string tension {res.unconstrained_sigma:.4f} exceeds 0.05 threshold"
    )


def test_t2_2_1_benchmark_dataset_acid_provenance() -> None:
    """Verify existence, header metadata provenance, and numerical integrity of benchmark table."""
    benchmark_path = Path("data/processed/benchmarks/m2_2/vortex_potential_benchmark.dat")
    assert benchmark_path.exists(), f"Benchmark file {benchmark_path} does not exist"

    with open(benchmark_path, encoding="utf-8") as f:
        lines = f.readlines()

    header_lines = [line for line in lines if line.strip().startswith("#")]
    data_lines = [line for line in lines if not line.strip().startswith("#") and line.strip()]

    assert len(header_lines) >= 8
    assert any("Git Commit:" in line for line in header_lines)
    assert any("Generated:" in line for line in header_lines)
    assert any("Theoretical Exact Vison Gap:" in line for line in header_lines)
    assert len(data_lines) >= 20

    for dl in data_lines:
        tokens = dl.strip().split()
        assert len(tokens) == 13
        # Tokens: L, h, p1_1, p1_2, p2_1, p2_2, ncut, R, E_0, E_2v, V_R, w_p1, w_p2
        L = int(tokens[0])
        h = float(tokens[1])
        R = float(tokens[7])
        V_R = float(tokens[10])
        w_p1 = float(tokens[11])
        w_p2 = float(tokens[12])

        assert L in (5, 7, 8)
        assert h in (0.0, 0.05)
        assert R > 0.0
        assert V_R > 0.0
        assert np.isclose(w_p1, -1.0, atol=1e-5)
        assert np.isclose(w_p2, -1.0, atol=1e-5)
