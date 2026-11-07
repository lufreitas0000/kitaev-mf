"""Unit tests for Phase 3 Task T3.W1.3: Thermodynamic Spectral Sum Rules up to L=37.

Validates Hypothesis Contracts:
- H_SPEC_SUMRULE_01: Sum rule saturation W_tot >= 0.45 across all clusters up to L=37.
- H_SPEC_SUMRULE_02: Universal excitation energy positivity (omega >= -0.02 |K|).
- H_SPEC_SUMRULE_03: Persistent MZM vs. bulk demarcation gap Delta omega >= 0.15 |K| up to L=37.
- H_SPEC_SUMRULE_04: Continuous quadrature normalization error < 1e-3.
"""

import pytest
import numpy as np

from kitaev_mf.services.phase_3.thermodynamic_spectral_benchmarks import (
    benchmark_thermodynamic_spectral_response,
)


@pytest.fixture(scope="module")
def spectral_benchmark_data():
    """Compute benchmark metrics for representative clusters."""
    # Test L in {6, 12, 18, 24, 30, 36, 37}
    return benchmark_thermodynamic_spectral_response(L_values=[6, 12, 18, 24, 30, 36, 37])


def test_h_spec_sumrule_01_saturation(spectral_benchmark_data):
    """Verify H_SPEC_SUMRULE_01: sum rule saturation W_tot >= 0.45 up to L=37."""
    for d in spectral_benchmark_data:
        assert d.core_sum_rule >= 0.45, (
            f"L={d.L}: Expected core sum rule >= 0.45, got {d.core_sum_rule:.4f}"
        )
        assert d.bulk_sum_rule >= 0.60, (
            f"L={d.L}: Expected bulk sum rule >= 0.60, got {d.bulk_sum_rule:.4f}"
        )
        if not d.is_singular_core:
            assert d.core_sum_rule >= 0.70, (
                f"L={d.L}: Regular core quench should achieve >= 0.70, got {d.core_sum_rule:.4f}"
            )
        if not d.is_singular_bulk:
            assert d.bulk_sum_rule >= 0.70, (
                f"L={d.L}: Regular bulk quench should achieve >= 0.70, got {d.bulk_sum_rule:.4f}"
            )


def test_h_spec_sumrule_02_energy_positivity(spectral_benchmark_data):
    """Verify H_SPEC_SUMRULE_02: physical excitation energies satisfy physical bounds (>= -0.02 |K|)."""
    for d in spectral_benchmark_data:
        assert d.omega_MZM >= -0.02, (
            f"L={d.L}: Expected omega_MZM >= -0.02 |K|, got {d.omega_MZM:.4f}"
        )
        if not d.is_singular_core:
            assert d.omega_MZM >= 0.0, (
                f"L={d.L}: Regular core quench should satisfy omega_MZM >= 0, got {d.omega_MZM:.4f}"
            )
        assert d.bulk_min_energy >= 0.15, (
            f"L={d.L}: Expected bulk minimum excitation >= 0.15 |K|, got {d.bulk_min_energy:.4f}"
        )


def test_h_spec_sumrule_03_demarcation_gap(spectral_benchmark_data):
    """Verify H_SPEC_SUMRULE_03: persistent demarcation gap Delta omega >= 0.15 |K| up to L=37."""
    for d in spectral_benchmark_data:
        assert d.demarcation_gap >= 0.15, (
            f"L={d.L}: Expected demarcation gap >= 0.15 |K|, got {d.demarcation_gap:.4f}"
        )


def test_h_spec_sumrule_04_continuous_quadrature(spectral_benchmark_data):
    """Verify H_SPEC_SUMRULE_04: continuous quadrature error < 1e-3 across all sizes."""
    for d in spectral_benchmark_data:
        assert d.quadrature_error < 1e-3, (
            f"L={d.L}: Expected continuous quadrature error < 1e-3, got {d.quadrature_error:.4e}"
        )
