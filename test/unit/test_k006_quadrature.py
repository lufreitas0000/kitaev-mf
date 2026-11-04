"""Contracts for reproducible K-006 Brillouin-zone midpoint quadrature."""

import numpy as np
import pytest

from kitaev_mf.services.quadrature import (
    midpoint_convergence,
    midpoint_form_factor_average,
)


def _direct_midpoint_average(mesh_size: int) -> float:
    points = (np.arange(mesh_size, dtype=np.float64) + 0.5) / mesh_size
    total = 0.0
    for k1 in points:
        for k2 in points:
            form_factor = (
                1.0
                + np.exp(-2.0j * np.pi * k1)
                + np.exp(-2.0j * np.pi * k2)
            )
            total += abs(form_factor)
    return total / (3.0 * mesh_size**2)


def test_k006_midpoint_matches_direct_small_mesh_sum() -> None:
    assert midpoint_form_factor_average(12, row_chunk_size=3) == pytest.approx(
        _direct_midpoint_average(12), abs=2.0e-15
    )


def test_k006_midpoint_is_independent_of_row_chunk_size() -> None:
    small_chunk = midpoint_form_factor_average(96, row_chunk_size=1)
    large_chunk = midpoint_form_factor_average(96, row_chunk_size=32)
    assert small_chunk == pytest.approx(large_chunk, abs=2.0e-15)


@pytest.mark.parametrize("mesh_size", [0, -3, 4, 10, 97])
def test_k006_midpoint_rejects_invalid_mesh_sizes(mesh_size: int) -> None:
    with pytest.raises(ValueError):
        midpoint_form_factor_average(mesh_size)


def test_k006_convergence_reports_cubic_order_and_error_estimate() -> None:
    estimates = midpoint_convergence((96, 192, 384, 768))
    middle = estimates[1]
    assert middle.observed_order == pytest.approx(3.0, abs=2.0e-3)
    assert middle.estimated_error is not None
    assert middle.estimated_error > 0.0
    assert middle.richardson_value is not None
    assert middle.richardson_value == pytest.approx(
        0.52486574585, abs=2.0e-10
    )


def test_k006_convergence_requires_doubling_meshes() -> None:
    with pytest.raises(ValueError):
        midpoint_convergence((96, 192, 576))
