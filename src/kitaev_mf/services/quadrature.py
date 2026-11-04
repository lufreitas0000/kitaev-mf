"""Deterministic midpoint quadrature for the pure-Kitaev matter correlator."""

from collections.abc import Sequence
from dataclasses import dataclass
from math import log2

import numpy as np


@dataclass(frozen=True)
class MidpointEstimate:
    """Midpoint value and optional three-mesh Richardson diagnostics."""

    mesh_size: int
    value: float
    observed_order: float | None = None
    estimated_error: float | None = None
    richardson_value: float | None = None


def midpoint_form_factor_average(
    mesh_size: int,
    row_chunk_size: int = 128,
) -> float:
    """Return (1/3) times the midpoint BZ average of |1 + exp(-ik1) + exp(-ik2)|.

    Fractional momenta are measured in reciprocal primitive coordinates, so
    each phase contains 2*pi. Meshes divisible by three avoid sampling the
    Dirac points on the half-shifted midpoint grid.
    """
    if mesh_size <= 0 or mesh_size % 3 != 0:
        raise ValueError("mesh_size must be a positive multiple of three")
    if row_chunk_size <= 0:
        raise ValueError("row_chunk_size must be positive")

    coordinates = (np.arange(mesh_size, dtype=np.float64) + 0.5) / mesh_size
    phases = 2.0 * np.pi * coordinates
    phase_k2 = np.exp(-1.0j * phases)
    total = 0.0

    for start in range(0, mesh_size, row_chunk_size):
        stop = min(start + row_chunk_size, mesh_size)
        phase_k1 = np.exp(-1.0j * phases[start:stop])
        form_factor = 1.0 + phase_k1[:, None] + phase_k2[None, :]
        total += float(np.sum(np.abs(form_factor), dtype=np.float64))

    return total / (3.0 * mesh_size * mesh_size)


def midpoint_convergence(
    mesh_sizes: Sequence[int],
    row_chunk_size: int = 128,
) -> tuple[MidpointEstimate, ...]:
    """Evaluate a doubling sequence and attach empirical Richardson data.

    For three consecutive meshes N/2, N, 2N, estimate p from the ratio of
    successive differences. The error estimate and extrapolated value assume
    asymptotic O(N**-p) behavior; they are diagnostics, not rigorous bounds.
    """
    meshes = tuple(mesh_sizes)
    if len(meshes) < 3:
        raise ValueError("at least three mesh sizes are required")
    if any(size <= 0 or size % 3 != 0 for size in meshes):
        raise ValueError("all mesh sizes must be positive multiples of three")
    if any(right != 2 * left for left, right in zip(meshes, meshes[1:], strict=False)):
        raise ValueError("mesh sizes must form an exact doubling sequence")

    values = tuple(midpoint_form_factor_average(size, row_chunk_size) for size in meshes)
    estimates = [
        MidpointEstimate(mesh_size=size, value=value)
        for size, value in zip(meshes, values, strict=True)
    ]

    for index in range(1, len(meshes) - 1):
        coarse_difference = values[index] - values[index - 1]
        fine_difference = values[index + 1] - values[index]
        if coarse_difference == 0.0 or fine_difference == 0.0:
            continue
        ratio = abs(coarse_difference / fine_difference)
        if ratio <= 1.0:
            continue
        order = log2(ratio)
        denominator = 2.0**order - 1.0
        error = abs(coarse_difference) / denominator
        extrapolated = values[index] + coarse_difference / denominator
        estimates[index] = MidpointEstimate(
            mesh_size=meshes[index],
            value=values[index],
            observed_order=order,
            estimated_error=error,
            richardson_value=extrapolated,
        )

    return tuple(estimates)
