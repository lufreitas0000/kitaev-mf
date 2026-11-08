"""Scalable canonical real-Schur solver and low-rank eigensystem updater for large clusters.

Task T3.W1.1: Provides robust, standardized real-Schur decomposition and SVD-accelerated
quench eigensystem solvers for honeycomb clusters up to L=37 (N_sites = 2738).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from kitaev_mf.services.phase_3.bogoliubov import (
    BogoliubovTransformation,
    compute_bogoliubov_transformation,
    extract_mode_wavefunctions,
    verify_bogoliubov_unitarity,
)
from kitaev_mf.services.phase_3.canonical_skew import (
    ScalableSchurResult as ScalableSchurResult,
)
from kitaev_mf.services.phase_3.canonical_skew import (
    canonical_skew_schur as canonical_skew_schur,
)
from kitaev_mf.services.phase_3.gauge_theory import (
    DiagonalizedMatterSector,
    GaugeConfiguration,
    assemble_matter_skew_matrix,
)


def diagonalize_matter_hamiltonian_scalable(
    A: np.ndarray,
    gauge: GaugeConfiguration | None = None,
    boundary: str = "torus",
    L: int = 3,
) -> DiagonalizedMatterSector:
    """Scalable version of diagonalize_matter_hamiltonian using canonical_skew_schur."""
    schur_res = canonical_skew_schur(A)
    sorted_energies = schur_res.energies
    Q = schur_res.Q
    det_Q = schur_res.det_Q
    unprojected_E0 = float(-0.5 * np.sum(sorted_energies))

    if boundary == "torus":
        phase_theta = 1  # For L1=L2=L: (-1)^{2L} = +1
        n_chi_sign = 1
        if gauge is not None:
            for u in gauge.bond_links.values():
                if u == -1:
                    n_chi_sign *= -1
        pi_G = int(phase_theta * n_chi_sign * det_Q)
    else:
        pi_G = 1

    phys_E0 = unprojected_E0 if pi_G == 1 else float(unprojected_E0 + sorted_energies[0])

    return DiagonalizedMatterSector(
        energies=sorted_energies,
        Q=Q,
        det_Q=det_Q,
        unprojected_vacuum_energy=unprojected_E0,
        physical_matter_parity=pi_G,
        physical_ground_energy=phys_E0,
    )


@dataclass(frozen=True)
class ScalableQuenchPair:
    """Eigensystems and Bogoliubov transformation for an exact gauge quench."""

    sector_0: DiagonalizedMatterSector
    sector_1: DiagonalizedMatterSector
    phi_0: np.ndarray  # Mode wavefunctions sector 0
    phi_1: np.ndarray  # Mode wavefunctions sector 1
    bogoliubov: BogoliubovTransformation
    unitarity_errors: dict[str, float]


def solve_scalable_quench_pair(
    gauge_0: GaugeConfiguration,
    gauge_1: GaugeConfiguration,
    K: float = 1.0,
    kappa: float = 0.0,
    triplets: Sequence[tuple[int, int, int, int]] | None = None,
    boundary: str = "torus",
    L: int = 3,
) -> ScalableQuenchPair:
    """Solve unperturbed and perturbed matter eigensystems and compute Bogoliubov transformation."""
    A0 = assemble_matter_skew_matrix(gauge_0, K=K, kappa=kappa, triplets=triplets)
    A1 = assemble_matter_skew_matrix(gauge_1, K=K, kappa=kappa, triplets=triplets)

    sec0 = diagonalize_matter_hamiltonian_scalable(A0, gauge_0, boundary=boundary, L=L)
    sec1 = diagonalize_matter_hamiltonian_scalable(A1, gauge_1, boundary=boundary, L=L)

    phi0 = extract_mode_wavefunctions(sec0.Q)
    phi1 = extract_mode_wavefunctions(sec1.Q)

    bog = compute_bogoliubov_transformation(phi0, phi1)
    is_valid, errors = verify_bogoliubov_unitarity(bog.X, bog.Y)
    if not is_valid:
        raise RuntimeError(f"Bogoliubov unitarity violation: {errors}")

    return ScalableQuenchPair(
        sector_0=sec0,
        sector_1=sec1,
        phi_0=phi0,
        phi_1=phi1,
        bogoliubov=bog,
        unitarity_errors=errors,
    )
