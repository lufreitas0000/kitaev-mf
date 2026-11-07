"""Cylinder geometry, chiral Majorana edge modes, and finite-size scaling.

Implements Milestone M3.3 features:
- Construction and analysis of honeycomb cylinders (periodic along a1, open along a2).
- Boundary localization and dispersion of chiral Majorana edge modes (nu = 1).
- Local dynamical spin correlation comparison between edge and bulk sites.
- Finite-size scaling of the Majorana Zero Mode hybridization doublet on tori (L in {4, 6, 8}).

Author: Lu Freitas
Date: 2026-11-07
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from kitaev_mf.services.phase_3.bogoliubov import (
    extract_mode_wavefunctions,
)
from kitaev_mf.services.phase_3.gauge_theory import (
    DiagonalizedMatterSector,
    DirectedBond,
    GaugeConfiguration,
    assemble_matter_skew_matrix,
    build_honeycomb_cylinder_cluster,
    build_uniform_gauge_configuration,
    diagonalize_matter_hamiltonian,
)
from kitaev_mf.services.phase_3.mzm_spectroscopy import (
    build_two_vortex_mzm_system,
)
from kitaev_mf.services.phase_3.spectral_engine import (
    ContinuousSpectrum,
    compute_local_spectral_poles,
    evaluate_spectral_function,
)
from kitaev_mf.services.phase_3.vortex_engine import (
    HoneycombGeometry,
)


@dataclass(frozen=True)
class CylinderSystem:
    """A honeycomb cylinder system with periodic boundary along a1 and open boundary along a2."""

    Lx: int
    Ly: int
    num_sites: int
    bonds: list[DirectedBond]
    triplets: list[tuple[int, int, int, int]]
    gauge: GaugeConfiguration
    diag: DiagonalizedMatterSector
    edge_sites_bottom: list[int]  # Sites with m2 = 0
    edge_sites_top: list[int]  # Sites with m2 = Ly - 1
    bulk_sites: list[int]  # Sites located in the middle ribbon of the cylinder
    bulk_gap_estimate: float  # 6 * sqrt(3) * kappa


def build_honeycomb_cylinder_system(
    Lx: int = 6,
    Ly: int = 12,
    K: float = 1.0,
    kappa: float = 0.15,
) -> CylinderSystem:
    """Construct and diagonalize a honeycomb cylinder cluster with chiral 3-spin coupling."""
    num_sites, bonds, triplets = build_honeycomb_cylinder_cluster(Lx=Lx, Ly=Ly)
    gauge_0 = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    A = assemble_matter_skew_matrix(gauge_0, K=K, kappa=kappa, triplets=triplets)
    diag = diagonalize_matter_hamiltonian(A, gauge=gauge_0, boundary="obc", L=Lx)

    def site_idx(m1: int, m2: int, sub: int) -> int:
        return 2 * ((m1 % Lx) + Lx * m2) + sub

    edge_bottom = [site_idx(m1, 0, sub) for m1 in range(Lx) for sub in (0, 1)]
    edge_top = [site_idx(m1, Ly - 1, sub) for m1 in range(Lx) for sub in (0, 1)]

    # Bulk ribbon in middle third
    m2_min = Ly // 3
    m2_max = 2 * Ly // 3
    bulk = [
        site_idx(m1, m2, sub) for m2 in range(m2_min, m2_max) for m1 in range(Lx) for sub in (0, 1)
    ]

    bulk_gap = float(6.0 * np.sqrt(3.0) * kappa)

    return CylinderSystem(
        Lx=Lx,
        Ly=Ly,
        num_sites=num_sites,
        bonds=bonds,
        triplets=triplets,
        gauge=gauge_0,
        diag=diag,
        edge_sites_bottom=edge_bottom,
        edge_sites_top=edge_top,
        bulk_sites=bulk,
        bulk_gap_estimate=bulk_gap,
    )


def analyze_cylinder_edge_modes(
    cylinder: CylinderSystem,
    max_modes: int = 4,
) -> dict[str, Any]:
    """Analyze the energy and spatial wavefunctions of in-gap chiral edge modes on the cylinder.

    Returns:
        Dictionary containing:
        - energies: list of lowest positive excitation energies
        - in_gap_count: number of states below bulk gap estimate
        - layer_densities: array of shape (max_modes, Ly) showing P_n(m2)
        - boundary_fractions: array of shape (max_modes,) showing fraction in outer 2 layers
    """
    diag = cylinder.diag
    phi = extract_mode_wavefunctions(diag.Q)
    Lx = cylinder.Lx
    Ly = cylinder.Ly

    def site_idx(m1: int, m2: int, sub: int) -> int:
        return 2 * ((m1 % Lx) + Lx * m2) + sub

    energies = [float(e) for e in diag.energies[:max_modes]]
    in_gap_count = int(np.sum(diag.energies < cylinder.bulk_gap_estimate))

    layer_densities = np.zeros((max_modes, Ly), dtype=float)
    boundary_fractions = np.zeros(max_modes, dtype=float)

    for n in range(max_modes):
        prob = np.abs(phi[n, :]) ** 2
        for m2 in range(Ly):
            layer_sum = 0.0
            for m1 in range(Lx):
                for sub in (0, 1):
                    layer_sum += prob[site_idx(m1, m2, sub)]
            layer_densities[n, m2] = layer_sum

        # Fraction in outer 2 layers at both boundaries: m2 in {0, 1, Ly-2, Ly-1}
        outer_sum = (
            layer_densities[n, 0]
            + layer_densities[n, 1]
            + layer_densities[n, Ly - 2]
            + layer_densities[n, Ly - 1]
        )
        total_norm = float(np.sum(prob))
        boundary_fractions[n] = float(outer_sum / total_norm) if total_norm > 1e-12 else 0.0

    return {
        "energies": energies,
        "in_gap_count": in_gap_count,
        "layer_densities": layer_densities,
        "boundary_fractions": boundary_fractions,
        "bulk_gap_estimate": cylinder.bulk_gap_estimate,
    }


def compute_cylinder_edge_vs_bulk_spectra(
    cylinder: CylinderSystem,
    omega_grid: np.ndarray | None = None,
    sigma: float = 0.06,
    K: float = 1.0,
    kappa: float = 0.15,
) -> tuple[ContinuousSpectrum, ContinuousSpectrum]:
    """Compute local dynamical spin correlation spectra comparing an edge site vs a bulk site.

    Parameters:
        cylinder: The initialized CylinderSystem.
        omega_grid: Frequency grid points.
        sigma: Gaussian broadening standard deviation.

    Returns:
        (edge_spectrum, bulk_spectrum) as ContinuousSpectrum objects.
    """
    if omega_grid is None:
        omega_grid = np.linspace(0.0, 4.0, 400)

    # Edge site: bottom boundary (m1=0, m2=0, A)
    edge_site = cylinder.edge_sites_bottom[0]
    # Bulk site: middle ribbon (m1=0, m2=Ly//2, A)
    bulk_site = cylinder.bulk_sites[0]

    # HoneycombGeometry dummy for bond lookups
    geom = HoneycombGeometry(L=cylinder.Lx)

    spec_edge = compute_local_spectral_poles(
        gauge=cylinder.gauge,
        site_j=edge_site,
        bond_type=3,  # z-bond
        geom=geom,
        K=K,
        kappa=kappa,
        boundary="obc",
        bonds=cylinder.bonds,
        triplets=cylinder.triplets,
        diag_initial=cylinder.diag,
    )

    spec_bulk = compute_local_spectral_poles(
        gauge=cylinder.gauge,
        site_j=bulk_site,
        bond_type=3,  # z-bond
        geom=geom,
        K=K,
        kappa=kappa,
        boundary="obc",
        bonds=cylinder.bonds,
        triplets=cylinder.triplets,
        diag_initial=cylinder.diag,
    )

    cont_edge = evaluate_spectral_function(
        spec_edge, omega_grid=omega_grid, broadening_type="gaussian", broadening_param=sigma
    )
    cont_bulk = evaluate_spectral_function(
        spec_bulk, omega_grid=omega_grid, broadening_type="gaussian", broadening_param=sigma
    )

    return cont_edge, cont_bulk


def sweep_finite_size_mzm_scaling(
    L_values: list[int] | None = None,
    K: float = 1.0,
    kappa: float = 0.15,
) -> list[dict[str, Any]]:
    """Perform finite-size scaling sweep across torus clusters L in {4, 6, 8}.

    Evaluates the exponential suppression of MZM doublet splitting delta_eps = |eps_1 - eps_0|
    with increasing vortex separation R_12.

    Returns:
        List of result dictionaries for each L.
    """
    if L_values is None:
        L_values = [4, 6, 8]

    records: list[dict[str, Any]] = []

    for L in L_values:
        p1 = (1, 1)
        p2 = (L // 2 + 1, L // 2 + 1)
        sys = build_two_vortex_mzm_system(L=L, p1=p1, p2=p2, K=K, kappa=kappa)
        eps_0, eps_1 = sys.mzm_energies
        delta_eps = float(abs(eps_1 - eps_0))

        records.append(
            {
                "L": L,
                "num_sites": 2 * L * L,
                "vortex_distance": float(sys.vortex_distance),
                "eps_0": float(eps_0),
                "eps_1": float(eps_1),
                "delta_eps": delta_eps,
                "in_gap_ratio": float(sys.in_gap_ratio),
                "bulk_gap_estimate": float(sys.bulk_gap_estimate),
            }
        )

    return records
