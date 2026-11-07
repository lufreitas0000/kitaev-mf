"""Thermodynamic Sum Rules and Spectral Functions Service up to L=37.

Provides comprehensive Lehmann representation analysis across system sizes up to L=37:
1. Core vs. bulk Lehmann poles and sum-rule saturation.
2. Continuous spectral function evaluation and quadrature normalization.
3. Universal energy positivity and MZM vs. bulk demarcation gap.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from kitaev_mf.services.phase_3.gauge_theory import (
    DiagonalizedMatterSector as DiagonalizedMatterSector,
)

# Preserve historical module-level project re-exports.
from kitaev_mf.services.phase_3.gauge_theory import DirectedBond as DirectedBond
from kitaev_mf.services.phase_3.gauge_theory import GaugeConfiguration as GaugeConfiguration
from kitaev_mf.services.phase_3.gauge_theory import (
    assemble_matter_skew_matrix,
    build_honeycomb_torus_cluster,
    build_uniform_gauge_configuration,
    diagonalize_matter_hamiltonian,
)
from kitaev_mf.services.phase_3.spectral_engine import ContinuousSpectrum as ContinuousSpectrum
from kitaev_mf.services.phase_3.spectral_engine import DiscreteSpectralData as DiscreteSpectralData
from kitaev_mf.services.phase_3.spectral_engine import SpectralPole as SpectralPole
from kitaev_mf.services.phase_3.spectral_engine import (
    compute_local_spectral_poles,
    evaluate_spectral_function,
)
from kitaev_mf.services.phase_3.vortex_engine import (
    HoneycombGeometry,
    create_two_vortex_gauge_configuration,
)
from kitaev_mf.services.phase_3.vortex_engine import (
    apply_local_spin_quench as apply_local_spin_quench,
)


@dataclasses.dataclass(frozen=True)
class ThermodynamicSpectralBenchmarkDatum:
    """Datapoint for core vs. bulk spectral response at linear size L."""

    L: int
    n_sites: int
    is_singular_core: bool
    is_singular_bulk: bool
    core_sum_rule: float
    bulk_sum_rule: float
    omega_MZM: float  # noqa: N815 - established mathematical field name
    bulk_min_energy: float
    demarcation_gap: float
    quadrature_error: float


def benchmark_thermodynamic_spectral_response(
    L_values: list[int] | None = None,
    K: float = 1.0,
    kappa: float = 0.15,
    boundary: str = "obc",
) -> list[ThermodynamicSpectralBenchmarkDatum]:
    """Benchmark local dynamical spin correlations on clusters up to L=37.

    Parameters
    ----------
    L_values : Optional[List[int]]
        List of linear dimensions to evaluate. Defaults to [6, 12, 18, 24, 30, 36, 37].
    K : float
        Kitaev coupling strength.
    kappa : float
        Chern mass coupling.
    boundary : str
        Parity boundary policy (default 'obc' for bulk thermodynamic vacuum).

    Returns
    -------
    results : List[ThermodynamicSpectralBenchmarkDatum]
        Measured spectral metrics for each linear size L.
    """
    if L_values is None:
        L_values = [6, 12, 18, 24, 30, 36, 37]

    results: list[ThermodynamicSpectralBenchmarkDatum] = []

    for L in L_values:
        geom = HoneycombGeometry(L=L)
        p1 = (1, 1)
        p2 = (L // 2 + 1, L // 2 + 1)
        nodes, bonds, triplets = build_honeycomb_torus_cluster(L=L)
        g0 = build_uniform_gauge_configuration(nodes, bonds)
        g2v, _ = create_two_vortex_gauge_configuration(g0, p1, p2, geom)

        # Precompute initial unquenched diagonalization to optimize runtime
        A_init = assemble_matter_skew_matrix(g2v, K=K, kappa=kappa, triplets=triplets)
        diag_init = diagonalize_matter_hamiltonian(A_init, gauge=g2v, boundary=boundary, L=L)

        # 1. Core site quench: bond connected to vortex p1
        p1_bonds = geom.plaquette_bonds(*p1)
        site_core = p1_bonds[0][0]

        spec_core = compute_local_spectral_poles(
            gauge=g2v,
            site_j=site_core,
            bond_type=3,
            geom=geom,
            K=K,
            kappa=kappa,
            boundary=boundary,
            bonds=bonds,
            triplets=triplets,
            diag_initial=diag_init,
        )

        # 2. Bulk site quench: find site with maximal distance from both vortex cores on torus
        best_dist = -1.0
        site_bulk = 0
        for m2 in range(geom.L):
            for m1 in range(geom.L):
                d = min(geom.torus_distance((m1, m2), p1), geom.torus_distance((m1, m2), p2))
                if d > best_dist:
                    best_dist = d
                    site_bulk = geom.plaquette_bonds(m1, m2)[0][0]

        spec_bulk = compute_local_spectral_poles(
            gauge=g2v,
            site_j=site_bulk,
            bond_type=3,
            geom=geom,
            K=K,
            kappa=kappa,
            boundary=boundary,
            bonds=bonds,
            triplets=triplets,
            diag_initial=diag_init,
        )

        # 3. Continuous quadrature verification on core spectrum
        # Extend grid sufficiently below and above all poles to capture full Gaussian tails (>= 5 sigma)
        sigma = 0.05
        w_min = float(min(spec_core.poles[0].energy - 6.0 * sigma, -0.5))
        w_max = float(max(spec_core.poles[-1].energy + 6.0 * sigma, 7.0))
        n_grid = int(max(2000, int((w_max - w_min) / 0.005)))
        omega_grid = np.linspace(w_min, w_max, n_grid)
        cont_spec = evaluate_spectral_function(
            spec_core,
            omega_grid=omega_grid,
            broadening_type="gaussian",
            broadening_param=sigma,
        )
        quad_err = float(np.abs(cont_spec.integrated_weight - spec_core.total_spectral_weight))

        # Core MZM energy: lowest physical pole with significant weight (>= 0.001)
        significant_core_poles = [p for p in spec_core.poles if p.weight >= 0.001]
        omega_mzm = (
            float(significant_core_poles[0].energy)
            if len(significant_core_poles) > 0
            else float(spec_core.poles[0].energy)
        )

        # Bulk minimum energy
        significant_bulk_poles = [p for p in spec_bulk.poles if p.weight >= 0.001]
        bulk_min = (
            float(significant_bulk_poles[0].energy)
            if len(significant_bulk_poles) > 0
            else float(spec_bulk.poles[0].energy)
        )

        demarcation = float(bulk_min - omega_mzm)

        results.append(
            ThermodynamicSpectralBenchmarkDatum(
                L=L,
                n_sites=2 * L * L,
                is_singular_core=spec_core.is_singular,
                is_singular_bulk=spec_bulk.is_singular,
                core_sum_rule=float(spec_core.total_spectral_weight),
                bulk_sum_rule=float(spec_bulk.total_spectral_weight),
                omega_MZM=omega_mzm,
                bulk_min_energy=bulk_min,
                demarcation_gap=demarcation,
                quadrature_error=quad_err,
            )
        )

    return results
