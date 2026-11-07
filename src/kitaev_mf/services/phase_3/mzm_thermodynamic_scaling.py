"""Thermodynamic Limit Scaling of MZM Hybridization Suppression and Zero-Bias Sharpness.

Services for quantifying:
1. MZM hybridization splitting suppression: delta_epsilon_0(R) ~ exp(-R / xi_M).
2. Zero-bias MZM peak convergence: omega_MZM(L) -> 0 as L -> infty.
3. Zero-bias spectral sharpness saturation: I_subgap / I_total >= 0.15 across all L.
4. Bulk-site subgap vanishing & vison gap protection: I_subgap_bulk < 10^-5 and omega_bulk >= Delta E_2v.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import scipy.stats

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
from kitaev_mf.services.phase_3.gauge_theory import (
    build_honeycomb_cylinder_cluster as build_honeycomb_cylinder_cluster,
)
from kitaev_mf.services.phase_3.gauge_theory import flip_bond_gauge_link as flip_bond_gauge_link
from kitaev_mf.services.phase_3.spectral_engine import DiscreteSpectralData as DiscreteSpectralData
from kitaev_mf.services.phase_3.spectral_engine import SpectralPole as SpectralPole
from kitaev_mf.services.phase_3.spectral_engine import (
    compute_local_spectral_poles,
)
from kitaev_mf.services.phase_3.vortex_engine import (
    HoneycombGeometry,
    create_two_vortex_gauge_configuration,
)


@dataclasses.dataclass(frozen=True)
class MZMHybridizationDatum:
    """Datapoint for MZM hybridization splitting at separation R."""

    L: int
    dy: int
    R_units_a: float
    delta_eps_0: float


@dataclasses.dataclass(frozen=True)
class MZMExponentialFitResult:
    """Results of log-linear regression for MZM hybridization splitting."""

    xi_M: float  # noqa: N815 - established mathematical field name
    xi_M_stderr: float  # noqa: N815 - established mathematical field name
    intercept_alpha: float
    r_squared: float
    p_value: float
    n_points: int


@dataclasses.dataclass(frozen=True)
class ThermodynamicSharpnessDatum:
    """Datapoint for dynamical quench spectral sharpness at linear dimension L."""

    L: int
    boundary: str
    omega_MZM: float  # noqa: N815 - established mathematical field name
    total_spectral_weight: float
    subgap_weight: float
    sharpness_ratio: float
    bulk_subgap_weight: float
    bulk_min_energy: float


def measure_mzm_hybridization_scaling(
    L_values: list[int] | None = None,
    K: float = 1.0,
    kappa: float = 0.15,
) -> tuple[list[MZMHybridizationDatum], MZMExponentialFitResult]:
    """Measure MZM hybridization splitting delta_eps_0(R_max) across commensurate sizes 3 | L.

    Parameters
    ----------
    L_values : Optional[List[int]]
        Commensurate linear sizes (default [6, 12, 18, 24]).
    K : float
        Kitaev coupling strength.
    kappa : float
        Time-reversal breaking Chern mass coupling.

    Returns
    -------
    data : List[MZMHybridizationDatum]
        Measured hybridization energies at antipodal separation R_max.
    fit : MZMExponentialFitResult
        Fitted exponential decay length xi_M and regression statistics.
    """
    if L_values is None:
        L_values = [6, 12, 18, 24]

    data: list[MZMHybridizationDatum] = []

    for L in L_values:
        geom = HoneycombGeometry(L=L)
        p1 = (1, 1)
        p2 = (L // 2 + 1, L // 2 + 1)
        nodes, bonds, triplets = build_honeycomb_torus_cluster(L=L)
        g0 = build_uniform_gauge_configuration(nodes, bonds)
        g2v, _ = create_two_vortex_gauge_configuration(g0, p1, p2, geom)

        A2v = assemble_matter_skew_matrix(g2v, K=K, kappa=kappa, triplets=triplets)
        d2v = diagonalize_matter_hamiltonian(A2v, gauge=g2v, boundary="torus", L=L)

        delta_eps0 = float(d2v.energies[0])
        R_max = float((L / 2) * np.sqrt(3.0))
        data.append(MZMHybridizationDatum(L=L, dy=L // 2, R_units_a=R_max, delta_eps_0=delta_eps0))

    # Perform log-linear regression: ln(delta_eps_0) = alpha - R / xi_M
    Rs = np.array([d.R_units_a for d in data])
    y_vals = np.log(np.array([d.delta_eps_0 for d in data]))

    slope, intercept, r_val, p_val, stderr = scipy.stats.linregress(Rs, y_vals)
    xi_M = float(-1.0 / slope) if slope < 0 else float("inf")
    xi_M_err = float(stderr / (slope**2)) if slope < 0 else float("inf")

    fit = MZMExponentialFitResult(
        xi_M=xi_M,
        xi_M_stderr=xi_M_err,
        intercept_alpha=float(intercept),
        r_squared=float(r_val**2),
        p_value=float(p_val),
        n_points=len(data),
    )
    return data, fit


def compute_thermodynamic_zero_bias_sharpness(
    L_values: list[int],
    K: float = 1.0,
    kappa: float = 0.15,
    boundary: str = "obc",
    delta_subgap_cutoff: float = 0.18,
) -> list[ThermodynamicSharpnessDatum]:
    """Compute local dynamical quench response and spectral sharpness across L values.

    Parameters
    ----------
    L_values : List[int]
        Linear sizes to benchmark.
    K : float
        Kitaev coupling strength.
    kappa : float
        Chern mass coupling.
    boundary : str
        Parity boundary policy ('obc' or 'torus'). Default 'obc' (thermodynamic bulk parity).
    delta_subgap_cutoff : float
        Subgap threshold. Defaults to 0.18 |K|, strictly below the two-vison creation gap
        Delta E_2v >= 0.209 |K| across all sizes, while including the hybridized MZM peak
        for small clusters (omega_MZM = 0.1564 at L=6) and converging to zero in the
        thermodynamic limit (omega_MZM = 0.00317 at L=20).

    Returns
    -------
    results : List[ThermodynamicSharpnessDatum]
        Sharpness and peak energy metrics for each size L.
    """
    subgap_window = delta_subgap_cutoff
    results: list[ThermodynamicSharpnessDatum] = []

    for L in L_values:
        geom = HoneycombGeometry(L=L)
        p1 = (1, 1)
        p2 = (L // 2 + 1, L // 2 + 1)

        nodes, bonds, triplets = build_honeycomb_torus_cluster(L=L)
        g0 = build_uniform_gauge_configuration(nodes, bonds)
        g2v, _ = create_two_vortex_gauge_configuration(g0, p1, p2, geom)

        # 1. Core site quench: bond connected to vortex core p1
        p1_bonds = geom.plaquette_bonds(*p1)
        site_core = p1_bonds[0][0]

        spec_core = compute_local_spectral_poles(
            gauge=g2v,
            site_j=site_core,
            bond_type=3,  # z-bond
            geom=geom,
            K=K,
            kappa=kappa,
            boundary=boundary,
            bonds=bonds,
            triplets=triplets,
        )

        tot_weight_core = spec_core.total_spectral_weight
        subgap_poles_core = [p for p in spec_core.poles if p.energy <= subgap_window]
        subgap_weight_core = float(sum(p.weight for p in subgap_poles_core))
        sharpness = subgap_weight_core / tot_weight_core if tot_weight_core > 0 else 0.0
        omega_mzm = (
            float(subgap_poles_core[0].energy)
            if len(subgap_poles_core) > 0
            else float(spec_core.poles[0].energy)
        )

        # 2. Bulk site quench: find site with maximal distance from both vortex cores on torus
        best_dist = -1.0
        bulk_site = 0
        for m2 in range(geom.L):
            for m1 in range(geom.L):
                d = min(geom.torus_distance((m1, m2), p1), geom.torus_distance((m1, m2), p2))
                if d > best_dist:
                    best_dist = d
                    bulk_site = geom.plaquette_bonds(m1, m2)[0][0]

        spec_bulk = compute_local_spectral_poles(
            gauge=g2v,
            site_j=bulk_site,
            bond_type=3,
            geom=geom,
            K=K,
            kappa=kappa,
            boundary=boundary,
            bonds=bonds,
            triplets=triplets,
        )

        subgap_poles_bulk = [p for p in spec_bulk.poles if p.energy <= subgap_window]
        bulk_subgap_weight = float(sum(p.weight for p in subgap_poles_bulk))
        bulk_min_energy = float(spec_bulk.poles[0].energy) if len(spec_bulk.poles) > 0 else 0.0

        results.append(
            ThermodynamicSharpnessDatum(
                L=L,
                boundary=boundary,
                omega_MZM=omega_mzm,
                total_spectral_weight=tot_weight_core,
                subgap_weight=subgap_weight_core,
                sharpness_ratio=sharpness,
                bulk_subgap_weight=bulk_subgap_weight,
                bulk_min_energy=bulk_min_energy,
            )
        )

    return results
