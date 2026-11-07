"""Majorana Zero Mode Spectroscopy, Real-Space Local Density of States, and arXiv:2209.12930 Reproduction.

Implements the MZM detection and spatial spectroscopy engine for arXiv:2209.12930:
- Two-vortex background construction and mid-gap MZM state isolation (eps_MZM << Delta_bulk).
- Real-space MZM wavefunction density profile rho_MZM(r) and localization length xi_M.
- Quantitative reproduction of Figure 3 (3-tier STS local spectra: Core, Ring, Bulk).
- Quantitative reproduction of Figure 4 (Spatial radial decay of zero-bias conductance).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kitaev_mf.services.phase_3.bogoliubov import (
    extract_mode_wavefunctions,
)
from kitaev_mf.services.phase_3.gauge_theory import (
    DiagonalizedMatterSector,
    DirectedBond,
    GaugeConfiguration,
    assemble_matter_skew_matrix,
    build_honeycomb_torus_cluster,
    build_uniform_gauge_configuration,
    diagonalize_matter_hamiltonian,
)
from kitaev_mf.services.phase_3.spectral_engine import (
    ContinuousSpectrum,
    SpectralPole,
    compute_local_spectral_poles,
    evaluate_spectral_function,
)
from kitaev_mf.services.phase_3.vortex_engine import (
    HoneycombGeometry,
    create_two_vortex_gauge_configuration,
)


@dataclass(frozen=True)
class TwoVortexMZMSystem:
    """An isolated two-vortex configuration binding a pair of in-gap Majorana Zero Modes."""

    L: int
    p1: tuple[int, int]  # Coordinate of first vortex plaquette
    p2: tuple[int, int]  # Coordinate of second vortex plaquette
    vortex_distance: float  # Shortest Euclidean distance between p1 and p2
    geom: HoneycombGeometry
    bonds: list[DirectedBond]
    triplets: list[tuple[int, int, int, int]]
    gauge: GaugeConfiguration
    diag: DiagonalizedMatterSector
    mzm_energies: tuple[float, float]  # (eps_0, eps_1) mid-gap excitation energies
    bulk_gap_estimate: float  # 6 * sqrt(3) * kappa
    in_gap_ratio: float  # eps_1 / eps_2 (ratio of highest MZM to lowest bulk mode)
    mzm_density: np.ndarray  # rho_MZM(j) = |phi_0(j)|^2 + |phi_1(j)|^2, shape (num_sites,)


@dataclass(frozen=True)
class Figure3Spectrogram:
    """Quantitative reproduction of arXiv:2209.12930 Figure 3 (Core vs Ring vs Bulk spectra)."""

    omega: np.ndarray  # Frequency grid
    core_spectrum: ContinuousSpectrum  # Tip on vortex core boundary (dominant zero-bias peak)
    ring_spectrum: ContinuousSpectrum  # Tip on first neighbor ring (suppressed zero-bias peak)
    bulk_spectrum: (
        ContinuousSpectrum  # Tip in bulk far from vortices (zero subgap weight, vison onset)
    )
    core_site: int
    ring_site: int
    bulk_site: int
    core_poles: list[SpectralPole]
    ring_poles: list[SpectralPole]
    bulk_poles: list[SpectralPole]
    core_lowest_energy: float
    bulk_lowest_energy: float


@dataclass(frozen=True)
class Figure4SpatialDecay:
    """Quantitative reproduction of arXiv:2209.12930 Figure 4 (Spatial conductance decay)."""

    radial_distances: np.ndarray  # Binned Euclidean distances r from nearest vortex core
    low_energy_weights: np.ndarray  # Mean low-bias spectral weight I_low(r)
    mzm_densities: np.ndarray  # Mean wavefunction density rho_MZM(r)
    fitted_xi_M: float  # noqa: N815  # Fitted Majorana localization length xi_M in units of lattice constant a
    fit_slope: float  # Linear regression slope d(ln rho) / dr
    fit_r_squared: float  # Coefficient of determination R^2
    omega_cutoff: float  # Energy threshold separating MZM from bulk continuum


def build_two_vortex_mzm_system(
    L: int = 6,
    p1: tuple[int, int] = (1, 1),
    p2: tuple[int, int] = (4, 4),
    K: float = 1.0,
    kappa: float = 0.15,
) -> TwoVortexMZMSystem:
    """Construct and diagonalize a two-vortex cluster binding isolated mid-gap MZMs.

    Parameters:
        L: Linear cluster dimension (default L=6, N=72 sites).
        p1, p2: Plaquette coordinates of the two trapped vortices.
        K: Kitaev exchange coupling (K > 0 for ferromagnetic).
        kappa: Chiral 3-spin coupling breaking time-reversal symmetry.
    """
    num_sites, bonds, triplets = build_honeycomb_torus_cluster(L=L)
    geom = HoneycombGeometry(L=L)
    gauge_0 = build_uniform_gauge_configuration(num_sites=num_sites, bonds=bonds)

    gauge_2v, _ = create_two_vortex_gauge_configuration(base_gauge=gauge_0, p1=p1, p2=p2, geom=geom)

    A_2v = assemble_matter_skew_matrix(gauge_2v, K=K, kappa=kappa, triplets=triplets)
    diag_2v = diagonalize_matter_hamiltonian(A_2v, gauge=gauge_2v, boundary="obc", L=L)

    phi = extract_mode_wavefunctions(diag_2v.Q)
    # Total MZM probability density from the two mid-gap modes
    rho_mzm = np.abs(phi[0, :]) ** 2 + np.abs(phi[1, :]) ** 2

    eps_0 = float(diag_2v.energies[0])
    eps_1 = float(diag_2v.energies[1])
    eps_2 = float(diag_2v.energies[2])
    bulk_gap = float(6.0 * np.sqrt(3.0) * kappa)
    in_gap_ratio = float(eps_1 / eps_2) if eps_2 > 1e-12 else 1.0
    vortex_dist = geom.torus_distance(p1, p2)

    return TwoVortexMZMSystem(
        L=L,
        p1=p1,
        p2=p2,
        vortex_distance=vortex_dist,
        geom=geom,
        bonds=bonds,
        triplets=triplets,
        gauge=gauge_2v,
        diag=diag_2v,
        mzm_energies=(eps_0, eps_1),
        bulk_gap_estimate=bulk_gap,
        in_gap_ratio=in_gap_ratio,
        mzm_density=rho_mzm,
    )


def reproduce_figure_3_spectrogram(
    system: TwoVortexMZMSystem,
    omega_grid: np.ndarray | None = None,
    sigma: float = 0.06,
    bond_type: int = 3,
    K: float = 1.0,
    kappa: float = 0.15,
) -> Figure3Spectrogram:
    """Reproduce published Figure 3 of arXiv:2209.12930.

    Extracts local dynamical spin correlation spectra S_{jj}^{gamma gamma}(omega) for:
    1. Core site (on the boundary of vortex plaquette p1).
    2. Ring site (on the first ring around p1, not touching p1).
    3. Bulk site (far from both p1 and p2).
    """
    geom = system.geom
    p1 = system.p1
    p2 = system.p2

    # 1. Identify Core site: first site bounding p1
    p1_bonds = geom.plaquette_bonds(p1[0], p1[1])
    core_site = p1_bonds[0][0]

    # 2. Identify Ring site: on adjacent plaquette (p1[0], p1[1] + 1) not touching p1
    p_adj = (p1[0], (p1[1] + 1) % geom.L)
    p1_sites = {s for b in p1_bonds for s in b}
    p_adj_bonds = geom.plaquette_bonds(p_adj[0], p_adj[1])
    ring_candidates = [b[0] for b in p_adj_bonds if b[0] not in p1_sites]
    ring_site = ring_candidates[0] if ring_candidates else p_adj_bonds[0][0]

    # 3. Identify Bulk site: plaquette maximizing distance from both p1 and p2
    best_dist = -1.0
    bulk_site = 0
    for m2 in range(geom.L):
        for m1 in range(geom.L):
            d = min(geom.torus_distance((m1, m2), p1), geom.torus_distance((m1, m2), p2))
            if d > best_dist:
                best_dist = d
                bulk_site = geom.plaquette_bonds(m1, m2)[0][0]

    if omega_grid is None:
        omega_grid = np.linspace(0.0, 8.0, 500)

    # Compute discrete poles
    spec_core = compute_local_spectral_poles(
        gauge=system.gauge,
        site_j=core_site,
        bond_type=bond_type,
        geom=geom,
        K=K,
        kappa=kappa,
        boundary="obc",
        bonds=system.bonds,
        triplets=system.triplets,
        diag_initial=system.diag,
    )
    spec_ring = compute_local_spectral_poles(
        gauge=system.gauge,
        site_j=ring_site,
        bond_type=bond_type,
        geom=geom,
        K=K,
        kappa=kappa,
        boundary="obc",
        bonds=system.bonds,
        triplets=system.triplets,
        diag_initial=system.diag,
    )
    spec_bulk = compute_local_spectral_poles(
        gauge=system.gauge,
        site_j=bulk_site,
        bond_type=bond_type,
        geom=geom,
        K=K,
        kappa=kappa,
        boundary="obc",
        bonds=system.bonds,
        triplets=system.triplets,
        diag_initial=system.diag,
    )

    # Convolve with Gaussian broadening
    cont_core = evaluate_spectral_function(
        spec_core, omega_grid=omega_grid, broadening_type="gaussian", broadening_param=sigma
    )
    cont_ring = evaluate_spectral_function(
        spec_ring, omega_grid=omega_grid, broadening_type="gaussian", broadening_param=sigma
    )
    cont_bulk = evaluate_spectral_function(
        spec_bulk, omega_grid=omega_grid, broadening_type="gaussian", broadening_param=sigma
    )

    return Figure3Spectrogram(
        omega=omega_grid,
        core_spectrum=cont_core,
        ring_spectrum=cont_ring,
        bulk_spectrum=cont_bulk,
        core_site=core_site,
        ring_site=ring_site,
        bulk_site=bulk_site,
        core_poles=spec_core.poles,
        ring_poles=spec_ring.poles,
        bulk_poles=spec_bulk.poles,
        core_lowest_energy=spec_core.poles[0].energy,
        bulk_lowest_energy=spec_bulk.poles[0].energy,
    )


def reproduce_figure_4_spatial_decay(
    system: TwoVortexMZMSystem,
    omega_cutoff: float = 0.35,
    max_fit_distance: float = 3.5,
    bond_type: int = 3,
    K: float = 1.0,
    kappa: float = 0.15,
) -> Figure4SpatialDecay:
    """Reproduce published Figure 4 of arXiv:2209.12930.

    Evaluates the spatial decay of the zero-bias / low-energy conductance weight
    I_low(r) = sum_{lambda: omega_lambda <= omega_cutoff} w_lambda
    as a function of distance r from the vortex core, fitting the Majorana localization length xi_M.

    Note on radial binning:
        Distance r is measured from unit-cell center (m1, m2) to the nearest vortex plaquette center,
        providing a natural discrete coarse-graining across hexagonal rings on the lattice.
    """
    geom = system.geom
    p1 = system.p1
    p2 = system.p2
    num_sites = geom.num_sites

    site_distances: list[float] = []
    low_weights: list[float] = []
    densities: list[float] = []

    for s in range(num_sites):
        m1, m2, sub = geom.site_coordinates(s)
        d = min(geom.torus_distance((m1, m2), p1), geom.torus_distance((m1, m2), p2))
        site_distances.append(d)
        densities.append(float(system.mzm_density[s]))

        # Compute low-energy spectral weight
        spec = compute_local_spectral_poles(
            gauge=system.gauge,
            site_j=s,
            bond_type=bond_type,
            geom=geom,
            K=K,
            kappa=kappa,
            boundary="obc",
            bonds=system.bonds,
            triplets=system.triplets,
            diag_initial=system.diag,
        )
        w_low = float(sum(p.weight for p in spec.poles if p.energy <= omega_cutoff))
        low_weights.append(w_low)

    # Bin by rounded distance
    unique_dists = sorted({round(d, 3) for d in site_distances})
    binned_r: list[float] = []
    binned_w: list[float] = []
    binned_rho: list[float] = []

    for d in unique_dists:
        idxs = [i for i, dist in enumerate(site_distances) if round(dist, 3) == d]
        binned_r.append(d)
        binned_w.append(float(np.mean([low_weights[i] for i in idxs])))
        binned_rho.append(float(np.mean([densities[i] for i in idxs])))

    r_arr = np.array(binned_r)
    rho_arr = np.array(binned_rho)

    # Fit exponential decay ln(rho) = slope * r + intercept for r <= max_fit_distance
    fit_mask = r_arr <= max_fit_distance
    if np.sum(fit_mask) >= 2:
        r_fit = r_arr[fit_mask]
        log_rho = np.log(rho_arr[fit_mask])
        slope, intercept = np.polyfit(r_fit, log_rho, 1)

        # Coefficient of determination R^2
        y_pred = slope * r_fit + intercept
        ss_res = float(np.sum((log_rho - y_pred) ** 2))
        ss_tot = float(np.sum((log_rho - np.mean(log_rho)) ** 2))
        r_squared = float(1.0 - ss_res / ss_tot) if ss_tot > 1e-12 else 1.0

        xi_M = float(-1.0 / slope) if slope < 0 else float("inf")
    else:
        slope = 0.0
        r_squared = 0.0
        xi_M = float("inf")

    return Figure4SpatialDecay(
        radial_distances=r_arr,
        low_energy_weights=np.array(binned_w),
        mzm_densities=rho_arr,
        fitted_xi_M=xi_M,
        fit_slope=float(slope),
        fit_r_squared=float(r_squared),
        omega_cutoff=omega_cutoff,
    )
