"""Unit tests for Wave 1 HC-03 inhomogeneous energy, covariant NNN chiral bridge, and two-vison scaling."""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import InhomogeneousMeanFieldState
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field
from kitaev_mf.services.vortex_pair_scaling import (
    bond_signs_from_state_gauge,
    build_covariant_nnn_matrix,
    compare_exp_vs_powerlaw_scaling,
    compute_inhomogeneous_gaussian_energy,
    measure_two_vison_point,
)


def test_inhomogeneous_energy_matches_homogeneous_to_machine_precision() -> None:
    cell = make_honeycomb_2site()
    cfg = SolverConfig(k_grid_shape=(18, 18), tol_correlator=1e-11, tol_gauge=1e-11)
    params = ModelParameters(K=-1.0, h=(0.12, 0.12, 0.12))
    st_h, met = solve_self_consistent_mean_field(cell, params, make_qsl_seed(cell), cfg)
    assert met.converged

    L = 6
    e_inh = compute_inhomogeneous_gaussian_energy(
        InhomogeneousMeanFieldState.from_homogeneous(st_h, L * L), L, params
    )
    e_hom = compute_gaussian_energy(st_h, cell, params)
    assert abs(e_inh.E_total - (L * L) * e_hom.E_total) < 1e-12


def test_covariant_nnn_gauge_invariance_and_kwg_bridge_match() -> None:
    L = 12
    cluster = RealSpaceTorusCluster(L)
    cell = make_honeycomb_2site()
    cfg = SolverConfig(k_grid_shape=(24, 24), tol_correlator=1e-10, tol_gauge=1e-10)
    st0, _ = solve_self_consistent_mean_field(
        cell, ModelParameters(K=-1.0), make_qsl_seed(cell), cfg
    )

    # Adjacent vison pair across bond (0, 0, 0) corresponds to kappa_static = kappa_KWG / 2 = 0.04
    kappa_static = 0.04
    cut1 = [(0, 0, 0)]
    cut2 = [(0, 0, 1), (0, 0, 2)]  # local Z2 gauge transform at A(0, 0) flips all 3 bonds at site

    specs = []
    for cut in (cut1, cut2):
        signs = np.ones((L, L, 3), dtype=np.float64)
        for m1, m2, b in cut:
            signs[m1, m2, b] = -1.0
        a_tot = cluster.assemble_skew_hamiltonian(
            ModelParameters(K=-1.0), st0, flipped_bonds=cut
        ) + build_covariant_nnn_matrix(cluster, kappa_static, signs)
        w2 = np.linalg.eigvalsh(0.5 * ((-a_tot @ a_tot) + (-a_tot @ a_tot).T))
        specs.append(np.sqrt(np.maximum(w2, 0.0)))

    assert np.max(np.abs(specs[0] - specs[1])) < 1e-12
    # Single-particle E0 / Delta_K = eps_min / (12 * sqrt(3) * kappa_static)
    e0_over_delta = float(specs[0][0]) / (12.0 * np.sqrt(3.0) * kappa_static)
    # Cross-repo exact match against kitaev-walk-geometry chiral_vison_tmatrix_scaling.dat
    assert e0_over_delta == pytest.approx(0.28256246, abs=1e-6)


def test_self_consistent_two_vison_hybridization_and_exp_decay() -> None:
    cache: dict[tuple[int, float, float], tuple[float, float, InhomogeneousMeanFieldState]] = {}
    sample_grid = ((8, 2), (8, 4), (10, 3), (10, 5), (12, 2), (12, 4), (12, 6))
    rows = [
        measure_two_vison_point(
            L,
            r,
            kappa_static=0.06,
            h_field=0.0,
            e0_cache=cache,
            use_eigh=True,
        )
        for L, r in sample_grid
    ]
    assert all(row.converged for row in rows)
    assert all(row.defect_count == 2 for row in rows)
    assert all(row.max_gauge_residual < 1e-11 for row in rows)
    assert all(row.min_abs_pfaffian_wp > 0.99 for row in rows)
    # Strict monotonic decay with separation on L=12
    l12 = [row.omega_over_delta for row in rows if row.L == 12]
    assert l12 == sorted(l12, reverse=True)
    comp = compare_exp_vs_powerlaw_scaling(rows)
    assert comp.rss_exp < comp.rss_pow
    assert comp.delta_bic_exp_over_pow > 0.0
    _ = bond_signs_from_state_gauge
