"""Unit tests for Wave 1 HC-01 cubic Feshbach-Schur decomposition and HC-02 Gutzwiller projection."""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.services.cubic_gap_perturbation import analyze_cubic_gap_at_point
from kitaev_mf.services.gutzwiller_projection import (
    bca_confidence_interval,
    evaluate_exact_fock_equivalence,
    evaluate_exact_gutzwiller,
    evaluate_vmc_gutzwiller,
)
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def test_hc01_bare_zeeman_vanishes_while_kkt_dressed_feshbach_generates_cubic_gap() -> None:
    decomp, _ = analyze_cubic_gap_at_point(0.01, 24)
    assert decomp.converged
    assert decomp.delta_schur_bare_zeeman < 1e-13
    assert decomp.norm_dv_pp < 1e-13
    assert abs(decomp.norm_lambda_over_h - float(np.sqrt(5.0))) < 0.02
    assert decomp.delta_can_direct / (0.01**3) > 9.5
    rel_schur = abs(decomp.delta_schur_full - decomp.delta_can_direct) / decomp.delta_can_direct
    assert rel_schur < 5e-3
    rel_3rd = abs(decomp.delta_3rd_dressed - decomp.delta_can_direct) / decomp.delta_can_direct
    assert rel_3rd < 5e-3
    assert decomp.delta_3rd_zeeman_uv_no_lambda < 1e-13
    rel_qq_lam = (
        abs(decomp.delta_3rd_qq_lambda_only - decomp.delta_3rd_dressed) / decomp.delta_3rd_dressed
    )
    assert rel_qq_lam < 1e-3
    assert decomp.delta_3rd_lambda_only > 0.0


def test_hc02_exact_fock_subspace_equivalence_n2() -> None:
    rng = np.random.default_rng(7)
    for _ in range(10):
        mat = rng.standard_normal((8, 8))
        skew = mat - mat.T
        try:
            err = evaluate_exact_fock_equivalence(skew, n_sites=2)
            assert err < 1e-12
            return
        except RuntimeError:
            continue
    raise AssertionError("Failed to sample an even-parity N=2 Majorana matrix")


def test_hc02_exact_and_vmc_gutzwiller_n8_and_zero_variance_n12() -> None:
    cell = make_honeycomb_2site()
    cfg = SolverConfig(k_grid_shape=(24, 24), tol_correlator=1e-10, tol_gauge=1e-10)

    st0, met0 = solve_self_consistent_mean_field(
        cell, ModelParameters(K=-1.0), make_qsl_seed(cell), cfg
    )
    assert met0.converged
    ex8 = evaluate_exact_gutzwiller(2, 2, 0.0, st0)
    assert ex8.e_proj_per_site == pytest.approx(-(1.0 + np.sqrt(5.0)) / 4.0, abs=1e-9)

    ex12 = evaluate_exact_gutzwiller(2, 3, 0.0, st0)
    assert ex12.e_proj_per_site == pytest.approx(-0.8166902250, abs=1e-8)
    vmc12 = evaluate_vmc_gutzwiller(
        2,
        3,
        0.0,
        st0,
        n_chains=12,
        burn_in_sweeps=8,
        samples_per_chain=20,
        preferred_bc=ex12.bc_twist,
    )
    assert vmc12.e_loc_std < 1e-10
    assert vmc12.e_proj_mean == pytest.approx(ex12.e_proj_per_site, abs=1e-9)

    h = 0.25
    comp = h / float(np.sqrt(3.0))
    p25 = ModelParameters(K=-1.0, h=(comp, comp, comp))
    st25, met25 = solve_self_consistent_mean_field(cell, p25, st0, cfg)
    assert met25.converged
    ex8_h = evaluate_exact_gutzwiller(2, 2, h, st25)
    assert ex8_h.e_proj_per_site >= -0.935247756233 - 1e-9  # Bounded below by ED(N=8, h=0.25)
    vmc8_h = evaluate_vmc_gutzwiller(
        2,
        2,
        h,
        st25,
        n_chains=50,
        burn_in_sweeps=15,
        samples_per_chain=40,
        preferred_bc=ex8_h.bc_twist,
    )
    half_width = max(0.01, 1.5 * (vmc8_h.e_proj_ci_high - vmc8_h.e_proj_ci_low))
    assert abs(vmc8_h.e_proj_mean - ex8_h.e_proj_per_site) < half_width

    lo, hi = bca_confidence_interval(np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0]))
    assert lo < hi
    const_lo, const_hi = bca_confidence_interval(np.ones(12))
    assert const_lo == const_hi == 1.0
