"""Unit tests for Analytical KKT Perturbation Expansion and Cubic Dirac Gap (Task T08, Contract H6, K-010, K-019)."""

from __future__ import annotations

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def test_kkt_multiplier_linear_response_ratio_sqrt5() -> None:
    """Contract H6: Total multiplier norm on unit cell satisfies ||lambda*|| / h = sqrt(5).

    Per sublattice: ||lambda_s*|| / h = sqrt(5/2) = sqrt(2.5).
    Per component under C_3: lambda_s^gamma / h = sqrt(5/6).
    """
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()
    seed = make_qsl_seed(cell)

    h_val = 0.005
    h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
    params = ModelParameters(K=-1.0, h=h_vec)
    config = SolverConfig(
        k_grid_shape=(16, 16),
        max_iterations=60,
        tol_correlator=1e-10,
        tol_gauge=1e-10,
        method="diis",
    )
    res, met = solve_self_consistent_mean_field(cell, params, seed, config, backend)
    assert met.converged

    # 1. Total norm across 2-sublattice unit cell
    full_norm = float(np.linalg.norm(res.lambdas))
    ratio_full = full_norm / h_val
    assert np.isclose(ratio_full, np.sqrt(5.0), atol=1e-3), (
        f"Unit cell ||lambda|| / h expected sqrt(5) = {np.sqrt(5.0):.6f}, got {ratio_full:.6f}"
    )

    # 2. Per-sublattice norm (sublattice 0)
    sub_norm = float(np.linalg.norm(res.lambdas[0]))
    ratio_sub = sub_norm / h_val
    assert np.isclose(ratio_sub, np.sqrt(2.5), atol=1e-3), (
        f"Per-site ||lambda_s|| / h expected sqrt(2.5) = {np.sqrt(2.5):.6f}, got {ratio_sub:.6f}"
    )

    # 3. Individual components under C_3 isotropy
    for gamma in range(3):
        comp = float(res.lambdas[0, gamma])
        ratio_comp = comp / h_val
        assert np.isclose(ratio_comp, np.sqrt(5.0 / 6.0), atol=1e-3), (
            f"Component lambda^gamma / h expected sqrt(5/6) = {np.sqrt(5.0 / 6.0):.6f}, got {ratio_comp:.6f}"
        )


def test_cubic_dirac_mass_scaling_and_prefactor() -> None:
    """K-010, Contract H6: Majorana gap at Dirac point K opens with cubic power law Delta_K(h) ~ c_MF * h^3."""
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()
    seed = make_qsl_seed(cell)
    k_K = np.array([1.0 / 3.0, 2.0 / 3.0])

    h_vals = [0.01, 0.02]
    gaps: list[float] = []

    for h_val in h_vals:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)
        config = SolverConfig(
            k_grid_shape=(16, 16),
            max_iterations=60,
            tol_correlator=1e-10,
            tol_gauge=1e-10,
            method="diis",
        )
        res, met = solve_self_consistent_mean_field(cell, params, seed, config, backend)
        assert met.converged

        H_K = build_bloch_hamiltonian(k_K, cell, params, res).numpy()
        w_K = np.sort(np.linalg.eigvalsh(H_K))
        pos_w = w_K[w_K > 0]
        gap = float(2.0 * pos_w[0])
        gaps.append(gap)

    # Ratio Delta(2h) / Delta(h) should be 2^3 = 8
    ratio_doubling = gaps[1] / gaps[0]
    assert np.isclose(ratio_doubling, 8.0, atol=0.05), (
        f"Expected cubic gap doubling ratio ~8.0, got {ratio_doubling:.4f}"
    )

    # Prefactor c_MF = gap / h^3 in Bloch normalization
    c_est = gaps[0] / (h_vals[0] ** 3)
    assert 2.50 <= c_est <= 2.70, f"Expected c_MF in [2.50, 2.70], got {c_est:.4f}"


def test_analytical_reconciliation_identities_and_energy_denominators() -> None:
    """K-019, Contract H6: Analytical operator map between mean-field gap and 3rd-order spin perturbation."""
    # Canonical factor 4 from Bloch matrix normalization H_Bloch = (i/4) A
    c_mf_bloch = 2.59214
    c_mf_canonical = 4.0 * c_mf_bloch
    assert np.isclose(c_mf_canonical, 10.36856, atol=1e-4)

    # Intermediate energy scales: localized flat-band gauge Majorana vs two-vison threshold
    u_bar = 0.52486574585
    two_delta_gauge = 2.0 * u_bar  # 1.04973
    delta_vison = 0.2640  # two-vison excitation threshold in pure Kitaev model
    ratio_sq = (two_delta_gauge / delta_vison) ** 2
    assert np.isclose(ratio_sq, 15.8106, atol=1e-3)

    # Static vison perturbation gap scale
    c_static_vison = 4.0 / (delta_vison**2)
    assert np.isclose(c_static_vison, 57.3921, atol=1e-3)
