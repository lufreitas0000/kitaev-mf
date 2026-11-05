"""Integration test suite for variational response functions and Hellmann-Feynman scaling (T3).

Contracts:
    TST-VAR-01: Variational Hellmann-Feynman magnetization derivative and O(delta^2) error scaling
    TST-VAR-02: Thermodynamic magnetic concavity and non-negative susceptibility
"""

from math import log2

import pytest

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.correlators import compute_spectral_projector
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.response import (
    compute_riemannian_gradient_norm,
    evaluate_finite_difference_magnetization,
    evaluate_thermodynamic_susceptibility,
)


def test_tst_var_01_hellmann_feynman_quadratic_scaling(cpu_backend: TorchBackend) -> None:
    """TST-VAR-01: Variational Hellmann-Feynman derivative and quadratic truncation scaling.

    Asserts:
        1. At delta = 0.002, discrepancy |(Delta E / Delta h^z) - (-<sigma_tot^z>)| < 5e-7.
        2. Observed convergence order p in [1.90, 2.10] across delta in {0.004, 0.002}.
        3. Riemannian KKT commutator norm ||[H(k), P(k)]||_F < 1e-13 at self-consistency.
    """
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, J=0.04, h=(0.02, 0.02, 0.03))
    config = SolverConfig(
        k_grid_shape=(16, 16),
        tol_correlator=1e-9,
        tol_gauge=1e-9,
        max_iterations=60,
    )

    delta_coarse = 0.004
    delta_fine = 0.002

    res_coarse = evaluate_finite_difference_magnetization(
        cell=cell,
        params=params,
        h_dir=2,
        delta=delta_coarse,
        config=config,
        backend=cpu_backend,
    )
    res_fine = evaluate_finite_difference_magnetization(
        cell=cell,
        params=params,
        h_dir=2,
        delta=delta_fine,
        config=config,
        backend=cpu_backend,
    )

    # 1. Error magnitude at delta = 0.002
    assert res_fine.discrepancy < 5e-7

    # 2. Asymptotic quadratic scaling order p = log2(Delta_coarse / Delta_fine) ~ 2.0
    order_p = log2(res_coarse.discrepancy / res_fine.discrepancy)
    assert 1.90 <= order_p <= 2.10

    # 3. KKT Riemannian gradient norm on the base converged state
    k_grid = cell.make_k_grid(16, 16)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, res_fine.state_base, cpu_backend)
    p_batch = compute_spectral_projector(h_batch, cpu_backend)
    kkt_norm = compute_riemannian_gradient_norm(h_batch, p_batch)
    assert kkt_norm < 1e-13


def test_tst_var_02_thermodynamic_susceptibility_positive(cpu_backend: TorchBackend) -> None:
    """TST-VAR-02: Thermodynamic magnetic concavity and positive susceptibility.

    Asserts:
        chi^zz = d<sigma_tot^z>/dh^z > 0.5 under finite field along the z-axis.
    """
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, J=0.03, h=(0.0, 0.0, 0.03))
    config = SolverConfig(
        k_grid_shape=(16, 16),
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        max_iterations=50,
    )

    chi_zz = evaluate_thermodynamic_susceptibility(
        cell=cell,
        params=params,
        h_dir=2,
        delta=1e-3,
        config=config,
        backend=cpu_backend,
    )

    assert chi_zz > 0.5


def test_response_service_input_validation(cpu_backend: TorchBackend) -> None:
    """Verify error contracts on invalid Cartesian indices and non-positive deltas."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    config = SolverConfig(k_grid_shape=(8, 8), max_iterations=5)

    with pytest.raises(ValueError, match="h_dir must be in"):
        evaluate_finite_difference_magnetization(
            cell=cell, params=params, h_dir=3, delta=1e-3, config=config, backend=cpu_backend
        )

    with pytest.raises(ValueError, match="delta must be strictly positive"):
        evaluate_finite_difference_magnetization(
            cell=cell, params=params, h_dir=0, delta=0.0, config=config, backend=cpu_backend
        )

    with pytest.raises(ValueError, match="h_dir must be in"):
        evaluate_thermodynamic_susceptibility(
            cell=cell, params=params, h_dir=-1, delta=1e-3, config=config, backend=cpu_backend
        )

    with pytest.raises(ValueError, match="delta must be strictly positive"):
        evaluate_thermodynamic_susceptibility(
            cell=cell, params=params, h_dir=1, delta=-1e-4, config=config, backend=cpu_backend
        )
