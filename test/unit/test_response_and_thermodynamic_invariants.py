"""Unit tests for response functions, thermodynamic susceptibility, Klein duality, and Riemannian gradient (K-010, K-014).

Verifies:
1. Hellmann-Feynman theorem: dE/dh^gamma = - <sum_i sigma_i^gamma> via central differences.
2. Thermodynamic linear susceptibility: chi^{gamma,gamma} = d<sigma^gamma>/dh^gamma > 0.
3. Four-sublattice Klein duality energy symmetry: |E(p) - E(dual)| < tol.
4. Riemannian gradient stationarity commutator: ||[H(k), P(k)]||_F < 1e-10.
5. Input validations on direction indices and finite-difference step sizes.
"""

import numpy as np
import pytest

from kitaev_mf.domain.lattice import UnitCell, make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.correlators import compute_spectral_projector
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.response import (
    compute_riemannian_gradient_norm,
    evaluate_finite_difference_magnetization,
    evaluate_klein_dual_pair,
    evaluate_thermodynamic_susceptibility,
)


@pytest.fixture
def primitive_honeycomb_cell() -> UnitCell:
    """Two-sublattice primitive honeycomb unit cell."""
    return make_honeycomb_2site()


@pytest.fixture
def fast_solver_config() -> SolverConfig:
    """Fast solver configuration for response evaluation."""
    return SolverConfig(
        max_iterations=30,
        tol_correlator=1e-7,
        tol_gauge=1e-7,
        k_grid_shape=(6, 6),
        method="diis",
        diis_history_depth=4,
    )


def test_hellmann_feynman_central_difference(
    primitive_honeycomb_cell: UnitCell, fast_solver_config: SolverConfig
) -> None:
    """Verify Hellmann-Feynman theorem dE/dh^z = - <sigma_tot^z>."""
    params = ModelParameters(
        K=-1.0,
        h=(0.04, 0.04, 0.04),
    )
    res = evaluate_finite_difference_magnetization(
        cell=primitive_honeycomb_cell,
        params=params,
        h_dir=2,
        delta=0.005,
        config=fast_solver_config,
    )

    assert abs(res.discrepancy) < 5e-3, f"Hellmann-Feynman discrepancy too large: {res.discrepancy}"
    assert np.isclose(res.de_numerical, res.neg_sigma_tot, atol=5e-3)


def test_thermodynamic_susceptibility_positivity(
    primitive_honeycomb_cell: UnitCell, fast_solver_config: SolverConfig
) -> None:
    """Verify that thermodynamic linear susceptibility chi^{z,z} > 0."""
    params = ModelParameters(
        K=-1.0,
        h=(0.03, 0.03, 0.03),
    )
    chi_zz = evaluate_thermodynamic_susceptibility(
        cell=primitive_honeycomb_cell,
        params=params,
        h_dir=2,
        delta=0.005,
        config=fast_solver_config,
    )

    assert chi_zz > 0.0, f"Susceptibility must be positive (paramagnetic), got {chi_zz}"


def test_klein_duality_energy_invariance(
    primitive_honeycomb_cell: UnitCell, fast_solver_config: SolverConfig
) -> None:
    """Verify Klein duality energy symmetry across symmetric perturbation pair."""
    res = evaluate_klein_dual_pair(
        cell=primitive_honeycomb_cell,
        K0=-1.0,
        eps_J=0.02,
        eps_Gamma=0.01,
        eps_Gamma_prime=0.0,
        config=fast_solver_config,
    )

    assert res.energy_diff < 1e-4, f"Klein dual energy difference too large: {res.energy_diff}"


def test_riemannian_gradient_stationarity_at_solution(
    primitive_honeycomb_cell: UnitCell, fast_solver_config: SolverConfig
) -> None:
    """Verify that [H(k), P(k)] vanishes for converged mean-field fixed point."""
    k_pts = np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float64)
    # Assemble single-particle Hamiltonian
    U_test = np.zeros((3, 4, 4), dtype=np.float64)
    U_test[0, 3, 3] = 1.0
    U_test[1, 1, 1] = 1.0
    U_test[2, 2, 2] = 1.0
    V_test = np.zeros((2, 4, 4), dtype=np.float64)
    lam_test = np.zeros((2, 3), dtype=np.float64)

    state_test = MeanFieldState(U=U_test, V=V_test, lambdas=lam_test)
    H_torch = build_bloch_hamiltonian(
        k_pts=k_pts,
        cell=primitive_honeycomb_cell,
        params=ModelParameters(K=-1.0),
        state=state_test,
    )
    P_torch = compute_spectral_projector(H_torch)

    grad_norm = compute_riemannian_gradient_norm(H_torch, P_torch)
    assert grad_norm < 1e-12, f"Commutator norm must vanish for eigenprojector, got {grad_norm}"

    # Also test 2D single k-point input
    grad_norm_2d = compute_riemannian_gradient_norm(H_torch[0], P_torch[0])
    assert grad_norm_2d < 1e-12


def test_response_service_input_validations(
    primitive_honeycomb_cell: UnitCell, fast_solver_config: SolverConfig
) -> None:
    """Verify defensive exceptions on invalid direction or step size."""
    params = ModelParameters(K=-1.0)
    with pytest.raises(ValueError, match="h_dir must be in"):
        evaluate_finite_difference_magnetization(
            primitive_honeycomb_cell, params, h_dir=3, delta=0.01, config=fast_solver_config
        )

    with pytest.raises(ValueError, match="delta must be strictly positive"):
        evaluate_finite_difference_magnetization(
            primitive_honeycomb_cell, params, h_dir=0, delta=-0.01, config=fast_solver_config
        )

    with pytest.raises(ValueError, match="h_dir must be in"):
        evaluate_thermodynamic_susceptibility(
            primitive_honeycomb_cell, params, h_dir=-1, delta=0.01, config=fast_solver_config
        )

    with pytest.raises(ValueError, match="delta must be strictly positive"):
        evaluate_thermodynamic_susceptibility(
            primitive_honeycomb_cell, params, h_dir=0, delta=0.0, config=fast_solver_config
        )
