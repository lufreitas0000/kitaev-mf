"""Integration test suite for solvable dualities and boundary asymptotics (T3).

Contracts:
    TST-RED-03: Bipartite sign-inversion duality (K = -1 vs K = +1)
    TST-DUA-01: Symmetrically parameterized four-sublattice Klein duality
"""

import numpy as np

from kitaev_mf.domain.lattice import (
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.response import evaluate_klein_dual_pair
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def test_tst_red_03_bipartite_sign_inversion(cpu_backend: TorchBackend) -> None:
    """TST-RED-03: Bipartite sign-inversion duality in the pure Kitaev limit.

    The bipartite gauge transformation c_A^0 -> +c_A^0, c_B^0 -> -c_B^0 maps K -> -K.
    In the canonical zero-flux sector (U_b^{gg} = +1.0 for directed A -> B bonds):
        1. Single-particle Bloch spectra spec(H(k)) are strictly identical.
        2. Ground-state energies per cell and site are strictly identical.
        3. Matter hoppings satisfy U_b^{00}(K=+1) = -U_b^{00}(K=-1).
    """
    cell = make_honeycomb_2site()
    params_neg = ModelParameters(K=-1.0)
    params_pos = ModelParameters(K=1.0)

    # Crucial gauge sign: matter hopping is negative for K = -1, positive for K = +1
    init_neg = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5)
    init_pos = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=0.5)

    config = SolverConfig(
        k_grid_shape=(24, 24),
        max_iterations=10,
        tol_correlator=1e-10,
        tol_gauge=1e-10,
        mixing_alpha=0.8,
    )

    state_neg, metrics_neg = solve_self_consistent_mean_field(
        cell, params_neg, init_neg, config, cpu_backend
    )
    state_pos, metrics_pos = solve_self_consistent_mean_field(
        cell, params_pos, init_pos, config, cpu_backend
    )

    assert metrics_neg.converged is True
    assert metrics_pos.converged is True

    # 1. Energy equality
    e_neg = compute_gaussian_energy(state_neg, cell, params_neg)
    e_pos = compute_gaussian_energy(state_pos, cell, params_pos)
    assert np.isclose(e_neg.per_cell, e_pos.per_cell, atol=1e-13)
    assert np.isclose(e_neg.per_site, e_pos.per_site, atol=1e-13)

    # 2. Matter bond correlator sign inversion and gauge link invariance
    for b_idx, bond in enumerate(cell.bonds):
        gamma = bond.bond_type
        # U^{00}(+1) == -U^{00}(-1)
        assert np.isclose(state_pos.U[b_idx, 0, 0], -state_neg.U[b_idx, 0, 0], atol=1e-13)
        # U^{gg}(+1) == U^{gg}(-1) == 1.0
        assert np.isclose(state_pos.U[b_idx, gamma, gamma], 1.0, atol=1e-12)
        assert np.isclose(state_neg.U[b_idx, gamma, gamma], 1.0, atol=1e-12)

    # 3. Spectral identity across BZ
    k_grid = cell.make_k_grid(12, 12)
    h_neg = build_bloch_hamiltonian(k_grid, cell, params_neg, state_neg, cpu_backend)
    h_pos = build_bloch_hamiltonian(k_grid, cell, params_pos, state_pos, cpu_backend)

    evals_neg, _ = cpu_backend.batched_eigh(h_neg)
    evals_pos, _ = cpu_backend.batched_eigh(h_pos)

    assert np.allclose(
        cpu_backend.to_numpy(evals_neg),
        cpu_backend.to_numpy(evals_pos),
        atol=1e-13,
    )


def test_tst_dua_01_four_sublattice_klein_duality(cpu_backend: TorchBackend) -> None:
    """TST-DUA-01: Symmetrically parameterized four-sublattice Klein duality.

    In the 4-sublattice rectangular cell, Klein duality maps:
        (J, K, Gamma, Gamma') -> (-J, K + 2J, -Gamma, Gamma')
    Fixing K + J = -1.0 with eps_J = eps_Gamma = 0.02 and eps_Gamma' = -0.01:
        Finite-cell boundary phase splitting satisfies |delta E_site| < 1e-5.
    """
    cell = make_honeycomb_4site_rectangular()
    config = SolverConfig(
        k_grid_shape=(12, 12),
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        max_iterations=50,
    )

    res = evaluate_klein_dual_pair(
        cell=cell,
        K0=-1.0,
        eps_J=0.02,
        eps_Gamma=0.02,
        eps_Gamma_prime=-0.01,
        config=config,
        backend=cpu_backend,
    )

    assert res.energy_diff < 1e-5
