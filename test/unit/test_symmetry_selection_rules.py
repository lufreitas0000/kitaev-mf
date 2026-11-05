"""Unit tests for Tier-2 discrete symmetry selection rules.

Contracts:
    TST-SYM-01: Time-reversal-odd correlator vanishing
    TST-SYM-02: Sublattice inversion parity
    TST-SYM-03: C3 ray permutation invariance
"""

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.correlators import compute_correlators_from_bloch
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian


def test_tst_sym_01_time_reversal_odd_correlators_vanish(cpu_backend: TorchBackend) -> None:
    """TST-SYM-01: Time-reversal-odd correlator vanishing.

    Under time-reversal symmetry Theta (h = 0, D = 0), matter Majoranas c^0
    are even and gauge Majoranas c^gamma are odd. All hybrid matter-gauge
    two-point correlators U_b^{0,gamma} and U_b^{gamma,0} must vanish identically.
    """
    cell = make_honeycomb_2site()
    # Non-zero J, K, Gamma, Gamma' preserving time reversal
    params = ModelParameters(
        J=0.06,
        K=-1.0,
        Gamma=0.15,
        Gamma_prime=-0.04,
        D=0.0,
        h=(0.0, 0.0, 0.0),
    )
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.35)

    k_grid = cell.make_k_grid(14, 14)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)
    u_out, _ = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)

    for b_idx in range(len(cell.bonds)):
        u_b = u_out[b_idx]
        for gamma in (1, 2, 3):
            # Hybrid matter-gauge cross hoppings must vanish
            assert np.isclose(u_b[0, gamma], 0.0, atol=1e-14)
            assert np.isclose(u_b[gamma, 0], 0.0, atol=1e-14)


def test_tst_sym_02_sublattice_inversion_parity(cpu_backend: TorchBackend) -> None:
    """TST-SYM-02: Sublattice inversion parity.

    Inversion symmetry maps sublattice A to B. For spatially homogeneous
    parameters with zero staggered fields, on-site correlation matrices
    on sublattices A and B must be strictly identical: V_A == V_B.
    """
    cell = make_honeycomb_2site()
    params = ModelParameters(
        J=0.05,
        K=-1.0,
        Gamma=0.1,
        h=(0.03, 0.03, 0.03),
    )
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.35)

    k_grid = cell.make_k_grid(16, 16)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)
    _, v_out = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)

    v_a = v_out[0]
    v_b = v_out[1]

    # V_A must equal V_B to machine precision
    diff_inv = float(np.max(np.abs(v_a - v_b)))
    assert diff_inv < 1e-14


def test_tst_sym_03_c3_ray_permutation_invariance(cpu_backend: TorchBackend) -> None:
    """TST-SYM-03: C3 ray permutation invariance.

    Under C3 rotation about a plaquette center, x-bond -> y-bond -> z-bond.
    For isotropic couplings and [1, 1, 1] field:
    1. Matter bond hoppings are identical: U_x^{00} == U_y^{00} == U_z^{00}.
    2. Cyclic gauge link hoppings match: U_x^{11} == U_y^{22} == U_z^{33}.
    """
    cell = make_honeycomb_2site()
    params = ModelParameters(
        J=0.04,
        K=-1.0,
        Gamma=0.08,
        h=(0.02, 0.02, 0.02),
    )
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.35)

    # Use 36x36 grid to preserve C3 symmetry on the torus
    k_grid = cell.make_k_grid(36, 36)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)
    u_out, _ = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)

    # In make_honeycomb_2site(): bond 0 is z, bond 1 is x, bond 2 is y
    u_z = u_out[0]
    u_x = u_out[1]
    u_y = u_out[2]

    # Matter hoppings must be invariant under C3
    assert np.isclose(u_x[0, 0], u_y[0, 0], atol=1e-12)
    assert np.isclose(u_y[0, 0], u_z[0, 0], atol=1e-12)

    # Cyclic gauge hoppings must match: U_x^{11} == U_y^{22} == U_z^{33}
    assert np.isclose(u_x[1, 1], u_y[2, 2], atol=1e-12)
    assert np.isclose(u_y[2, 2], u_z[3, 3], atol=1e-12)
