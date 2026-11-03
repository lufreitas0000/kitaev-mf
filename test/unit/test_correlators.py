"""Test reciprocal spectral projection, real-space BZ correlators, and gauge residuals."""

import numpy as np

from kitaev_mf.domain.lattice import (
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.correlators import (
    compute_correlators_from_bloch,
    compute_gauge_residuals,
    compute_spectral_projector,
)
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian


def test_correlator_reality_and_skew_symmetry(cpu_backend: TorchBackend) -> None:
    """Verify that computed real-space U and V are float64 real and V is skew-symmetric."""
    cell = make_honeycomb_2site()
    params = ModelParameters(
        J=0.08,
        K=(-1.0, -1.05, -0.95),
        Gamma=0.15,
        Gamma_prime=-0.03,
        D=0.025,
        h=(0.04, 0.02, 0.06),
    )
    rng = np.random.default_rng(101)
    u_init = rng.standard_normal((len(cell.bonds), 4, 4))
    v_raw = rng.standard_normal((cell.num_sublattices, 4, 4))
    v_init = 0.5 * (v_raw - np.swapaxes(v_raw, -1, -2))
    state = MeanFieldState(U=u_init, V=v_init, lambdas=np.zeros((cell.num_sublattices, 3)))

    k_grid = cell.make_k_grid(8, 8)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)

    U_out, V_out = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)

    assert U_out.shape == (3, 4, 4)
    assert V_out.shape == (2, 4, 4)
    assert np.issubdtype(U_out.dtype, np.floating)
    assert np.issubdtype(V_out.dtype, np.floating)

    skew_error = np.max(np.abs(V_out + np.swapaxes(V_out, -1, -2)))
    assert skew_error < 1e-14


def test_car_trace_normalization_sum_rule(cpu_backend: TorchBackend) -> None:
    """Verify CAR sum rule: (1/4) * tr(<c_{s,R} c_{s,R}^T>) = 1.0 on all sublattices."""
    cell = make_honeycomb_2site()
    params = ModelParameters(J=0.1, K=-1.0, Gamma=0.2, h=(0.05, 0.05, 0.05))
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.3)

    k_grid = cell.make_k_grid(12, 12)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)

    p_tensor = compute_spectral_projector(h_batch, cpu_backend)
    p_np = cpu_backend.to_numpy(p_tensor)

    for s in range(cell.num_sublattices):
        diag_block = p_np[:, 4 * s : 4 * s + 4, 4 * s : 4 * s + 4].real
        avg_trace = np.mean(np.trace(diag_block, axis1=-2, axis2=-1)) / 4.0
        assert np.isclose(avg_trace, 1.0, atol=1e-14)


def test_pure_kitaev_correlator_invariant_k006(cpu_backend: TorchBackend) -> None:
    """Verify invariant K-006: U_b^{00} approx -0.52486 and U_b^{gg} = +1.0 for pure Kitaev."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    # Physical ground state has negative matter hopping for K = -1, u_g = +1
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.3)

    k_grid = cell.make_k_grid(36, 36)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)

    U_out, V_out = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)

    for b_idx, bond in enumerate(cell.bonds):
        gamma = bond.bond_type
        assert np.isclose(U_out[b_idx, 0, 0], -0.5248645, atol=3e-3)
        assert np.isclose(U_out[b_idx, gamma, gamma], 1.0, atol=1e-12)

    assert np.allclose(V_out, 0.0, atol=1e-14)
    residuals = compute_gauge_residuals(V_out)
    assert np.allclose(residuals, 0.0, atol=1e-14)


def test_zero_mode_exact_gauge_neutrality(cpu_backend: TorchBackend) -> None:
    """Verify that flat gauge bands at zero energy (u_matter=0) yield zero gauge residuals."""
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0)
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=0.0)

    k_grid = cell.make_k_grid(16, 16)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)

    _, V_out = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)
    residuals = compute_gauge_residuals(V_out)
    assert np.allclose(residuals, 0.0, atol=1e-14)


def test_brillouin_zone_folding_correlator_parity(cpu_backend: TorchBackend) -> None:
    """Verify 2-site folded correlators match 4-site rectangular correlators to Riemann mesh accuracy."""
    cell_2 = make_honeycomb_2site()
    cell_4 = make_honeycomb_4site_rectangular()
    params = ModelParameters(K=-1.0, Gamma=0.1, h=(0.02, 0.02, 0.02))

    state_2 = MeanFieldState.from_pure_kitaev(cell_2, u_gauge=1.0, u_matter=-0.25)
    state_4 = MeanFieldState.from_pure_kitaev(cell_4, u_gauge=1.0, u_matter=-0.25)

    k2_grid = cell_2.make_k_grid(36, 36)
    k4_grid = cell_4.make_k_grid(18, 18)

    h2_batch = build_bloch_hamiltonian(k2_grid, cell_2, params, state_2, cpu_backend)
    h4_batch = build_bloch_hamiltonian(k4_grid, cell_4, params, state_4, cpu_backend)

    U2, V2 = compute_correlators_from_bloch(h2_batch, cell_2, k2_grid, cpu_backend)
    U4, V4 = compute_correlators_from_bloch(h4_batch, cell_4, k4_grid, cpu_backend)

    # Site 0 in 2-site corresponds to Sublattice 0 (A1) in 4-site (Riemann mesh error < 2e-3)
    assert np.allclose(V2[0], V4[0], atol=2e-3)
    # Bond 0 (z-bond) in 2-site corresponds to Bond 0 in 4-site
    assert np.allclose(U2[0], U4[0], atol=2e-3)
