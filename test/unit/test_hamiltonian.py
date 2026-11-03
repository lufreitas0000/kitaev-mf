"""Test momentum-space Bloch-Majorana Hamiltonian builders across unit cells."""

import numpy as np
import pytest
import torch

from kitaev_mf.domain.lattice import (
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.hamiltonian import (
    build_bloch_hamiltonian,
    build_bloch_hamiltonian_2site,
    build_bloch_hamiltonian_4site,
    build_on_site_matrices,
)


def test_mean_field_state_validation() -> None:
    cell_2 = make_honeycomb_2site()
    state_valid = MeanFieldState.from_pure_kitaev(cell_2, u_gauge=1.0)
    assert state_valid.U.shape == (3, 4, 4)
    assert state_valid.V.shape == (2, 4, 4)
    assert state_valid.lambdas.shape == (2, 3)

    # Verify skew-symmetry check fails on non-skew V
    invalid_v = np.zeros((2, 4, 4), dtype=np.float64)
    invalid_v[0, 1, 2] = 1.0
    with pytest.raises(ValueError, match="skew-symmetric"):
        MeanFieldState(U=state_valid.U, V=invalid_v, lambdas=state_valid.lambdas)


def test_hermiticity_arbitrary_couplings(cpu_backend: TorchBackend) -> None:
    cell = make_honeycomb_2site()
    params = ModelParameters(
        J=0.1, K=(-1.0, -1.0, -1.2), Gamma=0.2, Gamma_prime=-0.05, D=0.03, h=(0.05, -0.02, 0.08)
    )

    rng = np.random.default_rng(42)
    u_rand = rng.standard_normal((3, 4, 4))
    v_raw = rng.standard_normal((2, 4, 4))
    v_skew = 0.5 * (v_raw - np.swapaxes(v_raw, -1, -2))
    lam_rand = rng.standard_normal((2, 3))
    state = MeanFieldState(U=u_rand, V=v_skew, lambdas=lam_rand)

    k_grid = cell.make_k_grid(6, 6)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)

    # Check H(k) = H^\dagger(k) to machine precision
    h_batch_adj = torch.conj(h_batch.transpose(-1, -2))
    max_hermitian_error = torch.max(torch.abs(h_batch - h_batch_adj)).item()
    assert max_hermitian_error < 1e-14


def test_particle_hole_skew_symmetry(cpu_backend: TorchBackend) -> None:
    """Verify Majorana basis particle-hole condition H(-k) = -H^*(k)."""
    cell = make_honeycomb_2site()
    params = ModelParameters(J=0.12, K=-1.0, Gamma=0.18, Gamma_prime=0.04, D=0.02, h=(0.04, 0.04, 0.04))
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=0.3)

    rng = np.random.default_rng(123)
    k_pts = rng.uniform(0.0, 1.0, size=(10, 2))
    k_pts_neg = -k_pts

    h_pos = build_bloch_hamiltonian(k_pts, cell, params, state, cpu_backend)
    h_neg = build_bloch_hamiltonian(k_pts_neg, cell, params, state, cpu_backend)

    # H(-k) + H^*(k) == 0
    ph_error = torch.max(torch.abs(h_neg + torch.conj(h_pos))).item()
    assert ph_error < 1e-14


def test_pure_kitaev_analytical_dispersion(cpu_backend: TorchBackend) -> None:
    """Verify matter fermion dispersion matches Freitas et al. (2024) and gauge bands are flat."""
    k_val = -1.0
    u_gauge = 1.0
    u_matter = 0.42
    params = ModelParameters(K=k_val)
    cell = make_honeycomb_2site()
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=u_gauge, u_matter=u_matter)

    k_grid = cell.make_k_grid(12, 12)
    h_batch = build_bloch_hamiltonian_2site(k_grid, params, state, cpu_backend)
    eigenvalues, _ = cpu_backend.batched_eigh(h_batch)
    e_np = cpu_backend.to_numpy(eigenvalues)

    # Analytical dispersion for matter fermion c^0:
    # E(k) = +/- (|K| / 2) * |1 + exp(-i * 2pi * k1) + exp(-i * 2pi * k2)|
    k1 = k_grid[:, 0]
    k2 = k_grid[:, 1]
    f_k = 1.0 + np.exp(-1j * 2.0 * np.pi * k1) + np.exp(-1j * 2.0 * np.pi * k2)
    e_matter_expected = 0.5 * abs(k_val) * np.abs(f_k)

    # Analytical flat band for gauge fermions c^gamma:
    # E_gauge = +/- (|K| / 2) * |u_matter|
    e_gauge_expected = 0.5 * abs(k_val) * abs(u_matter)

    # At Gamma point k = (0, 0), |f(0)| = 3.0 => E(0) = 1.5 * |K|
    idx_gamma = 0
    assert np.isclose(e_np[idx_gamma, 7], 1.5 * abs(k_val), atol=1e-14)
    assert np.isclose(e_np[idx_gamma, 0], -1.5 * abs(k_val), atol=1e-14)

    for i in range(len(k_grid)):
        vals = np.sort(e_np[i])
        # Outer bands are matter fermions (or flat gauge if matter dispersion drops below)
        expected_spectrum = np.array(
            [
                -e_matter_expected[i],
                -e_gauge_expected,
                -e_gauge_expected,
                -e_gauge_expected,
                e_gauge_expected,
                e_gauge_expected,
                e_gauge_expected,
                e_matter_expected[i],
            ]
        )
        expected_spectrum = np.sort(expected_spectrum)
        assert np.allclose(vals, expected_spectrum, atol=1e-14)


def test_brillouin_zone_folding_parity_pure_kitaev(cpu_backend: TorchBackend) -> None:
    """Verify 2-site folded spectrum matches 4-site rectangular spectrum to machine precision."""
    params = ModelParameters(K=-1.0)
    cell_2 = make_honeycomb_2site()
    cell_4 = make_honeycomb_4site_rectangular()

    state_2 = MeanFieldState.from_pure_kitaev(cell_2, u_gauge=1.0, u_matter=0.35)
    state_4 = MeanFieldState.from_pure_kitaev(cell_4, u_gauge=1.0, u_matter=0.35)

    # Test random sampling of rectangular BZ momenta (K1, K2)
    rng = np.random.default_rng(777)
    k4_pts = rng.uniform(0.0, 1.0, size=(10, 2))

    for K1, K2 in k4_pts:
        # Folding relations:
        # k1 = (K1 + K2) / 2, k2 = (K2 - K1) / 2
        # Second folded point is shifted by (1/2, 1/2)
        k2_pt_a = np.array([(K1 + K2) / 2.0 % 1.0, (K2 - K1) / 2.0 % 1.0])
        k2_pt_b = np.array([(K1 + K2 + 1.0) / 2.0 % 1.0, (K2 - K1 + 1.0) / 2.0 % 1.0])

        h4 = build_bloch_hamiltonian_4site(np.array([K1, K2]), params, state_4, cpu_backend)
        h2_a = build_bloch_hamiltonian_2site(k2_pt_a, params, state_2, cpu_backend)
        h2_b = build_bloch_hamiltonian_2site(k2_pt_b, params, state_2, cpu_backend)

        e4, _ = cpu_backend.eigh(h4)
        e2_a, _ = cpu_backend.eigh(h2_a)
        e2_b, _ = cpu_backend.eigh(h2_b)

        vals_4 = np.sort(cpu_backend.to_numpy(e4))
        vals_2_folded = np.sort(np.concatenate([cpu_backend.to_numpy(e2_a), cpu_backend.to_numpy(e2_b)]))

        assert np.allclose(vals_4, vals_2_folded, atol=1e-14)


def test_brillouin_zone_folding_parity_general_couplings(cpu_backend: TorchBackend) -> None:
    """Verify folding parity holds under general anisotropic J-K-Gamma-Gamma'-D-h parameters."""
    params = ModelParameters(
        J=0.08,
        K=(-1.0, -1.05, -0.95),
        Gamma=0.15,
        Gamma_prime=-0.03,
        D=0.025,
        h=(0.04, 0.02, 0.06),
    )
    cell_2 = make_honeycomb_2site()
    cell_4 = make_honeycomb_4site_rectangular()

    # Use uniform states across cells
    state_2 = MeanFieldState.from_pure_kitaev(cell_2, u_gauge=1.0, u_matter=0.2)
    state_4 = MeanFieldState.from_pure_kitaev(cell_4, u_gauge=1.0, u_matter=0.2)

    rng = np.random.default_rng(888)
    k4_pts = rng.uniform(0.0, 1.0, size=(5, 2))

    for K1, K2 in k4_pts:
        k2_pt_a = np.array([(K1 + K2) / 2.0 % 1.0, (K2 - K1) / 2.0 % 1.0])
        k2_pt_b = np.array([(K1 + K2 + 1.0) / 2.0 % 1.0, (K2 - K1 + 1.0) / 2.0 % 1.0])

        h4 = build_bloch_hamiltonian_4site(np.array([K1, K2]), params, state_4, cpu_backend)
        h2_a = build_bloch_hamiltonian_2site(k2_pt_a, params, state_2, cpu_backend)
        h2_b = build_bloch_hamiltonian_2site(k2_pt_b, params, state_2, cpu_backend)

        e4, _ = cpu_backend.eigh(h4)
        e2_a, _ = cpu_backend.eigh(h2_a)
        e2_b, _ = cpu_backend.eigh(h2_b)

        vals_4 = np.sort(cpu_backend.to_numpy(e4))
        vals_2_folded = np.sort(np.concatenate([cpu_backend.to_numpy(e2_a), cpu_backend.to_numpy(e2_b)]))

        assert np.allclose(vals_4, vals_2_folded, atol=1e-14)


def test_on_site_b_matrix_hartree_and_dm_transposition() -> None:
    """Verify on-site B_i matrix correctly transposes J_{ji} = J_{ij}^T on reverse bonds."""
    cell = make_honeycomb_2site()
    # Pure Dzyaloshinskii-Moriya coupling makes exchange matrix purely antisymmetric off-diagonal
    params = ModelParameters(J=0.0, K=0.0, Gamma=0.0, D=0.5, h=(0.1, -0.2, 0.3))

    num_sub = 2
    v_arr = np.zeros((num_sub, 4, 4), dtype=np.float64)
    # Give site B (index 1) an unpolarized magnetization in the z-direction: sigma^3 = i c^0 c^3
    # N^3 element (0, 3) = 1, (3, 0) = -1
    v_arr[1, 0, 3] = 0.5
    v_arr[1, 3, 0] = -0.5

    lam_arr = np.array([[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]], dtype=np.float64)
    state = MeanFieldState(U=np.zeros((3, 4, 4)), V=v_arr, lambdas=lam_arr)

    b_mats = build_on_site_matrices(cell, params, state, torch.device("cpu"))
    b_np = b_mats.detach().cpu().numpy()

    # B_i must be strictly real skew-symmetric
    for s in range(num_sub):
        assert np.allclose(b_np[s], -b_np[s].T, atol=1e-15)
