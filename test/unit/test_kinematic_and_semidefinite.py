"""Unit tests for Tier-1 kinematic sum rules and semidefinite representability gates.

Contracts:
    TST-KIN-01: Parseval-Frobenius spectral moment identity
    TST-KIN-02: Pointwise Bloch Hamiltonian and spectral projector commutation
    TST-KIN-03: Spectral projector idempotency, trace filling, and reflection sum rule
    TST-SDP-01: Local two-site Majorana LMI positivity
    TST-SDP-02: Covariance operator norm and singular value bounds
    TST-SDP-03: Dipole Bloch-ball and fermion parity envelope bounds
"""

from collections.abc import Callable

import numpy as np
import pytest
import torch

from kitaev_mf.domain.lattice import (
    UnitCell,
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.correlators import (
    compute_correlators_from_bloch,
    compute_spectral_projector,
)
from kitaev_mf.services.hamiltonian import (
    build_bloch_hamiltonian,
    build_bond_hopping_matrices,
    build_on_site_matrices,
)


@pytest.fixture
def sample_perturbed_parameters() -> ModelParameters:
    """Generic non-integrable model couplings with broken TR and DM exchange."""
    return ModelParameters(
        J=0.08,
        K=(-1.0, -1.05, -0.95),
        Gamma=0.12,
        Gamma_prime=-0.04,
        D=0.03,
        h=(0.04, 0.02, 0.05),
    )


def _make_nonstationary_state(cell: UnitCell, seed: int = 42) -> MeanFieldState:
    """Construct an un-converged state with active cross-terms U^{0,gamma} and V != 0."""
    rng = np.random.default_rng(seed)
    num_bonds = len(cell.bonds)
    num_sub = cell.num_sublattices

    u_arr = 0.25 * rng.standard_normal((num_bonds, 4, 4))
    v_raw = 0.15 * rng.standard_normal((num_sub, 4, 4))
    v_arr = 0.5 * (v_raw - np.swapaxes(v_raw, -1, -2))
    lam_arr = 0.1 * rng.standard_normal((num_sub, 3))
    return MeanFieldState(U=u_arr, V=v_arr, lambdas=lam_arr)


# ---------------------------------------------------------------------------
# TST-KIN-01: Parseval-Frobenius Spectral Moment Identity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cell_factory", [make_honeycomb_2site, make_honeycomb_4site_rectangular])
def test_tst_kin_01_parseval_frobenius_spectral_moment(
    cpu_backend: TorchBackend,
    sample_perturbed_parameters: ModelParameters,
    cell_factory: Callable[[], UnitCell],
) -> None:
    """Verify BZ sum of squared eigenvalues matches sum of squared block norms.

    Condition of validity: N1 >= 2, N2 >= 2 to prevent Fourier aliasing
    between parallel bonds sharing the same sublattice endpoints.
    """
    cell = cell_factory()
    state = _make_nonstationary_state(cell, seed=101)
    k_grid = cell.make_k_grid(12, 12)
    device = cpu_backend.device

    # Single-particle real on-site and hopping blocks
    b_tensors = build_on_site_matrices(cell, sample_perturbed_parameters, state, device)
    a_tensors = build_bond_hopping_matrices(cell, sample_perturbed_parameters, state, device)

    b_np = b_tensors.detach().cpu().numpy()
    a_np = a_tensors.detach().cpu().numpy()

    # Sum of squared Frobenius norms: (1/16) * sum ||B_s||_F^2 + (1/8) * sum ||A_b||_F^2
    sum_b_sq = float(np.sum([np.trace(b @ b.T) for b in b_np]))
    sum_a_sq = float(np.sum([np.trace(a @ a.T) for a in a_np]))
    expected_moment = (1.0 / 16.0) * sum_b_sq + (1.0 / 8.0) * sum_a_sq

    # Numerical BZ average of sum_n E_n(k)^2
    h_batch = build_bloch_hamiltonian(k_grid, cell, sample_perturbed_parameters, state, cpu_backend)
    eigenvalues, _ = cpu_backend.batched_eigh(h_batch)
    e_np = eigenvalues.detach().cpu().numpy()

    actual_moment = float(np.mean(np.sum(e_np**2, axis=-1)))

    assert np.isclose(actual_moment, expected_moment, atol=1e-12, rtol=0.0)


# ---------------------------------------------------------------------------
# TST-KIN-02: Pointwise Commutation [H(k), P(k)] = 0
# ---------------------------------------------------------------------------


def test_tst_kin_02_pointwise_hamiltonian_projector_commutation(
    cpu_backend: TorchBackend,
    sample_perturbed_parameters: ModelParameters,
) -> None:
    """Verify [H(k), P(k)] == 0 pointwise across the entire Brillouin zone."""
    cell = make_honeycomb_2site()
    state = _make_nonstationary_state(cell, seed=202)
    k_grid = cell.make_k_grid(10, 10)

    h_batch = build_bloch_hamiltonian(k_grid, cell, sample_perturbed_parameters, state, cpu_backend)
    p_batch = compute_spectral_projector(h_batch, cpu_backend)

    # Commutator [H, P] = H @ P - P @ H
    hp = torch.matmul(h_batch, p_batch)
    ph = torch.matmul(p_batch, h_batch)
    comm = hp - ph

    max_comm = float(torch.max(torch.abs(comm)).item())
    assert max_comm < 1e-14


# ---------------------------------------------------------------------------
# TST-KIN-03: Projector Idempotency, Trace Filling, and Reflection Sum Rule
# ---------------------------------------------------------------------------


def test_tst_kin_03_projector_algebraic_properties(
    cpu_backend: TorchBackend,
    sample_perturbed_parameters: ModelParameters,
) -> None:
    """Verify trace filling, reflection sum rule, and idempotency for gapped states."""
    cell = make_honeycomb_2site()
    # High field guarantees a finite spectral gap across the entire BZ
    gapped_params = ModelParameters(K=-1.0, h=(0.25, 0.25, 0.25))
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.4)

    # 1. Pointwise idempotency P(k)^2 == 2 P(k) on a gapped dispersion
    k_grid = cell.make_k_grid(8, 8)
    h_batch = build_bloch_hamiltonian(k_grid, cell, gapped_params, state, cpu_backend)
    p_batch = compute_spectral_projector(h_batch, cpu_backend)

    p_sq = torch.matmul(p_batch, p_batch)
    idempotency_res = float(torch.max(torch.abs(p_sq - 2.0 * p_batch)).item())
    assert idempotency_res < 1e-13

    # 2. Brillouin-zone trace filling: (1/N_k) * sum_k tr(P(k)) == 4 * N_sub
    p_np = p_batch.detach().cpu().numpy()
    mean_trace = float(np.mean(np.trace(p_np.real, axis1=-2, axis2=-1)))
    expected_trace = 4.0 * cell.num_sublattices
    assert np.isclose(mean_trace, expected_trace, atol=1e-13, rtol=0.0)

    # 3. Reflection sum rule: P(-k) + P^T(k) == 2 * I
    k_single = np.array([[0.17, 0.39]], dtype=np.float64)
    k_neg = np.array([[-0.17, -0.39]], dtype=np.float64)

    h_pos = build_bloch_hamiltonian(k_single, cell, sample_perturbed_parameters, state, cpu_backend)
    h_neg = build_bloch_hamiltonian(k_neg, cell, sample_perturbed_parameters, state, cpu_backend)

    p_pos = compute_spectral_projector(h_pos, cpu_backend)[0]
    p_neg = compute_spectral_projector(h_neg, cpu_backend)[0]

    dim = 4 * cell.num_sublattices
    identity_dim = torch.eye(dim, dtype=torch.complex128, device=cpu_backend.device)

    reflection_matrix = p_neg + p_pos.transpose(-1, -2)
    diff_reflection = float(torch.max(torch.abs(reflection_matrix - 2.0 * identity_dim)).item())
    assert diff_reflection < 1e-14


# ---------------------------------------------------------------------------
# TST-SDP-01: Local Two-Site Majorana LMI Positivity
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cell_factory", [make_honeycomb_2site, make_honeycomb_4site_rectangular])
def test_tst_sdp_01_local_two_site_lmi_positivity(
    cpu_backend: TorchBackend,
    sample_perturbed_parameters: ModelParameters,
    cell_factory: Callable[[], UnitCell],
) -> None:
    """Verify C_ij = <c_ij c_ij^T> is Hermitian and positive semi-definite (lambda_min >= 0)."""
    cell = cell_factory()
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.3)
    k_grid = cell.make_k_grid(16, 16)

    h_batch = build_bloch_hamiltonian(k_grid, cell, sample_perturbed_parameters, state, cpu_backend)
    u_out, v_out = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)

    identity_4 = np.eye(4, dtype=np.complex128)

    for b_idx, bond in enumerate(cell.bonds):
        s_i = bond.sublattice_i
        s_j = bond.sublattice_j

        v_i = v_out[s_i]
        v_j = v_out[s_j]
        u_b = u_out[b_idx]

        # C_ij = [[I - i*V_i, -i*U_b], [i*U_b^T, I - i*V_j]]
        c_mat = np.zeros((8, 8), dtype=np.complex128)
        c_mat[:4, :4] = identity_4 - 1j * v_i
        c_mat[:4, 4:] = -1j * u_b
        c_mat[4:, :4] = 1j * u_b.T
        c_mat[4:, 4:] = identity_4 - 1j * v_j

        # Hermiticity
        herm_err = np.max(np.abs(c_mat - c_mat.conj().T))
        assert herm_err < 1e-14

        # Positive semi-definiteness: lambda_min >= 0
        eigenvalues = np.linalg.eigvalsh(c_mat)
        assert np.min(eigenvalues) >= -1e-12


# ---------------------------------------------------------------------------
# TST-SDP-02: Covariance Operator Norm and Singular Value Bounds
# ---------------------------------------------------------------------------


def test_tst_sdp_02_operator_norm_and_singular_value_bounds(
    cpu_backend: TorchBackend,
    sample_perturbed_parameters: ModelParameters,
) -> None:
    """Verify sigma_max(U_b) <= 1.0 and |nu_k(V_s)| <= 1.0."""
    cell = make_honeycomb_2site()
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.4)
    k_grid = cell.make_k_grid(16, 16)

    h_batch = build_bloch_hamiltonian(k_grid, cell, sample_perturbed_parameters, state, cpu_backend)
    u_out, v_out = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)

    # 1. Bond singular values: sigma_max(U_b) <= 1.0
    for b_idx in range(len(cell.bonds)):
        s_vals = np.linalg.svd(u_out[b_idx], compute_uv=False)
        assert np.max(s_vals) <= 1.0 + 1e-12

    # 2. On-site singular values: max |nu_k(V_s)| <= 1.0
    for s in range(cell.num_sublattices):
        v_s_vals = np.linalg.svd(v_out[s], compute_uv=False)
        assert np.max(v_s_vals) <= 1.0 + 1e-12


# ---------------------------------------------------------------------------
# TST-SDP-03: Dipole Bloch-Ball and Fermion Parity Envelopes
# ---------------------------------------------------------------------------


def test_tst_sdp_03_dipole_and_parity_envelopes(
    cpu_backend: TorchBackend,
    sample_perturbed_parameters: ModelParameters,
) -> None:
    """Verify ||<sigma_s>||_2 <= 1.0, Pf(V_s) in [-1, 1], and <T_s^2> in [0, 0.75]."""
    cell = make_honeycomb_2site()
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.4)
    k_grid = cell.make_k_grid(16, 16)

    h_batch = build_bloch_hamiltonian(k_grid, cell, sample_perturbed_parameters, state, cpu_backend)
    _, v_out = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)

    for s in range(cell.num_sublattices):
        v_s = v_out[s]

        # Dipole moment: <sigma_s^gamma> = V_s^{0, gamma}
        dipole_norm = float(np.sqrt(v_s[0, 1] ** 2 + v_s[0, 2] ** 2 + v_s[0, 3] ** 2))
        assert dipole_norm <= 1.0 + 1e-12

        # Fermion parity: Pf(V_s) = V_01*V_23 - V_02*V_13 + V_03*V_12
        pfaffian = float(v_s[0, 1] * v_s[2, 3] - v_s[0, 2] * v_s[1, 3] + v_s[0, 3] * v_s[1, 2])
        assert -1.0 - 1e-12 <= pfaffian <= 1.0 + 1e-12

        # Gauge Casimir expectation: <T_s^2> = (3/8) * (1 - Pf(V_s))
        gauge_casimir = 0.375 * (1.0 - pfaffian)
        assert 0.0 - 1e-12 <= gauge_casimir <= 0.75 + 1e-12
