"""Pure physics service assembling momentum-space Bloch-Majorana Hamiltonian matrices."""

import numpy as np
import torch

from kitaev_mf.domain.constants import G_MATRICES, N_MATRICES
from kitaev_mf.domain.lattice import (
    UnitCell,
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.ports import LinalgBackendPort
from kitaev_mf.domain.state import MeanFieldState


def build_on_site_matrices(
    cell: UnitCell,
    params: ModelParameters,
    state: MeanFieldState,
    device: torch.device,
) -> torch.Tensor:
    """Construct real skew-symmetric on-site matrices B_i of shape (num_sublattices, 4, 4).

    B_i = sum_{j in V_i} sum_{alpha,beta} J_{ij}^{alpha,beta} N^alpha tr(V_j^T N^beta)
          + sum_gamma (lambda_i^gamma G^gamma - 2 h^gamma N^gamma)
    """
    num_sub = cell.num_sublattices
    b_tensors = torch.zeros((num_sub, 4, 4), dtype=torch.float64, device=device)

    # Pre-convert algebraic matrices to target device
    n_torch = {
        a: torch.from_numpy(N_MATRICES[a]).to(device=device, dtype=torch.float64) for a in (1, 2, 3)
    }
    g_torch = {
        a: torch.from_numpy(G_MATRICES[a]).to(device=device, dtype=torch.float64) for a in (1, 2, 3)
    }

    v_torch = torch.from_numpy(state.V).to(device=device, dtype=torch.float64)
    lam_torch = torch.from_numpy(state.lambdas).to(device=device, dtype=torch.float64)
    h_vec = torch.from_numpy(params.h_vec).to(device=device, dtype=torch.float64)

    # Pre-calculate local spin expectations tr(V_j^T N^beta)
    tr_v_n = torch.zeros((num_sub, 3), dtype=torch.float64, device=device)
    for site in range(num_sub):
        for beta in (1, 2, 3):
            tr_v_n[site, beta - 1] = torch.trace(v_torch[site].T @ n_torch[beta])

    # Accumulate Hartree exchange fields from nearest neighbors
    for bond in cell.bonds:
        s_i = bond.sublattice_i
        s_j = bond.sublattice_j
        j_mat = torch.from_numpy(params.exchange_matrix(bond.bond_type)).to(
            device=device, dtype=torch.float64
        )

        # Contribution to site i from neighbor j (exchange matrix J_ij)
        for alpha in (1, 2, 3):
            for beta in (1, 2, 3):
                coeff_i = j_mat[alpha - 1, beta - 1] * tr_v_n[s_j, beta - 1]
                if coeff_i != 0.0:
                    b_tensors[s_i] += coeff_i * n_torch[alpha]

        # Contribution to site j from neighbor i (transposed exchange matrix J_ji = J_ij^T)
        for alpha in (1, 2, 3):
            for beta in (1, 2, 3):
                coeff_j = j_mat[beta - 1, alpha - 1] * tr_v_n[s_i, beta - 1]
                if coeff_j != 0.0:
                    b_tensors[s_j] += coeff_j * n_torch[alpha]

    # Add gauge Lagrange multipliers and external Zeeman fields
    for site in range(num_sub):
        for gamma in (1, 2, 3):
            lam_val = lam_torch[site, gamma - 1]
            if lam_val != 0.0:
                b_tensors[site] += lam_val * g_torch[gamma]

            h_val = h_vec[gamma - 1]
            if h_val != 0.0:
                b_tensors[site] -= 2.0 * h_val * n_torch[gamma]

    return b_tensors


def build_bond_hopping_matrices(
    cell: UnitCell,
    params: ModelParameters,
    state: MeanFieldState,
    device: torch.device,
) -> torch.Tensor:
    """Construct real 4x4 hopping matrices A_b for each directed bond b in E.

    A_b = 2 * sum_{alpha,beta} J_b^{alpha,beta} N^alpha U_b N^beta
    """
    num_bonds = len(cell.bonds)
    a_tensors = torch.zeros((num_bonds, 4, 4), dtype=torch.float64, device=device)

    n_torch = {
        a: torch.from_numpy(N_MATRICES[a]).to(device=device, dtype=torch.float64) for a in (1, 2, 3)
    }
    u_torch = torch.from_numpy(state.U).to(device=device, dtype=torch.float64)

    for b_idx, bond in enumerate(cell.bonds):
        j_mat = torch.from_numpy(params.exchange_matrix(bond.bond_type)).to(
            device=device, dtype=torch.float64
        )
        u_b = u_torch[b_idx]

        for alpha in (1, 2, 3):
            for beta in (1, 2, 3):
                j_val = j_mat[alpha - 1, beta - 1]
                if j_val != 0.0:
                    a_tensors[b_idx] += 2.0 * j_val * (n_torch[alpha] @ u_b @ n_torch[beta])

    return a_tensors


def build_bloch_hamiltonian(
    k_pts: np.ndarray | torch.Tensor,
    cell: UnitCell,
    params: ModelParameters,
    state: MeanFieldState,
    backend: LinalgBackendPort | None = None,
) -> torch.Tensor:
    """Assemble the momentum-space Bloch-Majorana Hamiltonian H(k) = (i / 4) * M(k).

    Parameters:
        k_pts: Fractional reciprocal coordinates (k1, k2) in [0, 1)^2.
               Shape: (2,) for a single point or (N_k, 2) for a batch.
        cell: Lattice unit cell geometry Lambda = (V, E).
        params: Model exchange couplings and magnetic field.
        state: Mean-field correlation matrices (U, V, lambdas).
        backend: Optional hardware adapter. Defaults to CPU if None.

    Returns:
        torch.Tensor of shape (N_k, 4*N_sub, 4*N_sub) or (4*N_sub, 4*N_sub), complex128.
    """
    device = backend.device if backend is not None else torch.device("cpu")

    single_k = False
    if isinstance(k_pts, np.ndarray):
        k_tensor = torch.from_numpy(k_pts).to(device=device, dtype=torch.float64)
    else:
        k_tensor = k_pts.to(device=device, dtype=torch.float64)

    if k_tensor.ndim == 1 and k_tensor.shape[0] == 2:
        k_tensor = k_tensor.unsqueeze(0)
        single_k = True
    elif k_tensor.ndim != 2 or k_tensor.shape[1] != 2:
        raise ValueError(f"k_pts must have shape (2,) or (N_k, 2), got {k_tensor.shape}")

    num_k = k_tensor.shape[0]
    num_sub = cell.num_sublattices
    dim = 4 * num_sub

    # Construct single-particle real on-site and hopping blocks
    b_mat = build_on_site_matrices(cell, params, state, device)
    a_mat = build_bond_hopping_matrices(cell, params, state, device)

    # Allocate complex Hermitian Hamiltonian tensor
    h_bloch = torch.zeros((num_k, dim, dim), dtype=torch.complex128, device=device)

    # 1. On-site skew-symmetric terms: H_ii = (i / 4) * B_i
    for s in range(num_sub):
        b_complex = b_mat[s].to(dtype=torch.complex128)
        h_bloch[:, 4 * s : 4 * s + 4, 4 * s : 4 * s + 4] = 0.25j * b_complex

    # 2. Directed bond terms with explicit 2*pi cell-offset phases
    k1 = k_tensor[:, 0]
    k2 = k_tensor[:, 1]

    for b_idx, bond in enumerate(cell.bonds):
        s_i = bond.sublattice_i
        s_j = bond.sublattice_j
        n1, n2 = bond.cell_offset

        # Phase factor: phi_b(k) = 2*pi*(k1*n1 + k2*n2)
        phase = 2.0 * np.pi * (k1 * n1 + k2 * n2)
        exp_phase = torch.exp(1j * phase)[:, None, None]
        exp_phase_conj = torch.conj(exp_phase)

        a_complex = a_mat[b_idx].to(dtype=torch.complex128)

        # Block (s_i, s_j): +(i / 4) * A_b * exp(i * phi_b)
        h_bloch[:, 4 * s_i : 4 * s_i + 4, 4 * s_j : 4 * s_j + 4] += 0.25j * a_complex * exp_phase

        # Block (s_j, s_i): -(i / 4) * A_b^T * exp(-i * phi_b)
        h_bloch[:, 4 * s_j : 4 * s_j + 4, 4 * s_i : 4 * s_i + 4] -= (
            0.25j * a_complex.T * exp_phase_conj
        )

    return h_bloch.squeeze(0) if single_k else h_bloch


def build_bloch_hamiltonian_2site(
    k_pts: np.ndarray | torch.Tensor,
    params: ModelParameters,
    state: MeanFieldState,
    backend: LinalgBackendPort | None = None,
) -> torch.Tensor:
    """Thin facade assembling H(k) for the standard primitive 2-sublattice honeycomb cell."""
    cell = make_honeycomb_2site()
    return build_bloch_hamiltonian(k_pts, cell, params, state, backend)


def build_bloch_hamiltonian_4site(
    k_pts: np.ndarray | torch.Tensor,
    params: ModelParameters,
    state: MeanFieldState,
    backend: LinalgBackendPort | None = None,
) -> torch.Tensor:
    """Thin facade assembling H(k) for the rectangular 4-sublattice honeycomb cell."""
    cell = make_honeycomb_4site_rectangular()
    return build_bloch_hamiltonian(k_pts, cell, params, state, backend)
