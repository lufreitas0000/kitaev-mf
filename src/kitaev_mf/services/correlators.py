"""Pure physics service for Majorana spectral projection, BZ integration, and gauge residuals."""

import numpy as np
import torch

from kitaev_mf.domain.constants import G_MATRICES
from kitaev_mf.domain.lattice import UnitCell
from kitaev_mf.domain.ports import LinalgBackendPort


def compute_spectral_projector(
    h_k_batch: torch.Tensor,
    backend: LinalgBackendPort | None = None,
    zero_tol: float = 1e-12,
) -> torch.Tensor:
    """Compute reciprocal two-point correlation matrix P(k) = 2 * Pi^-(k).

    Parameters:
        h_k_batch: Complex Bloch Hamiltonian tensor of shape (N_k, dim, dim).
        backend: Optional hardware adapter.
        zero_tol: Threshold below which eigenvalues are treated as exact zero-modes.

    Returns:
        torch.Tensor of shape (N_k, dim, dim), complex128.
    """
    if backend is not None:
        eigenvalues, eigenvectors = backend.batched_eigh(h_k_batch)
    else:
        with torch.no_grad():
            eigenvalues, eigenvectors = torch.linalg.eigh(h_k_batch)

    # Weights: 1.0 for occupied (E < -zero_tol), 0.0 for empty (E > zero_tol), 0.5 for zero-modes
    weights = torch.where(
        eigenvalues < -zero_tol,
        torch.tensor(1.0, dtype=torch.float64, device=h_k_batch.device),
        torch.where(
            eigenvalues > zero_tol,
            torch.tensor(0.0, dtype=torch.float64, device=h_k_batch.device),
            torch.tensor(0.5, dtype=torch.float64, device=h_k_batch.device),
        ),
    )

    # Vectorized projector: P_{ab}(k) = 2 * sum_n weights_n(k) * v_{an}(k) * conj(v_{bn}(k))
    weights_c = weights.to(dtype=torch.complex128).unsqueeze(-2)  # (N_k, 1, dim)
    scaled_v = eigenvectors * weights_c  # (N_k, dim, dim)
    p_k = 2.0 * torch.matmul(scaled_v, eigenvectors.conj().transpose(-1, -2))
    return p_k


def compute_correlators_from_bloch(
    h_k_batch: torch.Tensor,
    cell: UnitCell,
    k_pts: np.ndarray | torch.Tensor,
    backend: LinalgBackendPort | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate spectral projector P(k) over the BZ returning real-space (U, V).

    Parameters:
        h_k_batch: Complex Bloch Hamiltonian tensor of shape (N_k, dim, dim).
        cell: Lattice unit cell geometry Lambda = (V, E).
        k_pts: Fractional coordinates of the BZ mesh, shape (N_k, 2).
        backend: Optional hardware adapter.

    Returns:
        Tuple of (U, V):
            U: Bond correlation matrices of shape (num_bonds, 4, 4), float64.
            V: On-site skew-symmetric correlation matrices of shape (num_sublattices, 4, 4), float64.
    """
    p_k = compute_spectral_projector(h_k_batch, backend)
    p_np = p_k.detach().cpu().numpy()

    if isinstance(k_pts, torch.Tensor):
        k_np = k_pts.detach().cpu().numpy()
    else:
        k_np = np.asarray(k_pts, dtype=np.float64)

    num_k = k_np.shape[0]
    num_bonds = len(cell.bonds)
    num_sub = cell.num_sublattices

    U_out = np.zeros((num_bonds, 4, 4), dtype=np.float64)
    V_out = np.zeros((num_sub, 4, 4), dtype=np.float64)

    # 1. Real-space bond correlators: U_b = (1 / N_k) * sum_k [ Im(P) cos(phi_b) - Re(P) sin(phi_b) ]
    k1 = k_np[:, 0]
    k2 = k_np[:, 1]

    for b_idx, bond in enumerate(cell.bonds):
        s_i = bond.sublattice_i
        s_j = bond.sublattice_j
        n1, n2 = bond.cell_offset

        phi_b = 2.0 * np.pi * (k1 * n1 + k2 * n2)  # (N_k,)
        cos_phi = np.cos(phi_b)[:, None, None]
        sin_phi = np.sin(phi_b)[:, None, None]

        p_block = p_np[:, 4 * s_i : 4 * s_i + 4, 4 * s_j : 4 * s_j + 4]
        u_b = np.sum(p_block.imag * cos_phi - p_block.real * sin_phi, axis=0) / float(num_k)
        U_out[b_idx] = u_b

    # 2. Real-space on-site correlators: V_s = (1 / N_k) * sum_k Im(P_{ss}(k))
    for s in range(num_sub):
        p_diag = p_np[:, 4 * s : 4 * s + 4, 4 * s : 4 * s + 4]
        v_s = np.sum(p_diag.imag, axis=0) / float(num_k)
        # Enforce exact floating-point skew-symmetry
        V_out[s] = 0.5 * (v_s - v_s.T)

    return U_out, V_out


def compute_gauge_residuals(V: np.ndarray) -> np.ndarray:
    """Evaluate local gauge constraint residuals g_s^gamma = tr(G^gamma V_s).

    Parameters:
        V: Real skew-symmetric on-site correlators of shape (num_sublattices, 4, 4).

    Returns:
        np.ndarray of shape (num_sublattices, 3), float64.
    """
    num_sub = V.shape[0]
    residuals = np.zeros((num_sub, 3), dtype=np.float64)

    for s in range(num_sub):
        v_s = V[s]
        for gamma in (1, 2, 3):
            residuals[s, gamma - 1] = float(np.trace(G_MATRICES[gamma] @ v_s))

    return residuals
