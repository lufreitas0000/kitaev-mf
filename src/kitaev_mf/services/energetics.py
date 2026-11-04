"""Pure physics service for variational Gaussian energetics and spectral consistency."""

import numpy as np
import torch

from kitaev_mf.domain.constants import G_MATRICES, N_MATRICES
from kitaev_mf.domain.lattice import UnitCell
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.ports import LinalgBackendPort
from kitaev_mf.domain.state import EnergyBreakdown, MeanFieldState


def compute_gaussian_energy(
    state: MeanFieldState,
    cell: UnitCell,
    params: ModelParameters,
) -> EnergyBreakdown:
    """Evaluate variational Hartree-Fock energy E_MF^G[U, V, lambda].

    Parameters:
        state: Variational mean-field state (U, V, lambdas).
        cell: Lattice unit cell geometry Lambda = (V, E).
        params: Model exchange couplings and magnetic field.

    Returns:
        EnergyBreakdown containing extensive components and per-bond/cell/site densities.
    """
    num_sub = cell.num_sublattices
    num_bonds = len(cell.bonds)

    moments = np.zeros((num_sub, 3), dtype=np.float64)
    gauge_res = np.zeros((num_sub, 3), dtype=np.float64)

    for s in range(num_sub):
        v_s = state.V[s]
        for a in (1, 2, 3):
            moments[s, a - 1] = float(np.trace(v_s.T @ N_MATRICES[a]))
            gauge_res[s, a - 1] = float(np.trace(G_MATRICES[a] @ v_s))

    e_exchange = 0.0
    for b_idx, bond in enumerate(cell.bonds):
        s_i = bond.sublattice_i
        s_j = bond.sublattice_j
        j_mat = params.exchange_matrix(bond.bond_type)
        u_b = state.U[b_idx]

        for a in (1, 2, 3):
            for b in (1, 2, 3):
                j_val = j_mat[a - 1, b - 1]
                if j_val == 0.0:
                    continue
                contract = float(np.trace(N_MATRICES[a].T @ u_b @ N_MATRICES[b] @ u_b.T))
                e_exchange += j_val * (
                    0.25 * moments[s_i, a - 1] * moments[s_j, b - 1] - 0.5 * contract
                )

    e_field = 0.0
    h_vec = params.h_vec
    for s in range(num_sub):
        for a in (1, 2, 3):
            if h_vec[a - 1] != 0.0:
                e_field -= 0.5 * h_vec[a - 1] * moments[s, a - 1]

    e_multiplier = 0.0
    for s in range(num_sub):
        for gamma in (1, 2, 3):
            lam_val = state.lambdas[s, gamma - 1]
            if lam_val != 0.0:
                e_multiplier -= 0.25 * lam_val * gauge_res[s, gamma - 1]

    e_total = e_exchange + e_field + e_multiplier

    return EnergyBreakdown(
        E_exchange=float(e_exchange),
        E_field=float(e_field),
        E_multiplier=float(e_multiplier),
        E_total=float(e_total),
        num_bonds=num_bonds,
        num_cells=1,
        num_sites=num_sub,
    )


def compute_spectral_ground_state_energy(
    h_k_batch: torch.Tensor,
    cell: UnitCell,
    state: MeanFieldState,
    params: ModelParameters,
    backend: LinalgBackendPort | None = None,
) -> float:
    """Compute E_quad,0 - C from occupied negative Bloch eigenvalues.

    Parameters:
        h_k_batch: Complex Bloch Hamiltonian tensor of shape (N_k, dim, dim).
        cell: Lattice unit cell geometry Lambda = (V, E).
        state: Variational mean-field state (U, V, lambdas).
        params: Model exchange couplings and magnetic field.
        backend: Optional hardware adapter.

    Returns:
        Total energy per unit cell, float.
    """
    if backend is not None:
        eigenvalues, _ = backend.batched_eigh(h_k_batch)
        e_np = backend.to_numpy(eigenvalues)
    else:
        with torch.no_grad():
            eigenvalues = torch.linalg.eigvalsh(h_k_batch)
        e_np = eigenvalues.detach().cpu().numpy()

    # E_{quad,0} / N_cell = -2.0 * mean_k( sum_{E_n > 0} E_n(k) )
    pos_e = np.where(e_np > 1e-12, e_np, 0.0)
    e_quad_0 = -2.0 * float(np.mean(np.sum(pos_e, axis=-1)))

    # Subtract C = E_exchange (scalar double-counting subtraction)
    breakdown = compute_gaussian_energy(state, cell, params)
    c_scalar = breakdown.E_exchange

    e_spectral = e_quad_0 - c_scalar
    return float(e_spectral)


def check_energy_consistency(
    state: MeanFieldState,
    cell: UnitCell,
    params: ModelParameters,
    h_k_batch: torch.Tensor,
    backend: LinalgBackendPort | None = None,
    tolerance: float = 1e-10,
) -> float:
    """Assert |(E_quad,0 - C) - E_total^G| < tolerance, returning the residual.

    Parameters:
        state: Variational mean-field state.
        cell: Lattice unit cell geometry.
        params: Model parameters.
        h_k_batch: Bloch Hamiltonian batch over the BZ.
        backend: Optional hardware adapter.
        tolerance: Maximum allowed absolute difference between spectral and Wick energies.

    Returns:
        float absolute energy difference.
    """
    e_spectral = compute_spectral_ground_state_energy(h_k_batch, cell, state, params, backend)
    breakdown = compute_gaussian_energy(state, cell, params)
    residual = abs(e_spectral - breakdown.per_cell)

    if residual > tolerance:
        raise AssertionError(
            f"Energy consistency check failed: |E_spectral - E_total^G| = {residual} "
            f"> tolerance={tolerance}. (E_spectral={e_spectral}, E_total={breakdown.per_cell})"
        )
    return float(residual)
