"""Pure physics service for Gaussian spin observables, correlations, and Wilson loops."""

import numpy as np

from kitaev_mf.domain.constants import N_MATRICES
from kitaev_mf.domain.lattice import UnitCell
from kitaev_mf.domain.state import MeanFieldState


def compute_sublattice_magnetizations(
    state: MeanFieldState,
    cell: UnitCell,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute local magnetization vectors on each sublattice.

    Parameters:
        state: Active variational mean-field state (U, V, lambdas).
        cell: Lattice unit cell geometry Lambda = (V, E).

    Returns:
        Tuple of (sigma_mag, spin_mag):
            sigma_mag: Shape (num_sublattices, 3), float64, <sigma_s^a>_G = V_s^{0a}.
            spin_mag:  Shape (num_sublattices, 3), float64, <S_s^a>_G = 0.5 * V_s^{0a}.
    """
    num_sub = cell.num_sublattices
    sigma_mag = np.zeros((num_sub, 3), dtype=np.float64)
    spin_mag = np.zeros((num_sub, 3), dtype=np.float64)

    for s in range(num_sub):
        v_s = state.V[s]
        for a in (1, 2, 3):
            m_val = float(np.trace(v_s.T @ N_MATRICES[a]))
            sigma_mag[s, a - 1] = 0.5 * m_val
            spin_mag[s, a - 1] = 0.25 * m_val

    return sigma_mag, spin_mag


def compute_bond_spin_correlators(
    state: MeanFieldState,
    cell: UnitCell,
) -> np.ndarray:
    """Compute 3x3 spin-spin correlation tensor for every directed bond b in E.

    Parameters:
        state: Variational mean-field state (U, V, lambdas).
        cell: Lattice unit cell geometry Lambda = (V, E).

    Returns:
        np.ndarray of shape (num_bonds, 3, 3), float64, where
        spin_corr[b, a-1, b'-1] = <sigma_i^a sigma_j^{b'}>_G.
    """
    num_bonds = len(cell.bonds)
    corrs = np.zeros((num_bonds, 3, 3), dtype=np.float64)

    for b_idx, bond in enumerate(cell.bonds):
        s_i = bond.sublattice_i
        s_j = bond.sublattice_j
        u_b = state.U[b_idx]
        v_i = state.V[s_i]
        v_j = state.V[s_j]

        m_i = np.array(
            [float(np.trace(v_i.T @ N_MATRICES[a])) for a in (1, 2, 3)],
            dtype=np.float64,
        )
        m_j = np.array(
            [float(np.trace(v_j.T @ N_MATRICES[b])) for b in (1, 2, 3)],
            dtype=np.float64,
        )

        for a in (1, 2, 3):
            for b in (1, 2, 3):
                contract = float(np.trace(N_MATRICES[a].T @ u_b @ N_MATRICES[b] @ u_b.T))
                corrs[b_idx, a - 1, b - 1] = 0.25 * m_i[a - 1] * m_j[b - 1] - 0.5 * contract

    return corrs


def compute_plaquette_fluxes(
    state: MeanFieldState,
    cell: UnitCell,
) -> np.ndarray:
    """Compute local Z_2 Wilson loop expectation values W_p = prod_{b in d p} u_b.

    Parameters:
        state: Variational mean-field state (U, V, lambdas).
        cell: Lattice unit cell geometry Lambda = (V, E).

    Returns:
        np.ndarray of shape (len(cell.plaquettes),), float64 in [-1.0, +1.0].
        In the canonical zero-flux ground state, all entries equal +1.0.
    """
    num_plaq = len(cell.plaquettes)
    fluxes = np.ones(num_plaq, dtype=np.float64)

    for p, plaq in enumerate(cell.plaquettes):
        if plaq.oriented_bonds:
            b_indices = [b_idx for b_idx, _ in plaq.oriented_bonds]
        else:
            b_indices = list(plaq.bond_indices)

        prod_val = 1.0
        for b_idx in b_indices:
            gamma = cell.bonds[b_idx].bond_type
            prod_val *= float(state.U[b_idx, gamma, gamma])
        fluxes[p] = prod_val

    return fluxes
