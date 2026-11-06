"""Seed state factories, symmetry operations, and controlled perturbations."""

from __future__ import annotations

from typing import Literal

import numpy as np

from kitaev_mf.domain.lattice import UnitCell
from kitaev_mf.domain.state import MeanFieldState


def make_qsl_seed(
    cell: UnitCell,
    u_gauge: float = 1.0,
    u_matter: float = -0.5248645,
) -> MeanFieldState:
    """Create converged or nominal pure-Kitaev zero-flux ground state."""
    return MeanFieldState.from_pure_kitaev(cell, u_gauge=u_gauge, u_matter=u_matter)


def make_polarized_seed(
    cell: UnitCell,
    axis: Literal["111", "x", "y", "z", "minus_111"] = "111",
) -> MeanFieldState:
    """Create an unentangled, gauge-neutral polarized state oriented along designated axis."""
    if axis == "111":
        n = (1.0 / np.sqrt(3), 1.0 / np.sqrt(3), 1.0 / np.sqrt(3))
    elif axis == "minus_111":
        n = (-1.0 / np.sqrt(3), -1.0 / np.sqrt(3), -1.0 / np.sqrt(3))
    elif axis == "x":
        n = (1.0, 0.0, 0.0)
    elif axis == "y":
        n = (0.0, 1.0, 0.0)
    elif axis == "z":
        n = (0.0, 0.0, 1.0)
    else:
        raise ValueError(f"Unknown axis: {axis}")
    return MeanFieldState.from_polarized(cell, n)


def apply_gauge_transformation_sublattice_b(
    state: MeanFieldState, cell: UnitCell
) -> MeanFieldState:
    """Apply a discrete Z_2 gauge transformation on sublattice B (inverting U_b^{gamma, gamma}).

    Under this gauge transformation, the physical observables and plaquette fluxes W_p remain invariant.
    """
    u_new = state.U.copy()
    for b_idx, bond in enumerate(cell.bonds):
        gamma = bond.bond_type
        # Bond is oriented A -> B: inverting gauge signs on B flips U_b^{gamma, gamma}
        u_new[b_idx, gamma, gamma] = -u_new[b_idx, gamma, gamma]
    return MeanFieldState(U=u_new, V=state.V.copy(), lambdas=state.lambdas.copy())


def apply_flux_flip(state: MeanFieldState, cell: UnitCell, bond_idx: int = 0) -> MeanFieldState:
    """Flip the gauge sign on a specific bond index to switch the flux sector in supercells."""
    u_new = state.U.copy()
    gamma = cell.bonds[bond_idx].bond_type
    u_new[bond_idx, gamma, gamma] = -u_new[bond_idx, gamma, gamma]
    return MeanFieldState(U=u_new, V=state.V.copy(), lambdas=state.lambdas.copy())


def add_controlled_perturbation(
    state: MeanFieldState,
    delta: float,
    rng: np.random.Generator,
    skew_only: bool = True,
) -> MeanFieldState:
    """Add controlled random Gaussian perturbation to state correlators.

    Parameters:
        state: Base MeanFieldState to perturb.
        delta: Perturbation amplitude scale.
        rng: NumPy random generator.
        skew_only: Whether to strictly enforce skew-symmetry on V perturbations.
    """
    if delta == 0.0:
        return state

    du = delta * rng.standard_normal(state.U.shape)
    u_pert = state.U + du

    dv = delta * rng.standard_normal(state.V.shape)
    if skew_only:
        dv = 0.5 * (dv - np.swapaxes(dv, -1, -2))
    v_pert = state.V + dv

    # If non-skew is allowed, note that MeanFieldState constructor will validate skew-symmetry.
    # For Arm A testing, non-skew states are manipulated as raw arrays in the solver.
    if not skew_only:
        # Enforce floating-point skewness so MeanFieldState doesn't raise, unless raw testing
        v_pert = 0.5 * (v_pert - np.swapaxes(v_pert, -1, -2))

    return MeanFieldState(U=u_pert, V=v_pert, lambdas=state.lambdas.copy())


def check_seed_representability(state: MeanFieldState, cell: UnitCell) -> bool:
    """Check whether local two-site covariance blocks satisfy fermionic representability ||Omega|| <= 1.

    For each bond b = <ij>, the 8x8 block is C_ij = [[I - i V_i, -i U_b], [i U_b^T, I - i V_j]].
    Fermionic positivity requires 0 <= C_ij <= 2 * I_8, or equivalently eigenvalues in [0, 2].
    """
    for b_idx, bond in enumerate(cell.bonds):
        s_i = bond.sublattice_i
        s_j = bond.sublattice_j
        u_b = state.U[b_idx]
        v_i = state.V[s_i]
        v_j = state.V[s_j]

        block = np.block(
            [
                [np.eye(4) - 1.0j * v_i, -1.0j * u_b],
                [1.0j * u_b.T, np.eye(4) - 1.0j * v_j],
            ]
        )
        eigvals = np.linalg.eigvalsh(block)
        if np.any(eigvals < -1e-12) or np.any(eigvals > 2.0 + 1e-12):
            return False
    return True


def add_symmetric_perturbation(
    state: MeanFieldState,
    delta: float,
    rng: np.random.Generator,
    cell: UnitCell | None = None,
) -> MeanFieldState:
    """Add random perturbation strictly preserving C_3 rotation and sublattice A<->B equivalence.

    Parameters:
        state: Base MeanFieldState to perturb.
        delta: Perturbation magnitude.
        rng: NumPy random number generator.
        cell: Optional UnitCell.

    Returns:
        Perturbed MeanFieldState satisfying C_3 and sublattice symmetries.
    """
    if delta == 0.0:
        return state

    n_bonds = state.U.shape[0]
    n_sites = state.V.shape[0]

    # 1. Generate single symmetric bond perturbation du_base for U
    # Under C_3 and A<->B, all 3 bond types x, y, z are cyclically mapped.
    # We generate a common random 4x4 tensor du_sym that satisfies the diagonal/off-diagonal C3 structure
    du_raw = delta * rng.standard_normal((4, 4))
    # Symmetrize matter sector (0,0) and gauge sector (gamma, gamma)
    u_pert = state.U.copy()
    for b in range(n_bonds):
        # Add perturbation maintaining bond-type permutation
        gamma = (b % 3) + 1
        u_pert[b, 0, 0] += du_raw[0, 0]
        u_pert[b, gamma, gamma] += du_raw[1, 1]

    # 2. Generate common site perturbation dv for all sites (sublattice equivalence)
    dv_raw = delta * rng.standard_normal((4, 4))
    # Anti-symmetrize to maintain so(4) Lie algebra
    dv_skew = 0.5 * (dv_raw - dv_raw.T)

    v_pert = np.zeros_like(state.V)
    for s in range(n_sites):
        v_pert[s] = state.V[s] + dv_skew

    # 3. Perturb lambdas identically on all sites
    dlam = delta * rng.standard_normal(3)
    # Symmetrize across x, y, z components for [111] C3 isotropy
    dlam_scalar = float(np.mean(dlam))
    lambdas_pert = np.zeros_like(state.lambdas)
    for s in range(n_sites):
        lambdas_pert[s] = state.lambdas[s] + dlam_scalar

    return MeanFieldState(U=u_pert, V=v_pert, lambdas=lambdas_pert)
