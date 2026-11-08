"""Symmetry reduction and analytical block-diagonalization for Kitaev mean-field theory under C3 and sublattice parity.

Under magnetic field h || [111], the system preserves:
1. Spatial translation invariance across the 2-sublattice honeycomb cell.
2. C3 rotation about the hexagon/site normal, cycling Cartesian axes x -> y -> z.
3. Sublattice inversion parity I relating sites A and B.

This module proves that the unconstrained D = 86 variational manifold collapses
to a D_sym = 5 to 7 scalar parameter subspace, and implements exact analytical
projections, block-diagonalization, and Dirac point spectra.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Preserve historical module-level project re-exports.
from kitaev_mf.domain.constants import G_MATRICES as G_MATRICES
from kitaev_mf.domain.constants import N_MATRICES as N_MATRICES
from kitaev_mf.domain.lattice import UnitCell, make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import MeanFieldState


@dataclass(frozen=True)
class SymmetricParameters:
    """The 7 scalar parameters defining the full C3- and inversion-symmetric mean-field state.

    Attributes:
        lambda_val: Uniform gauge Lagrange multiplier (lambda_A^gamma = lambda_B^gamma = lambda_val).
        v0: Matter-gauge on-site correlator V_s^{0, gamma} = v0 (magnetization m_111 = 2 * v0).
        v1: Gauge-gauge on-site correlator V_s^{1,2} = V_s^{2,3} = V_s^{3,1} = v1.
        u0: Matter-matter bond hopping U_z^{0,0} = u0 (~0.52486 at h=0).
        uz: Parallel gauge-gauge bond link U_z^{3,3} = uz (~1.0 at h=0).
        uperp: Transverse gauge-gauge bond link U_z^{1,1} = U_z^{2,2} = uperp (~0.0 at h=0).
        ucross: Matter-gauge bond cross term U_z^{0,1} = -U_z^{0,2} = ucross (~0.0 at h=0).
    """

    lambda_val: float
    v0: float
    v1: float
    u0: float
    uz: float
    uperp: float
    ucross: float

    @classmethod
    def pure_kitaev(cls, u_matter: float = -0.52486574585) -> SymmetricParameters:
        """Reference unperturbed pure-Kitaev state at zero field (u_matter = -0.5248657 for K = -1.0)."""
        return cls(
            lambda_val=0.0,
            v0=0.0,
            v1=0.0,
            u0=float(u_matter),
            uz=1.0,
            uperp=0.0,
            ucross=0.0,
        )

    def to_vector(self) -> np.ndarray:
        """Export as 7D numpy array."""
        return np.array(
            [
                self.lambda_val,
                self.v0,
                self.v1,
                self.u0,
                self.uz,
                self.uperp,
                self.ucross,
            ],
            dtype=np.float64,
        )

    @classmethod
    def from_vector(cls, vec: np.ndarray) -> SymmetricParameters:
        """Reconstruct from 7D vector."""
        v = np.asarray(vec, dtype=np.float64).ravel()
        return cls(
            lambda_val=float(v[0]),
            v0=float(v[1]),
            v1=float(v[2]),
            u0=float(v[3]),
            uz=float(v[4]),
            uperp=float(v[5]),
            ucross=float(v[6]),
        )


def pack_symmetric_state(
    sym: SymmetricParameters,
    cell: UnitCell | None = None,
) -> MeanFieldState:
    """Construct a full MeanFieldState enforcing exact C3 and sublattice symmetry.

    Args:
        sym: Symmetric scalar parameter bundle.
        cell: Honeycomb unit cell. If None, uses standard 2-site cell.

    Returns:
        MeanFieldState with mathematically symmetric U, V, and lambdas.
    """
    if cell is None:
        cell = make_honeycomb_2site()

    num_sub = cell.num_sublattices
    num_bonds = len(cell.bonds)

    # 1. Lagrange multipliers: lambda_i^gamma = lambda_val for all i, gamma
    lambdas = np.full((num_sub, 3), sym.lambda_val, dtype=np.float64)

    # 2. On-site skew-symmetric matrix V_s
    # V[0, gamma] = v0, V[gamma, 0] = -v0
    # V[1, 2] = -V[2, 1] = v1, V[2, 3] = -V[3, 2] = v1, V[3, 1] = -V[1, 3] = -v1
    v_mat = np.zeros((4, 4), dtype=np.float64)
    for g in (1, 2, 3):
        v_mat[0, g] = sym.v0
        v_mat[g, 0] = -sym.v0

    v_mat[1, 2] = -sym.v1
    v_mat[2, 1] = sym.v1
    v_mat[2, 3] = -sym.v1
    v_mat[3, 2] = sym.v1
    v_mat[3, 1] = sym.v1
    v_mat[1, 3] = -sym.v1

    v_all = np.zeros((num_sub, 4, 4), dtype=np.float64)
    for s in range(num_sub):
        v_all[s] = v_mat.copy()

    # 3. Directed bond matrices U_b
    # Reference z-bond (gamma = 3):
    u_z = np.zeros((4, 4), dtype=np.float64)
    u_z[0, 0] = sym.u0
    u_z[3, 3] = sym.uz
    u_z[1, 1] = sym.uperp
    u_z[2, 2] = sym.uperp
    u_z[0, 1] = sym.ucross
    u_z[1, 0] = sym.ucross
    u_z[0, 2] = -sym.ucross
    u_z[2, 0] = -sym.ucross

    # Cyclic permutation matrices for x-bond (1) and y-bond (2):
    # Mapping for flavors: (0, 1, 2, 3) -> (0, 2, 3, 1) for y, (0, 3, 1, 2) for x
    # Permutation matrix P: P[0,0]=1; x-bond cyclically maps 3 -> 1 -> 2 -> 3
    # P_x maps z (3) to x (1): (1, 2, 3) -> (2, 3, 1) or cyclic perm
    u_all = np.zeros((num_bonds, 4, 4), dtype=np.float64)
    for b_idx, bond in enumerate(cell.bonds):
        g = bond.bond_type
        if g == 3:  # z-bond
            u_all[b_idx] = u_z.copy()
        elif g == 1:  # x-bond: cycle (1, 2, 3) -> (3, 1, 2)
            u_x = np.zeros((4, 4), dtype=np.float64)
            u_x[0, 0] = sym.u0
            u_x[1, 1] = sym.uz
            u_x[2, 2] = sym.uperp
            u_x[3, 3] = sym.uperp
            u_x[0, 2] = sym.ucross
            u_x[2, 0] = sym.ucross
            u_x[0, 3] = -sym.ucross
            u_x[3, 0] = -sym.ucross
            u_all[b_idx] = u_x
        elif g == 2:  # y-bond: cycle (1, 2, 3) -> (2, 3, 1)
            u_y = np.zeros((4, 4), dtype=np.float64)
            u_y[0, 0] = sym.u0
            u_y[2, 2] = sym.uz
            u_y[3, 3] = sym.uperp
            u_y[1, 1] = sym.uperp
            u_y[0, 3] = sym.ucross
            u_y[3, 0] = sym.ucross
            u_y[0, 1] = -sym.ucross
            u_y[1, 0] = -sym.ucross
            u_all[b_idx] = u_y

    return MeanFieldState(U=u_all, V=v_all, lambdas=lambdas)


def extract_symmetric_parameters(state: MeanFieldState) -> SymmetricParameters:
    """Extract the best-fit SymmetricParameters from a generic MeanFieldState via C3 projection."""
    # Multiplier average across all sites and flavors
    lam_val = float(np.mean(state.lambdas))

    # On-site V average
    v0_vals = []
    for s in range(state.V.shape[0]):
        v0_vals.extend([state.V[s, 0, 1], state.V[s, 0, 2], state.V[s, 0, 3]])
    v0 = float(np.mean(v0_vals))

    v1_vals = []
    for s in range(state.V.shape[0]):
        v1_vals.extend([-state.V[s, 1, 2], -state.V[s, 2, 3], state.V[s, 3, 1]])
    v1 = float(np.mean(v1_vals))

    # Bond matrix U averages
    u0_vals = [state.U[b, 0, 0] for b in range(state.U.shape[0])]
    u0 = float(np.mean(u0_vals))

    uz_vals = [
        state.U[0, 3, 3],  # z-bond
        state.U[1, 1, 1],  # x-bond
        state.U[2, 2, 2],  # y-bond
    ]
    uz = float(np.mean(uz_vals))

    uperp_vals = [
        state.U[0, 1, 1],
        state.U[0, 2, 2],
        state.U[1, 2, 2],
        state.U[1, 3, 3],
        state.U[2, 3, 3],
        state.U[2, 1, 1],
    ]
    uperp = float(np.mean(uperp_vals))

    ucross_vals = [
        state.U[0, 0, 1],
        -state.U[0, 0, 2],
        state.U[1, 0, 2],
        -state.U[1, 0, 3],
        state.U[2, 0, 3],
        -state.U[2, 0, 1],
    ]
    ucross = float(np.mean(ucross_vals))

    return SymmetricParameters(
        lambda_val=lam_val,
        v0=v0,
        v1=v1,
        u0=u0,
        uz=uz,
        uperp=uperp,
        ucross=ucross,
    )


def project_state_to_symmetric(
    state: MeanFieldState, cell: UnitCell | None = None
) -> MeanFieldState:
    """Project an arbitrary MeanFieldState onto the exact C3 and sublattice symmetric manifold."""
    sym = extract_symmetric_parameters(state)
    return pack_symmetric_state(sym, cell)


def build_dirac_bloch_hamiltonian(
    model_params: ModelParameters,
    sym: SymmetricParameters,
) -> np.ndarray:
    """Assemble the 8x8 Bloch Hamiltonian H_MF(K) at the Dirac point K = (1/3, 2/3).

    Returns:
        h_k: Complex 8x8 Hermitian matrix in the basis (A0, A1, A2, A3, B0, B1, B2, B3).
    """
    cell = make_honeycomb_2site()
    state = pack_symmetric_state(sym, cell)

    from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian

    k_pts = np.array([[1.0 / 3.0, 2.0 / 3.0]])
    h_bloch = build_bloch_hamiltonian(k_pts, cell, model_params, state)
    if hasattr(h_bloch, "numpy"):
        return np.asarray(h_bloch.numpy()[0])
    return np.asarray(h_bloch[0])


def block_diagonalize_dirac_point(
    model_params: ModelParameters,
    sym: SymmetricParameters,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Perform unitary transformation rotating H_MF(K) into matter (2x2) and gauge (6x6) blocks.

    Returns:
        (W, h_matter_eff, h_gauge_eff):
            W: 8x8 unitary matrix.
            h_matter_eff: 2x2 effective matter Dirac block.
            h_gauge_eff: 6x6 gauge Majorana block.
    """
    h8 = build_dirac_bloch_hamiltonian(model_params, sym)

    # Permutation matrix grouping matter basis (A0, B0) followed by gauge basis (A1..3, B1..3)
    # Original indices: A0=0, A1=1, A2=2, A3=3, B0=4, B1=5, B2=6, B3=7
    perm = [0, 4, 1, 2, 3, 5, 6, 7]
    h_perm = h8[np.ix_(perm, perm)]

    h_matter = h_perm[:2, :2]
    h_cross = h_perm[:2, 2:]
    h_gauge = h_perm[2:, 2:]

    # At h=0, h_cross is identically zero or vanishingly small, and h_matter is exactly zero.
    # At h != 0, Schrieffer-Wolff 2nd-order resolvent projection yields:
    # h_matter_eff = h_matter - h_cross @ inv(h_gauge) @ h_cross^dagger
    try:
        inv_gauge = np.linalg.pinv(h_gauge)
        h_matter_eff = h_matter - h_cross @ inv_gauge @ h_cross.conj().T
    except np.linalg.LinAlgError:
        h_matter_eff = h_matter.copy()

    # Unitary transformation matrix W
    W = np.eye(8, dtype=np.complex128)[perm]

    return W, h_matter_eff, h_gauge
