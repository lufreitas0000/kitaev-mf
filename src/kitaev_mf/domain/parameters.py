"""Immutable parameter dataclasses and bond exchange tensors for the extended Kitaev model."""

from dataclasses import dataclass

import numpy as np

BondCoupling = float | tuple[float, float, float]


def _resolve_bond_coupling(val: BondCoupling, bond_type: int) -> float:
    """Extract coupling constant for bond direction gamma in {1, 2, 3}."""
    if bond_type not in (1, 2, 3):
        raise ValueError(f"bond_type gamma must be in {{1, 2, 3}}, got {bond_type}")
    if isinstance(val, tuple):
        return float(val[bond_type - 1])
    return float(val)


@dataclass(frozen=True)
class ModelParameters:
    """Coupling constants for the extended J-K-Gamma-Gamma'-D-h Kitaev model.

    Each exchange parameter (J, K, Gamma, Gamma_prime, D) may be specified either as an
    isotropic float or as a directional 3-tuple (val_x, val_y, val_z) for gamma in {1, 2, 3}.
    """

    J: BondCoupling = 0.0
    K: BondCoupling = -1.0
    Gamma: BondCoupling = 0.0
    Gamma_prime: BondCoupling = 0.0
    D: BondCoupling = 0.0
    h: tuple[float, float, float] = (0.0, 0.0, 0.0)
    V0: float = 0.0

    @property
    def h_vec(self) -> np.ndarray:
        """External Zeeman field vector h = (hx, hy, hz)^T."""
        return np.array(self.h, dtype=np.float64)

    def exchange_matrix(self, bond_type: int) -> np.ndarray:
        """Construct the 3x3 exchange matrix J_{<ij>_gamma} for a directed bond <ij>_gamma (i in A, j in B).

        For gamma in {1, 2, 3} with (alpha, beta, gamma) a cyclic permutation of (1, 2, 3):
          J^{mu,nu} = J * delta_{mu,nu} + K * delta_{mu,gamma}*delta_{nu,gamma}
                      + Gamma * (delta_{mu,alpha}*delta_{nu,beta} + delta_{mu,beta}*delta_{nu,alpha})
                      + Gamma' * (delta_{mu,alpha}*delta_{nu,gamma} + delta_{mu,gamma}*delta_{nu,alpha}
                                  + delta_{mu,beta}*delta_{nu,gamma} + delta_{mu,gamma}*delta_{nu,beta})
                      + D * (delta_{mu,alpha}*delta_{nu,beta} - delta_{mu,beta}*delta_{nu,alpha})
        """
        j_val = _resolve_bond_coupling(self.J, bond_type)
        k_val = _resolve_bond_coupling(self.K, bond_type)
        g_val = _resolve_bond_coupling(self.Gamma, bond_type)
        gp_val = _resolve_bond_coupling(self.Gamma_prime, bond_type)
        d_val = _resolve_bond_coupling(self.D, bond_type)

        # 0-based cyclic indices (alpha, beta, gamma) in {0, 1, 2}
        gamma_idx = bond_type - 1
        alpha_idx = (gamma_idx + 1) % 3
        beta_idx = (gamma_idx + 2) % 3

        j_mat = np.eye(3, dtype=np.float64) * j_val
        j_mat[gamma_idx, gamma_idx] += k_val

        # Symmetric off-diagonal Gamma and antisymmetric DM coupling D
        j_mat[alpha_idx, beta_idx] += g_val + d_val
        j_mat[beta_idx, alpha_idx] += g_val - d_val

        # Trigonal distortion off-diagonal Gamma'
        j_mat[alpha_idx, gamma_idx] += gp_val
        j_mat[gamma_idx, alpha_idx] += gp_val
        j_mat[beta_idx, gamma_idx] += gp_val
        j_mat[gamma_idx, beta_idx] += gp_val

        return j_mat
