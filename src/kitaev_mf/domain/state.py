"""Immutable mean-field state entity storing ultralocal bond and on-site correlators."""

from dataclasses import dataclass

import numpy as np

from kitaev_mf.domain.lattice import UnitCell


@dataclass(frozen=True)
class MeanFieldState:
    """Ultralocal variational mean-field state (U, V, lambdas).

    Attributes:
        U: Bond correlation matrices U_b^{mu,nu} = <i c_{i,0}^mu c_{j,R_b}^nu>
           for directed bonds b in E. Shape: (num_bonds, 4, 4), float64.
        V: On-site skew-symmetric correlation matrices V_i^{mu,nu} = <i c_i^mu c_i^nu>.
           Shape: (num_sublattices, 4, 4), float64, V_i^T = -V_i.
        lambdas: Local gauge Lagrange multipliers enforcing <c_i^T G^gamma c_i> = 0.
           Shape: (num_sublattices, 3), float64.
    """

    U: np.ndarray
    V: np.ndarray
    lambdas: np.ndarray

    def __post_init__(self) -> None:
        if self.U.ndim != 3 or self.U.shape[1:] != (4, 4):
            raise ValueError(f"U must have shape (num_bonds, 4, 4), got {self.U.shape}")
        if self.V.ndim != 3 or self.V.shape[1:] != (4, 4):
            raise ValueError(f"V must have shape (num_sublattices, 4, 4), got {self.V.shape}")
        if self.lambdas.ndim != 2 or self.lambdas.shape[1] != 3:
            raise ValueError(
                f"lambdas must have shape (num_sublattices, 3), got {self.lambdas.shape}"
            )
        if self.V.shape[0] != self.lambdas.shape[0]:
            raise ValueError(
                f"Number of sublattices in V ({self.V.shape[0]}) does not match "
                f"lambdas ({self.lambdas.shape[0]})"
            )

        # Enforce real skew-symmetry of on-site correlator V_i^T = -V_i
        skew_error = np.max(np.abs(self.V + np.swapaxes(self.V, -1, -2)))
        if skew_error > 1e-12:
            raise ValueError(
                f"On-site correlation matrix V must be skew-symmetric, max error: {skew_error}"
            )

    @classmethod
    def from_pure_kitaev(
        cls,
        cell: UnitCell,
        u_gauge: float = 1.0,
        u_matter: float = 0.0,
    ) -> "MeanFieldState":
        """Initialize the exact zero-flux pure Kitaev spin liquid state.

        On directed bond b = <ij>_gamma oriented from A to B:
          U_b^{gamma, gamma} = u_gauge (canonical convention: +1.0 for A -> B)
          U_b^{0, 0} = u_matter (0.0 for un-iterated seed, ~0.52486 at self-consistency)
          V_i = 0, lambdas_i = 0
        """
        num_bonds = len(cell.bonds)
        num_sub = cell.num_sublattices

        u_arr = np.zeros((num_bonds, 4, 4), dtype=np.float64)
        for b_idx, bond in enumerate(cell.bonds):
            gamma = bond.bond_type
            u_arr[b_idx, gamma, gamma] = float(u_gauge)
            u_arr[b_idx, 0, 0] = float(u_matter)

        v_arr = np.zeros((num_sub, 4, 4), dtype=np.float64)
        lam_arr = np.zeros((num_sub, 3), dtype=np.float64)

        return cls(U=u_arr, V=v_arr, lambdas=lam_arr)
