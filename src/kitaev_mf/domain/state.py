"""Immutable mean-field state and energy breakdown value objects."""

from dataclasses import dataclass

import numpy as np

from kitaev_mf.domain.lattice import UnitCell


@dataclass(frozen=True)
class EnergyBreakdown:
    """Variational Gaussian energy components and normalized densities.

    Attributes:
        E_exchange: Extensive exchange interaction energy (summed over directed bonds).
        E_field: Extensive external Zeeman energy.
        E_multiplier: Extensive gauge constraint penalty energy.
        E_total: Extensive total energy E_exchange + E_field + E_multiplier.
        num_bonds: Number of physical directed bonds in the cell.
        num_cells: Number of unit cells represented (typically 1 for local cell).
        num_sites: Number of physical lattice sites (sublattices) in the cell.
    """

    E_exchange: float
    E_field: float
    E_multiplier: float
    E_total: float
    num_bonds: int
    num_cells: int
    num_sites: int

    @property
    def per_bond(self) -> float:
        """Energy normalized per directed bond."""
        return self.E_total / float(self.num_bonds)

    @property
    def per_cell(self) -> float:
        """Energy normalized per unit cell."""
        return self.E_total / float(self.num_cells)

    @property
    def per_site(self) -> float:
        """Energy normalized per physical site."""
        return self.E_total / float(self.num_sites)


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

    @classmethod
    def from_polarized(
        cls,
        cell: UnitCell,
        n_vec: tuple[float, float, float] | np.ndarray | None = None,
    ) -> "MeanFieldState":
        """Initialize an unentangled, gauge-neutral polarized state with magnetization n_vec.

        Using the SO(4) spin generators M^a, the on-site matrix V_s = sum_a n_a M^a satisfies:
          <sigma_s^a> = n_a,
          tr(G^gamma V_s) = 0 for all gamma in {1, 2, 3} (strictly gauge neutral),
          V_s^2 = -||n||^2 * I_4 (pure state when ||n|| = 1),
          U_b = 0, lambdas_s = 0.
        """
        from kitaev_mf.domain.constants import M_MATRICES

        if n_vec is None:
            inv_sqrt3 = 1.0 / float(np.sqrt(3.0))
            n_arr = np.array([inv_sqrt3, inv_sqrt3, inv_sqrt3], dtype=np.float64)
        else:
            n_arr = np.asarray(n_vec, dtype=np.float64)
        if n_arr.shape != (3,):
            raise ValueError(f"Polarization vector must have length 3, got {n_arr.shape}")

        num_bonds = len(cell.bonds)
        num_sub = cell.num_sublattices

        u_arr = np.zeros((num_bonds, 4, 4), dtype=np.float64)
        v_single = n_arr[0] * M_MATRICES[1] + n_arr[1] * M_MATRICES[2] + n_arr[2] * M_MATRICES[3]
        v_arr = np.stack([v_single.copy() for _ in range(num_sub)], axis=0)
        lam_arr = np.zeros((num_sub, 3), dtype=np.float64)

        return cls(U=u_arr, V=v_arr, lambdas=lam_arr)


@dataclass(frozen=True)
class InhomogeneousMeanFieldState:
    """Real-space inhomogeneous variational mean-field state (U, V, lambdas).

    Indexed to preserve spatial locality and avoid artificial translational averaging:
    - U: Bond correlation matrices U_{(m1, m2), b}^{mu, nu} = <i c_{A(m1, m2)}^mu c_{B(m1+d1, m2+d2)}^nu>
         for the 3 directed bonds emanating from site A in cell (m1, m2).
         Shape: (num_cells, 3, 4, 4), float64, where num_cells = L^2.
         The second dimension b in {0, 1, 2} indexes bonds:
           b=0: z-bond (gamma=3), displacement (0, 0)
           b=1: x-bond (gamma=1), displacement (0, -1)
           b=2: y-bond (gamma=2), displacement (-1, 0)
    - V: On-site skew-symmetric correlation matrices V_{(m1, m2), sub}^{mu, nu} = <i c_s^mu c_s^nu>
         for sublattices sub in {0, 1} (0 for A, 1 for B).
         Shape: (num_cells, 2, 4, 4), float64, satisfying V^T = -V.
    - lambdas: Local gauge Lagrange multipliers enforcing <c_s^T G^gamma c_s> = 0 on each site.
         Shape: (num_cells, 2, 3), float64.
    """

    U: np.ndarray
    V: np.ndarray
    lambdas: np.ndarray

    def __post_init__(self) -> None:
        if self.U.ndim != 4 or self.U.shape[1:] != (3, 4, 4):
            raise ValueError(f"U must have shape (num_cells, 3, 4, 4), got {self.U.shape}")
        if self.V.ndim != 4 or self.V.shape[1:] != (2, 4, 4):
            raise ValueError(f"V must have shape (num_cells, 2, 4, 4), got {self.V.shape}")
        if self.lambdas.ndim != 3 or self.lambdas.shape[1:] != (2, 3):
            raise ValueError(f"lambdas must have shape (num_cells, 2, 3), got {self.lambdas.shape}")
        num_cells = self.U.shape[0]
        if self.V.shape[0] != num_cells or self.lambdas.shape[0] != num_cells:
            raise ValueError(
                f"Cell count mismatch: U ({num_cells}), V ({self.V.shape[0]}), "
                f"lambdas ({self.lambdas.shape[0]})"
            )

        skew_error = float(np.max(np.abs(self.V + np.swapaxes(self.V, -1, -2))))
        if skew_error > 1e-12:
            raise ValueError(
                f"On-site correlation tensor V must be skew-symmetric, max error: {skew_error}"
            )

    @property
    def num_cells(self) -> int:
        """Number of unit cells L^2."""
        return int(self.U.shape[0])

    @property
    def num_sites(self) -> int:
        """Total physical lattice sites 2 * L^2."""
        return 2 * self.num_cells

    @property
    def num_bonds(self) -> int:
        """Total directed bonds 3 * L^2."""
        return 3 * self.num_cells

    @property
    def dim(self) -> int:
        """Total variational dimension D = 86 * L^2."""
        return 86 * self.num_cells

    def flatten(self) -> np.ndarray:
        """Pack all variational parameters into a 1D vector of length 86 * num_cells."""
        return np.concatenate(
            [self.U.ravel(), self.V.ravel(), self.lambdas.ravel()],
            dtype=np.float64,
        )

    @classmethod
    def from_flattened(cls, x: np.ndarray, num_cells: int) -> "InhomogeneousMeanFieldState":
        """Reconstruct state from 1D parameter vector, projecting V onto skew-symmetric subspace."""
        len_u = num_cells * 3 * 16
        len_v = num_cells * 2 * 16
        len_lam = num_cells * 2 * 3
        expected_dim = len_u + len_v + len_lam
        if x.size != expected_dim:
            raise ValueError(f"Expected flattened vector of size {expected_dim}, got {x.size}")

        u_arr = x[:len_u].reshape((num_cells, 3, 4, 4)).copy()
        v_arr = x[len_u : len_u + len_v].reshape((num_cells, 2, 4, 4)).copy()
        v_arr = 0.5 * (v_arr - np.swapaxes(v_arr, -1, -2))
        lam_arr = x[len_u + len_v :].reshape((num_cells, 2, 3)).copy()

        return cls(U=u_arr, V=v_arr, lambdas=lam_arr)

    def to_homogeneous_mean(self) -> MeanFieldState:
        """Compute the translationally averaged MeanFieldState across all unit cells."""
        u_mean = np.mean(self.U, axis=0)
        v_mean = np.mean(self.V, axis=0)
        v_mean = 0.5 * (v_mean - np.swapaxes(v_mean, -1, -2))
        lam_mean = np.mean(self.lambdas, axis=0)
        return MeanFieldState(U=u_mean, V=v_mean, lambdas=lam_mean)

    @classmethod
    def from_homogeneous(
        cls, state: MeanFieldState, num_cells: int
    ) -> "InhomogeneousMeanFieldState":
        """Replicate a translationally invariant MeanFieldState across num_cells."""
        u_arr = np.repeat(state.U[np.newaxis, ...], num_cells, axis=0)
        v_arr = np.repeat(state.V[np.newaxis, ...], num_cells, axis=0)
        lam_arr = np.repeat(state.lambdas[np.newaxis, ...], num_cells, axis=0)
        return cls(U=u_arr.copy(), V=v_arr.copy(), lambdas=lam_arr.copy())

    @classmethod
    def from_pure_kitaev(
        cls,
        num_cells: int,
        u_gauge: float = 1.0,
        u_matter: float = 0.0,
    ) -> "InhomogeneousMeanFieldState":
        """Initialize the exact zero-flux pure Kitaev state on num_cells."""
        u_arr = np.zeros((num_cells, 3, 4, 4), dtype=np.float64)
        # b=0: z-bond (gamma=3)
        u_arr[:, 0, 3, 3] = float(u_gauge)
        u_arr[:, 0, 0, 0] = float(u_matter)
        # b=1: x-bond (gamma=1)
        u_arr[:, 1, 1, 1] = float(u_gauge)
        u_arr[:, 1, 0, 0] = float(u_matter)
        # b=2: y-bond (gamma=2)
        u_arr[:, 2, 2, 2] = float(u_gauge)
        u_arr[:, 2, 0, 0] = float(u_matter)

        v_arr = np.zeros((num_cells, 2, 4, 4), dtype=np.float64)
        lam_arr = np.zeros((num_cells, 2, 3), dtype=np.float64)

        return cls(U=u_arr, V=v_arr, lambdas=lam_arr)
