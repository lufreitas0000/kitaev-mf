"""Real-space torus cluster Hamiltonian assembly and eigensolver benchmarks.

Represents an L x L periodic honeycomb torus cluster with:
- N_sites = 2 * L^2 physical sites
- M = 8 * L^2 Majorana fermion modes
- Exact real skew-symmetric matrix representation A in R^{M x M} (A^T = -A).
"""

from dataclasses import dataclass

import numpy as np
import torch

from kitaev_mf.domain.lattice import DEFAULT_D_NN, make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.hamiltonian import (
    build_bond_hopping_matrices,
    build_on_site_matrices,
)


@dataclass(frozen=True)
class RealSpaceTorusCluster:
    """Honeycomb torus cluster of L x L unit cells with periodic boundary conditions.

    Attributes:
        L: Linear dimension of the torus (L >= 2).
        d_nn: Nearest-neighbor bond distance.
    """

    L: int
    d_nn: float = DEFAULT_D_NN

    @property
    def num_cells(self) -> int:
        """Total unit cells N_cells = L^2."""
        return self.L * self.L

    @property
    def num_sites(self) -> int:
        """Total physical lattice sites N_sites = 2 * L^2."""
        return 2 * self.num_cells

    @property
    def num_majoranas(self) -> int:
        """Total Majorana fermion modes M = 4 * N_sites = 8 * L^2."""
        return 4 * self.num_sites

    def site_index(self, m1: int, m2: int, sublattice: int) -> int:
        """Map unit cell coordinate (m1, m2) in [0, L)^2 and sublattice in {0, 1} to site index."""
        r1 = m1 % self.L
        r2 = m2 % self.L
        return 2 * (r1 + self.L * r2) + sublattice

    def majorana_slice(self, site: int) -> slice:
        """Return index slice of the 4 Majoranas [c_{s}^0, c_{s}^1, c_{s}^2, c_{s}^3] at physical site s."""
        return slice(4 * site, 4 * site + 4)

    def assemble_skew_hamiltonian(
        self,
        params: ModelParameters,
        state: MeanFieldState,
    ) -> np.ndarray:
        """Assemble the exact real skew-symmetric M x M matrix A where H = (i/4) c^T A c."""
        cell2 = make_honeycomb_2site(self.d_nn)
        device = torch.device("cpu")

        # Use canonical production services for on-site B_s and hopping A_b
        b_mat = build_on_site_matrices(cell2, params, state, device).cpu().numpy()
        a_mat = build_bond_hopping_matrices(cell2, params, state, device).cpu().numpy()

        M = self.num_majoranas
        A = np.zeros((M, M), dtype=np.float64)

        # 1. On-site skew-symmetric blocks
        for s in range(self.num_sites):
            sub = s % cell2.num_sublattices
            sl = self.majorana_slice(s)
            A[sl, sl] = b_mat[sub]

        # 2. Directed bond hopping blocks
        disp = [(0, 0, 3), (0, -1, 1), (-1, 0, 2)]

        for m2 in range(self.L):
            for m1 in range(self.L):
                s_A = self.site_index(m1, m2, 0)
                sl_A = self.majorana_slice(s_A)

                for b_idx, (d1, d2, _) in enumerate(disp):
                    s_B = self.site_index(m1 + d1, m2 + d2, 1)
                    sl_B = self.majorana_slice(s_B)

                    # In Bloch Hamiltonian, block (s_A, s_B) is +(i/4) * A_b
                    # In real-space Majorana form H = (i/4) c^T A c:
                    # c_A^T A_b c_B => in skew-symmetric A:
                    # A[A, B] = A_b, A[B, A] = -A_b^T
                    A[sl_A, sl_B] += a_mat[b_idx]
                    A[sl_B, sl_A] -= a_mat[b_idx].T

        # Enforce exact skew-symmetry
        return np.asarray(0.5 * (A - A.T), dtype=np.float64)

    def solve_covariance_complex_eigh(self, A: np.ndarray) -> tuple[np.ndarray, float]:
        """Compute ground-state covariance matrix Omega = i sgn(i A) via complex Hermitian eigensolver zheevd.

        Returns:
            Tuple of (Omega, ground_state_energy).
        """
        # H_c = i * A is Hermitian: H_c^\dagger = -i A^T = i A = H_c
        H_c = 1j * A.astype(np.complex128)
        w, v = np.linalg.eigh(H_c)

        # Negative energy projector: P^- = sum_{eps_k < 0} |v_k><v_k|
        # Omega = i * (P^+ - P^-) = i * (I - 2 P^-)
        # Or equivalently: Omega = i * v @ diag(sgn(w)) @ v^\dagger
        sgn_w = np.sign(w)
        # Handle zero modes if any
        sgn_w[np.abs(w) < 1e-12] = 0.0

        omega_c = 1j * (v @ (sgn_w[:, None] * v.conj().T))
        omega_real = np.real(omega_c)

        # Ground state energy: E = -0.5 * sum_{eps_k > 0} eps_k = 0.5 * sum_{w < 0} w
        gs_energy = 0.5 * float(np.sum(w[w < 0]))
        return omega_real, gs_energy

    def solve_covariance_real_eigh(self, A: np.ndarray) -> tuple[np.ndarray, float]:
        """Compute ground-state covariance matrix Omega via pure real symmetric eigensolver dsyevd on A^2.

        Since A is real skew-symmetric, S = -A^2 = A A^T is real symmetric positive semidefinite.
        Its eigenvalues are eps_k^2 with eigenvectors identical to A.
        Omega = - A @ S^{-1/2} = i sgn(i A).
        """
        S = -A @ A
        # S is real symmetric
        S = 0.5 * (S + S.T)
        w2, v = np.linalg.eigh(S)
        w2 = np.maximum(w2, 0.0)
        w = np.sqrt(w2)

        # Omega = - A @ (S)^{-1/2} = - A @ v @ diag(1 / sqrt(w2)) @ v^T
        # Regularize near zero
        inv_w = np.zeros_like(w)
        mask = w > 1e-11
        inv_w[mask] = 1.0 / w[mask]

        omega = -A @ (v @ (inv_w[:, None] * v.T))
        omega = 0.5 * (omega - omega.T)

        # S = A A^T has doubled eigenvalues eps_k^2, so sum(w) = 2 * sum(eps_k > 0)
        # Therefore E_GS = -0.5 * sum(eps_k > 0) = -0.25 * sum(w)
        gs_energy = -0.25 * float(np.sum(w))
        return omega, gs_energy
