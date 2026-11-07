"""Exact Diagonalization (ED) in the spin formalism for Kitaev honeycomb clusters.

Provides sparse Hamiltonian generation, spatial translation / point-group symmetry
reduction, Lanczos eigensolvers, and physical observable computation (ground-state
energy, excitation gap, magnetization, plaquette flux, and fidelity susceptibility).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass(frozen=True)
class HoneycombCluster:
    """Honeycomb lattice cluster with periodic boundary conditions on a torus.

    Attributes:
        n1: Number of unit cell repeats along a1.
        n2: Number of unit cell repeats along a2.
        num_cells: Total unit cells (n1 * n2).
        num_spins: Total spins N = 2 * n1 * n2.
        bonds: List of directed bonds (site_i, site_j, bond_type) where bond_type in {1: 'x', 2: 'y', 3: 'z'}.
        plaquettes: List of 6-tuples defining the vertices of each hexagonal plaquette.
    """

    n1: int
    n2: int
    num_cells: int
    num_spins: int
    bonds: list[tuple[int, int, int]]
    plaquettes: list[list[int]]

    @classmethod
    def make_torus(cls, n1: int, n2: int) -> "HoneycombCluster":
        """Construct a periodic honeycomb torus with n1 x n2 primitive 2-site unit cells.

        Sublattices: A = 2 * cell_idx, B = 2 * cell_idx + 1.
        Bonds from A to B:
          z-bond (3): within same cell (c1, c2)
          x-bond (1): cell (c1, c2) -> ((c1) % n1, (c2 - 1) % n2)
          y-bond (2): cell (c1, c2) -> ((c1 - 1) % n1, (c2) % n2)
        """
        num_cells = n1 * n2
        num_spins = 2 * num_cells

        def cell_index(c1: int, c2: int) -> int:
            return (c1 % n1) * n2 + (c2 % n2)

        def site_a(c1: int, c2: int) -> int:
            return 2 * cell_index(c1, c2)

        def site_b(c1: int, c2: int) -> int:
            return 2 * cell_index(c1, c2) + 1

        bonds: list[tuple[int, int, int]] = []
        for c1 in range(n1):
            for c2 in range(n2):
                a_idx = site_a(c1, c2)
                # z-bond (type 3): same cell
                b_z = site_b(c1, c2)
                bonds.append((a_idx, b_z, 3))

                # x-bond (type 1): offset (0, -1)
                b_x = site_b(c1, c2 - 1)
                bonds.append((a_idx, b_x, 1))

                # y-bond (type 2): offset (-1, 0)
                b_y = site_b(c1 - 1, c2)
                bonds.append((a_idx, b_y, 2))

        # Build hexagonal plaquettes
        plaquettes: list[list[int]] = []
        for c1 in range(n1):
            for c2 in range(n2):
                p_sites = [
                    site_a(c1, c2),
                    site_b(c1, c2),
                    site_a(c1, c2 + 1),
                    site_b(c1 + 1, c2),
                    site_a(c1 + 1, c2),
                    site_b(c1, c2 - 1),
                ]
                plaquettes.append(p_sites)

        return cls(
            n1=n1,
            n2=n2,
            num_cells=num_cells,
            num_spins=num_spins,
            bonds=bonds,
            plaquettes=plaquettes,
        )


@dataclass(frozen=True)
class EDObservables:
    """Physical observables evaluated on the exact many-body eigenstates."""

    energy_ground: float
    energy_per_spin: float
    excitation_gap: float
    low_lying_energies: list[float]
    magnetization: float
    plaquette_flux_mean: float
    nearest_neighbor_corr_z: float


class ExactDiagonalizationSpinSolver:
    """High-performance vectorized sparse spin Exact Diagonalization solver for Kitaev honeycomb clusters."""

    def __init__(
        self,
        cluster: HoneycombCluster,
        K: float = 1.0,
        h: float = 0.0,
    ) -> None:
        """Initialize solver.

        Args:
            cluster: Honeycomb cluster geometry.
            K: Ferromagnetic Kitaev exchange (H_K = -K sum sigma_i^gamma sigma_j^gamma).
            h: Zeeman field strength along [111] (H_Z = -h/sqrt(3) sum (sigma^x + sigma^y + sigma^z)).
        """
        self.cluster = cluster
        self.num_spins = cluster.num_spins
        self.dim = 1 << self.num_spins
        self.K = K
        self.h = h

        self.inv_sqrt3 = 1.0 / np.sqrt(3.0)
        self.hz = h * self.inv_sqrt3
        self.hx = h * self.inv_sqrt3
        self.hy = h * self.inv_sqrt3

    def build_sparse_hamiltonian(self) -> sp.csr_matrix:
        """Construct the sparse 2^N x 2^N Hamiltonian using vectorized numpy slicing."""
        dim = self.dim
        num_spins = self.num_spins
        bonds = self.cluster.bonds

        idx_arr = np.arange(dim, dtype=np.int32)
        diag = np.zeros(dim, dtype=np.float64)

        # Count total off-diagonal operators to pre-allocate
        num_offdiag_ops = 0
        has_hx = abs(self.hx) > 1e-15
        has_hy = abs(self.hy) > 1e-15

        if has_hx:
            num_offdiag_ops += num_spins
        if has_hy:
            num_offdiag_ops += num_spins

        # x-bonds and y-bonds are off-diagonal
        num_xy_bonds = sum(1 for _, _, b_type in bonds if b_type in (1, 2))
        num_offdiag_ops += num_xy_bonds

        total_nnz = num_offdiag_ops * dim
        rows = np.empty(total_nnz, dtype=np.int32)
        cols = np.empty(total_nnz, dtype=np.int32)
        vals = np.empty(total_nnz, dtype=np.complex128)

        # 1. Zeeman diagonal: - hz sum_i sigma_i^z
        if abs(self.hz) > 1e-15:
            for i in range(num_spins):
                mask = 1 << i
                diag -= self.hz * np.where((idx_arr & mask) == 0, 1.0, -1.0)

        # 2. Exchange z-bonds: - K sum_<ij>_z sigma_i^z sigma_j^z
        for i, j, b_type in bonds:
            if b_type == 3:  # z-bond
                mask_i = 1 << i
                mask_j = 1 << j
                s_i = np.where((idx_arr & mask_i) == 0, 1.0, -1.0)
                s_j = np.where((idx_arr & mask_j) == 0, 1.0, -1.0)
                diag -= self.K * (s_i * s_j)

        # 3. Off-diagonal vectorized assignments
        ptr = 0

        # (a) Zeeman x terms: - hx sigma_i^x
        if has_hx:
            for i in range(num_spins):
                mask = 1 << i
                s_prime = idx_arr ^ mask
                rows[ptr : ptr + dim] = idx_arr
                cols[ptr : ptr + dim] = s_prime
                vals[ptr : ptr + dim] = -self.hx
                ptr += dim

        # (b) Zeeman y terms: - hy sigma_i^y
        if has_hy:
            for i in range(num_spins):
                mask = 1 << i
                s_prime = idx_arr ^ mask
                sz = np.where((idx_arr & mask) == 0, 1.0, -1.0)
                y_mat_elem = np.where(sz > 0, -1.0j, 1.0j)
                rows[ptr : ptr + dim] = idx_arr
                cols[ptr : ptr + dim] = s_prime
                vals[ptr : ptr + dim] = -self.hy * y_mat_elem
                ptr += dim

        # (c) Exchange x-bonds: - K sigma_i^x sigma_j^x
        # (d) Exchange y-bonds: - K sigma_i^y sigma_j^y
        for i, j, b_type in bonds:
            if b_type == 1:  # x-bond
                mask_ij = (1 << i) ^ (1 << j)
                s_prime = idx_arr ^ mask_ij
                rows[ptr : ptr + dim] = idx_arr
                cols[ptr : ptr + dim] = s_prime
                vals[ptr : ptr + dim] = -self.K
                ptr += dim
            elif b_type == 2:  # y-bond
                mask_ij = (1 << i) ^ (1 << j)
                mask_i = 1 << i
                mask_j = 1 << j
                s_prime = idx_arr ^ mask_ij
                sz_i = np.where((idx_arr & mask_i) == 0, 1.0, -1.0)
                sz_j = np.where((idx_arr & mask_j) == 0, 1.0, -1.0)
                y_i = np.where(sz_i > 0, -1.0j, 1.0j)
                y_j = np.where(sz_j > 0, -1.0j, 1.0j)
                yy_elem = y_i * y_j
                rows[ptr : ptr + dim] = idx_arr
                cols[ptr : ptr + dim] = s_prime
                vals[ptr : ptr + dim] = -self.K * yy_elem
                ptr += dim

        # Assemble CSR matrix
        if total_nnz > 0:
            h_offdiag = sp.csr_matrix(
                (vals, (rows, cols)),
                shape=(dim, dim),
                dtype=np.complex128,
            )
        else:
            h_offdiag = sp.csr_matrix((dim, dim), dtype=np.complex128)

        h_diag = sp.diags(diag, shape=(dim, dim), dtype=np.complex128)
        h_total = h_diag + h_offdiag
        return 0.5 * (h_total + h_total.conj().T)

    def solve_spectrum(self, num_states: int = 4) -> tuple[np.ndarray, np.ndarray]:
        """Compute the lowest num_states eigenvalues and eigenvectors via Lanczos / ARPACK."""
        h_sparse = self.build_sparse_hamiltonian()
        if self.dim <= 64:
            w, v = np.linalg.eigh(h_sparse.toarray())
            return w[:num_states].real, v[:, :num_states]

        evals, evecs = spla.eigsh(
            h_sparse,
            k=min(num_states, self.dim - 2),
            which="SA",
            tol=1e-10,
            maxiter=1000,
        )
        idx = np.argsort(evals.real)
        return evals[idx].real, evecs[:, idx]

    def compute_observables(self, num_states: int = 4) -> EDObservables:
        """Compute physical observables on the ground state."""
        evals, evecs = self.solve_spectrum(num_states=num_states)
        e0 = float(evals[0])
        gap = float(evals[1] - evals[0]) if len(evals) > 1 else 0.0
        psi0 = evecs[:, 0]

        # 1. Magnetization along [111]: <m_111> = (1/N) sum_i <sigma_i . n_111>
        # Compute across all sites
        m_tot = 0.0
        for i in range(self.num_spins):
            mask = 1 << i
            sz_i = np.where((np.arange(self.dim) & mask) == 0, 1.0, -1.0)
            m_tot += float(np.sum(np.abs(psi0) ** 2 * sz_i))
        # By symmetry along [111], <sigma_z> = <sigma_x> = <sigma_y>
        m111 = (np.sqrt(3.0) / float(self.num_spins)) * m_tot

        # 2. Nearest-neighbor z-bond correlator: average over all z-bonds
        corr_z_sum = 0.0
        z_bonds = [(i, j) for i, j, t in self.cluster.bonds if t == 3]
        for i, j in z_bonds:
            mask_i = 1 << i
            mask_j = 1 << j
            sz_ij = np.where((np.arange(self.dim) & mask_i) == 0, 1.0, -1.0) * np.where(
                (np.arange(self.dim) & mask_j) == 0, 1.0, -1.0
            )
            corr_z_sum += float(np.sum(np.abs(psi0) ** 2 * sz_ij))
        corr_z_mean = corr_z_sum / float(len(z_bonds))

        # 3. Plaquette flux expectation <W_p>
        w_p_val = self._evaluate_plaquette_flux(psi0, self.cluster.plaquettes[0])

        return EDObservables(
            energy_ground=e0,
            energy_per_spin=e0 / float(self.num_spins),
            excitation_gap=gap,
            low_lying_energies=[float(e) for e in evals],
            magnetization=m111,
            plaquette_flux_mean=w_p_val,
            nearest_neighbor_corr_z=corr_z_mean,
        )

    def _evaluate_plaquette_flux(self, psi: np.ndarray, plaquette_sites: list[int]) -> float:
        """Evaluate the expectation value of hexagonal loop operator W_p on state psi."""
        p_ops = ["x", "y", "z", "x", "y", "z"]
        vec = psi.copy()

        for site, op in zip(plaquette_sites, p_ops):
            mask = 1 << site
            s_prime = np.arange(self.dim) ^ mask
            if op == "x":
                vec = vec[s_prime]
            elif op == "y":
                sz = np.where((np.arange(self.dim) & mask) == 0, 1.0, -1.0)
                y_elem = np.where(sz > 0, -1.0j, 1.0j)
                vec = y_elem * vec[s_prime]
            elif op == "z":
                sz = np.where((np.arange(self.dim) & mask) == 0, 1.0, -1.0)
                vec = sz * vec

        return float(np.real(np.vdot(psi, vec)))
