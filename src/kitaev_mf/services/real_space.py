"""Real-space torus cluster Hamiltonian assembly and eigensolver benchmarks.

Represents an L x L periodic honeycomb torus cluster with:
- N_sites = 2 * L^2 physical sites
- M = 8 * L^2 Majorana fermion modes
- Exact real skew-symmetric matrix representation A in R^{M x M} (A^T = -A).
"""

from collections import deque
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, Literal

import numpy as np
import torch

from kitaev_mf.domain.constants import G_MATRICES, N_MATRICES
from kitaev_mf.domain.lattice import DEFAULT_D_NN, make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import InhomogeneousMeanFieldState, MeanFieldState
from kitaev_mf.services.fixed_point import JointAndersonSolver, LimitedMemoryBroydenSolver
from kitaev_mf.services.hamiltonian import (
    build_bond_hopping_matrices,
    build_on_site_matrices,
)
from kitaev_mf.services.purification import purify_newton_schulz


def compute_skew_pfaffian(A: np.ndarray, tol: float = 1e-14) -> float:
    """Compute the exact Pfaffian of a real 2n x 2n skew-symmetric matrix via Parlett-Reid elimination.

    Args:
        A: Real skew-symmetric matrix of shape (2n, 2n).
        tol: Numerical tolerance for zero pivots.

    Returns:
        Exact real Pfaffian Pf(A) satisfying Pf(A)^2 = det(A).
    """
    n = A.shape[0]
    if n % 2 != 0:
        return 0.0
    M = A.copy().astype(np.float64)

    sign = 1.0
    pf = 1.0

    for k in range(0, n, 2):
        sub = np.abs(M[k:, k:])
        max_idx = np.unravel_index(np.argmax(sub), sub.shape)
        p, q = int(max_idx[0] + k), int(max_idx[1] + k)

        if sub[max_idx] < tol:
            return 0.0

        if p != k:
            M[[k, p], :] = M[[p, k], :]
            M[:, [k, p]] = M[:, [p, k]]
            sign *= -1.0

        if q == k:
            q = p
        if q != k + 1:
            M[[k + 1, q], :] = M[[q, k + 1], :]
            M[:, [k + 1, q]] = M[:, [q, k + 1]]
            sign *= -1.0

        a = M[k, k + 1]
        pf *= a

        if k + 2 < n:
            b0 = M[k, (k + 2) :]
            b1 = M[k + 1, (k + 2) :]
            M[(k + 2) :, (k + 2) :] += (np.outer(b1, b0) - np.outer(b0, b1)) / a

    return float(sign * pf)


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
        state: MeanFieldState | InhomogeneousMeanFieldState,
        flipped_bonds: Sequence[tuple[int, int, int]] | set[tuple[int, int, int]] | None = None,
        bond_disorder: dict[tuple[int, int, int], float] | None = None,
        site_vacancies: set[int] | Sequence[int] | None = None,
    ) -> np.ndarray:
        """Assemble the exact real skew-symmetric M x M matrix A where H = (i/4) c^T A c.

        Supports both translationally invariant MeanFieldState and spatially resolved
        InhomogeneousMeanFieldState without loss of defect localization.

        Args:
            params: Exchange couplings and external field.
            state: Mean-field correlation state (homogeneous or inhomogeneous).
            flipped_bonds: Optional collection of directed bond keys (m1, m2, b_idx)
                           whose gauge link signs are inverted (u_b -> -u_b).
            bond_disorder: Optional mapping from bond key (m1, m2, b_idx) to relative
                           exchange modulation factor delta_b such that K_b = K * (1 + delta_b).
            site_vacancies: Optional collection of physical site indices to remove
                            (dilute non-magnetic vacancies with zeroed couplings).
        """
        cell2 = make_honeycomb_2site(self.d_nn)
        flipped_set = set(flipped_bonds) if flipped_bonds is not None else set()
        vacancies_set = set(site_vacancies) if site_vacancies is not None else set()
        disorder_map = bond_disorder if bond_disorder is not None else {}
        M = self.num_majoranas
        A = np.zeros((M, M), dtype=np.float64)
        disp = [(0, 0, 3), (0, -1, 1), (-1, 0, 2)]

        if isinstance(state, InhomogeneousMeanFieldState):
            if state.num_cells != self.num_cells:
                raise ValueError(
                    f"State cell count {state.num_cells} does not match cluster {self.num_cells}"
                )

            # 1. On-site skew-symmetric blocks B_s
            # Precompute local spin expectations: tr(V_s^T N^beta)
            tr_v_n = np.zeros((self.num_sites, 3), dtype=np.float64)
            for m2 in range(self.L):
                for m1 in range(self.L):
                    c = m1 + self.L * m2
                    for sub in (0, 1):
                        s = self.site_index(m1, m2, sub)
                        v_s = state.V[c, sub]
                        for beta in (1, 2, 3):
                            tr_v_n[s, beta - 1] = float(np.trace(v_s.T @ N_MATRICES[beta]))

            # On-site accumulation: Zeeman and Lagrange multipliers
            for m2 in range(self.L):
                for m1 in range(self.L):
                    c = m1 + self.L * m2
                    for sub in (0, 1):
                        s = self.site_index(m1, m2, sub)
                        if s in vacancies_set:
                            continue
                        sl = self.majorana_slice(s)
                        b_s = np.zeros((4, 4), dtype=np.float64)

                        lam = state.lambdas[c, sub]
                        for gamma in (1, 2, 3):
                            if lam[gamma - 1] != 0.0:
                                b_s += lam[gamma - 1] * G_MATRICES[gamma]
                            h_val = params.h_vec[gamma - 1]
                            if h_val != 0.0:
                                b_s -= 2.0 * h_val * N_MATRICES[gamma]

                        A[sl, sl] = b_s

            # Neighbor Hartree contributions to on-site B_s:
            for m2 in range(self.L):
                for m1 in range(self.L):
                    c = m1 + self.L * m2
                    s_A = self.site_index(m1, m2, 0)
                    sl_A = self.majorana_slice(s_A)

                    for _b_idx, (d1, d2, gamma) in enumerate(disp):
                        s_B = self.site_index(m1 + d1, m2 + d2, 1)
                        sl_B = self.majorana_slice(s_B)

                        if s_A in vacancies_set or s_B in vacancies_set:
                            continue

                        j_mat = params.exchange_matrix(gamma)
                        for alpha in (1, 2, 3):
                            for beta in (1, 2, 3):
                                c_A = j_mat[alpha - 1, beta - 1] * tr_v_n[s_B, beta - 1]
                                if c_A != 0.0:
                                    A[sl_A, sl_A] += c_A * N_MATRICES[alpha]
                                c_B = j_mat[beta - 1, alpha - 1] * tr_v_n[s_A, beta - 1]
                                if c_B != 0.0:
                                    A[sl_B, sl_B] += c_B * N_MATRICES[alpha]

            # 2. Directed bond hopping blocks
            for m2 in range(self.L):
                for m1 in range(self.L):
                    c = m1 + self.L * m2
                    s_A = self.site_index(m1, m2, 0)
                    sl_A = self.majorana_slice(s_A)

                    for b_idx, (d1, d2, gamma) in enumerate(disp):
                        s_B = self.site_index(m1 + d1, m2 + d2, 1)
                        sl_B = self.majorana_slice(s_B)

                        if s_A in vacancies_set or s_B in vacancies_set:
                            continue

                        bond_key = (m1, m2, b_idx)
                        sign = -1.0 if bond_key in flipped_set else 1.0
                        disorder_factor = 1.0 + disorder_map.get(bond_key, 0.0)

                        j_mat = params.exchange_matrix(gamma)
                        u_b = state.U[c, b_idx]
                        a_b = np.zeros((4, 4), dtype=np.float64)
                        for alpha in (1, 2, 3):
                            for beta in (1, 2, 3):
                                j_ab = j_mat[alpha - 1, beta - 1]
                                if j_ab != 0.0:
                                    a_b += 2.0 * j_ab * (N_MATRICES[alpha] @ u_b @ N_MATRICES[beta])

                        A[sl_A, sl_B] += sign * disorder_factor * a_b
                        A[sl_B, sl_A] -= sign * disorder_factor * a_b.T

            return np.asarray(0.5 * (A - A.T), dtype=np.float64)

        # Homogeneous fallback
        device = torch.device("cpu")
        b_mat = build_on_site_matrices(cell2, params, state, device).cpu().numpy()
        a_mat = build_bond_hopping_matrices(cell2, params, state, device).cpu().numpy()

        # 1. On-site skew-symmetric blocks
        for s in range(self.num_sites):
            if s in vacancies_set:
                continue
            sub = s % cell2.num_sublattices
            sl = self.majorana_slice(s)
            A[sl, sl] = b_mat[sub]

        # 2. Directed bond hopping blocks
        for m2 in range(self.L):
            for m1 in range(self.L):
                s_A = self.site_index(m1, m2, 0)
                sl_A = self.majorana_slice(s_A)

                for b_idx, (d1, d2, _) in enumerate(disp):
                    s_B = self.site_index(m1 + d1, m2 + d2, 1)
                    sl_B = self.majorana_slice(s_B)

                    if s_A in vacancies_set or s_B in vacancies_set:
                        continue

                    bond_key = (m1, m2, b_idx)
                    sign = -1.0 if bond_key in flipped_set else 1.0
                    disorder_factor = 1.0 + disorder_map.get(bond_key, 0.0)

                    A[sl_A, sl_B] += sign * disorder_factor * a_mat[b_idx]
                    A[sl_B, sl_A] -= sign * disorder_factor * a_mat[b_idx].T

        return np.asarray(0.5 * (A - A.T), dtype=np.float64)

    def build_dual_graph(
        self,
    ) -> dict[tuple[int, int], list[tuple[tuple[int, int], tuple[int, int, int]]]]:
        """Construct the dual triangular lattice graph connecting hexagonal plaquettes.

        Returns:
            Dictionary mapping each plaquette coordinate (m1, m2) to a list of
            tuples: (neighbor_plaquette, crossing_bond_key), where crossing_bond_key
            is (m1, m2, b_idx).
        """
        L = self.L
        bond_to_plaquettes: dict[tuple[int, int, int], list[tuple[int, int]]] = {}

        for m2 in range(L):
            for m1 in range(L):
                p = (m1, m2)
                bonds = [
                    (m1, m2, 0),
                    (m1, (m2 + 1) % L, 1),
                    (m1, (m2 + 1) % L, 2),
                    ((m1 - 1) % L, (m2 + 1) % L, 0),
                    ((m1 - 1) % L, (m2 + 1) % L, 1),
                    (m1, m2, 2),
                ]
                for b in bonds:
                    bond_to_plaquettes.setdefault(b, []).append(p)

        adj: dict[tuple[int, int], list[tuple[tuple[int, int], tuple[int, int, int]]]] = {
            (m1, m2): [] for m2 in range(L) for m1 in range(L)
        }
        for b, pls in bond_to_plaquettes.items():
            pA, pB = pls[0], pls[1]
            adj[pA].append((pB, b))
            adj[pB].append((pA, b))

        return adj

    def generate_vortex_pair_cut(
        self,
        p1: tuple[int, int],
        p2: tuple[int, int],
    ) -> list[tuple[int, int, int]]:
        """Generate a minimal-length branch cut (string of bonds) connecting plaquettes p1 and p2.

        Inverting the gauge links along this cut produces an isolated pair of vortices
        with W_{p1} = -1 and W_{p2} = -1, leaving all other plaquettes flux-free (W = +1).

        Args:
            p1: Source plaquette (m1, m2).
            p2: Target plaquette (m1, m2).

        Returns:
            List of directed bond keys (m1, m2, b_idx) forming the dual branch cut.
        """
        adj = self.build_dual_graph()
        visited = {p1}
        queue: deque[tuple[tuple[int, int], list[tuple[int, int, int]]]] = deque([(p1, [])])

        while queue:
            curr, path_bonds = queue.popleft()
            if curr == p2:
                return path_bonds
            for nbr, b in adj[curr]:
                if nbr not in visited:
                    visited.add(nbr)
                    queue.append((nbr, path_bonds + [b]))

        return []

    def create_unpinned_vortex_pair_seed(
        self,
        base_state: InhomogeneousMeanFieldState,
        p1: tuple[int, int],
        p2: tuple[int, int],
    ) -> tuple[InhomogeneousMeanFieldState, list[tuple[int, int, int]]]:
        """Create an unpinned two-vortex mean-field seed on the torus.

        Finds the minimal-length dual branch cut connecting plaquettes p1 and p2,
        and inverts the corresponding gauge link correlators U_{c, b}^{gamma, gamma} -> -1.0.

        Args:
            base_state: Inhomogeneous base state (e.g. from_pure_kitaev).
            p1: Source plaquette coordinate (m1, m2).
            p2: Target plaquette coordinate (m1, m2).

        Returns:
            Tuple of (unpinned_vortex_state, branch_cut_bonds).
        """
        cut = self.generate_vortex_pair_cut(p1, p2)
        disp_types = [3, 1, 2]

        u_new = base_state.U.copy()
        for m1, m2, b_idx in cut:
            c = m1 + self.L * m2
            gamma = disp_types[b_idx]
            curr_val = u_new[c, b_idx, gamma, gamma]
            u_new[c, b_idx, gamma, gamma] = -abs(curr_val) if curr_val != 0.0 else -1.0

        vortex_state = InhomogeneousMeanFieldState(
            U=u_new,
            V=base_state.V.copy(),
            lambdas=base_state.lambdas.copy(),
        )
        return vortex_state, cut

    def compute_plaquette_fluxes(
        self,
        omega: np.ndarray,
    ) -> dict[tuple[int, int], float]:
        """Compute all L^2 hexagonal Wilson loops W_p from the Majorana covariance matrix Omega.

        For each plaquette p, W_p = prod_{b in partial p} u_b in [-1, 1], where
        u_b = Omega_{4 s_A + gamma, 4 s_B + gamma} is the gauge link expectation.

        Args:
            omega: Skew-symmetric Majorana covariance matrix in R^{M x M}.

        Returns:
            Dictionary mapping plaquette coordinate (m1, m2) to Wilson loop W_p.
        """
        L = self.L
        disp = [(0, 0, 3), (0, -1, 1), (-1, 0, 2)]

        def get_bond_u(m1: int, m2: int, b_idx: int) -> float:
            s_A = self.site_index(m1, m2, 0)
            d1, d2, gamma = disp[b_idx]
            s_B = self.site_index(m1 + d1, m2 + d2, 1)
            idx_A = 4 * s_A + gamma
            idx_B = 4 * s_B + gamma
            return float(omega[idx_A, idx_B])

        fluxes: dict[tuple[int, int], float] = {}
        for m2 in range(L):
            for m1 in range(L):
                ub1 = get_bond_u(m1, m2, 0)
                ub2 = get_bond_u(m1, (m2 + 1) % L, 1)
                ub3 = get_bond_u(m1, (m2 + 1) % L, 2)
                ub4 = get_bond_u((m1 - 1) % L, (m2 + 1) % L, 0)
                ub5 = get_bond_u((m1 - 1) % L, (m2 + 1) % L, 1)
                ub6 = get_bond_u(m1, m2, 2)
                fluxes[(m1, m2)] = ub1 * ub2 * ub3 * ub4 * ub5 * ub6

        return fluxes

    def compute_plaquette_gauge_majoranas(self, p: tuple[int, int]) -> list[int]:
        """Return the 12 gauge Majorana indices on the perimeter of plaquette p = (m1, m2).

        The indices are grouped into 6 bond pairs (A_k, B_k) around the hexagonal perimeter:
        b1 (z-bond): (m1, m2) A <-> (m1, m2) B
        b2 (x-bond): (m1, m2+1) A <-> (m1, m2) B
        b3 (y-bond): (m1, m2+1) A <-> (m1-1, m2+1) B
        b4 (z-bond): (m1-1, m2+1) A <-> (m1-1, m2+1) B
        b5 (x-bond): (m1-1, m2+1) A <-> (m1-1, m2) B
        b6 (y-bond): (m1, m2) A <-> (m1-1, m2) B

        Args:
            p: Plaquette coordinate (m1, m2).

        Returns:
            List of 12 integer Majorana mode indices in {0, ..., 8*L^2 - 1}.
        """
        m1, m2 = p
        L = self.L
        bonds = [
            (m1, m2, 3, 0, 0),
            (m1, (m2 + 1) % L, 1, 0, -1),
            (m1, (m2 + 1) % L, 2, -1, 0),
            ((m1 - 1) % L, (m2 + 1) % L, 3, 0, 0),
            ((m1 - 1) % L, (m2 + 1) % L, 1, 0, -1),
            (m1, m2, 2, -1, 0),
        ]
        indices = []
        for m1_a, m2_a, gamma, d1, d2 in bonds:
            s_A = self.site_index(m1_a, m2_a, 0)
            s_B = self.site_index(m1_a + d1, m2_a + d2, 1)
            indices.extend([4 * s_A + gamma, 4 * s_B + gamma])
        return indices

    def compute_plaquette_pfaffian_flux(
        self,
        omega: np.ndarray,
        p: tuple[int, int],
    ) -> float:
        """Compute the exact 12-Majorana Pfaffian Wilson loop W_p = Pf(Omega|_{partial p}).

        Args:
            omega: Skew-symmetric Majorana covariance matrix in R^{M x M}.
            p: Plaquette coordinate (m1, m2).

        Returns:
            Real-valued Pfaffian expectation in [-1, 1].
        """
        indices = self.compute_plaquette_gauge_majoranas(p)
        sub_omega = omega[np.ix_(indices, indices)]
        return compute_skew_pfaffian(sub_omega)

    def compute_all_plaquette_pfaffian_fluxes(
        self,
        omega: np.ndarray,
    ) -> dict[tuple[int, int], float]:
        """Compute exact 12-Majorana Pfaffian fluxes W_p across all L^2 plaquettes.

        Args:
            omega: Skew-symmetric Majorana covariance matrix in R^{M x M}.

        Returns:
            Dictionary mapping plaquette coordinate (m1, m2) to Pfaffian flux.
        """
        return {
            (m1, m2): self.compute_plaquette_pfaffian_flux(omega, (m1, m2))
            for m2 in range(self.L)
            for m1 in range(self.L)
        }

    def compute_plaquette_boundary_sites(self, p: tuple[int, int]) -> set[int]:
        """Compute the 6 physical site indices bounding hexagonal plaquette p = (m1, m2).

        Args:
            p: Plaquette coordinate (m1, m2).

        Returns:
            Set of 6 physical site indices in {0, ..., 2*L^2 - 1}.
        """
        m1, m2 = p
        L = self.L
        return {
            self.site_index(m1, m2, 0),
            self.site_index(m1, m2, 1),
            self.site_index(m1, (m2 + 1) % L, 0),
            self.site_index((m1 - 1) % L, (m2 + 1) % L, 1),
            self.site_index((m1 - 1) % L, (m2 + 1) % L, 0),
            self.site_index((m1 - 1) % L, m2, 1),
        }

    def compute_vortex_mode_core_weight(
        self,
        psi: np.ndarray,
        p1: tuple[int, int],
        p2: tuple[int, int],
    ) -> tuple[float, float, float]:
        """Compute the core localization weight W_core of mode psi around two vortices at p1 and p2.

        Args:
            psi: Single-particle eigenvector of dimension M = 8*L^2 (or 2*L^2).
            p1: Plaquette coordinate of first vortex.
            p2: Plaquette coordinate of second vortex.

        Returns:
            Tuple of (W_core, W_uniform, ratio):
                - W_core: Sum of probability density on the boundary sites of p1 and p2.
                - W_uniform: Expected weight for a uniformly delocalized state (|C| / N_sites).
                - ratio: W_core / W_uniform enrichment factor.
        """
        core_sites = self.compute_plaquette_boundary_sites(
            p1
        ) | self.compute_plaquette_boundary_sites(p2)
        num_core_sites = len(core_sites)
        num_total_sites = self.num_sites  # 2 * L^2

        psi_norm_sq = np.abs(psi) ** 2
        total_norm = float(np.sum(psi_norm_sq))
        if total_norm < 1e-15:
            return 0.0, float(num_core_sites / num_total_sites), 0.0

        if len(psi) == self.num_majoranas:
            site_weights = np.sum(psi_norm_sq.reshape((num_total_sites, 4)), axis=1) / total_norm
        elif len(psi) == num_total_sites:
            site_weights = psi_norm_sq / total_norm
        else:
            raise ValueError(
                f"Eigenvector length {len(psi)} must be {self.num_majoranas} or {num_total_sites}"
            )

        w_core = float(np.sum([site_weights[s] for s in core_sites]))
        w_uniform = float(num_core_sites / num_total_sites)
        ratio = float(w_core / w_uniform) if w_uniform > 0.0 else 0.0
        return w_core, w_uniform, ratio

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

    def extract_mean_field_state(
        self,
        omega: np.ndarray,
        current_lambdas: np.ndarray,
    ) -> MeanFieldState:
        """Extract translationally averaged bond (U) and on-site (V) correlators from Omega.

        Args:
            omega: Skew-symmetric covariance matrix in R^{M x M}.
            current_lambdas: Current Lagrange multipliers of shape (num_sublattices, 3).

        Returns:
            MeanFieldState representing the translation-averaged correlation tensors.
        """
        disp = [(0, 0, 3), (0, -1, 1), (-1, 0, 2)]
        num_bonds = len(disp)
        num_cells = self.num_cells

        U_out = np.zeros((num_bonds, 4, 4), dtype=np.float64)
        for b_idx, (d1, d2, _) in enumerate(disp):
            u_acc = np.zeros((4, 4), dtype=np.float64)
            for m2 in range(self.L):
                for m1 in range(self.L):
                    s_A = self.site_index(m1, m2, 0)
                    s_B = self.site_index(m1 + d1, m2 + d2, 1)
                    u_acc += omega[self.majorana_slice(s_A), self.majorana_slice(s_B)]
            U_out[b_idx] = u_acc / float(num_cells)

        V_out = np.zeros((2, 4, 4), dtype=np.float64)
        for sub in (0, 1):
            v_acc = np.zeros((4, 4), dtype=np.float64)
            for m2 in range(self.L):
                for m1 in range(self.L):
                    s = self.site_index(m1, m2, sub)
                    v_acc += omega[self.majorana_slice(s), self.majorana_slice(s)]
            # Enforce exact skew-symmetry
            v_mean = v_acc / float(num_cells)
            V_out[sub] = 0.5 * (v_mean - v_mean.T)

        return MeanFieldState(
            U=U_out, V=V_out, lambdas=np.asarray(current_lambdas, dtype=np.float64)
        )

    def extract_inhomogeneous_state(
        self,
        omega: np.ndarray,
        current_lambdas: np.ndarray | None = None,
    ) -> InhomogeneousMeanFieldState:
        """Extract exact local bond (U) and on-site (V) correlators from Omega without spatial averaging.

        Args:
            omega: Skew-symmetric covariance matrix in R^{M x M}.
            current_lambdas: Optional local Lagrange multipliers of shape (num_cells, 2, 3) or (2, 3).
                             If None, defaults to zeros.

        Returns:
            InhomogeneousMeanFieldState retaining full point-wise defect and screening profiles.
        """
        disp = [(0, 0, 3), (0, -1, 1), (-1, 0, 2)]
        num_cells = self.num_cells

        U_out = np.zeros((num_cells, 3, 4, 4), dtype=np.float64)
        for m2 in range(self.L):
            for m1 in range(self.L):
                c = m1 + self.L * m2
                s_A = self.site_index(m1, m2, 0)
                sl_A = self.majorana_slice(s_A)
                for b_idx, (d1, d2, _) in enumerate(disp):
                    s_B = self.site_index(m1 + d1, m2 + d2, 1)
                    sl_B = self.majorana_slice(s_B)
                    U_out[c, b_idx] = omega[sl_A, sl_B]

        V_out = np.zeros((num_cells, 2, 4, 4), dtype=np.float64)
        for m2 in range(self.L):
            for m1 in range(self.L):
                c = m1 + self.L * m2
                for sub in (0, 1):
                    s = self.site_index(m1, m2, sub)
                    sl = self.majorana_slice(s)
                    v_local = omega[sl, sl]
                    V_out[c, sub] = 0.5 * (v_local - v_local.T)

        if current_lambdas is None:
            lam_out = np.zeros((num_cells, 2, 3), dtype=np.float64)
        else:
            lam_arr = np.asarray(current_lambdas, dtype=np.float64)
            if lam_arr.shape == (2, 3):
                lam_out = np.repeat(lam_arr[np.newaxis, ...], num_cells, axis=0)
            elif lam_arr.shape == (num_cells, 2, 3):
                lam_out = lam_arr.copy()
            else:
                raise ValueError(
                    f"lambdas must have shape ({num_cells}, 2, 3) or (2, 3), got {lam_arr.shape}"
                )

        return InhomogeneousMeanFieldState(U=U_out, V=V_out, lambdas=lam_out)

    def solve_real_space_purification(
        self,
        params: ModelParameters,
        initial_state: MeanFieldState,
        max_iter: int = 35,
        tol: float = 1e-8,
        accelerator: Literal["picard", "diis", "lbroyden"] = "diis",
        mixing_alpha: float = 0.5,
        lambda_step: float = 0.2,
        diis_history_depth: int = 6,
        lbroyden_history_size: int = 8,
        warm_start: bool = False,
        purify_tol: float = 1e-11,
        purify_max_iter: int = 35,
        flipped_bonds: Sequence[tuple[int, int, int]] | None = None,
        bond_disorder: dict[tuple[int, int, int], float] | None = None,
        site_vacancies: set[int] | None = None,
    ) -> tuple[np.ndarray, MeanFieldState, dict[str, Any]]:
        """Execute self-consistent mean-field iteration entirely via matrix sign purification.

        Computes the sequence of Hamiltonians and covariance matrices:
            A_k = assemble_skew_hamiltonian(params, state_k)
            Omega_{k+1} = purify_newton_schulz(A_k)
        without calculating any eigenvectors or invoking complex Hermitian routines.

        Args:
            params: Exchange couplings and Zeeman field.
            initial_state: Initial mean-field seed.
            max_iter: Maximum self-consistent outer iterations.
            tol: Outer convergence threshold on correlator updates and gauge residuals.
            accelerator: Fixed-point acceleration algorithm ('picard', 'diis', 'lbroyden').
            mixing_alpha: Linear mixing factor in (0, 1].
            lambda_step: Step length for dual gradient ascent on Lagrange multipliers.
            diis_history_depth: Subspace depth for DIIS.
            lbroyden_history_size: Memory history size for L-Broyden.
            warm_start: Whether to pass previous covariance into purification.
            purify_tol: Inner Newton-Schulz idempotency tolerance.
            purify_max_iter: Inner Newton-Schulz max iterations.
            flipped_bonds: Optional vortex string / branch cut.
            bond_disorder: Optional bond disorder map.
            site_vacancies: Optional site vacancy set.

        Returns:
            Tuple of (converged_omega, converged_state, metrics_dict).
        """
        curr_state = initial_state
        prev_omega: np.ndarray | None = None
        total_gemm = 0
        total_purify_iter = 0
        history_res: list[float] = []
        converged = False

        accel_engine: Any = None
        if accelerator == "diis":
            accel_engine = JointAndersonSolver(history_depth=diis_history_depth, beta=mixing_alpha)
        elif accelerator == "lbroyden":
            accel_engine = LimitedMemoryBroydenSolver(
                history_size=lbroyden_history_size, alpha=mixing_alpha
            )

        for _it in range(max_iter):
            # 1. Assemble real skew-symmetric Hamiltonian
            A = self.assemble_skew_hamiltonian(
                params,
                curr_state,
                flipped_bonds=flipped_bonds,
                bond_disorder=bond_disorder,
                site_vacancies=site_vacancies,
            )

            # 2. Evaluate ground-state covariance via cubic Newton-Schulz purification
            warm_arg = prev_omega if warm_start else None
            omega_next, p_met = purify_newton_schulz(
                A, warm_start=warm_arg, tol=purify_tol, max_iter=purify_max_iter
            )
            total_gemm += p_met.gemm_calls
            total_purify_iter += p_met.iterations
            prev_omega = omega_next

            # 3. Extract updated physical correlators
            state_next_raw = self.extract_mean_field_state(omega_next, curr_state.lambdas)

            # 4. Compute gauge residuals on on-site Majoranas
            g_res = np.zeros((2, 3), dtype=np.float64)
            for sub in (0, 1):
                v_s = state_next_raw.V[sub]
                for gamma in (1, 2, 3):
                    g_res[sub, gamma - 1] = float(np.trace(G_MATRICES[gamma] @ v_s))

            # 5. Dual gradient step on Lagrange multipliers
            lambda_target = curr_state.lambdas - lambda_step * g_res

            # 6. Check residuals
            diff_u = float(np.max(np.abs(state_next_raw.U - curr_state.U)))
            diff_v = float(np.max(np.abs(state_next_raw.V - curr_state.V)))
            max_g = float(np.max(np.abs(g_res)))
            total_res = max(diff_u, diff_v, max_g)
            history_res.append(total_res)

            if total_res < tol:
                converged = True
                curr_state = MeanFieldState(
                    U=state_next_raw.U,
                    V=state_next_raw.V,
                    lambdas=lambda_target,
                )
                break

            # 7. Advance state via chosen accelerator
            if accel_engine is not None:
                x_curr = np.concatenate(
                    [curr_state.U.ravel(), curr_state.V.ravel(), curr_state.lambdas.ravel()]
                )
                phi_curr = np.concatenate(
                    [state_next_raw.U.ravel(), state_next_raw.V.ravel(), lambda_target.ravel()]
                )
                x_next, _ = accel_engine.step(x_curr, phi_curr)

                len_u = curr_state.U.size
                len_v = curr_state.V.size
                u_n = x_next[:len_u].reshape(curr_state.U.shape)
                v_n = x_next[len_u : len_u + len_v].reshape(curr_state.V.shape)
                v_n = 0.5 * (v_n - np.swapaxes(v_n, -1, -2))
                lam_n = x_next[len_u + len_v :].reshape(curr_state.lambdas.shape)
                curr_state = MeanFieldState(U=u_n, V=v_n, lambdas=lam_n)
            else:
                # Damped Picard
                u_mixed = (1.0 - mixing_alpha) * curr_state.U + mixing_alpha * state_next_raw.U
                v_mixed = (1.0 - mixing_alpha) * curr_state.V + mixing_alpha * state_next_raw.V
                v_mixed = 0.5 * (v_mixed - np.swapaxes(v_mixed, -1, -2))
                curr_state = MeanFieldState(U=u_mixed, V=v_mixed, lambdas=lambda_target)

        assert prev_omega is not None
        metrics: dict[str, Any] = {
            "converged": converged,
            "iterations": len(history_res),
            "final_residual": history_res[-1] if history_res else 1.0,
            "total_gemm_calls": total_gemm,
            "total_purify_iterations": total_purify_iter,
            "residual_history": history_res,
        }
        return prev_omega, curr_state, metrics

    def solve_inhomogeneous_purification(
        self,
        params: ModelParameters,
        initial_state: InhomogeneousMeanFieldState,
        max_iter: int = 35,
        tol: float = 1e-8,
        accelerator: Literal["picard", "diis", "lbroyden"] = "diis",
        mixing_alpha: float = 0.5,
        lambda_step: float = 0.2,
        diis_history_depth: int = 6,
        lbroyden_history_size: int = 8,
        warm_start: bool = False,
        purify_tol: float = 1e-11,
        purify_max_iter: int = 35,
        flipped_bonds: Sequence[tuple[int, int, int]] | None = None,
        bond_disorder: dict[tuple[int, int, int], float] | None = None,
        site_vacancies: set[int] | None = None,
    ) -> tuple[np.ndarray, InhomogeneousMeanFieldState, dict[str, Any]]:
        """Execute self-consistent inhomogeneous mean-field iteration entirely via matrix sign purification.

        Maintains the full D = 86 * L^2 parameter space without spatial averaging, allowing
        unpinned vortices, branch cuts, and local disorder configurations to relax faithfully.

        Args:
            params: Exchange couplings and external Zeeman field.
            initial_state: Initial inhomogeneous mean-field seed.
            max_iter: Maximum self-consistent outer iterations.
            tol: Outer convergence threshold on correlator updates and gauge residuals.
            accelerator: Fixed-point acceleration algorithm ('picard', 'diis', 'lbroyden').
            mixing_alpha: Linear mixing factor in (0, 1].
            lambda_step: Step length for dual gradient ascent on Lagrange multipliers.
            diis_history_depth: Subspace depth for DIIS.
            lbroyden_history_size: Memory history size for L-Broyden.
            warm_start: Whether to pass previous covariance into purification.
            purify_tol: Inner Newton-Schulz idempotency tolerance.
            purify_max_iter: Inner Newton-Schulz max iterations.
            flipped_bonds: Optional vortex string / branch cut.
            bond_disorder: Optional bond disorder map.
            site_vacancies: Optional site vacancy set.

        Returns:
            Tuple of (converged_omega, converged_inhomogeneous_state, metrics_dict).
        """
        curr_state = initial_state
        prev_omega: np.ndarray | None = None
        total_gemm = 0
        total_purify_iter = 0
        history_res: list[float] = []
        converged = False

        accel_engine: Any = None
        if accelerator == "diis":
            accel_engine = JointAndersonSolver(history_depth=diis_history_depth, beta=mixing_alpha)
        elif accelerator == "lbroyden":
            accel_engine = LimitedMemoryBroydenSolver(
                history_size=lbroyden_history_size, alpha=mixing_alpha
            )

        for _it in range(max_iter):
            # 1. Assemble inhomogeneous real skew-symmetric Hamiltonian
            A = self.assemble_skew_hamiltonian(
                params,
                curr_state,
                flipped_bonds=flipped_bonds,
                bond_disorder=bond_disorder,
                site_vacancies=site_vacancies,
            )

            # 2. Evaluate ground-state covariance via cubic Newton-Schulz purification
            warm_arg = prev_omega if warm_start else None
            omega_next, p_met = purify_newton_schulz(
                A, warm_start=warm_arg, tol=purify_tol, max_iter=purify_max_iter
            )
            total_gemm += p_met.gemm_calls
            total_purify_iter += p_met.iterations
            prev_omega = omega_next

            # 3. Extract updated physical correlators without spatial averaging
            state_next_raw = self.extract_inhomogeneous_state(omega_next, curr_state.lambdas)

            # 4. Compute gauge residuals on on-site Majoranas across all cells
            g_res = np.zeros((self.num_cells, 2, 3), dtype=np.float64)
            for c in range(self.num_cells):
                for sub in (0, 1):
                    v_s = state_next_raw.V[c, sub]
                    for gamma in (1, 2, 3):
                        g_res[c, sub, gamma - 1] = float(np.trace(G_MATRICES[gamma] @ v_s))

            # 5. Dual gradient step on Lagrange multipliers
            lambda_target = curr_state.lambdas - lambda_step * g_res

            # 6. Check residuals
            diff_u = float(np.max(np.abs(state_next_raw.U - curr_state.U)))
            diff_v = float(np.max(np.abs(state_next_raw.V - curr_state.V)))
            max_g = float(np.max(np.abs(g_res)))
            total_res = max(diff_u, diff_v, max_g)
            history_res.append(total_res)

            if total_res < tol:
                converged = True
                curr_state = InhomogeneousMeanFieldState(
                    U=state_next_raw.U,
                    V=state_next_raw.V,
                    lambdas=lambda_target,
                )
                break

            # 7. Advance state via chosen accelerator
            if accel_engine is not None:
                x_curr = curr_state.flatten()
                phi_curr = InhomogeneousMeanFieldState(
                    U=state_next_raw.U,
                    V=state_next_raw.V,
                    lambdas=lambda_target,
                ).flatten()
                x_next, _ = accel_engine.step(x_curr, phi_curr)
                curr_state = InhomogeneousMeanFieldState.from_flattened(x_next, self.num_cells)
            else:
                # Damped Picard
                u_mixed = (1.0 - mixing_alpha) * curr_state.U + mixing_alpha * state_next_raw.U
                v_mixed = (1.0 - mixing_alpha) * curr_state.V + mixing_alpha * state_next_raw.V
                v_mixed = 0.5 * (v_mixed - np.swapaxes(v_mixed, -1, -2))
                curr_state = InhomogeneousMeanFieldState(
                    U=u_mixed, V=v_mixed, lambdas=lambda_target
                )

        assert prev_omega is not None
        metrics: dict[str, Any] = {
            "converged": converged,
            "iterations": len(history_res),
            "final_residual": history_res[-1] if history_res else 1.0,
            "total_gemm_calls": total_gemm,
            "total_purify_iterations": total_purify_iter,
            "residual_history": history_res,
        }
        return prev_omega, curr_state, metrics
