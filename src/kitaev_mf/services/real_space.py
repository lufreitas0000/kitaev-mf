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

from kitaev_mf.domain.constants import G_MATRICES
from kitaev_mf.domain.lattice import DEFAULT_D_NN, make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.fixed_point import JointAndersonSolver, LimitedMemoryBroydenSolver
from kitaev_mf.services.hamiltonian import (
    build_bond_hopping_matrices,
    build_on_site_matrices,
)
from kitaev_mf.services.purification import purify_newton_schulz


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
        flipped_bonds: Sequence[tuple[int, int, int]] | set[tuple[int, int, int]] | None = None,
        bond_disorder: dict[tuple[int, int, int], float] | None = None,
        site_vacancies: set[int] | Sequence[int] | None = None,
    ) -> np.ndarray:
        """Assemble the exact real skew-symmetric M x M matrix A where H = (i/4) c^T A c.

        Args:
            params: Exchange couplings and external field.
            state: Mean-field correlation state.
            flipped_bonds: Optional collection of directed bond keys (m1, m2, b_idx)
                           whose gauge link signs are inverted (u_b -> -u_b).
            bond_disorder: Optional mapping from bond key (m1, m2, b_idx) to relative
                           exchange modulation factor delta_b such that K_b = K * (1 + delta_b).
            site_vacancies: Optional collection of physical site indices to remove
                            (dilute non-magnetic vacancies with zeroed couplings).
        """
        cell2 = make_honeycomb_2site(self.d_nn)
        device = torch.device("cpu")

        # Use canonical production services for on-site B_s and hopping A_b
        b_mat = build_on_site_matrices(cell2, params, state, device).cpu().numpy()
        a_mat = build_bond_hopping_matrices(cell2, params, state, device).cpu().numpy()

        flipped_set = set(flipped_bonds) if flipped_bonds is not None else set()
        vacancies_set = set(site_vacancies) if site_vacancies is not None else set()
        disorder_map = bond_disorder if bond_disorder is not None else {}

        M = self.num_majoranas
        A = np.zeros((M, M), dtype=np.float64)

        # 1. On-site skew-symmetric blocks
        for s in range(self.num_sites):
            if s in vacancies_set:
                continue
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

                    if s_A in vacancies_set or s_B in vacancies_set:
                        continue

                    bond_key = (m1, m2, b_idx)
                    sign = -1.0 if bond_key in flipped_set else 1.0
                    disorder_factor = 1.0 + disorder_map.get(bond_key, 0.0)

                    # In Bloch Hamiltonian, block (s_A, s_B) is +(i/4) * A_b
                    # In real-space Majorana form H = (i/4) c^T A c:
                    # c_A^T A_b c_B => in skew-symmetric A:
                    # A[A, B] = A_b, A[B, A] = -A_b^T
                    A[sl_A, sl_B] += sign * disorder_factor * a_mat[b_idx]
                    A[sl_B, sl_A] -= sign * disorder_factor * a_mat[b_idx].T

        # Enforce exact skew-symmetry
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
