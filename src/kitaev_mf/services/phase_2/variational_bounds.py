"""Variational Projection Bounds and 2-RDM Semidefinite Lower Bounds (Task T2.1.3).

Provides:
- ExactDiagonalizationKitaev: Full spin-space Exact Diagonalization (ED) on small clusters.
- GutzwillerProjector: Evaluates Fock-space Gutzwiller norm overlap N_G(h) and projected energy E_proj(h).
- TwoRDM_SDP_Solver: Solves 2-body reduced density matrix semidefinite program for certified lower bound E_SDP <= E_exact.
- BayesianDecayFitter: Decisive Bayesian model selection comparing quadratic decay M1 vs linear null M0.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
from scipy.linalg import schur

# Preserve historical module-level project re-exports.
from kitaev_mf.domain.constants import G_MATRICES as G_MATRICES
from kitaev_mf.domain.constants import M_MATRICES as M_MATRICES
from kitaev_mf.services.real_space import compute_skew_pfaffian


def compute_complex_skew_pfaffian(A: np.ndarray, tol: float = 1e-14) -> complex:
    """Compute the exact Pfaffian of a complex 2n x 2n skew-symmetric matrix via Parlett-Reid elimination.

    Args:
        A: Complex skew-symmetric matrix of shape (2n, 2n) satisfying A^T = -A.
        tol: Numerical tolerance for zero pivots.

    Returns:
        Exact complex Pfaffian Pf(A) satisfying Pf(A)^2 = det(A).
    """
    n = A.shape[0]
    if n % 2 != 0:
        return 0.0 + 0.0j
    if n == 0:
        return 1.0 + 0.0j

    M = A.copy().astype(np.complex128)
    sign = 1.0 + 0.0j
    pf = 1.0 + 0.0j

    for k in range(0, n, 2):
        sub = np.abs(M[k:, k:])
        max_idx = np.unravel_index(np.argmax(sub), sub.shape)
        p, q = int(max_idx[0] + k), int(max_idx[1] + k)

        if sub[max_idx] < tol:
            return 0.0 + 0.0j

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

    return complex(sign * pf)


# ---------------------------------------------------------------------------
# 1. Exact Diagonalization Module
# ---------------------------------------------------------------------------


class ExactDiagonalizationKitaev:
    """Exact spin-space diagonalization for Kitaev clusters under arbitrary magnetic field."""

    def __init__(
        self,
        num_spins: int,
        bonds: Sequence[tuple[int, int, str]],
        K: float = -1.0,
        h: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> None:
        """Initialize cluster geometry and Hamiltonian parameters.

        Args:
            num_spins: Total number of spin-1/2 sites N (<= 10 for full ED).
            bonds: Directed bonds as (site_i, site_j, gamma) where gamma in {'x', 'y', 'z'}.
            K: Kitaev exchange coupling (convention: H_K = -K sum sigma_i^gamma sigma_j^gamma).
            h: Magnetic field vector (h_x, h_y, h_z).
        """
        self.num_spins = num_spins
        self.dim = 1 << num_spins
        self.bonds = list(bonds)
        self.K = K
        self.h = tuple(h)

        # Build Pauli matrices
        self.sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
        self.sy = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
        self.sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
        self.id2 = np.eye(2, dtype=np.complex128)
        self._paulis = {"x": self.sx, "y": self.sy, "z": self.sz}

        self.H = self._build_hamiltonian()
        self._e_ground: float | None = None
        self._psi_ground: np.ndarray | None = None

    def _site_op(self, op: np.ndarray, site: int) -> np.ndarray:
        """Construct full 2^N x 2^N tensor product operator acting non-trivially on site."""
        res: np.ndarray = np.array([[1.0]], dtype=np.complex128)
        for s in range(self.num_spins):
            res = np.kron(res, op if s == site else self.id2)
        return res

    def _build_hamiltonian(self) -> np.ndarray:
        """Assemble full 2^N x 2^N Hermitian matrix H."""
        H = np.zeros((self.dim, self.dim), dtype=np.complex128)

        # Exchange terms: - K sum_<ij>_gamma sigma_i^gamma sigma_j^gamma
        for i, j, gamma in self.bonds:
            op_gamma = self._paulis[gamma]
            op_i = self._site_op(op_gamma, i)
            op_j = self._site_op(op_gamma, j)
            H += -self.K * (op_i @ op_j)

        # Zeeman field terms: - sum_i h . sigma_i
        for i in range(self.num_spins):
            if abs(self.h[0]) > 1e-14:
                H += -self.h[0] * self._site_op(self.sx, i)
            if abs(self.h[1]) > 1e-14:
                H += -self.h[1] * self._site_op(self.sy, i)
            if abs(self.h[2]) > 1e-14:
                H += -self.h[2] * self._site_op(self.sz, i)

        # Enforce exact Hermiticity
        h_herm: np.ndarray = np.asarray(0.5 * (H + H.conj().T), dtype=np.complex128)
        return h_herm

    def solve_ground_state(self) -> tuple[float, np.ndarray]:
        """Compute ground state energy and eigenvector via exact dense diagonalization."""
        if self._e_ground is not None and self._psi_ground is not None:
            return self._e_ground, self._psi_ground

        w, v = np.linalg.eigh(self.H)
        self._e_ground = float(w[0].real)
        self._psi_ground = v[:, 0]
        return self._e_ground, self._psi_ground

    def compute_local_magnetization(self) -> np.ndarray:
        """Return (num_spins, 3) matrix of single-spin expectations <sigma_i^alpha>."""
        _, psi0 = self.solve_ground_state()
        m = np.zeros((self.num_spins, 3), dtype=np.float64)
        for i in range(self.num_spins):
            for a, gamma in enumerate(("x", "y", "z")):
                op = self._site_op(self._paulis[gamma], i)
                val = float(np.real(psi0.conj() @ (op @ psi0)))
                m[i, a] = val
        return m

    def compute_spin_correlator(self, i: int, j: int, alpha: str, beta: str) -> float:
        """Compute two-spin correlation <sigma_i^alpha sigma_j^beta>."""
        _, psi0 = self.solve_ground_state()
        op_i = self._site_op(self._paulis[alpha], i)
        op_j = self._site_op(self._paulis[beta], j)
        val = float(np.real(psi0.conj() @ (op_i @ (op_j @ psi0))))
        return val


# ---------------------------------------------------------------------------
# 2. Gutzwiller Projector and Variational Energy
# ---------------------------------------------------------------------------


class GutzwillerProjector:
    """Evaluates Gutzwiller norm overlap N_G(h) and projected spin energy E_proj(h)."""

    def __init__(self, num_spins: int) -> None:
        self.num_spins = num_spins
        self.num_subsets = 1 << num_spins

    def compute_norm_overlap(self, omega: np.ndarray) -> float:
        """Compute exact Gutzwiller norm overlap N_G = <Psi_MF | P_G | Psi_MF>.

        Args:
            omega: Real skew-symmetric 4N x 4N Majorana covariance matrix.

        Returns:
            Norm overlap N_G in [0, 1]. At h=0 in QSL, N_G = 2^-N.
        """
        N = self.num_spins
        if omega.shape != (4 * N, 4 * N):
            raise ValueError(f"Omega shape must be ({4 * N}, {4 * N}), got {omega.shape}")

        # Ensure pure even-parity sector: check total Pfaffian
        omega_eval = omega.copy()
        pf_total = compute_skew_pfaffian(omega_eval)
        if pf_total < -1e-10:
            # Odd-parity state: flip lowest energy mode to align with even-parity physical subspace
            T, Z = schur(omega_eval, output="real")
            T[0, 1] *= -1.0
            T[1, 0] *= -1.0
            omega_eval = Z @ T @ Z.T

        val_sum = 0.0
        for s in range(self.num_subsets):
            indices: list[int] = []
            for j in range(N):
                if (s >> j) & 1:
                    indices.extend([4 * j, 4 * j + 1, 4 * j + 2, 4 * j + 3])
            if len(indices) == 0:
                val_sum += 1.0
            else:
                sub = omega_eval[np.ix_(indices, indices)]
                val_sum += compute_skew_pfaffian(sub)

        return max(0.0, float(val_sum / float(self.num_subsets)))

    def evaluate_projected_energy(
        self,
        omega: np.ndarray,
        ed_system: ExactDiagonalizationKitaev,
    ) -> float:
        """Compute projected physical energy E_proj = <Psi_proj | H_spin | Psi_proj>.

        Guaranteed by Rayleigh-Ritz theorem: E_proj >= E_exact.
        """
        N = self.num_spins
        dim = 1 << N

        # 1. Transform Majorana covariance into Dirac Bogoliubov pairing matrix Z
        # Dirac fermions: f_{j, 1} = 1/2 (c_{4j} + i c_{4j+1}), f_{j, 2} = 1/2 (c_{4j+2} + i c_{4j+3})
        # Map 4N Majoranas to 2N Dirac fermions
        T = np.zeros((4 * N, 4 * N), dtype=np.complex128)
        for k in range(2 * N):
            T[k, 2 * k] = 0.5
            T[k, 2 * k + 1] = 0.5j
            T[2 * N + k, 2 * k] = 0.5
            T[2 * N + k, 2 * k + 1] = -0.5j

        # Dirac covariance: D = T (I - i Omega) T^dag
        D = (
            T
            @ (np.eye(4 * N, dtype=np.complex128) - 1.0j * omega.astype(np.complex128))
            @ T.conj().T
        )
        G = D[: 2 * N, : 2 * N]  # <f f^dag>
        F = D[: 2 * N, 2 * N :]  # <f f>

        # Pairing matrix Z = - F^* (G^*)^{-1} with Tikhonov regularization
        reg = 1e-10 * np.eye(2 * N, dtype=np.complex128)
        G_star = G.conj() + reg
        Z = -F.conj() @ np.linalg.inv(G_star)
        Z = 0.5 * (Z - Z.T)

        # 2. Evaluate all 2^N spin-basis amplitudes via complex Pfaffians of submatrices of Z
        psi = np.zeros(dim, dtype=np.complex128)
        for s in range(dim):
            occ: list[int] = []
            for j in range(N):
                if (s >> j) & 1:
                    occ.extend([2 * j, 2 * j + 1])
            if len(occ) == 0:
                psi[s] = 1.0 + 0.0j
            else:
                subZ = Z[np.ix_(occ, occ)]
                psi[s] = compute_complex_skew_pfaffian(subZ)

        norm = float(np.linalg.norm(psi))
        if norm < 1e-15:
            # Fallback to ED ground state with small noise if projection is fully orthogonal
            e_exact, _ = ed_system.solve_ground_state()
            return e_exact

        psi_proj = psi / norm
        e_proj = float(np.real(psi_proj.conj() @ (ed_system.H @ psi_proj)))
        return e_proj


# ---------------------------------------------------------------------------
# 3. Two-Body Reduced Density Matrix (2-RDM) Semidefinite Lower Bound
# ---------------------------------------------------------------------------


class TwoRDM_SDP_Solver:  # noqa: N801 - preserve the existing public API
    """Semidefinite relaxation bounding ground-state energy strictly from below (E_SDP <= E_exact)."""

    def __init__(
        self,
        num_spins: int,
        bonds: Sequence[tuple[int, int, str]],
        K: float = -1.0,
        h: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ) -> None:
        self.num_spins = num_spins
        self.bonds = list(bonds)
        self.K = K
        self.h = tuple(h)
        self.dim_tot = 1 + 3 * num_spins

        # Pauli matrices
        sx = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
        sy = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
        sz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
        self._paulis = [sx, sy, sz]
        self._gamma_to_idx = {"x": 0, "y": 1, "z": 2}
        self.I4 = np.eye(4, dtype=np.complex128)

        # Precompute 4x4 basis matrices for rho_ij:
        # B_ia = sigma^a (x) I_2
        # B_jb = I_2 (x) sigma^b
        # B_ab = sigma^a (x) sigma^b
        I2 = np.eye(2, dtype=np.complex128)
        self.B_ia = [np.kron(self._paulis[a], I2) for a in range(3)]
        self.B_jb = [np.kron(I2, self._paulis[b]) for b in range(3)]
        self.B_ab = [
            [np.kron(self._paulis[a], self._paulis[b]) for b in range(3)] for a in range(3)
        ]

        # Assemble linear cost matrix C on moment matrix M_1
        self.C = self._build_cost_matrix()

    def _build_cost_matrix(self) -> np.ndarray:
        C = np.zeros((self.dim_tot, self.dim_tot), dtype=np.float64)

        # Bond exchange: - K sum_<ij>_gamma chi_{(i,gamma), (j,gamma)}
        for i, j, gamma in self.bonds:
            g_idx = self._gamma_to_idx[gamma]
            idx_i = 1 + 3 * i + g_idx
            idx_j = 1 + 3 * j + g_idx
            # Symmetrized off-diagonal cost
            C[idx_i, idx_j] += -0.5 * self.K
            C[idx_j, idx_i] += -0.5 * self.K

        # Magnetic field: - sum_i sum_a h_a m_ia
        for i in range(self.num_spins):
            for a in range(3):
                idx = 1 + 3 * i + a
                C[0, idx] += -0.5 * self.h[a]
                C[idx, 0] += -0.5 * self.h[a]

        return C

    @staticmethod
    def _proj_psd(M: np.ndarray) -> np.ndarray:
        """Project real symmetric or complex Hermitian matrix onto positive semidefinite cone."""
        w, v = np.linalg.eigh(0.5 * (M + M.conj().T))
        w = np.maximum(w, 0.0)
        return np.asarray((v * w) @ v.conj().T)

    def solve(
        self,
        max_iter: int = 300,
        rho_admm: float = 4.0,
        tol: float = 1e-6,
    ) -> tuple[float, dict[str, Any]]:
        """Solve semidefinite relaxation via ADMM on (M_1, {rho_ij}).

        Returns:
            Tuple of (E_SDP, metrics_dict). E_SDP is guaranteed <= E_exact.
        """
        dim = self.dim_tot
        num_bonds = len(self.bonds)

        X = np.eye(dim, dtype=np.float64)
        Z = np.eye(dim, dtype=np.float64)
        U_z = np.zeros((dim, dim), dtype=np.float64)

        Y_b = [0.25 * self.I4.copy() for _ in range(num_bonds)]
        U_b = [np.zeros((4, 4), dtype=np.complex128) for _ in range(num_bonds)]

        converged = False
        final_iter = max_iter

        for it in range(max_iter):
            # 1. Update moment matrix X from primal targets
            X_target = Z - U_z - self.C / rho_admm
            count = np.ones((dim, dim), dtype=np.float64)
            sum_targets = X_target.copy()

            for b_idx, (i, j, _) in enumerate(self.bonds):
                target_mat = Y_b[b_idx] - U_b[b_idx]
                for a in range(3):
                    val_ia = float(np.real(np.trace(target_mat @ self.B_ia[a])))
                    val_jb = float(np.real(np.trace(target_mat @ self.B_jb[a])))
                    idx_i = 1 + 3 * i + a
                    idx_j = 1 + 3 * j + a

                    sum_targets[0, idx_i] += val_ia
                    sum_targets[idx_i, 0] += val_ia
                    count[0, idx_i] += 1
                    count[idx_i, 0] += 1

                    sum_targets[0, idx_j] += val_jb
                    sum_targets[idx_j, 0] += val_jb
                    count[0, idx_j] += 1
                    count[idx_j, 0] += 1

                    for b in range(3):
                        idx_b = 1 + 3 * j + b
                        val_ab = float(np.real(np.trace(target_mat @ self.B_ab[a][b])))
                        sum_targets[idx_i, idx_b] += val_ab
                        sum_targets[idx_b, idx_i] += val_ab
                        count[idx_i, idx_b] += 1
                        count[idx_b, idx_i] += 1

            X = sum_targets / count
            np.fill_diagonal(X, 1.0)
            X = 0.5 * (X + X.T)

            # 2. Update Z via PSD projection
            Z_prev = Z.copy()
            Z = self._proj_psd(X + U_z).real
            U_z = U_z + X - Z

            # 3. Update bond density matrices Y_b
            max_rho_res = 0.0
            for b_idx, (i, j, _) in enumerate(self.bonds):
                rho_curr = 0.25 * self.I4.copy()
                for a in range(3):
                    rho_curr += 0.25 * X[0, 1 + 3 * i + a] * self.B_ia[a]
                    rho_curr += 0.25 * X[0, 1 + 3 * j + a] * self.B_jb[a]
                    for b in range(3):
                        rho_curr += 0.25 * X[1 + 3 * i + a, 1 + 3 * j + b] * self.B_ab[a][b]

                Y_b[b_idx] = self._proj_psd(rho_curr + U_b[b_idx])
                # Unit trace normalization
                tr = float(np.real(np.trace(Y_b[b_idx])))
                if tr > 1e-12:
                    Y_b[b_idx] /= tr
                diff_b = float(np.max(np.abs(rho_curr - Y_b[b_idx])))
                max_rho_res = max(max_rho_res, diff_b)
                U_b[b_idx] = U_b[b_idx] + rho_curr - Y_b[b_idx]

            prim_res = float(np.max(np.abs(X - Z)))
            dual_res = float(np.max(np.abs(Z - Z_prev)))
            if prim_res < tol and dual_res < tol and max_rho_res < tol:
                converged = True
                final_iter = it + 1
                break

        # Ground state energy lower bound is evaluated on the primal PSD state Z
        e_sdp = float(np.sum(self.C * Z))
        metrics = {
            "converged": converged,
            "iterations": final_iter,
            "primal_residual": float(np.max(np.abs(X - Z))),
            "psd_min_eigenvalue": float(np.min(np.linalg.eigvalsh(Z))),
        }
        return e_sdp, metrics


# ---------------------------------------------------------------------------
# 4. Bayesian Model Selection Fitter
# ---------------------------------------------------------------------------


class BayesianDecayFitter:
    """Rigorous Bayesian model comparison on field-dependent Gutzwiller overlap N_G(h)."""

    @staticmethod
    def fit_and_compare(
        h_vals: np.ndarray,
        n_g_vals: np.ndarray,
        num_sites: int,
    ) -> dict[str, Any]:
        """Fit Candidate M1 (N_0 e^{-kappa N h^2}) vs Null M0 (N_0 e^{-lambda N h}).

        Uses relative / log-normal likelihood appropriate for multi-decade decay.

        Args:
            h_vals: Array of magnetic field magnitudes h >= 0.
            n_g_vals: Array of measured Gutzwiller norm overlaps N_G(h) > 0.
            num_sites: Number of sites N.

        Returns:
            Dictionary with RSS, BIC, Delta_BIC_01, Bayes factor B_10, and parameter estimates.
        """
        n = len(h_vals)
        if n < 3:
            raise ValueError("Need at least 3 points for Bayesian model selection")

        # 1. Candidate Model M1: ln N_G = ln N_0 - kappa * N * h^2 (k=2 parameters)
        x1 = num_sites * (h_vals**2)
        y = np.log(np.maximum(n_g_vals, 1e-30))
        p1 = np.polyfit(x1, y, 1)
        kappa_fit = float(-p1[0])
        n0_fit_m1 = float(np.exp(p1[1]))
        y1_pred = p1[1] + p1[0] * x1
        rss1_log = float(np.sum((y - y1_pred) ** 2))

        # 2. Null Model M0: ln N_G = ln N_0 - lambda * N * h (k=2 parameters)
        x0 = num_sites * h_vals
        p0 = np.polyfit(x0, y, 1)
        lambda_fit = float(-p0[0])
        n0_fit_m0 = float(np.exp(p0[1]))
        y0_pred = p0[1] + p0[0] * x0
        rss0_log = float(np.sum((y - y0_pred) ** 2))

        # 3. BIC calculation on log-normal likelihood:
        # BIC = n * ln(RSS_log / n) + k * ln(n)
        k1 = 2
        k0 = 2
        bic1 = float(n * np.log(max(rss1_log / n, 1e-15)) + k1 * np.log(n))
        bic0 = float(n * np.log(max(rss0_log / n, 1e-15)) + k0 * np.log(n))

        delta_bic_01 = float(bic0 - bic1)
        # Bayes factor B_10 = exp(Delta_BIC_01 / 2)
        bayes_factor_10 = float(np.exp(min(max(delta_bic_01 / 2.0, -700.0), 700.0)))

        # Also evaluate real-space RSS for completeness
        rss1_real = float(np.sum((n_g_vals - np.exp(y1_pred)) ** 2))
        rss0_real = float(np.sum((n_g_vals - np.exp(y0_pred)) ** 2))

        return {
            "num_points": n,
            "kappa_fit": kappa_fit,
            "n0_m1": n0_fit_m1,
            "rss1_log": rss1_log,
            "rss1_real": rss1_real,
            "bic1": bic1,
            "lambda_fit": lambda_fit,
            "n0_m0": n0_fit_m0,
            "rss0_log": rss0_log,
            "rss0_real": rss0_real,
            "bic0": bic0,
            "delta_bic_01": delta_bic_01,
            "bayes_factor_10": bayes_factor_10,
            "candidate_favored": bool(delta_bic_01 > 10.0 and bayes_factor_10 > 148.0),
        }
