"""Model-agnostic fixed-point accelerators and quasi-Newton root-finding algorithms.

This module operates strictly on Euclidean state vectors x in R^D and fixed-point maps
Phi: R^D -> R^D, with residual G(x) = Phi(x) - x. It is independent of specific lattice,
Majorana, or Hamiltonian abstractions.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Protocol

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres


class FixedPointAccelerator(Protocol):
    """Protocol defining the interface for fixed-point iteration accelerators."""

    def reset(self) -> None:
        """Clear historical subspace and state memory."""
        ...

    def step(
        self,
        x: np.ndarray,
        phi_x: np.ndarray,
        residual_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    ) -> tuple[np.ndarray, int]:
        """Propose the next state iterate x_{k+1}.

        Parameters:
            x: Current state vector x_k of shape (D,), float64.
            phi_x: Evaluated map output Phi(x_k) of shape (D,), float64.
            residual_fn: Optional callable computing G(x') = Phi(x') - x' for extra map
                evaluations (used by Jacobian-Free Newton-Krylov).

        Returns:
            Tuple of (x_next, extra_phi_evaluations):
                x_next: Proposed next iterate x_{k+1} of shape (D,).
                extra_phi_evaluations: Number of additional Phi evaluations incurred (0 for FP/DIIS/Broyden).
        """
        ...


class PlainFixedPoint(FixedPointAccelerator):
    """Unaccelerated standard fixed-point iteration: x_{k+1} = Phi(x_k)."""

    def reset(self) -> None:
        """No internal state to clear."""

    def step(
        self,
        x: np.ndarray,
        phi_x: np.ndarray,
        residual_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    ) -> tuple[np.ndarray, int]:
        """Advance one step via direct fixed-point mapping."""
        return np.asarray(phi_x, dtype=np.float64).copy(), 0


class DampedPicard(FixedPointAccelerator):
    """Damped linear mixing fixed-point iteration: x_{k+1} = (1 - alpha) x_k + alpha Phi(x_k)."""

    def __init__(self, alpha: float = 0.5) -> None:
        if not (0.0 < alpha <= 1.0):
            raise ValueError(f"Damping alpha must be in (0, 1], got {alpha}")
        self.alpha = float(alpha)

    def reset(self) -> None:
        """No internal state to clear."""

    def step(
        self,
        x: np.ndarray,
        phi_x: np.ndarray,
        residual_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    ) -> tuple[np.ndarray, int]:
        """Advance one step via convex combination."""
        x_next = (1.0 - self.alpha) * x + self.alpha * phi_x
        return np.asarray(x_next, dtype=np.float64), 0


class JointAndersonSolver(FixedPointAccelerator):
    """Joint Anderson Acceleration (DIIS) minimizing residual norm over historical subspace."""

    def __init__(
        self,
        history_depth: int = 6,
        reg: float = 1e-12,
        beta: float = 1.0,
    ) -> None:
        if history_depth < 1:
            raise ValueError(f"history_depth must be at least 1, got {history_depth}")
        if reg < 0.0:
            raise ValueError(f"reg must be non-negative, got {reg}")
        self.history_depth = history_depth
        self.reg = reg
        self.beta = beta
        self.reset()

    def reset(self) -> None:
        """Clear historical subspace buffers."""
        self._history_x: list[np.ndarray] = []
        self._history_r: list[np.ndarray] = []

    def step(
        self,
        x: np.ndarray,
        phi_x: np.ndarray,
        residual_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    ) -> tuple[np.ndarray, int]:
        """Perform DIIS extrapolation over retained subspace."""
        residual = phi_x - x

        self._history_x.append(np.asarray(x, dtype=np.float64).copy())
        self._history_r.append(np.asarray(residual, dtype=np.float64).copy())

        if len(self._history_x) > self.history_depth:
            self._history_x.pop(0)
            self._history_r.pop(0)

        m = len(self._history_x)
        if m == 1:
            x_next = x + self.beta * residual
            return np.asarray(x_next, dtype=np.float64), 0

        # Build Gram matrix Omega_{pq} = <R_p, R_q>
        omega = np.zeros((m, m), dtype=np.float64)
        for p in range(m):
            for q in range(p, m):
                dot_val = float(np.dot(self._history_r[p], self._history_r[q]))
                omega[p, q] = dot_val
                omega[q, p] = dot_val

        diag_scale = np.diag(omega)
        for p in range(m):
            omega[p, p] += self.reg * max(float(diag_scale[p]), 1e-14)

        # Assemble (m + 1) x (m + 1) KKT system: [[Omega, -1], [-1^T, 0]] [c, nu]^T = [0, -1]^T
        kkt = np.zeros((m + 1, m + 1), dtype=np.float64)
        kkt[:m, :m] = omega
        kkt[:m, m] = -1.0
        kkt[m, :m] = -1.0

        rhs = np.zeros(m + 1, dtype=np.float64)
        rhs[m] = -1.0

        try:
            sol = np.linalg.solve(kkt, rhs)
            weights = sol[:m]
        except np.linalg.LinAlgError:
            # Fallback to standard step if KKT matrix is singular
            x_fallback = x + self.beta * residual
            return np.asarray(x_fallback, dtype=np.float64), 0

        # Construct accelerated state: X_next = sum_j c_j * (X_j + beta * R_j)
        x_next = np.zeros_like(x, dtype=np.float64)
        for j in range(m):
            x_next += weights[j] * (self._history_x[j] + self.beta * self._history_r[j])

        return np.asarray(x_next, dtype=np.float64), 0


class BroydenSolver(FixedPointAccelerator):
    """Good Broyden rank-1 quasi-Newton update on the inverse Jacobian approximation H_k."""

    def __init__(self, alpha: float = 0.5, tol_singular: float = 1e-12) -> None:
        if not (0.0 < alpha <= 1.0):
            raise ValueError(f"alpha must be in (0, 1], got {alpha}")
        self.alpha = float(alpha)
        self.tol_singular = float(tol_singular)
        self.reset()

    def reset(self) -> None:
        """Clear inverse Jacobian approximation and trajectory memory."""
        self._H: np.ndarray | None = None
        self._prev_x: np.ndarray | None = None
        self._prev_g: np.ndarray | None = None

    def step(
        self,
        x: np.ndarray,
        phi_x: np.ndarray,
        residual_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    ) -> tuple[np.ndarray, int]:
        """Update inverse Jacobian H_k via Sherman-Morrison secant formula and propose x_{k+1}."""
        g = phi_x - x
        dim = x.size

        if self._H is None:
            # Initialize inverse Jacobian approximation as scaled negative identity: H_0 = -alpha * I
            # Because G(x) = Phi(x) - x, at a root G(x*) = 0, Newton's step is x_{k+1} = x_k - J^{-1} G_k
            # If J approx -I, then -J^{-1} G_k approx G_k, yielding x_{k+1} = x_k + alpha * G_k (damped Picard)
            self._H = -self.alpha * np.eye(dim, dtype=np.float64)
            x_next = x - self._H @ g
            self._prev_x = np.asarray(x, dtype=np.float64).copy()
            self._prev_g = np.asarray(g, dtype=np.float64).copy()
            return np.asarray(x_next, dtype=np.float64), 0

        assert self._prev_x is not None
        assert self._prev_g is not None

        s = x - self._prev_x
        y = g - self._prev_g

        # Sherman-Morrison rank-1 update to inverse Jacobian H:
        # H_{k+1} = H_k + ((s - H_k y) s^T H_k) / (s^T H_k y)
        s_th = s @ self._H
        denom = float(np.dot(s_th, y))
        norm_s = float(np.linalg.norm(s))
        norm_hy = float(np.linalg.norm(self._H @ y))

        if abs(denom) > self.tol_singular * max(norm_s * norm_hy, 1e-14):
            diff = s - self._H @ y
            self._H += np.outer(diff, s_th) / denom
        else:
            # Singular secant safeguard: soft reset to isotropic scaling
            self._H = -self.alpha * np.eye(dim, dtype=np.float64)

        x_next = x - self._H @ g
        self._prev_x = np.asarray(x, dtype=np.float64).copy()
        self._prev_g = np.asarray(g, dtype=np.float64).copy()
        return np.asarray(x_next, dtype=np.float64), 0


class JFNKSolver(FixedPointAccelerator):
    """Jacobian-Free Newton-Krylov solver evaluating Newton steps via matrix-free GMRES."""

    def __init__(
        self,
        gmres_tol: float = 1e-4,
        gmres_maxiter: int = 15,
        eps_fd: float = 1e-7,
        line_search: bool = True,
    ) -> None:
        self.gmres_tol = float(gmres_tol)
        self.gmres_maxiter = int(gmres_maxiter)
        self.eps_fd = float(eps_fd)
        self.line_search = bool(line_search)
        self.reset()

    def reset(self) -> None:
        """No multi-step memory to clear."""

    def step(
        self,
        x: np.ndarray,
        phi_x: np.ndarray,
        residual_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    ) -> tuple[np.ndarray, int]:
        """Solve J(x) delta = -G(x) via matrix-free GMRES and advance state."""
        g = phi_x - x
        dim = x.size

        if residual_fn is None:
            # If no residual function is provided for directional derivatives, fallback to Picard
            return np.asarray(phi_x, dtype=np.float64).copy(), 0

        phi_eval_count = 0
        norm_x = float(np.linalg.norm(x))

        def matvec(v: np.ndarray) -> np.ndarray:
            nonlocal phi_eval_count
            norm_v = float(np.linalg.norm(v))
            if norm_v < 1e-14:
                return np.zeros_like(v)
            # Finite difference perturbation parameter epsilon
            eps = self.eps_fd * (1.0 + norm_x) / norm_v
            x_pert = x + eps * v
            g_pert = residual_fn(x_pert)
            phi_eval_count += 1
            diff = (g_pert - g) / eps
            return np.asarray(diff, dtype=np.float64)

        j_op = LinearOperator((dim, dim), matvec=matvec, dtype=np.float64)

        # Solve J delta = -G(x)
        delta, info = gmres(
            j_op,
            -g,
            rtol=self.gmres_tol,
            maxiter=self.gmres_maxiter,
            atol=1e-12,
        )

        norm_g = float(np.linalg.norm(g))
        t = 1.0

        if self.line_search:
            # Armijo backtracking line search
            c1 = 1e-4
            max_backtracks = 4
            for _ in range(max_backtracks):
                x_trial = x + t * delta
                g_trial = residual_fn(x_trial)
                phi_eval_count += 1
                norm_trial = float(np.linalg.norm(g_trial))
                if norm_trial <= (1.0 - c1 * t) * norm_g:
                    break
                t *= 0.5

        x_next = x + t * delta
        return np.asarray(x_next, dtype=np.float64), phi_eval_count
