"""Matrix sign purification algorithms for evaluating the mean-field map without eigenvectors.

Computes the single-particle Majorana covariance matrix:
    Omega = i * sgn(i * A) in R^{M x M} (Omega^T = -Omega, Omega^2 = -I)
in pure real arithmetic using high-performance matrix multiplications (dgemm).
"""

import time
from dataclasses import dataclass
from typing import Literal

import numpy as np


@dataclass(frozen=True)
class PurificationMetrics:
    """Performance metrics for matrix sign purification."""

    iterations: int
    wall_time_ms: float
    gemm_calls: int
    final_idempotency_error: float
    ground_state_energy: float


def estimate_spectral_bound(A: np.ndarray) -> float:
    """Compute an upper bound alpha >= ||A||_2 for real skew-symmetric matrix A."""
    # For any matrix, ||A||_2 <= min(||A||_1, ||A||_inf, ||A||_F)
    # Since A^T = -A, ||A||_1 = ||A||_inf
    norm_inf = float(np.max(np.sum(np.abs(A), axis=1)))
    norm_fro = float(np.linalg.norm(A))
    # Add small safety buffer to guarantee ||Y_0||_2 < 1.0
    bound = min(norm_inf, norm_fro) * 1.05
    return max(bound, 1e-12)


def _initialize_skew_iterate(
    A: np.ndarray,
    warm_start: np.ndarray | None,
    alpha: float,
) -> tuple[np.ndarray, int]:
    """Initialize skew-symmetric iterate Y_0 with spectral norm strictly below 1."""
    M = A.shape[0]
    if warm_start is not None and warm_start.shape == (M, M):
        # Warm start: refine previous projector with the new Hamiltonian A.
        # Direct Y_0 = warm_start freezes because Omega^2 = -I is a fixed point.
        # Injecting grad_comm = A + omega_0 A omega_0 rotates the occupied-virtual eigenbasis.
        # For any valid Majorana covariance (-I <= omega_0^2 <= 0), ||omega_0||_2 <= 1 even when
        # zero modes exist (3 | L or vacancies). By the triangle inequality,
        # ||Y_init||_2 <= 1 + tau ||grad_comm||_2 <= 1 + delta_bound, so dividing by
        # (1 + delta_bound) * 1.02 guarantees ||Y_0||_2 < 1 in O(M^2) without any O(M^3) SVD.
        omega_0 = 0.5 * (warm_start - warm_start.T)
        grad_comm = A + omega_0 @ A @ omega_0
        tau = 0.25 / max(alpha, 1e-12)
        Y_init = omega_0 - tau * grad_comm
        delta_inf = float(np.max(np.sum(np.abs(grad_comm), axis=1)))
        delta_fro = float(np.linalg.norm(grad_comm))
        delta_bound = tau * min(delta_inf, delta_fro)
        scale = (1.0 + delta_bound) * 1.02
        return 0.5 * (Y_init - Y_init.T) / scale, 2

    # Cold start: scale Hamiltonian so spectral radius is strictly < 1
    # Since H = (i/4) c^T A c, the ground-state covariance satisfies Omega = i sgn(i A) = -A / |A|.
    return -(A / alpha), 0


def purify_newton_schulz(
    A: np.ndarray,
    warm_start: np.ndarray | None = None,
    tol: float = 1e-10,
    max_iter: int = 35,
) -> tuple[np.ndarray, PurificationMetrics]:
    """Compute Omega = i * sgn(i * A) via Newton-Schulz / McWeeny cubic purification.

    Iteration:
        Y_{k+1} = 0.5 * Y_k @ (3 * I + Y_k^2)

    For real skew-symmetric A, Y_k remains strictly real skew-symmetric at every step.

    Parameters:
        A: Real skew-symmetric Hamiltonian matrix of shape (M, M).
        warm_start: Optional previous covariance matrix Omega_{prev} for one first-order
            commutator-gradient basis rotation step (note: subsequent polynomial iterations
            commute with Y_0, so outer SCF loops should keep ``warm_start=False`` with
            outer DIIS/L-Broyden acceleration for exact fixed-point evaluations).
        tol: Convergence tolerance on relative Cauchy difference ||Y_{k+1} - Y_k||_F / sqrt(M).
        max_iter: Maximum purification iterations (must be >= 1).

    Returns:
        Tuple of (Omega, metrics). Note ``metrics.final_idempotency_error`` records the
        terminal relative Cauchy increment from the stopping criterion.
    """
    if max_iter < 1:
        raise ValueError(f"max_iter must be >= 1, got {max_iter}")
    M = A.shape[0]
    t0 = time.perf_counter()
    I_mat = np.eye(M, dtype=np.float64)

    alpha = estimate_spectral_bound(A)
    Y, gemm_count = _initialize_skew_iterate(A, warm_start, alpha)

    final_err = 1.0
    iter_count = 0

    for k in range(max_iter):
        # Compute Y^2 (real symmetric negative semidefinite)
        Y2 = Y @ Y
        gemm_count += 1

        # Step: Y_{k+1} = 0.5 * Y @ (3 * I + Y^2)
        step_factor = 3.0 * I_mat + Y2
        Y_next = 0.5 * (Y @ step_factor)
        gemm_count += 1

        # Maintain exact skew-symmetry
        Y_next = 0.5 * (Y_next - Y_next.T)

        # Cauchy increment: ||Y_{k+1} - Y_k||_F / sqrt(M)
        cauchy_diff = float(np.linalg.norm(Y_next - Y)) / np.sqrt(M)
        final_err = cauchy_diff
        Y = Y_next
        iter_count = k + 1

        if cauchy_diff < tol:
            break

    wall_ms = (time.perf_counter() - t0) * 1000.0

    # Ensure final skew-symmetry
    Omega = 0.5 * (Y - Y.T)

    # Ground-state energy: E = 0.25 * sum_{ij} A_{ij} Omega_{ij}
    gs_energy = 0.25 * float(np.sum(A * Omega))

    metrics = PurificationMetrics(
        iterations=iter_count,
        wall_time_ms=wall_ms,
        gemm_calls=gemm_count,
        final_idempotency_error=final_err,
        ground_state_energy=gs_energy,
    )

    return Omega, metrics


def purify_pade_sign(
    A: np.ndarray,
    warm_start: np.ndarray | None = None,
    tol: float = 1e-10,
    max_iter: int = 35,
    order: Literal["32", "21"] = "32",
) -> tuple[np.ndarray, PurificationMetrics]:
    """Compute Omega = i * sgn(i * A) via Kenney-Laub rational Pade purification.

    Uses the Pade [3/2] or [2/1] rational approximation to the matrix sign function:
    - Pade [3/2] (order="32"):
        Y_{k+1} = Y_k @ (5 * I - 10 * Y_k^2 + Y_k^4) @ (I - 10 * Y_k^2 + 5 * Y_k^4)^{-1}
      Local convergence rate is quintic (p=5), with origin slope f'(0) = 5.
    - Pade [2/1] (order="21", Halley's method):
        Y_{k+1} = Y_k @ (3 * I - Y_k^2) @ (I - 3 * Y_k^2)^{-1}
      Local convergence rate is cubic (p=3), with origin slope f'(0) = 3.

    For real skew-symmetric A with ||Y_k||_2 <= 1, Y_k remains strictly real skew-symmetric
    and the SPD denominator D satisfies cond(D) <= 16 for [3/2] and cond(D) <= 4 for [2/1].

    Parameters:
        A: Real skew-symmetric Hamiltonian matrix of shape (M, M).
        warm_start: Optional previous covariance matrix Omega_{prev}.
        tol: Convergence tolerance on relative Cauchy difference ||Y_{k+1} - Y_k||_F / sqrt(M).
        max_iter: Maximum purification iterations (must be >= 1).
        order: Pade order, either "32" (default, Pade [3/2]) or "21" (Pade [2/1] / Halley).

    Returns:
        Tuple of (Omega, metrics). Note ``metrics.final_idempotency_error`` records the
        terminal relative Cauchy increment from the stopping criterion.
    """
    if max_iter < 1:
        raise ValueError(f"max_iter must be >= 1, got {max_iter}")
    if order not in ("32", "21"):
        raise ValueError(f"order must be '32' or '21', got {order!r}")
    M = A.shape[0]
    t0 = time.perf_counter()
    I_mat = np.eye(M, dtype=np.float64)

    alpha = estimate_spectral_bound(A)
    Y, gemm_count = _initialize_skew_iterate(A, warm_start, alpha)

    final_err = 1.0
    iter_count = 0

    for k in range(max_iter):
        Y2 = Y @ Y  # symmetric negative semidefinite (Y2 <= 0)
        gemm_count += 1

        if order == "32":
            # Pade [3/2]: N = 5 I - 10 Y^2 + Y^4, D = I - 10 Y^2 + 5 Y^4
            Y4 = Y2 @ Y2  # Y4 >= 0
            gemm_count += 1
            N_mat = 5.0 * I_mat - 10.0 * Y2 + Y4
            D_mat = I_mat - 10.0 * Y2 + 5.0 * Y4
        else:
            # Pade [2/1] (Halley): N = 3 I - Y^2, D = I - 3 Y^2
            N_mat = 3.0 * I_mat - Y2
            D_mat = I_mat - 3.0 * Y2

        Y_step = Y @ N_mat
        gemm_count += 1
        Y_next = np.linalg.solve(D_mat, Y_step.T).T

        # Maintain exact skew-symmetry
        Y_next = 0.5 * (Y_next - Y_next.T)

        cauchy_diff = float(np.linalg.norm(Y_next - Y)) / np.sqrt(M)
        final_err = cauchy_diff
        Y = Y_next
        iter_count = k + 1

        if cauchy_diff < tol:
            break

    wall_ms = (time.perf_counter() - t0) * 1000.0

    # Ensure final skew-symmetry
    Omega = 0.5 * (Y - Y.T)

    # Ground-state energy: E = 0.25 * sum_{ij} A_{ij} Omega_{ij}
    gs_energy = 0.25 * float(np.sum(A * Omega))

    metrics = PurificationMetrics(
        iterations=iter_count,
        wall_time_ms=wall_ms,
        gemm_calls=gemm_count,
        final_idempotency_error=final_err,
        ground_state_energy=gs_energy,
    )

    return Omega, metrics
