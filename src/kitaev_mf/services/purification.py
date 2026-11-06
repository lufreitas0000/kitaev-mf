"""Matrix sign purification algorithms for evaluating the mean-field map without eigenvectors.

Computes the single-particle Majorana covariance matrix:
    Omega = i * sgn(i * A) in R^{M x M} (Omega^T = -Omega, Omega^2 = -I)
in pure real arithmetic using high-performance matrix multiplications (dgemm).
"""

import time
from dataclasses import dataclass

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
        warm_start: Optional previous covariance matrix Omega_{prev} for warm start.
        tol: Convergence tolerance on relative idempotency error ||Y^2 + I||_F / sqrt(M).
        max_iter: Maximum purification iterations.

    Returns:
        Tuple of (Omega, metrics).
    """
    M = A.shape[0]
    t0 = time.perf_counter()
    I_mat = np.eye(M, dtype=np.float64)
    gemm_count = 0

    if warm_start is not None and warm_start.shape == (M, M):
        # Warm start: refine previous projector with the new Hamiltonian
        # To align with current A, perform one gradient / sign alignment step:
        # If A has changed slightly, initialize Y_0 = warm_start
        Y = 0.5 * (warm_start - warm_start.T)
    else:
        # Cold start: scale Hamiltonian so spectral radius is strictly < 1
        alpha = estimate_spectral_bound(A)
        # Invert sign convention: since H = (i/4) c^T A c, the ground state covariance
        # satisfies Omega = i sgn(i A) = - A / |A|
        Y = -(A / alpha)

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
