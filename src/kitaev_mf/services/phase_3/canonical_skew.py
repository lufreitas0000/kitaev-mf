"""Canonical real-skew blocks, including separated Schur null columns."""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
from scipy.linalg import schur


@dataclass(frozen=True)
class ScalableSchurResult:
    """Finite canonical decomposition; null-pair orientation is a state convention."""

    energies: np.ndarray
    Q: np.ndarray
    det_Q: float  # noqa: N815 - preserved public field; Pfaffian sign only if nonsingular
    orthogonality_error: float
    reconstruction_error: float
    elapsed_seconds: float


def canonical_skew_schur(A: np.ndarray, tol: float = 1e-12) -> ScalableSchurResult:
    """Return A = Q (direct sum eps_k J) Q^T to scale-relative precision.

    Real Schur forms have genuine 2x2 rotation blocks and possibly separated
    1x1 zero blocks. Collect the null columns before pairing them. ``tol`` is a
    relative threshold against max(abs(A)), not a claim of physical degeneracy.
    Below-threshold rotations are represented as zero and their reconstruction
    error is reported. No unique null-mode vacuum or parity is inferred.
    """
    start = time.perf_counter()
    if A.ndim != 2 or A.shape[0] != A.shape[1] or A.shape[0] == 0 or A.shape[0] % 2:
        raise ValueError("A must be a nonempty even-dimensional square matrix")
    if not np.isrealobj(A) or not np.all(np.isfinite(A)):
        raise ValueError("A must be real and finite")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be finite and positive")
    size = A.shape[0]
    scale = float(np.max(np.abs(A)))
    threshold = tol * scale
    if np.max(np.abs(A + A.T)) > threshold:
        raise ValueError("A must be skew-symmetric to the declared relative tolerance")
    triangular, basis = schur(A, output="real")
    pairs: list[tuple[float, int, int]] = []
    null_columns: list[int] = []
    index = 0
    while index < size:
        if index + 1 < size and abs(triangular[index + 1, index]) > threshold:
            omega = float(0.5 * (triangular[index, index + 1] - triangular[index + 1, index]))
            first, second = (index, index + 1) if omega >= 0 else (index + 1, index)
            pairs.append((abs(omega), first, second))
            index += 2
        else:
            if abs(triangular[index, index]) > threshold:
                raise ValueError("Unexpected nonzero real Schur eigenvalue for a skew matrix")
            null_columns.append(index)
            index += 1
    if len(null_columns) % 2:
        raise ValueError("Unpaired numerical null column; revise the tolerance")
    pairs.extend(
        (0.0, null_columns[i], null_columns[i + 1]) for i in range(0, len(null_columns), 2)
    )
    pairs.sort(key=lambda pair: pair[0])
    energies = np.array([pair[0] for pair in pairs])
    columns = [column for _, first, second in pairs for column in (first, second)]
    canonical = basis[:, columns]
    sigma = np.zeros_like(A, dtype=float)
    for i, energy in enumerate(energies):
        sigma[2 * i, 2 * i + 1] = energy
        sigma[2 * i + 1, 2 * i] = -energy
    orthogonality = float(np.max(np.abs(canonical.T @ canonical - np.eye(size))))
    reconstruction = float(np.max(np.abs(canonical @ sigma @ canonical.T - A)))
    # Allow accumulated backward error and the explicitly truncated tiny blocks.
    # A badly parsed Schur form must not be returned as a valid eigensystem.
    bound = 10 * size * max(threshold, np.finfo(float).eps * scale)
    if reconstruction > bound:
        raise ValueError("Canonical Schur reconstruction exceeds the declared precision")
    return ScalableSchurResult(
        energies=energies,
        Q=canonical,
        det_Q=float(np.round(np.linalg.det(canonical))),
        orthogonality_error=orthogonality,
        reconstruction_error=reconstruction,
        elapsed_seconds=time.perf_counter() - start,
    )
