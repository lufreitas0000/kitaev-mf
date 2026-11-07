"""Bogoliubov transformation, Thouless state representation, and many-body matrix elements.

Implements the exact overlap machinery for arXiv:2209.12930:
- Bogoliubov blocks X, Y from real Schur mode wavefunctions.
- Canonical Bogoliubov identities verification.
- Regular Thouless pairing matrix Z_b and vacuum overlap.
- Particle-hole interchange for singular transformations (det(X) ~ 0).
- Analytical matrix elements of Majorana operators c_j between distinct gauge sectors.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


def extract_mode_wavefunctions(Q: np.ndarray) -> np.ndarray:
    """Extract complex Dirac mode wavefunctions phi in C^{N x 2N} from real Schur matrix Q in O(2N).

    For canonical 2x2 blocks [[0, eps_k], [-eps_k, 0]], the annihilation operator is:
    a_k = 1/2 (c'_{2k} - i c'_{2k+1}) = sum_m phi_{k, m} c_m
    where phi_{k, m} = 1/2 (Q_{m, 2k} - i Q_{m, 2k+1}).
    """
    N2 = Q.shape[0]
    N = N2 // 2
    phi = np.zeros((N, N2), dtype=np.complex128)
    for k in range(N):
        phi[k, :] = 0.5 * (Q[:, 2 * k] - 1j * Q[:, 2 * k + 1])
    return phi


@dataclass(frozen=True)
class BogoliubovTransformation:
    """Bogoliubov transformation connecting mode bases a (sector G) and b (sector G').

    b_m = sum_k (X_{mk} a_k + Y_{mk} a_k^dagger)
    b_m^dagger = sum_k (Y^*_{mk} a_k + X^*_{mk} a_k^dagger)
    """

    X: np.ndarray  # Shape (N, N), complex128
    Y: np.ndarray  # Shape (N, N), complex128
    Z_b: np.ndarray | None  # Shape (N, N), pairing matrix -(X^dagger)^{-1} Y^T if regular
    det_X: complex  # noqa: N815
    vacuum_overlap_magnitude: float  # sqrt(|det(X)|) if regular, 0.0 if singular
    is_singular: bool
    singular_mode: int | None
    singular_overlap_magnitude: float | None  # sqrt(|det(X')|) after particle-hole flip


def verify_bogoliubov_unitarity(
    X: np.ndarray,
    Y: np.ndarray,
    tol: float = 1e-12,
) -> tuple[bool, dict[str, float]]:
    """Verify the three canonical Bogoliubov relations:

    1. X X^dagger + Y Y^dagger = I
    2. X^dagger X + Y^T Y^* = I
    3. X Y^T + Y X^T = 0
    """
    N = X.shape[0]
    eye_N = np.eye(N, dtype=np.complex128)

    rel1 = X @ X.conj().T + Y @ Y.conj().T
    rel2 = X.conj().T @ X + Y.T @ Y.conj()
    rel3 = X @ Y.T + Y @ X.T

    err1 = float(np.max(np.abs(rel1 - eye_N)))
    err2 = float(np.max(np.abs(rel2 - eye_N)))
    err3 = float(np.max(np.abs(rel3)))

    errors = {
        "max_err_XX_YY": err1,
        "max_err_XdagX_YtYstar": err2,
        "max_err_XYt_YXt": err3,
    }
    is_valid = bool(max(err1, err2, err3) < tol)
    return is_valid, errors


def apply_particle_hole_interchange(
    X: np.ndarray,
    Y: np.ndarray,
    mode_idx: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply particle-hole interchange on mode_idx per Eq. (29) of arXiv:2209.12930.

    Transforms: tilde{b}_mu = b_mu^dagger, tilde{b}_mu^dagger = b_mu.
    Swaps row mode_idx between X and Y:
      X'_{mu, k} = Y^*_{mu, k}
      Y'_{mu, k} = X^*_{mu, k}
    """
    X_prime = X.copy()
    Y_prime = Y.copy()
    X_prime[mode_idx, :] = Y[mode_idx, :].conj()
    Y_prime[mode_idx, :] = X[mode_idx, :].conj()
    return X_prime, Y_prime


def compute_bogoliubov_transformation(
    phi_a: np.ndarray,
    phi_b: np.ndarray,
    tol: float = 1e-6,
) -> BogoliubovTransformation:
    """Compute Bogoliubov matrices X = 2 phi_b phi_a^* and Y = 2 phi_b phi_a^T."""
    X = 2.0 * phi_b @ phi_a.conj().T
    Y = 2.0 * phi_b @ phi_a.T

    det_X = complex(np.linalg.det(X))
    det_mag = float(np.abs(det_X))
    is_singular = det_mag < tol

    singular_mode: int | None = None
    singular_overlap: float | None = None
    Z_b: np.ndarray | None = None

    if is_singular:
        # Transformation between opposite fermion-parity sectors: vacuum overlap vanishes
        vacuum_overlap = 0.0
        # Identify mode whose particle-hole flip maximizes det(X')
        N = X.shape[0]
        best_mode = 0
        best_det = 0.0
        for m in range(N):
            X_m, _ = apply_particle_hole_interchange(X, Y, mode_idx=m)
            d_m = float(np.abs(np.linalg.det(X_m)))
            if d_m > best_det:
                best_det = d_m
                best_mode = m
        singular_mode = best_mode
        singular_overlap = float(np.sqrt(best_det))
    else:
        vacuum_overlap = float(np.sqrt(det_mag))
        # Regular Thouless pairing matrix: Z_b = - (X^dagger)^{-1} Y^T
        X_inv_H = np.linalg.inv(X.conj().T)
        Z_b = -X_inv_H @ Y.T
        # Symmetrization cleanup: ensure exact skew-symmetry
        Z_b = 0.5 * (Z_b - Z_b.T)

    return BogoliubovTransformation(
        X=X,
        Y=Y,
        Z_b=Z_b,
        det_X=det_X,
        vacuum_overlap_magnitude=vacuum_overlap,
        is_singular=is_singular,
        singular_mode=singular_mode,
        singular_overlap_magnitude=singular_overlap,
    )


def compute_majorana_matrix_elements(
    bog: BogoliubovTransformation,
    phi_a: np.ndarray,
    phi_b: np.ndarray,
    site_j: int,
    initial_parity: int = -1,  # -1 for one-fermion ground state (torus), +1 for vacuum
    occupied_mode: int | None = None,  # Mode occupied when initial_parity == -1
    overlap_phase: complex = 1.0 + 0.0j,
) -> dict[str, Any]:
    """Compute analytical transition matrix elements <lambda(G') | c_j | phi_0(G)>.

    Uses exact Wick contractions with the Thouless pairing matrix Z_b and inverse Bogoliubov
    blocks (X^dagger)^{-1}.
    """
    N = phi_a.shape[0]
    X = bog.X
    Y = bog.Y

    if not bog.is_singular and bog.Z_b is not None:
        Z_b = bog.Z_b
        X_inv_H = np.linalg.inv(X.conj().T)
        overlap_0 = overlap_phase * bog.vacuum_overlap_magnitude
        phi_b_eff = phi_b
    else:
        # Singular case: particle-hole interchange on singular mode
        mu = bog.singular_mode if bog.singular_mode is not None else 0
        X_prime, Y_prime = apply_particle_hole_interchange(X, Y, mode_idx=mu)
        X_inv_H = np.linalg.inv(X_prime.conj().T)
        Z_b = -X_inv_H @ Y_prime.T
        Z_b = 0.5 * (Z_b - Z_b.T)
        sing_mag = (
            bog.singular_overlap_magnitude
            if bog.singular_overlap_magnitude is not None
            else float(np.sqrt(np.abs(np.linalg.det(X_prime))))
        )
        overlap_0 = overlap_phase * sing_mag
        phi_b_eff = phi_b.copy()
        phi_b_eff[mu, :] = phi_b[mu, :].conj()

    amplitudes: dict[tuple[int, ...], complex] = {}

    if initial_parity == -1:
        # Initial state is one-fermion state: |phi_0> = a_p^dagger |0_a>
        # c_j |phi_0> has EVEN fermion parity, connecting to even quasiparticle states in b-basis:
        p = occupied_mode

        # 1. Vacuum transition (0 quasiparticles in b-basis):
        # <0_b | c_j | a_p^dagger 0_a> = 2 overlap_0 * sum_m phi^*_{b, mj} [(X^dagger)^{-1}]_{mp}
        elem_0 = 2.0 * overlap_0 * complex(np.dot(phi_b_eff[:, site_j].conj(), X_inv_H[:, p]))
        amplitudes[()] = elem_0

        # 2. Two-quasiparticle transitions (m1, m2): <m1, m2 | c_j | a_p^dagger 0_a>
        # where <m1, m2| = <0_b| b_{m2} b_{m1} = - <0_b| b_{m1} b_{m2}
        for m1 in range(N):
            for m2 in range(m1 + 1, N):
                t1 = 2.0 * phi_b_eff[m2, site_j] * overlap_0 * X_inv_H[m1, p]
                t2 = -2.0 * phi_b_eff[m1, site_j] * overlap_0 * X_inv_H[m2, p]
                t3 = 0.0j
                for k in range(N):
                    wick_4 = overlap_0 * (
                        X_inv_H[k, p] * Z_b[m1, m2]
                        - X_inv_H[m2, p] * Z_b[m1, k]
                        + X_inv_H[m1, p] * Z_b[m2, k]
                    )
                    t3 += 2.0 * phi_b_eff[k, site_j].conj() * wick_4
                amplitudes[(m1, m2)] = -(t1 + t2 + t3)

    else:
        # Initial state is vacuum: |phi_0> = |0_a> (even parity)
        # c_j |0_a> has ODD fermion parity, connecting to odd quasiparticle states:
        for m in range(N):
            val = 2.0 * overlap_0 * complex(np.dot(phi_a[:, site_j], X_inv_H[m, :]))
            amplitudes[(m,)] = val

    sum_sq = float(sum(np.abs(amp) ** 2 for amp in amplitudes.values()))

    return {
        "amplitudes": amplitudes,
        "initial_parity": initial_parity,
        "occupied_mode": occupied_mode,
        "sum_of_squares": sum_sq,
        "is_singular": bog.is_singular,
    }
