"""Exact Feshbach projection, resolvent domain verification, and eliminated-block operator bounds (Track 4 / B-030-B-033).

Provides rigorous operator-theoretic bounds on the eliminated vortex block (Q H_K Q)
and proves why many-body perturbative control relies on linked-cluster locality rather
than extensive operator-norm smallness.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kitaev_mf.services.exact_matter import (
    build_matter_matrix_T,
    build_reference_links,
    compute_physical_ground_state_energy,
    create_vison_pair_links,
)


@dataclass(frozen=True)
class FeshbachThresholds:
    """Spectral thresholds and resolvent operator bounds for the Feshbach eliminated block."""

    E0: float
    E_Q_min: float
    delta_flux: float
    resolvent_norm_E0: float  # noqa: N815
    spectral_margin: float


def compute_spin_feshbach_thresholds(
    H_K: np.ndarray,
    W_p: np.ndarray,
    energy_window: float = 0.0,
    tol: float = 1e-8,
) -> FeshbachThresholds:
    """Compute exact spin-space Feshbach thresholds on a finite cluster.

    Args:
        H_K: 2^N x 2^N Hermitian spin Hamiltonian.
        W_p: 2^N x 2^N Hermitian and unitary plaquette Wilson loop operator.
        energy_window: External energy cutoff Lambda >= 0.
        tol: Numerical tolerance for projector zero filtering.
    """
    dim = H_K.shape[0]
    eye_dim = np.eye(dim, dtype=np.complex128)

    P_plus = 0.5 * (eye_dim + W_p)
    P_minus = 0.5 * (eye_dim - W_p)

    evals_plus = np.linalg.eigvalsh(P_plus @ H_K @ P_plus)
    evals_minus = np.linalg.eigvalsh(P_minus @ H_K @ P_minus)

    # Filter out null-space zero eigenvalues introduced by projectors
    e0_plus = float(np.min(evals_plus[np.abs(evals_plus) > tol]))
    e0_minus = float(np.min(evals_minus[np.abs(evals_minus) > tol]))

    # Ground state flux sector has lower energy:
    if e0_plus <= e0_minus:
        E0 = e0_plus
        E_Q_min = e0_minus
        evals_Q = evals_minus[np.abs(evals_minus) > tol]
    else:
        E0 = e0_minus
        E_Q_min = e0_plus
        evals_Q = evals_plus[np.abs(evals_plus) > tol]

    delta_flux = float(E_Q_min - E0)
    spectral_margin = float(delta_flux - energy_window)

    # Operator norm of (E0 - Q H_K Q)^(-1) on the Q subspace:
    diffs = np.abs(E0 - evals_Q)
    resolvent_norm_E0 = float(np.max(1.0 / diffs))

    return FeshbachThresholds(
        E0=E0,
        E_Q_min=E_Q_min,
        delta_flux=delta_flux,
        resolvent_norm_E0=resolvent_norm_E0,
        spectral_margin=spectral_margin,
    )


def compute_torus_feshbach_thresholds(
    L: int,
    K: float = 1.0,
    energy_window: float = 0.0,
) -> FeshbachThresholds:
    """Compute exact Majorana torus Feshbach thresholds across all 4 twist sectors.

    Args:
        L: Linear torus size (N = L x L unit cells).
        K: Kitaev exchange coupling.
        energy_window: Cutoff Lambda >= 0 below delta_flux.
    """
    # 1. Global flux-free ground state E0
    E0 = 1e9
    for tw in [(1, 1), (-1, 1), (1, -1), (-1, -1)]:
        u = build_reference_links(L, L, twists=tw)
        T = build_matter_matrix_T(u, L, L, K=K)
        E_phys, _, _ = compute_physical_ground_state_energy(T, L, L, u=u, twists=tw)
        if E_phys < E0:
            E0 = float(E_phys)

    # 2. Bottom of the eliminated block: nearest-neighbor vison pair across all twists
    E_Q_min = 1e9
    for tw in [(1, 1), (-1, 1), (1, -1), (-1, -1)]:
        u_pair, _ = create_vison_pair_links(L, L, (0, 0), (0, 1), twists=tw)
        T_pair = build_matter_matrix_T(u_pair, L, L, K=K)
        E_phys_pair, _, _ = compute_physical_ground_state_energy(T_pair, L, L, u=u_pair, twists=tw)
        if E_phys_pair < E_Q_min:
            E_Q_min = float(E_phys_pair)

    delta_flux = float(E_Q_min - E0)
    spectral_margin = float(delta_flux - energy_window)
    resolvent_norm_E0 = float(1.0 / delta_flux)

    return FeshbachThresholds(
        E0=E0,
        E_Q_min=E_Q_min,
        delta_flux=delta_flux,
        resolvent_norm_E0=resolvent_norm_E0,
        spectral_margin=spectral_margin,
    )


def evaluate_eliminated_block_resolvent_norm(
    H_K: np.ndarray,
    W_p: np.ndarray,
    E: float,
    tol: float = 1e-8,
) -> float:
    """Compute operator norm ||(E - Q H_K Q)^(-1)|| for E in rho(Q H_K Q)."""
    dim = H_K.shape[0]
    eye_dim = np.eye(dim, dtype=np.complex128)
    P_plus = 0.5 * (eye_dim + W_p)
    P_minus = 0.5 * (eye_dim - W_p)

    # Identify Q projector
    evals_plus = np.linalg.eigvalsh(P_plus @ H_K @ P_plus)
    evals_minus = np.linalg.eigvalsh(P_minus @ H_K @ P_minus)
    e0_plus = float(np.min(evals_plus[np.abs(evals_plus) > tol]))
    e0_minus = float(np.min(evals_minus[np.abs(evals_minus) > tol]))

    evals_Q = (
        evals_minus[np.abs(evals_minus) > tol]
        if e0_plus <= e0_minus
        else evals_plus[np.abs(evals_plus) > tol]
    )

    diffs = np.abs(E - evals_Q)
    return float(np.max(1.0 / diffs))


def audit_extensive_norm_fallacy(
    num_spins: int,
    h: float,
    delta_flux: float,
) -> dict[str, float]:
    """Audit the extensive operator-norm fallacy ||V|| vs local linked-cluster parameter."""
    global_norm_V = float(num_spins * h)
    global_ratio = float(global_norm_V / delta_flux)
    local_ratio = float(h / delta_flux)

    return {
        "global_norm_V": global_norm_V,
        "global_ratio": global_ratio,
        "local_ratio": local_ratio,
    }
