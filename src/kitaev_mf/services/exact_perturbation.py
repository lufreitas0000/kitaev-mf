"""Exact third-order spin perturbation theory, virtual vison selection rules, and Dirac gap reconciliation (Track 4 / B-034-B-037).

Computes the exact third-order chiral Majorana coupling kappa_exact, canonical Dirac mass gap,
and reconciles 4-flavor mean-field scaling (c_MF = 2.59214) with exact perturbation theory
via intermediate two-vison vs flat-band gauge scales and the 1/4 Bloch normalization.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ExactPerturbationResult:
    """Exact third-order chiral perturbation parameters and mean-field reconciliation metrics."""

    h: float
    K: float
    delta_flux: float
    kappa_exact: float
    delta_canonical: float
    delta_mf_predicted: float
    delta_mf_canonical: float
    scale_ratio: float
    relative_fifth_order_bound: float


def compute_exact_third_order_coupling(
    h: float,
    K: float = 1.0,
    delta_flux: float = 0.2641,
    delta_gauge: float = 0.52486575,
) -> ExactPerturbationResult:
    """Compute exact third-order chiral coupling and reconcile with 4-flavor mean-field theory.

    Args:
        h: Zeeman field strength along [111] axis (h_x = h_y = h_z = h / sqrt(3)).
        K: Kitaev exchange coupling (default 1.0).
        delta_flux: Physical two-vison creation gap (default 0.2641 from Kitaev 2006).
        delta_gauge: Flat-band gauge Majorana excitation scale K * U_bar (default 0.52486575).

    Returns:
        ExactPerturbationResult containing kappa_exact, delta_canonical, delta_mf, and bounds.
    """
    hx = hy = hz = h / np.sqrt(3.0)
    # kappa_exact = hx * hy * hz / Delta_E^2 = h^3 / (3 sqrt(3) Delta_E^2)
    kappa = (hx * hy * hz) / (delta_flux**2)

    # Canonical Dirac excitation gap Delta_canonical = 12 sqrt(3) kappa = 4 h^3 / Delta_E^2
    delta_canonical = 12.0 * np.sqrt(3.0) * kappa

    # Mean-field cubic prefactor (Phase 1 certified invariant K-019)
    c_mf = 2.59214
    delta_mf_predicted = c_mf * ((h / K) ** 3)

    # Canonical mean-field gap accounting for 1/4 Bloch normalization (H_Bloch = i A / 4)
    delta_mf_canonical = 4.0 * delta_mf_predicted

    # Bare intermediate scale ratio: (2 * Delta_gauge / Delta_flux)^2
    two_delta_gauge = 2.0 * delta_gauge
    scale_ratio = float((two_delta_gauge / delta_flux) ** 2)

    # Relative 5th-order truncation bound: (h / Delta_flux)^2
    relative_fifth_order_bound = float((h / delta_flux) ** 2)

    return ExactPerturbationResult(
        h=h,
        K=K,
        delta_flux=delta_flux,
        kappa_exact=float(kappa),
        delta_canonical=float(delta_canonical),
        delta_mf_predicted=float(delta_mf_predicted),
        delta_mf_canonical=float(delta_mf_canonical),
        scale_ratio=scale_ratio,
        relative_fifth_order_bound=relative_fifth_order_bound,
    )


def build_spin_zeeman_perturbation(
    num_spins: int,
    h_vec: tuple[float, float, float],
) -> np.ndarray:
    """Build the exact 2^N x 2^N Zeeman perturbation V = - sum_i h . sigma_i."""
    dim = 1 << num_spins
    hx, hy, hz = h_vec
    V = np.zeros((dim, dim), dtype=np.complex128)

    for i in range(num_spins):
        for s in range(dim):
            si = (s >> i) & 1
            # sigma_i^z
            if hz != 0.0:
                V[s, s] -= hz * (1.0 if si == 0 else -1.0)
            # sigma_i^x
            if hx != 0.0:
                s_prime = s ^ (1 << i)
                V[s_prime, s] -= hx * 1.0
            # sigma_i^y
            if hy != 0.0:
                s_prime = s ^ (1 << i)
                val = -1.0j if si == 0 else 1.0j
                V[s_prime, s] -= hy * val

    return V


def compute_spin_resolvent_third_order(
    H_K: np.ndarray,
    P: np.ndarray,
    h_vec: tuple[float, float, float],
    tol: float = 1e-10,
) -> dict[str, float | np.ndarray]:
    """Evaluate exact Rayleigh-Schroedinger perturbation theory through 3rd order in spin space.

    Args:
        H_K: 2^N x 2^N unperturbed Kitaev Hamiltonian.
        P: Projector onto unperturbed physical ground subspace.
        h_vec: Zeeman field vector (hx, hy, hz).
        tol: Tolerance for null-space filtering.

    Returns:
        Dictionary with norms of H_eff^(1), H_eff^(2), and H_eff^(3).
    """
    dim = H_K.shape[0]
    num_spins = int(np.round(np.log2(dim)))
    V = build_spin_zeeman_perturbation(num_spins, h_vec)

    evals, evecs = np.linalg.eigh(H_K)
    E0 = float(np.min(evals))

    # Ground state subspace projector
    gs_mask = np.abs(evals - E0) < tol
    Q_mask = ~gs_mask

    # Resolvent R_Q(E0) = sum_{m in Q} |m><m| / (E0 - Em)
    R_Q = np.zeros((dim, dim), dtype=np.complex128)
    for m in range(dim):
        if Q_mask[m]:
            denom = E0 - evals[m]
            psi_m = evecs[:, m : m + 1]
            R_Q += (psi_m @ psi_m.conj().T) / denom

    # 1st-order: P V P
    H_eff_1 = P @ V @ P
    norm_1 = float(np.linalg.norm(H_eff_1))

    # 2nd-order: P V R_Q V P
    H_eff_2 = P @ V @ R_Q @ V @ P
    norm_2 = float(np.linalg.norm(H_eff_2))

    # 3rd-order: P V R_Q V R_Q V P
    H_eff_3 = P @ V @ R_Q @ V @ R_Q @ V @ P
    norm_3 = float(np.linalg.norm(H_eff_3))

    return {
        "norm_H_eff_1": norm_1,
        "norm_H_eff_2": norm_2,
        "norm_H_eff_3": norm_3,
        "H_eff_3": H_eff_3,
        "E0": E0,
    }


def audit_perturbation_validity_bounds(
    h: float,
    K: float = 1.0,
    delta_flux: float = 0.2641,
    L: int = 12,
) -> dict[str, float | bool]:
    """Audit mathematical validity and asymptotic convergence bounds for 3rd-order truncation.

    Args:
        h: Magnetic field strength along [111].
        K: Kitaev coupling.
        delta_flux: Vison creation energy scale.
        L: Linear torus dimension.

    Returns:
        Dictionary of validity indicators and bounds.
    """
    lambda_intensive = float(h / delta_flux)
    rel_fifth_order = float(lambda_intensive**2)
    finite_size_bound = float(2.0 / L)

    is_perturbative = lambda_intensive < 0.30
    is_fifth_order_controlled = rel_fifth_order < 0.10

    return {
        "lambda_intensive": lambda_intensive,
        "relative_fifth_order_error": rel_fifth_order,
        "finite_size_bound": finite_size_bound,
        "is_perturbative": is_perturbative,
        "is_fifth_order_controlled": is_fifth_order_controlled,
    }
