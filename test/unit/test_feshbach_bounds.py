"""Unit tests for exact Feshbach projectors, resolvent domain bounds, and eliminated-block operator norms.

Validates Backlog Items B-030, B-031, B-032, and B-033.
"""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_mf.services.feshbach import (
    audit_extensive_norm_fallacy,
    compute_spin_feshbach_thresholds,
    compute_torus_feshbach_thresholds,
    evaluate_eliminated_block_resolvent_norm,
)


def build_spin_kitaev_hamiltonian(
    num_spins: int,
    bonds: list[tuple[int, int, str]],
    K: float = 1.0,
) -> np.ndarray:
    """Assemble exact 2^N x 2^N spin-space Hamiltonian H_K = -K sum_{<ij>_a} sigma_i^a sigma_j^a."""
    dim = 1 << num_spins
    H = np.zeros((dim, dim), dtype=np.complex128)
    for i, j, alpha in bonds:
        for s in range(dim):
            si = (s >> i) & 1
            sj = (s >> j) & 1
            if alpha == "z":
                H[s, s] -= K * ((-1.0) ** (si + sj))
            elif alpha == "x":
                s_prime = s ^ (1 << i) ^ (1 << j)
                H[s_prime, s] -= K * 1.0
            elif alpha == "y":
                s_prime = s ^ (1 << i) ^ (1 << j)
                val = -((-1.0) ** (si + sj))
                H[s_prime, s] -= K * val
    return H


def test_isolated_hexagon_spin_feshbach_projectors_and_resolvent() -> None:
    """Verify exact Feshbach projector algebra and resolvent norm on an isolated 6-spin hexagon."""
    num_spins = 6
    dim = 1 << num_spins
    bonds = [
        (0, 1, "z"),
        (1, 2, "x"),
        (2, 3, "y"),
        (3, 4, "z"),
        (4, 5, "x"),
        (5, 0, "y"),
    ]
    H_K = build_spin_kitaev_hamiltonian(num_spins, bonds, K=1.0)

    # Hexagonal Wilson loop W_p = sigma_0^x sigma_1^y sigma_2^z sigma_3^x sigma_4^y sigma_5^z
    W = np.zeros((dim, dim), dtype=np.complex128)
    for s in range(dim):
        s_new = s ^ (1 << 0) ^ (1 << 3) ^ (1 << 1) ^ (1 << 4)
        s1 = (s >> 1) & 1
        s2 = (s >> 2) & 1
        s4 = (s >> 4) & 1
        s5 = (s >> 5) & 1
        phase = ((-1.0) ** (s2 + s5)) * (-((-1.0) ** (s1 + s4)))
        W[s_new, s] = phase

    # Projectors P_ff and Q_flux
    P = 0.5 * (np.eye(dim) + W)
    Q = 0.5 * (np.eye(dim) - W)

    # Projector algebra
    assert np.allclose(P @ P, P), "P must be idempotent"
    assert np.allclose(Q @ Q, Q), "Q must be idempotent"
    assert np.allclose(P @ Q, 0.0), "P and Q must be mutually orthogonal"
    assert np.allclose(P + Q, np.eye(dim)), "P + Q must equal identity"

    # Block-diagonality of unperturbed Hamiltonian
    assert np.allclose(P @ H_K @ Q, 0.0, atol=1e-12), "Off-diagonal block P H Q must vanish"
    assert np.allclose(Q @ H_K @ P, 0.0, atol=1e-12), "Off-diagonal block Q H P must vanish"

    # Compute Feshbach thresholds
    thresh = compute_spin_feshbach_thresholds(H_K, W, energy_window=0.0)
    assert thresh.delta_flux > 0.0, f"Expected positive flux gap, got {thresh.delta_flux}"
    assert np.isclose(thresh.resolvent_norm_E0, 1.0 / thresh.delta_flux, atol=1e-12)

    # Test resolvent norm for E below E0
    for offset in [0.0, 0.5, 1.0]:
        E = thresh.E0 - offset
        expected_norm = 1.0 / (thresh.E_Q_min - E)
        actual_norm = evaluate_eliminated_block_resolvent_norm(H_K, W, E=E)
        assert np.isclose(actual_norm, expected_norm, atol=1e-12), (
            f"Resolvent norm mismatch at E={E}: {actual_norm} vs {expected_norm}"
        )


def test_double_hexagon_spin_feshbach_thresholds() -> None:
    """Verify Feshbach eliminated-block resolvent on a 10-spin two-hexagon cluster."""
    num_spins = 10
    dim = 1 << num_spins
    bonds = [
        (0, 1, "x"),
        (1, 2, "y"),
        (2, 3, "x"),
        (3, 4, "z"),
        (4, 5, "y"),
        (5, 0, "z"),
        (4, 6, "x"),
        (6, 7, "y"),
        (7, 8, "z"),
        (8, 9, "x"),
        (9, 3, "y"),
    ]
    H_K = build_spin_kitaev_hamiltonian(num_spins, bonds, K=1.0)

    # Plaquette 1: W1 = sigma_0^y sigma_1^z sigma_2^z sigma_3^y sigma_4^x sigma_5^x
    W1 = np.zeros((dim, dim), dtype=np.complex128)
    for s in range(dim):
        s_new = s ^ (1 << 0) ^ (1 << 3) ^ (1 << 4) ^ (1 << 5)
        s0, s1, s2, s3 = (s >> 0) & 1, (s >> 1) & 1, (s >> 2) & 1, (s >> 3) & 1
        phase = ((-1.0) ** (s1 + s2)) * (-((-1.0) ** (s0 + s3)))
        W1[s_new, s] = phase

    thresh = compute_spin_feshbach_thresholds(H_K, W1, energy_window=0.0)
    assert thresh.delta_flux > 0.0, f"Expected positive flux gap, got {thresh.delta_flux}"
    assert np.isclose(thresh.delta_flux, 0.292373, atol=1e-4)
    assert np.isclose(thresh.resolvent_norm_E0, 1.0 / thresh.delta_flux, atol=1e-12)


@pytest.mark.parametrize("L", [4, 6, 8, 10, 12])
def test_torus_majorana_feshbach_thresholds_and_resolvent_norm(L: int) -> None:
    """Verify exact Feshbach thresholds and eliminated-block resolvent norms on finite tori."""
    thresh = compute_torus_feshbach_thresholds(L, K=1.0, energy_window=0.04)

    # Flux gap must be strictly positive and within physical bounds [0.04, 0.60] K
    assert 0.04 < thresh.delta_flux < 0.60, (
        f"Unphysical torus flux gap on L={L}: {thresh.delta_flux}"
    )

    # Resolvent operator norm equals 1 / Delta_flux
    assert thresh.resolvent_norm_E0 <= 25.0, (
        f"Excessive resolvent norm on L={L}: {thresh.resolvent_norm_E0}"
    )
    assert np.isclose(thresh.resolvent_norm_E0, 1.0 / thresh.delta_flux, atol=1e-12)

    # Spectral margin with cutoff Lambda = 0.04 must remain positive
    assert thresh.spectral_margin > 0.0, (
        f"Negative spectral margin on L={L}: {thresh.spectral_margin}"
    )


def test_extensive_operator_norm_fallacy_demarcation() -> None:
    """Prove why global ||V|| / Delta_flux diverges extensively while local parameter remains << 1."""
    h = 0.05
    delta_flux = 0.264  # Thermodynamic vison gap

    audit_10 = audit_extensive_norm_fallacy(num_spins=10, h=h, delta_flux=delta_flux)
    audit_100 = audit_extensive_norm_fallacy(num_spins=100, h=h, delta_flux=delta_flux)
    audit_1000 = audit_extensive_norm_fallacy(num_spins=1000, h=h, delta_flux=delta_flux)

    # Global ratio grows extensively with system size:
    assert audit_10["global_ratio"] < audit_100["global_ratio"] < audit_1000["global_ratio"]
    assert audit_100["global_ratio"] > 10.0, (
        "Global operator norm condition fails massively for N=100"
    )

    # Local linked-cluster ratio is intensive and remains strictly small:
    assert np.isclose(audit_10["local_ratio"], audit_100["local_ratio"])
    assert np.isclose(audit_100["local_ratio"], audit_1000["local_ratio"])
    assert audit_100["local_ratio"] < 0.25, "Local perturbation parameter must remain << 1"
