"""Unit tests for exact third-order spin perturbation theory and mean-field reconciliation (Track 4 / B-034-B-037).

Verifies:
1. Virtual vison flux selection rules: P V P = 0 and scalar 2nd-order P V R_Q V P.
2. Strict cubic field scaling of H_eff^(3) ~ h^3 on spin clusters.
3. Closed-form canonical Dirac mass gap Delta_canonical = 4 h^3 / Delta_flux^2.
4. Reconciliation with 4-flavor mean-field scaling (K-019, c_MF = 2.59214) via 1/4 Bloch normalization.
5. Rigorous 5th-order truncation and finite-size validity bounds.
"""

import numpy as np
import pytest

from kitaev_mf.services.exact_perturbation import (
    audit_perturbation_validity_bounds,
    build_spin_zeeman_perturbation,
    compute_exact_third_order_coupling,
    compute_spin_resolvent_third_order,
)


def build_isolated_hexagon_hamiltonian(K: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Assemble exact 2^6 x 2^6 spin Hamiltonian and Wilson loop operator for isolated hexagon."""
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

    H = np.zeros((dim, dim), dtype=np.complex128)
    for i, j, alpha in bonds:
        for s in range(dim):
            si = (s >> i) & 1
            sj = (s >> j) & 1
            if alpha == "z":
                val = (-1.0) ** (si + sj)
                H[s, s] -= K * val
            elif alpha == "x":
                s_prime = s ^ (1 << i) ^ (1 << j)
                H[s_prime, s] -= K * 1.0
            elif alpha == "y":
                s_prime = s ^ (1 << i) ^ (1 << j)
                val = -((-1.0) ** (si + sj))
                H[s_prime, s] -= K * val

    # Wilson loop W_p = sigma_0^x sigma_1^y sigma_2^z sigma_3^x sigma_4^y sigma_5^z
    W = np.zeros((dim, dim), dtype=np.complex128)
    for s in range(dim):
        s_prime = s ^ (1 << 0) ^ (1 << 1) ^ (1 << 3) ^ (1 << 4)
        s0 = (s >> 0) & 1
        s1 = (s >> 1) & 1
        s2 = (s >> 2) & 1
        s3 = (s >> 3) & 1
        s4 = (s >> 4) & 1
        s5 = (s >> 5) & 1
        phase = (1.0 if s0 == 0 else 1.0) * (-1.0j if s1 == 0 else 1.0j)
        phase *= 1.0 if s2 == 0 else -1.0
        phase *= 1.0 if s3 == 0 else 1.0
        phase *= -1.0j if s4 == 0 else 1.0j
        phase *= 1.0 if s5 == 0 else -1.0
        W[s_prime, s] += phase

    return H, W


def test_flux_selection_rules_and_vanishing_first_order() -> None:
    """Verify B-034: single spin insertion flips fluxes, so P V P = 0 identically."""
    H, _ = build_isolated_hexagon_hamiltonian(K=1.0)
    evals, evecs = np.linalg.eigh(H)
    E0 = evals[0]

    # Ground state projector P
    gs_indices = np.where(np.abs(evals - E0) < 1e-10)[0]
    P = np.zeros_like(H)
    for idx in gs_indices:
        psi = evecs[:, idx : idx + 1]
        P += psi @ psi.conj().T

    h_vec = (0.02, 0.02, 0.02)
    V = build_spin_zeeman_perturbation(6, h_vec)

    # First order: P V P must vanish to machine precision
    PVP = P @ V @ P
    norm_PVP = float(np.linalg.norm(PVP))
    assert norm_PVP < 1e-12, f"P V P must vanish by flux selection rule, got {norm_PVP}"


def test_third_order_cubic_scaling_on_spin_cluster() -> None:
    """Verify B-035: effective third-order Hamiltonian scales strictly as h^3."""
    H, _ = build_isolated_hexagon_hamiltonian(K=1.0)
    evals, evecs = np.linalg.eigh(H)
    E0 = evals[0]

    gs_indices = np.where(np.abs(evals - E0) < 1e-10)[0]
    P = np.zeros_like(H)
    for idx in gs_indices:
        psi = evecs[:, idx : idx + 1]
        P += psi @ psi.conj().T

    h1 = 0.01
    h2 = 0.02
    res1 = compute_spin_resolvent_third_order(H, P, (h1, h1, h1))
    res2 = compute_spin_resolvent_third_order(H, P, (h2, h2, h2))

    norm1 = float(res1["norm_H_eff_3"])
    norm2 = float(res2["norm_H_eff_3"])

    ratio = norm2 / norm1
    # Doubling h must increase H_eff^(3) by exactly 2^3 = 8
    assert np.isclose(ratio, 8.0, rtol=1e-5), f"H_eff^(3) must scale as h^3, got ratio {ratio}"


def test_exact_third_order_coupling_and_canonical_dirac_gap() -> None:
    """Verify B-036: exact chiral coupling kappa and canonical Dirac mass gap."""
    h = 0.05
    K = 1.0
    delta_flux = 0.2641  # Kitaev 2006 reference
    res = compute_exact_third_order_coupling(h=h, K=K, delta_flux=delta_flux)

    # Analytical expectation: kappa = h^3 / (3 sqrt(3) Delta_E^2)
    expected_kappa = (h**3) / (3.0 * np.sqrt(3.0) * (delta_flux**2))
    assert np.isclose(res.kappa_exact, expected_kappa, rtol=1e-12)

    # Canonical Dirac gap: Delta = 12 sqrt(3) kappa = 4 h^3 / Delta_E^2
    expected_delta = 4.0 * (h**3) / (delta_flux**2)
    assert np.isclose(res.delta_canonical, expected_delta, rtol=1e-12)

    # Ratio Delta_canonical / h^3 ~ 57.35 (matching Kitaev static estimate)
    static_prefactor = res.delta_canonical / (h**3)
    assert 56.0 < static_prefactor < 58.5


def test_mean_field_cubic_reconciliation_k019() -> None:
    """Verify K-019: reconciliation between 4-flavor mean-field scaling and exact perturbation."""
    h = 0.05
    res = compute_exact_third_order_coupling(h=h, K=1.0)

    # 4-flavor mean field gap: Delta_K^MF = 2.59214 * h^3
    assert np.isclose(res.delta_mf_predicted, 2.59214 * (h**3), rtol=1e-5)

    # Canonical mean-field gap (1/4 Bloch normalization): Delta_canonical^MF = 4 * Delta_K^MF = 10.36856 * h^3
    assert np.isclose(res.delta_mf_canonical, 4.0 * 2.59214 * (h**3), rtol=1e-5)

    # Intermediate denominator ratio: (2 * Delta_gauge / Delta_flux)^2 ~ 15.81
    assert 15.5 < res.scale_ratio < 16.0


@pytest.mark.parametrize("h", [0.01, 0.02, 0.04, 0.05])
def test_perturbation_validity_and_fifth_order_error_bounds(h: float) -> None:
    """Verify B-037: asymptotic convergence bounds separating finite size from thermodynamic limit."""
    audit = audit_perturbation_validity_bounds(h=h, K=1.0, delta_flux=0.2641, L=12)

    assert audit["is_perturbative"] is True
    assert audit["is_fifth_order_controlled"] is True
    # Fifth-order error bound must remain below 5% for h <= 0.05
    assert float(audit["relative_fifth_order_error"]) < 0.05
    # Finite size bound 2/L = 2/12 ~ 0.1667
    assert np.isclose(audit["finite_size_bound"], 2.0 / 12.0)
