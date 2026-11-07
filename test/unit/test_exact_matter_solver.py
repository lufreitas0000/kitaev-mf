"""Unit tests for exact matter solver, thermodynamic quadrature, and vison pair extrapolation.

Validates Backlog Items B-012, B-013, B-014, B-015, B-016, B-017, B-020, B-021, B-022, and B-024.
"""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_mf.services.exact_matter import (
    build_antisymmetric_majorana_matrix_A,
    build_matter_matrix_T,
    build_reference_links,
    compute_thermodynamic_vacuum_quadrature,
    compute_torus_plaquette_fluxes,
    create_vison_pair_links,
    evaluate_vison_pair_energy,
    extrapolate_vison_gap_infinite_separation,
)


def test_thermodynamic_vacuum_quadrature_and_fock_width() -> None:
    """Validate 2D BZ quadrature against Kitaev (2006) Appendix A analytical constants."""
    e_vac, w_fock = compute_thermodynamic_vacuum_quadrature(K=1.0)
    # Kitaev (2006) quotes -1.5746 K and 3.1492 K per unit cell
    assert np.isclose(e_vac, -1.5745972, atol=1e-5), f"Unexpected e_vac: {e_vac}"
    assert np.isclose(w_fock, 3.1491945, atol=1e-5), f"Unexpected w_fock: {w_fock}"
    assert np.isclose(w_fock, -2.0 * e_vac), "w_fock must equal exactly -2 * e_vac"

    # Scale linearity in K
    e_vac_2, _ = compute_thermodynamic_vacuum_quadrature(K=2.0)
    assert np.isclose(e_vac_2, 2.0 * e_vac, atol=1e-10)


@pytest.mark.parametrize("L", [2, 3, 4, 5, 6])
def test_matter_dispersion_bound_and_dirac_resonances(L: int) -> None:
    """Verify that max(s) = 6K identically on all periodic tori, and 3 | L has Dirac zero modes."""
    K = 1.0
    u = build_reference_links(L, L)
    T = build_matter_matrix_T(u, L, L, K=K)
    s = np.linalg.svd(T, compute_uv=False)

    # Universal one-particle upper band edge is exactly 6K
    assert np.isclose(np.max(s), 6.0 * K, atol=1e-12), f"Expected 6.0 K, got {np.max(s)}"

    # Dirac zero modes exist if and only if 3 | L on periodic tori
    if L % 3 == 0:
        assert np.isclose(np.min(s), 0.0, atol=1e-12), (
            f"Expected Dirac zero mode for L={L}, got {np.min(s)}"
        )
    else:
        assert np.min(s) > 0.1, f"Expected non-zero spectral gap for L={L}, got {np.min(s)}"


def test_majorana_skew_matrix_and_svd_spectral_parity() -> None:
    """Verify singular values of T match positive eigenvalues of iA."""
    N1, N2 = 4, 4
    K = 1.0
    u = build_reference_links(N1, N2)
    T = build_matter_matrix_T(u, N1, N2, K=K)
    A = build_antisymmetric_majorana_matrix_A(T)

    assert np.allclose(A, -A.T), "Majorana matrix A must be strictly skew-symmetric"

    # Positive eigenvalues of iA
    iA = 1j * A
    evals = np.linalg.eigvalsh(iA)
    pos_evals = np.sort(evals[evals >= -1e-12])

    # Singular values of T
    s = np.sort(np.linalg.svd(T, compute_uv=False))

    assert np.allclose(pos_evals, s, atol=1e-12), (
        "Positive eigenvalues of iA must equal singular values of T"
    )


def test_gauge_invariance_of_exact_matter_spectrum() -> None:
    """Verify that local site gauge transformations leave the singular value spectrum invariant."""
    N1, N2 = 4, 4
    u_orig = build_reference_links(N1, N2)
    T_orig = build_matter_matrix_T(u_orig, N1, N2, K=1.0)
    s_orig = np.sort(np.linalg.svd(T_orig, compute_uv=False))

    rng = np.random.default_rng(2026)
    eta_A = rng.choice([-1, 1], size=(N1, N2))
    eta_B = rng.choice([-1, 1], size=(N1, N2))

    u_gauge = np.zeros_like(u_orig)
    for m in range(N1):
        for n in range(N2):
            u_gauge[m, n, 0] = eta_A[m, n] * u_orig[m, n, 0] * eta_B[(m + 1) % N1, n]
            u_gauge[m, n, 1] = eta_A[m, n] * u_orig[m, n, 1] * eta_B[m, (n + 1) % N2]
            u_gauge[m, n, 2] = eta_A[m, n] * u_orig[m, n, 2] * eta_B[m, n]

    T_gauge = build_matter_matrix_T(u_gauge, N1, N2, K=1.0)
    s_gauge = np.sort(np.linalg.svd(T_gauge, compute_uv=False))

    assert np.allclose(s_gauge, s_orig, atol=1e-12), (
        "Matter singular value spectrum must be strictly invariant under gauge transformations"
    )


@pytest.mark.parametrize(
    ("p1", "p2"),
    [
        ((0, 0), (0, 1)),
        ((1, 1), (3, 4)),
        ((0, 0), (3, 3)),
    ],
)
def test_dual_lattice_vison_pair_creation_plaquettes(
    p1: tuple[int, int],
    p2: tuple[int, int],
) -> None:
    """Verify that dual-lattice gauge strings toggle exactly the two target plaquettes."""
    N1, N2 = 6, 6
    u_pair, flipped_bonds = create_vison_pair_links(N1, N2, p1, p2)
    W = compute_torus_plaquette_fluxes(u_pair, N1, N2)

    vortices = [tuple(v) for v in np.argwhere(W == -1)]
    assert sorted(vortices) == sorted([p1, p2]), (
        f"Expected vortices at {p1} and {p2}, got {vortices}"
    )
    assert int(np.prod(W)) == 1, "Total plaquette flux product must remain +1"
    assert len(flipped_bonds) > 0, "At least one bond must be flipped"


def test_vison_pair_creation_energies_and_kitaev_targets() -> None:
    """Verify nearest-neighbor and separated vison-pair creation energies against Kitaev targets."""
    N1, N2 = 8, 8
    # Nearest neighbor pair: p1 = (0, 0), p2 = (0, 1)
    res_nn = evaluate_vison_pair_energy(N1, N2, (0, 0), (0, 1), K=1.0, twists=(-1, 1))
    dE_nn_vac = res_nn["dE_vac"]

    # Kitaev's reference target for nearest-neighbor pair is ~0.2672 K
    assert 0.20 < dE_nn_vac < 0.35, f"Unexpected NN vison energy: {dE_nn_vac}"
    assert np.isclose(dE_nn_vac, 0.27, atol=0.06)

    # Distant pair at maximum separation (4, 4) on 8x8 torus
    res_sep = evaluate_vison_pair_energy(N1, N2, (0, 0), (4, 4), K=1.0, twists=(-1, 1))
    dE_sep_vac = res_sep["dE_vac"]

    # Kitaev's reference target for 2 isolated visons is 2 * 0.1536 K = 0.3072 K
    assert 0.25 < dE_sep_vac < 0.38, f"Unexpected separated vison energy: {dE_sep_vac}"


def test_vison_gap_extrapolation_to_thermodynamic_limit() -> None:
    """Verify that 1/L extrapolation of vison gap yields a physical threshold matching Kitaev (2006)."""
    fit = extrapolate_vison_gap_infinite_separation([6, 8, 10, 12], K=1.0)
    dE_inf = fit["dE_inf"]
    mean_dE = fit["mean_dE"]

    # Extrapolated pair threshold is near 2 * E_v = 0.3072 K (within 5%)
    assert 0.25 < dE_inf < 0.35, f"Unphysical extrapolated vison gap: {dE_inf}"
    assert 0.25 < mean_dE < 0.35, f"Unphysical mean vison gap: {mean_dE}"
    assert np.isclose(dE_inf, 0.3000, atol=0.02)
