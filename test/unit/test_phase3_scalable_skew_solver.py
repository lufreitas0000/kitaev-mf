"""Unit and integration test suite for Task T3.W1.1: Scalable Canonical Skew-Schur Solver.

Verifies:
- H_SCALE_01: Orthogonality ||Q^T Q - I|| < 1e-14 and canonical reconstruction ||Q Sigma Q^T - A|| < 1e-13.
- H_SCALE_02: SVD-based particle-hole mode selection matches exhaustive determinant search.
- H_SCALE_03: SVD mode detection agrees with exhaustive search for a declared zero-mode vacuum.
- H_SCALE_04: Large-cluster quench pairs on L in {6, 8} verify Bogoliubov unitarity and sum rules.
"""

from __future__ import annotations

import sys
from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import schur

from kitaev_mf.services.phase_3.bogoliubov import (
    apply_particle_hole_interchange,
    compute_bogoliubov_transformation,
    compute_majorana_matrix_elements,
    extract_mode_wavefunctions,
    find_optimal_particle_hole_mode_svd,
    verify_bogoliubov_unitarity,
)
from kitaev_mf.services.phase_3.gauge_theory import (
    assemble_matter_skew_matrix,
    build_honeycomb_cylinder_cluster,
    build_honeycomb_torus_cluster,
    build_uniform_gauge_configuration,
    flip_bond_gauge_link,
)
from kitaev_mf.services.phase_3.scalable_skew import (
    ScalableQuenchPair,
    canonical_skew_schur,
    solve_scalable_quench_pair,
)


def _get_canonical_mode_phi(A: np.ndarray) -> np.ndarray:
    """Helper to extract standardized mode wavefunctions phi from skew-symmetric A."""
    N = A.shape[0] // 2
    T, Q = schur(A, output="real")
    for k in range(N):
        if T[2 * k, 2 * k + 1] < 0:
            Q[:, [2 * k, 2 * k + 1]] = Q[:, [2 * k + 1, 2 * k]]
    return extract_mode_wavefunctions(Q)


def test_h_scale_01_canonical_skew_schur_invariants() -> None:
    """Verify Hypothesis H_SCALE_01: Canonical skew-Schur decomposition preserves exact invariants."""
    np.random.seed(42)
    for n in [16, 64, 128]:
        M = np.random.randn(n, n)
        A = M - M.T

        res = canonical_skew_schur(A)

        # 1. Non-negative sorted excitation energies
        assert np.all(res.energies >= 0.0)
        assert np.all(np.diff(res.energies) >= -1e-15)

        # 2. Orthogonality: ||Q^T Q - I||_infty < 1e-14
        assert res.orthogonality_error < 1e-14

        # 3. Determinant in O(2N)
        assert res.det_Q in (-1.0, 1.0)

        # 4. Reconstruction: ||Q Sigma Q^T - A||_infty < 1e-13
        assert res.reconstruction_error < 1e-13


def test_h_scale_02_svd_particle_hole_mode_exactness() -> None:
    """Verify Hypothesis H_SCALE_02: SVD left-singular projection identifies the exact optimal mode."""
    # Find singular transformations from physical skew pairs
    for N in [2, 4, 8]:
        n = 2 * N
        found = False
        for seed in range(200):
            np.random.seed(seed)
            A1 = np.random.randn(n, n)
            A1 = A1 - A1.T
            A2 = np.random.randn(n, n)
            A2 = A2 - A2.T

            phi_a = _get_canonical_mode_phi(A1)
            phi_b = _get_canonical_mode_phi(A2)

            X = 2.0 * phi_b @ phi_a.conj().T
            Y = 2.0 * phi_b @ phi_a.T

            det_mag = float(np.abs(np.linalg.det(X)))
            if det_mag < 1e-4:
                # 1. Exhaustive loop
                best_loop_mode = 0
                best_loop_det = 0.0
                for m in range(N):
                    X_m, _ = apply_particle_hole_interchange(X, Y, mode_idx=m)
                    d = float(np.abs(np.linalg.det(X_m)))
                    if d > best_loop_det:
                        best_loop_det = d
                        best_loop_mode = m

                # 2. SVD selector
                svd_mode, svd_det = find_optimal_particle_hole_mode_svd(X, Y)

                # Both mode and determinant must match
                assert svd_mode == best_loop_mode
                assert np.isclose(svd_det, best_loop_det, atol=1e-12)
                assert svd_det > 1e-4
                found = True
                break
        assert found, f"Expected to find singular Bogoliubov pair for N={N}"


def _declared_cylinder_quench() -> ScalableQuenchPair:
    """Choose an explicit vacuum in the quenched cylinder's exact null pair.

    Opposite initial/quenched canonical orientation exercises the singular
    overlap branch. The zero-mode swap leaves A unchanged within roundoff;
    singularity is a property of this declared state choice, not the boundary
    geometry alone. No assertion about a unique cylinder ground state follows.
    """
    n, bonds, triplets = build_honeycomb_cylinder_cluster(8, 8)
    gauge = build_uniform_gauge_configuration(n, bonds)
    bond = bonds[1]
    changed = flip_bond_gauge_link(gauge, bond.site_a, bond.site_b)
    pair = solve_scalable_quench_pair(
        gauge, changed, kappa=0.1, triplets=triplets, boundary="cylinder", L=8
    )
    matrix = assemble_matter_skew_matrix(changed, kappa=0.1, triplets=triplets)
    assert np.linalg.norm(matrix @ pair.sector_1.Q[:, :2]) < 1e-12
    if pair.sector_0.det_Q != pair.sector_1.det_Q:
        return pair
    q1 = pair.sector_1.Q.copy()
    q1[:, [0, 1]] = q1[:, [1, 0]]
    sector1 = replace(pair.sector_1, Q=q1, det_Q=-pair.sector_1.det_Q)
    phi1 = extract_mode_wavefunctions(q1)
    bog = compute_bogoliubov_transformation(pair.phi_0, phi1)
    valid, errors = verify_bogoliubov_unitarity(bog.X, bog.Y)
    assert valid
    return replace(pair, sector_1=sector1, phi_1=phi1, bogoliubov=bog, unitarity_errors=errors)


def test_declared_cylinder_vacuum_handles_both_null_orientations(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Exercise the alternate Schur null orientation that exposed the remote failure."""
    chosen = _declared_cylinder_quench()
    q1 = chosen.sector_1.Q.copy()
    q1[:, [0, 1]] = q1[:, [1, 0]]
    alternate = replace(
        chosen,
        sector_1=replace(chosen.sector_1, Q=q1, det_Q=-chosen.sector_1.det_Q),
        phi_1=extract_mode_wavefunctions(q1),
    )
    alternate = replace(
        alternate, bogoliubov=compute_bogoliubov_transformation(alternate.phi_0, alternate.phi_1)
    )
    monkeypatch.setattr(
        sys.modules[__name__], "solve_scalable_quench_pair", lambda *args, **kwargs: alternate
    )
    selected = _declared_cylinder_quench()
    assert selected.sector_0.det_Q == -selected.sector_1.det_Q
    assert selected.bogoliubov.is_singular
    assert abs(selected.bogoliubov.det_X) < 1e-12


def test_h_scale_03_svd_mode_selection_matches_exhaustive() -> None:
    """Check the optimal mode for a declared singular vacuum, without a wall-clock gate."""
    pair = _declared_cylinder_quench()
    X, Y = pair.bogoliubov.X, pair.bogoliubov.Y
    assert float(np.abs(np.linalg.det(X))) < 1e-4
    best_loop_mode = 0
    best_loop_det = 0.0
    for m in range(X.shape[0]):
        X_m, _ = apply_particle_hole_interchange(X, Y, mode_idx=m)
        d = float(np.abs(np.linalg.det(X_m)))
        if d > best_loop_det:
            best_loop_det = d
            best_loop_mode = m
    svd_mode, svd_det = find_optimal_particle_hole_mode_svd(X, Y)
    assert svd_mode == best_loop_mode
    assert np.isclose(svd_det, best_loop_det, atol=1e-12)
    assert svd_det > 0.4


def test_h_scale_04_large_cluster_quench_pair_and_sum_rule() -> None:
    """Verify Hypothesis H_SCALE_04: Large cluster quench pairs preserve unitarity and sum rule on L=8 cylinder (N=64) and torus."""
    # 1. Test singular boundary quench on L=8 cylinder (N_sites = 128, N = 64)
    L_cyl = 8
    _, bonds, _ = build_honeycomb_cylinder_cluster(Lx=L_cyl, Ly=L_cyl)
    b0 = bonds[1]

    quench_cyl = _declared_cylinder_quench()

    # Singular quench must be detected and resolved via SVD
    assert quench_cyl.bogoliubov.is_singular is True
    assert quench_cyl.bogoliubov.singular_mode is not None
    assert quench_cyl.bogoliubov.singular_overlap_magnitude is not None
    assert quench_cyl.bogoliubov.singular_overlap_magnitude > 0.5
    assert quench_cyl.unitarity_errors["max_err_XX_YY"] < 1e-12
    assert quench_cyl.unitarity_errors["max_err_XdagX_YtYstar"] < 1e-12
    assert quench_cyl.unitarity_errors["max_err_XYt_YXt"] < 1e-12

    mat_elem_cyl = compute_majorana_matrix_elements(
        bog=quench_cyl.bogoliubov,
        phi_a=quench_cyl.phi_0,
        phi_b=quench_cyl.phi_1,
        site_j=b0.site_a,
        initial_parity=quench_cyl.sector_0.physical_matter_parity,
    )
    # Cylinder boundary sum of squares saturates >= 0.85
    assert mat_elem_cyl["sum_of_squares"] >= 0.85
    assert mat_elem_cyl["sum_of_squares"] <= 1.0001

    # 2. Test torus geometry L=6 (N_sites = 72, N = 36)
    L_torus = 6
    num_sites_t, bonds_t, triplets_t = build_honeycomb_torus_cluster(L_torus)
    g0_t = build_uniform_gauge_configuration(num_sites_t, bonds_t)
    b0_t = bonds_t[0]
    g1_t = flip_bond_gauge_link(g0_t, b0_t.site_a, b0_t.site_b)

    quench_t = solve_scalable_quench_pair(
        gauge_0=g0_t,
        gauge_1=g1_t,
        K=1.0,
        kappa=0.1,
        triplets=triplets_t,
        boundary="torus",
        L=L_torus,
    )

    assert quench_t.unitarity_errors["max_err_XX_YY"] < 1e-12
    assert quench_t.unitarity_errors["max_err_XdagX_YtYstar"] < 1e-12
    assert quench_t.unitarity_errors["max_err_XYt_YXt"] < 1e-12

    mat_elem_t = compute_majorana_matrix_elements(
        bog=quench_t.bogoliubov,
        phi_a=quench_t.phi_0,
        phi_b=quench_t.phi_1,
        site_j=b0_t.site_a,
        initial_parity=quench_t.sector_0.physical_matter_parity,
    )
    assert mat_elem_t["sum_of_squares"] >= 0.95
    assert mat_elem_t["sum_of_squares"] <= 1.0001


def test_schur_separated_null_blocks_reconstruct(monkeypatch: pytest.MonkeyPatch) -> None:
    """A valid real Schur form may interleave 1x1 zeros with a 2x2 rotation."""
    matrix = np.zeros((4, 4))
    matrix[1, 2], matrix[2, 1] = 2.0, -2.0
    monkeypatch.setattr(
        sys.modules[canonical_skew_schur.__module__],
        "schur",
        lambda a, output: (matrix.copy(), np.eye(4)),
    )
    result = canonical_skew_schur(matrix)
    assert np.array_equal(result.energies, [0.0, 2.0])
    assert result.reconstruction_error < 1e-14


@pytest.mark.parametrize("null_pairs", [1, 2, 4])
def test_rank_deficient_schur_matches_independent_hermitian_spectrum(null_pairs: int) -> None:
    rng = np.random.default_rng(71 + null_pairs)
    orthogonal, _ = np.linalg.qr(rng.normal(size=(8, 8)))
    block = np.zeros((8, 8))
    for mode in range(4 - null_pairs):
        block[2 * mode, 2 * mode + 1] = mode + 1
        block[2 * mode + 1, 2 * mode] = -(mode + 1)
    matrix = orthogonal @ block @ orthogonal.T
    matrix = 0.5 * (matrix - matrix.T)
    result = canonical_skew_schur(matrix)
    expected = np.linalg.eigvalsh(1j * matrix)[4:]
    assert np.allclose(result.energies, expected, atol=1e-12, rtol=0)
    assert result.reconstruction_error < 1e-12
    assert result.orthogonality_error < 1e-12


def test_canonical_schur_rejects_invalid_matrices() -> None:
    for matrix in [
        np.eye(3),
        np.zeros((2, 3)),
        np.eye(2),
        np.full((2, 2), np.nan),
        np.eye(2, dtype=complex),
    ]:
        with pytest.raises(ValueError):
            canonical_skew_schur(matrix)
