"""Exact 16-state & 4096-state Fock parity tests for the cross-QSL convention bridge (kitaev-mf#25)."""

from __future__ import annotations

from math import sqrt
from typing import Any

import numpy as np
import pytest
from numpy.typing import NDArray
from scipy import sparse

from kitaev_mf.domain.constants import G_MATRICES, M_MATRICES, N_MATRICES
from kitaev_mf.domain.lattice import (
    Bond,
    Plaquette,
    UnitCell,
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
    make_honeycomb_supercell,
)
from kitaev_mf.domain.ports import QSLConventionBridgePort
from kitaev_mf.domain.qsl_bridge import (
    CARTESIAN_AXES,
    CONVENTIONS_VERSION,
    P0_PARITY_BRIDGE,
    SO4_L,
    SO4_L_STACK,
    SO4_R,
    SO4_R_STACK,
    DefaultQSLConventionBridge,
    QSLInhomogeneousPortState,
    QSLMeanFieldPortState,
    axis_to_bond_type,
    bond_type_to_axis,
    chiral_full_gap_canonical,
    chiral_kappa_kmf_static_to_kslsc,
    chiral_kappa_kmf_static_to_kwg,
    chiral_kappa_kwg_to_kmf_static,
    compute_plaquette_wilson_loops_from_u,
    conjugate_multisite_p0,
    conjugate_site_block_p0,
    convert_inhomogeneous_state_kmf_to_v2,
    convert_inhomogeneous_state_v2_to_kmf,
    convert_state_kmf_to_v2,
    convert_state_v2_to_kmf,
    extract_pure_kitaev_bond_scalars,
    extract_site_su2_expectations,
    kitaev_coupling_bilayer_to_pauli,
    kitaev_coupling_pauli_to_bilayer,
    kmf_gauge_generator,
    kmf_spin_generator,
    lambdas_kmf_to_v2,
    lambdas_v2_to_kmf,
    pauli_bilinear_matrix,
    reconstruct_pure_kitaev_bond_tensors,
)
from kitaev_mf.domain.state import InhomogeneousMeanFieldState, MeanFieldState
from kitaev_mf.services.observables import (
    compute_plaquette_fluxes,
    compute_sublattice_magnetizations,
)
from kitaev_mf.services.real_space import RealSpaceTorusCluster

EYE2 = np.eye(2, dtype=np.complex128)
PX = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PY = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
PZ = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
TOL = 1.0e-13


def _kron_dense(factors: list[NDArray[np.complex128]]) -> NDArray[np.complex128]:
    out = factors[0]
    for fac in factors[1:]:
        out = np.asarray(np.kron(out, fac), dtype=np.complex128)
    return out


def _jw_majoranas_dense(num_complex_modes: int) -> tuple[NDArray[np.complex128], ...]:
    """Independent Jordan-Wigner representation of 2*num_complex_modes Hermitian Majoranas."""
    ops: list[NDArray[np.complex128]] = []
    for mode in range(num_complex_modes):
        for local in (PX, PY):
            factors = [PZ if m < mode else EYE2 for m in range(num_complex_modes)]
            factors[mode] = local
            ops.append(_kron_dense(factors))
    return tuple(ops)


def _quad_fock(
    mat: NDArray[np.float64],
    ops: tuple[NDArray[np.complex128], ...] | list[NDArray[np.complex128]],
) -> NDArray[np.complex128]:
    """Evaluate H = (i/4) sum_{u,v} mat[u,v] ops[u] ops[v] in exact Fock space."""
    dim = ops[0].shape[0]
    out = np.zeros((dim, dim), dtype=np.complex128)
    n = len(ops)
    for u in range(n):
        for v in range(n):
            if mat[u, v] != 0.0:
                out += (0.25j * float(mat[u, v])) * (ops[u] @ ops[v])
    return out


def _sparse_jw_majoranas(num_complex_modes: int) -> list[sparse.csc_matrix]:
    """Sparse Jordan-Wigner Majorana operators for exact 6-site (12-qubit, dim=4096) plaquette tests."""
    s_i = sparse.eye(2, dtype=np.complex128, format="csc")
    s_x = sparse.csc_matrix(PX)
    s_y = sparse.csc_matrix(PY)
    s_z = sparse.csc_matrix(PZ)
    ops: list[sparse.csc_matrix] = []
    for mode in range(num_complex_modes):
        for local in (s_x, s_y):
            factors = [s_z if m < mode else s_i for m in range(num_complex_modes)]
            factors[mode] = local
            cur = factors[0]
            for fac in factors[1:]:
                cur = sparse.kron(cur, fac, format="csc")
            ops.append(cur)
    return ops


def test_so4_lr_lie_algebra_immutability_and_cartesian_keys() -> None:
    """Verify [L^a,L^b]=eps L^c, [R^a,R^b]=eps R^c, [L^a,R^b]=0, Frobenius norms, and strict keys."""
    assert CARTESIAN_AXES == ("x", "y", "z")
    assert SO4_L_STACK.shape == (3, 4, 4)
    assert SO4_R_STACK.shape == (3, 4, 4)
    assert not P0_PARITY_BRIDGE.flags.writeable
    assert not SO4_L_STACK.flags.writeable
    assert not SO4_R_STACK.flags.writeable

    for ax in CARTESIAN_AXES:
        assert not SO4_L[ax].flags.writeable
        assert not SO4_R[ax].flags.writeable
        assert np.allclose(SO4_L[ax].T, -SO4_L[ax], atol=TOL, rtol=0.0)
        assert np.allclose(SO4_R[ax].T, -SO4_R[ax], atol=TOL, rtol=0.0)
        assert np.isclose(float(np.trace(SO4_L[ax].T @ SO4_L[ax])), 1.0, atol=TOL, rtol=0.0)
        assert np.isclose(float(np.trace(SO4_R[ax].T @ SO4_R[ax])), 1.0, atol=TOL, rtol=0.0)

    for a, b, c in (("x", "y", "z"), ("y", "z", "x"), ("z", "x", "y")):
        assert np.allclose(SO4_L[a] @ SO4_L[b] - SO4_L[b] @ SO4_L[a], SO4_L[c], atol=TOL, rtol=0.0)
        assert np.allclose(SO4_R[a] @ SO4_R[b] - SO4_R[b] @ SO4_R[a], SO4_R[c], atol=TOL, rtol=0.0)

    for a in CARTESIAN_AXES:
        for b in CARTESIAN_AXES:
            assert np.allclose(SO4_L[a] @ SO4_R[b] - SO4_R[b] @ SO4_L[a], 0.0, atol=TOL, rtol=0.0)
            assert np.isclose(float(np.trace(SO4_L[a].T @ SO4_R[b])), 0.0, atol=TOL, rtol=0.0)

    # Strict Cartesian key policy (numeric keys 1,2,3 or '1' forbidden per conventions-v2.0)
    for bad_key in ("1", "2", "3", "0", "a"):
        with pytest.raises(KeyError):
            _ = SO4_L[bad_key]
        with pytest.raises(KeyError):
            _ = SO4_R[bad_key]
        with pytest.raises(KeyError):
            _ = axis_to_bond_type(bad_key)

    so4_any: Any = SO4_L
    with pytest.raises(KeyError):
        _ = so4_any[1]
    with pytest.raises(TypeError):
        so4_any["w"] = np.zeros((4, 4))


def test_p0_parity_bridge_swaps_left_right_and_matches_legacy_constants() -> None:
    """Verify exact dictionary to G_MATRICES/M_MATRICES/N_MATRICES and P_0 L^a P_0 = R^a."""
    assert np.allclose(P0_PARITY_BRIDGE, P0_PARITY_BRIDGE.T, atol=TOL, rtol=0.0)
    assert np.allclose(P0_PARITY_BRIDGE @ P0_PARITY_BRIDGE, np.eye(4), atol=TOL, rtol=0.0)
    assert np.isclose(float(np.linalg.det(P0_PARITY_BRIDGE)), -1.0, atol=TOL, rtol=0.0)

    for b_type, ax in ((1, "x"), (2, "y"), (3, "z")):
        assert axis_to_bond_type(ax) == b_type
        assert bond_type_to_axis(b_type) == ax

        # Direct dictionary locked in contracts/qsl_conventions.json
        assert np.allclose(0.5 * M_MATRICES[b_type], SO4_R[ax], atol=TOL, rtol=0.0)
        assert np.allclose(-0.5 * G_MATRICES[b_type], SO4_L[ax], atol=TOL, rtol=0.0)
        assert np.allclose(kmf_spin_generator(ax), 0.5 * M_MATRICES[b_type], atol=TOL, rtol=0.0)
        assert np.allclose(kmf_gauge_generator(ax), -0.5 * G_MATRICES[b_type], atol=TOL, rtol=0.0)
        assert np.allclose(N_MATRICES[b_type], SO4_R[ax] - SO4_L[ax], atol=TOL, rtol=0.0)

        # P_0 = diag(-1, 1, 1, 1) swaps Left <-> Right and negates N^a
        assert np.allclose(conjugate_site_block_p0(SO4_L[ax]), SO4_R[ax], atol=TOL, rtol=0.0)
        assert np.allclose(conjugate_site_block_p0(SO4_R[ax]), SO4_L[ax], atol=TOL, rtol=0.0)
        assert np.allclose(
            P0_PARITY_BRIDGE @ SO4_L[ax] @ P0_PARITY_BRIDGE, SO4_R[ax], atol=TOL, rtol=0.0
        )
        assert np.allclose(
            conjugate_site_block_p0(pauli_bilinear_matrix(ax, basis="kitaev_mf")),
            pauli_bilinear_matrix(ax, basis="conventions_v2"),
            atol=TOL,
            rtol=0.0,
        )
        assert np.allclose(
            pauli_bilinear_matrix(ax, basis="conventions_v2"),
            -N_MATRICES[b_type],
            atol=TOL,
            rtol=0.0,
        )


def test_exact_16state_two_site_fock_parity_spin_gauge_and_wick_equivalence() -> None:
    """Verify exact 16-state (2-site, 8-Majorana) Fock parity between (-c, b) and (+c, b) bases."""
    # Canonical v2.0 basis on two sites: gamma_s = (+c_s, b_s^x, b_s^y, b_s^z)
    gamma_all = _jw_majoranas_dense(4)
    assert len(gamma_all) == 8
    eye16 = np.eye(16, dtype=np.complex128)

    # kitaev-mf basis on two sites: c_s = P_0 gamma_s = (-c_s, b_s^x, b_s^y, b_s^z)
    kmf_all = tuple((-gamma_all[idx] if idx % 4 == 0 else gamma_all[idx]) for idx in range(8))

    projectors: list[NDArray[np.complex128]] = []
    for site in (0, 1):
        g_s = gamma_all[4 * site : 4 * site + 4]
        c_s = kmf_all[4 * site : 4 * site + 4]

        # Local Kitaev parity equality: +c^0 c^1 c^2 c^3 == +b^x b^y b^z c == -(-1)^{n1+n2}
        d_kmf = c_s[0] @ c_s[1] @ c_s[2] @ c_s[3]
        d_v2 = g_s[1] @ g_s[2] @ g_s[3] @ g_s[0]
        assert np.allclose(d_kmf, d_v2, atol=TOL, rtol=0.0)

        # Dirac occupation parity check in kitaev-mf basis
        d1 = 0.5 * (c_s[0] + 1.0j * c_s[1])
        d2 = 0.5 * (c_s[2] + 1.0j * c_s[3])
        occ_parity = (eye16 - 2.0 * (d1.conj().T @ d1)) @ (eye16 - 2.0 * (d2.conj().T @ d2))
        assert np.allclose(d_kmf, -occ_parity, atol=TOL, rtol=0.0)

        p_plus = 0.5 * (eye16 + d_v2)
        p_minus = 0.5 * (eye16 - d_v2)
        projectors.append(p_plus)

        s_ops_v2 = {ax: _quad_fock(SO4_L[ax], g_s) for ax in CARTESIAN_AXES}
        g_ops_v2 = {ax: _quad_fock(SO4_R[ax], g_s) for ax in CARTESIAN_AXES}
        s_ops_kmf = {ax: _quad_fock(kmf_spin_generator(ax), c_s) for ax in CARTESIAN_AXES}
        g_ops_kmf = {ax: _quad_fock(kmf_gauge_generator(ax), c_s) for ax in CARTESIAN_AXES}

        for b_type, ax in ((1, "x"), (2, "y"), (3, "z")):
            sigma_v2 = 1.0j * (g_s[b_type] @ g_s[0])
            sigma_kmf = 1.0j * (c_s[0] @ c_s[b_type])
            assert np.allclose(sigma_kmf, sigma_v2, atol=TOL, rtol=0.0)
            assert np.allclose(
                2.0 * _quad_fock(pauli_bilinear_matrix(ax, basis="kitaev_mf"), c_s),
                sigma_v2,
                atol=TOL,
                rtol=0.0,
            )
            assert np.allclose(
                2.0 * _quad_fock(pauli_bilinear_matrix(ax, basis="conventions_v2"), g_s),
                sigma_v2,
                atol=TOL,
                rtol=0.0,
            )

            # Exact operator equality across bases
            assert np.allclose(s_ops_kmf[ax], s_ops_v2[ax], atol=TOL, rtol=0.0)
            assert np.allclose(g_ops_kmf[ax], g_ops_v2[ax], atol=TOL, rtol=0.0)
            # Legacy T^a = (1/2) O[G^a] == -G_v2^a
            t_legacy = 0.5 * _quad_fock(G_MATRICES[b_type], c_s)
            assert np.allclose(t_legacy, -g_ops_v2[ax], atol=TOL, rtol=0.0)

            # Physical projection & unprojected decomposition sigma = 2(S - G_v2)
            assert np.allclose(2.0 * s_ops_v2[ax], p_plus @ sigma_v2, atol=TOL, rtol=0.0)
            assert np.allclose(2.0 * g_ops_v2[ax], -p_minus @ sigma_v2, atol=TOL, rtol=0.0)
            assert np.allclose(g_ops_v2[ax] @ p_plus, 0.0, atol=TOL, rtol=0.0)
            assert np.allclose(2.0 * (s_ops_v2[ax] - g_ops_v2[ax]), sigma_v2, atol=TOL, rtol=0.0)

        # Casimirs S^2 = (3/8)(1 + D), G^2 = (3/8)(1 - D)
        s_cas = sum(s_ops_v2[ax] @ s_ops_v2[ax] for ax in CARTESIAN_AXES)
        g_cas = sum(g_ops_v2[ax] @ g_ops_v2[ax] for ax in CARTESIAN_AXES)
        assert np.allclose(s_cas, 0.375 * (eye16 + d_v2), atol=TOL, rtol=0.0)
        assert np.allclose(g_cas, 0.375 * (eye16 - d_v2), atol=TOL, rtol=0.0)

    # Global 2-site physical subspace has dimension 2^2 = 4 inside 16-state Fock space
    p_global = projectors[0] @ projectors[1]
    assert np.isclose(np.trace(p_global).real, 4.0, atol=TOL, rtol=0.0)

    # Exact Lagrange multiplier Hamiltonian equivalence H_lam^kmf == H_lam^v2
    lam_kmf = np.array([[0.37, -0.21, 0.49], [-0.15, 0.62, -0.33]], dtype=np.float64)
    lam_v2 = lambdas_kmf_to_v2(lam_kmf)
    assert np.allclose(lambdas_v2_to_kmf(lam_v2), lam_kmf, atol=TOL, rtol=0.0)
    a_lam_kmf = np.zeros((8, 8), dtype=np.float64)
    a_lam_v2 = np.zeros((8, 8), dtype=np.float64)
    for site in (0, 1):
        sl = slice(4 * site, 4 * site + 4)
        for b_type, ax in ((1, "x"), (2, "y"), (3, "z")):
            a_lam_kmf[sl, sl] += lam_kmf[site, b_type - 1] * G_MATRICES[b_type]
            a_lam_v2[sl, sl] += lam_v2[site, b_type - 1] * SO4_R[ax]
    assert np.allclose(conjugate_multisite_p0(a_lam_kmf), a_lam_v2, atol=TOL, rtol=0.0)
    assert np.allclose(
        _quad_fock(a_lam_kmf, kmf_all), _quad_fock(a_lam_v2, gamma_all), atol=TOL, rtol=0.0
    )

    # Random non-degenerate 8-Majorana Gaussian state with nonzero off-diagonal U^{0,a} and <G> != 0
    rng = np.random.default_rng(20261124)
    raw = rng.normal(size=(8, 8))
    a_kmf = 0.5 * (raw - raw.T)
    a_v2 = conjugate_multisite_p0(a_kmf)
    h_kmf = _quad_fock(a_kmf, kmf_all)
    h_v2 = _quad_fock(a_v2, gamma_all)
    assert np.allclose(h_kmf, h_v2, atol=TOL, rtol=0.0)

    evals, evecs = np.linalg.eigh(h_kmf)
    assert evals[1] - evals[0] > 1.0e-3
    gs = evecs[:, 0]

    # Exact covariance matrices from direct Fock expectation values <i q_u q_v>_{u!=v}
    om_kmf = np.zeros((8, 8), dtype=np.float64)
    om_v2 = np.zeros((8, 8), dtype=np.float64)
    for u in range(8):
        for v in range(8):
            if u != v:
                om_kmf[u, v] = float(np.vdot(gs, (1.0j * kmf_all[u] @ kmf_all[v]) @ gs).real)
                om_v2[u, v] = float(np.vdot(gs, (1.0j * gamma_all[u] @ gamma_all[v]) @ gs).real)
    assert np.allclose(conjugate_multisite_p0(om_kmf), om_v2, atol=TOL, rtol=0.0)
    # Non-vacuity: verify off-diagonal matter-gauge blocks U_01^{0, 1..3} are genuinely nonzero and flip sign
    u_kmf_01 = om_kmf[:4, 4:]
    u_v2_01 = om_v2[:4, 4:]
    assert float(np.linalg.norm(u_kmf_01[0, 1:])) > 1.0e-2
    assert np.allclose(u_v2_01[0, 1:], -u_kmf_01[0, 1:], atol=TOL, rtol=0.0)
    assert np.allclose(u_v2_01[1:, 0], -u_kmf_01[1:, 0], atol=TOL, rtol=0.0)
    assert np.isclose(u_v2_01[0, 0], u_kmf_01[0, 0], atol=TOL, rtol=0.0)
    assert np.allclose(u_v2_01[1:, 1:], u_kmf_01[1:, 1:], atol=TOL, rtol=0.0)

    # Verify site SU(2)_L x SU(2)_R extraction matches exact Fock expectations & observables.py
    cell2 = make_honeycomb_2site()
    state_kmf = MeanFieldState(
        U=np.stack([u_kmf_01, u_kmf_01, u_kmf_01], axis=0),
        V=np.stack([om_kmf[:4, :4], om_kmf[4:, 4:]], axis=0),
        lambdas=lam_kmf,
    )
    sigma_obs, _ = compute_sublattice_magnetizations(state_kmf, cell2)

    for site in (0, 1):
        sl = slice(4 * site, 4 * site + 4)
        dec_kmf = extract_site_su2_expectations(om_kmf[sl, sl], basis="kitaev_mf")
        dec_v2 = extract_site_su2_expectations(om_v2[sl, sl], basis="conventions_v2")
        assert np.allclose(dec_kmf.spin_s, dec_v2.spin_s, atol=TOL, rtol=0.0)
        assert np.allclose(dec_kmf.gauge_g, dec_v2.gauge_g, atol=TOL, rtol=0.0)
        assert np.allclose(
            dec_kmf.pauli_sigma_unproj, dec_v2.pauli_sigma_unproj, atol=TOL, rtol=0.0
        )
        assert np.allclose(dec_kmf.pauli_sigma_unproj, sigma_obs[site], atol=TOL, rtol=0.0)
        # Non-vacuity: random Gaussian state has nonzero spin AND nonzero gauge residual
        assert float(np.linalg.norm(dec_kmf.spin_s)) > 1.0e-2
        assert float(np.linalg.norm(dec_kmf.gauge_g)) > 1.0e-2

        g_site = gamma_all[4 * site : 4 * site + 4]
        for idx, ax in enumerate(CARTESIAN_AXES):
            exp_s = float(np.vdot(gs, _quad_fock(SO4_L[ax], g_site) @ gs).real)
            exp_g = float(np.vdot(gs, _quad_fock(SO4_R[ax], g_site) @ gs).real)
            assert np.isclose(dec_kmf.spin_s[idx], exp_s, atol=TOL, rtol=0.0)
            assert np.isclose(dec_kmf.gauge_g[idx], exp_g, atol=TOL, rtol=0.0)

    # Exact 2-site Kitaev bond exchange parity: +K_Pauli sigma_0^z sigma_1^z == -K_bilayer sigma_0^z sigma_1^z
    k_pauli = -1.0
    k_bilayer = kitaev_coupling_pauli_to_bilayer(k_pauli)
    assert np.isclose(k_bilayer, +1.0, atol=TOL, rtol=0.0)
    assert np.isclose(kitaev_coupling_bilayer_to_pauli(k_bilayer), k_pauli, atol=TOL, rtol=0.0)
    sig0_z = 1.0j * (gamma_all[3] @ gamma_all[0])
    sig1_z = 1.0j * (gamma_all[7] @ gamma_all[4])
    u_op_z = 1.0j * (gamma_all[3] @ gamma_all[7])
    v_op_z = 1.0j * (gamma_all[0] @ gamma_all[4])
    h_bond_pauli = k_pauli * (sig0_z @ sig1_z)
    h_bond_bilayer = -k_bilayer * (sig0_z @ sig1_z)
    # Using sigma_0^z sigma_1^z = -(i b_0^z b_1^z)(i c_0 c_1) = -u_op_z v_op_z
    assert np.allclose(h_bond_pauli, h_bond_bilayer, atol=TOL, rtol=0.0)
    assert np.allclose(h_bond_bilayer, k_bilayer * (u_op_z @ v_op_z), atol=TOL, rtol=0.0)


def test_exact_6site_4096state_fock_plaquette_wp_sign_and_pfaffian_parity() -> None:
    """Prove P_+ W_p^{spin} P_+ == +P_+ prod_{b in partial p} (i b_{A_b}^{alpha_b} b_{B_b}^{alpha_b}) P_+ on 2^12 Fock space."""
    # 6 sites on a hexagon: s=0(A1), 1(B2), 2(A3), 3(B4), 4(A5), 5(B6)
    # Perimeter bonds (matching RealSpaceTorusCluster.compute_plaquette_gauge_majoranas):
    # b1: A1(0)->B2(1) z(3); b2: A3(2)->B2(1) x(1); b3: A3(2)->B4(3) y(2);
    # b4: A5(4)->B4(3) z(3); b5: A5(4)->B6(5) x(1); b6: A1(0)->B6(5) y(2).
    # Exterior Pauli legs around hexagon: s0:x, s1:y, s2:z, s3:x, s4:y, s5:z.
    gamma_24 = _sparse_jw_majoranas(12)  # 4096 x 4096 sparse CSC operators
    dim = 1 << 12
    eye_s = sparse.eye(dim, dtype=np.complex128, format="csc")

    # Project random 4096-state vector onto the exact 64-dim physical subspace D_s = +1 for all 6 sites
    rng = np.random.default_rng(7711)
    psi = rng.normal(size=dim) + 1.0j * rng.normal(size=dim)
    for s in range(6):
        c_s, bx_s, by_s, bz_s = gamma_24[4 * s : 4 * s + 4]
        d_s = bx_s @ by_s @ bz_s @ c_s
        psi = 0.5 * (psi + d_s @ psi)
    norm_psi = float(np.linalg.norm(psi))
    assert norm_psi > 1.0e-3
    psi /= norm_psi

    # Also verify exact kitaev-mf parity operator +c^0 c^1 c^2 c^3 leaves psi invariant (+1 eigenvalue)
    for s in range(6):
        c0_kmf = -gamma_24[4 * s]
        d_kmf = c0_kmf @ gamma_24[4 * s + 1] @ gamma_24[4 * s + 2] @ gamma_24[4 * s + 3]
        assert np.allclose(d_kmf @ psi, psi, atol=TOL, rtol=0.0)

    # Apply W_p^{spin} = sigma_0^x sigma_1^y sigma_2^z sigma_3^x sigma_4^y sigma_5^z
    ext_legs = ((0, 1), (1, 2), (2, 3), (3, 1), (4, 2), (5, 3))
    vec_spin = psi.copy()
    for site, axis_idx in reversed(ext_legs):
        # sigma_s^alpha = i b_s^alpha c_s = i c_s^0 c_s^alpha (both bases identical)
        sigma_op = 1.0j * (gamma_24[4 * site + axis_idx] @ gamma_24[4 * site])
        vec_spin = sigma_op @ vec_spin

    # Apply +prod_{b in partial p} u_{b, A->B} where u_{A->B}^alpha = i b_A^alpha b_B^alpha
    directed_bonds = (
        (0, 1, 3),  # b1: A1 -> B2 (z)
        (2, 1, 1),  # b2: A3 -> B2 (x)
        (2, 3, 2),  # b3: A3 -> B4 (y)
        (4, 3, 3),  # b4: A5 -> B4 (z)
        (4, 5, 1),  # b5: A5 -> B6 (x)
        (0, 5, 2),  # b6: A1 -> B6 (y)
    )
    vec_links = psi.copy()
    for s_a, s_b, axis_idx in reversed(directed_bonds):
        u_op = 1.0j * (gamma_24[4 * s_a + axis_idx] @ gamma_24[4 * s_b + axis_idx])
        vec_links = u_op @ vec_links

    assert np.allclose(vec_spin, vec_links, atol=TOL, rtol=0.0)
    # Non-vacuity: W_p^{spin} is a nontrivial involution (W_p^2 = I, W_p != I on physical subspace)
    vec_spin_sq = vec_spin.copy()
    for site, axis_idx in reversed(ext_legs):
        vec_spin_sq = (1.0j * (gamma_24[4 * site + axis_idx] @ gamma_24[4 * site])) @ vec_spin_sq
    assert np.allclose(vec_spin_sq, psi, atol=TOL, rtol=0.0)
    assert float(np.linalg.norm(vec_spin - psi)) > 1.0e-1
    assert eye_s.shape == (dim, dim)

    # Cross-check compute_plaquette_wilson_loops_from_u vs observables.compute_plaquette_fluxes & Pfaffian
    for cell in (
        make_honeycomb_2site(),
        make_honeycomb_4site_rectangular(),
        make_honeycomb_supercell(2, 2),
    ):
        state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5248657458506)
        # Flip one bond to create a non-trivial vison pair in 4-site and 2x2 cells
        u_mod = state.U.copy()
        u_mod[0, cell.bonds[0].bond_type, cell.bonds[0].bond_type] *= -1.0
        flipped_state = MeanFieldState(U=u_mod, V=state.V, lambdas=state.lambdas)

        scalars = extract_pure_kitaev_bond_scalars(flipped_state, cell)
        wp_bridge = compute_plaquette_wilson_loops_from_u(scalars.u_gauge, cell)
        wp_obs = compute_plaquette_fluxes(flipped_state, cell)
        assert np.allclose(wp_bridge, wp_obs, atol=TOL, rtol=0.0)

    # Verify RealSpaceTorusCluster Pfaffian flux matches W_p = prod u_b on L=4 torus with vison pair
    cluster = RealSpaceTorusCluster(L=4)
    base_inhom = InhomogeneousMeanFieldState.from_pure_kitaev(
        cluster.num_cells, u_gauge=1.0, u_matter=-0.5248657458506
    )
    vortex_inhom, _ = cluster.create_unpinned_vortex_pair_seed(base_inhom, (0, 0), (1, 0))
    # Build exact pure-gauge covariance block matching vortex_inhom links
    omega = np.zeros((cluster.num_majoranas, cluster.num_majoranas), dtype=np.float64)
    disp = [(0, 0, 3), (0, -1, 1), (-1, 0, 2)]
    for m2 in range(cluster.L):
        for m1 in range(cluster.L):
            c = m1 + cluster.L * m2
            s_a = cluster.site_index(m1, m2, 0)
            for b_idx, (d1, d2, gamma) in enumerate(disp):
                s_b = cluster.site_index(m1 + d1, m2 + d2, 1)
                u_val = float(vortex_inhom.U[c, b_idx, gamma, gamma])
                omega[4 * s_a + gamma, 4 * s_b + gamma] = u_val
                omega[4 * s_b + gamma, 4 * s_a + gamma] = -u_val
    pf_fluxes = cluster.compute_all_plaquette_pfaffian_fluxes(omega)
    prod_fluxes = cluster.compute_plaquette_fluxes(omega)
    for key in prod_fluxes:
        assert np.isclose(pf_fluxes[key], prod_fluxes[key], atol=TOL, rtol=0.0)
    assert np.isclose(pf_fluxes[(0, 0)], -1.0, atol=TOL, rtol=0.0)
    assert np.isclose(pf_fluxes[(1, 0)], -1.0, atol=TOL, rtol=0.0)
    assert sum(1 for val in pf_fluxes.values() if val < 0.0) == 2


def test_state_port_roundtrip_bond_scalars_and_qsl_coupling_contracts() -> None:
    """Verify lossless state round-trips, (u_ij, v_ij) extraction, protocol compliance, and contracts."""
    bridge = DefaultQSLConventionBridge()
    assert isinstance(bridge, QSLConventionBridgePort)
    assert bridge.conventions_version == CONVENTIONS_VERSION == "conventions-v2.0"

    cell = make_honeycomb_4site_rectangular()
    rng = np.random.default_rng(1337)
    u_rand = rng.normal(size=(len(cell.bonds), 4, 4))
    v_raw = rng.normal(size=(cell.num_sublattices, 4, 4))
    v_skew = 0.5 * (v_raw - np.swapaxes(v_raw, -1, -2))
    lam_rand = rng.normal(size=(cell.num_sublattices, 3))
    state = MeanFieldState(U=u_rand, V=v_skew, lambdas=lam_rand)

    port_state = bridge.to_v2_state(state)
    assert isinstance(port_state, QSLMeanFieldPortState)
    recovered = bridge.from_v2_state(port_state)
    recovered_fn = convert_state_v2_to_kmf(convert_state_kmf_to_v2(state))
    assert np.allclose(recovered.U, state.U, atol=TOL, rtol=0.0)
    assert np.allclose(recovered.V, state.V, atol=TOL, rtol=0.0)
    assert np.allclose(recovered.lambdas, state.lambdas, atol=TOL, rtol=0.0)
    assert np.allclose(recovered_fn.U, state.U, atol=TOL, rtol=0.0)
    assert np.allclose(recovered_fn.V, state.V, atol=TOL, rtol=0.0)
    assert np.allclose(recovered_fn.lambdas, state.lambdas, atol=TOL, rtol=0.0)
    site0_bridge = bridge.extract_site_su2(state.V[0], basis="kitaev_mf")
    site0_v2 = bridge.extract_site_su2(port_state.V_v2[0], basis="conventions_v2")
    assert np.allclose(site0_bridge.spin_s, site0_v2.spin_s, atol=TOL, rtol=0.0)
    assert np.allclose(site0_bridge.gauge_g, site0_v2.gauge_g, atol=TOL, rtol=0.0)

    # Pure-Kitaev (u_ij, v_ij) extraction identical from kmf and v2 states
    sc_kmf = bridge.extract_bond_scalars(state, cell)
    sc_v2 = bridge.extract_bond_scalars(port_state, cell)
    assert np.allclose(sc_kmf.u_gauge, sc_v2.u_gauge, atol=TOL, rtol=0.0)
    assert np.allclose(sc_kmf.v_matter, sc_v2.v_matter, atol=TOL, rtol=0.0)
    assert sc_kmf.bond_axes == sc_v2.bond_axes
    assert np.isclose(sc_kmf.off_diagonal_norm, sc_v2.off_diagonal_norm, atol=TOL, rtol=0.0)
    assert sc_kmf.off_diagonal_norm > 0.1

    pure_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5248657458506)
    pure_sc = bridge.extract_bond_scalars(pure_state, cell)
    assert np.allclose(pure_sc.u_gauge, 1.0, atol=TOL, rtol=0.0)
    assert np.allclose(pure_sc.v_matter, -0.5248657458506, atol=TOL, rtol=0.0)
    assert np.isclose(pure_sc.off_diagonal_norm, 0.0, atol=TOL, rtol=0.0)
    u_rebuilt = reconstruct_pure_kitaev_bond_tensors(pure_sc.u_gauge, pure_sc.v_matter, cell)
    assert np.allclose(u_rebuilt, pure_state.U, atol=TOL, rtol=0.0)
    assert np.allclose(
        bridge.compute_plaquette_wilson_loops(pure_sc.u_gauge, cell), 1.0, atol=TOL, rtol=0.0
    )

    # Inhomogeneous state lossless round-trip
    inhom = InhomogeneousMeanFieldState(
        U=rng.normal(size=(4, 3, 4, 4)),
        V=0.5 * ((v_inh := rng.normal(size=(4, 2, 4, 4))) - np.swapaxes(v_inh, -1, -2)),
        lambdas=rng.normal(size=(4, 2, 3)),
    )
    inhom_v2 = convert_inhomogeneous_state_kmf_to_v2(inhom)
    assert isinstance(inhom_v2, QSLInhomogeneousPortState)
    inhom_back = convert_inhomogeneous_state_v2_to_kmf(inhom_v2)
    assert np.allclose(inhom_back.U, inhom.U, atol=TOL, rtol=0.0)
    assert np.allclose(inhom_back.V, inhom.V, atol=TOL, rtol=0.0)
    assert np.allclose(inhom_back.lambdas, inhom.lambdas, atol=TOL, rtol=0.0)

    # Cross-QSL chiral kappa and gap contracts (contracts/qsl_conventions.json)
    for kappa_kwg in (0.05, 0.10, 0.20, -0.15):
        kappa_static = chiral_kappa_kwg_to_kmf_static(kappa_kwg)
        assert np.isclose(kappa_static, 0.5 * kappa_kwg, atol=TOL, rtol=0.0)
        assert np.isclose(
            chiral_kappa_kmf_static_to_kwg(kappa_static), kappa_kwg, atol=TOL, rtol=0.0
        )
        kappa_kslsc = chiral_kappa_kmf_static_to_kslsc(kappa_static)
        assert np.isclose(kappa_kslsc, kappa_static, atol=TOL, rtol=0.0)
        delta_can_mf = chiral_full_gap_canonical(kappa_static)
        delta_k_kwg = 6.0 * sqrt(3.0) * abs(kappa_kwg)
        delta_bdg_kslsc = 12.0 * sqrt(3.0) * abs(kappa_kslsc)
        assert np.isclose(delta_can_mf, delta_k_kwg, atol=TOL, rtol=0.0)
        assert np.isclose(delta_can_mf, delta_bdg_kslsc, atol=TOL, rtol=0.0)

    # Explicit port adapter validation guards
    with pytest.raises(ValueError):
        conjugate_site_block_p0(np.zeros((3, 3)))
    with pytest.raises(ValueError):
        conjugate_multisite_p0(np.zeros((6, 6)))
    with pytest.raises(ValueError):
        lambdas_kmf_to_v2(np.zeros((2, 4)))
    with pytest.raises(ValueError):
        extract_site_su2_expectations(np.eye(4))  # symmetric, not skew-symmetric
    with pytest.raises(ValueError):
        kitaev_coupling_pauli_to_bilayer(float("nan"))
    with pytest.raises(ValueError):
        chiral_kappa_kwg_to_kmf_static(float("inf"))
    with pytest.raises(ValueError):
        extract_pure_kitaev_bond_scalars(
            MeanFieldState.from_pure_kitaev(make_honeycomb_2site()), cell
        )


def test_port_state_finiteness_version_and_canonical_ab_bond_guards() -> None:
    """Verify UOW-07 F4(1) finiteness, conventions_version, and canonical A->B bond guards."""
    u_ok = np.zeros((3, 4, 4), dtype=np.float64)
    v_ok = np.zeros((2, 4, 4), dtype=np.float64)
    lam_ok = np.zeros((2, 3), dtype=np.float64)

    # Reject non-canonical conventions_version on both port states
    with pytest.raises(ValueError, match="conventions_version"):
        QSLMeanFieldPortState(
            U_v2=u_ok, V_v2=v_ok, lambdas_v2=lam_ok, conventions_version="conventions-v1.0"
        )
    u_inh = np.zeros((4, 3, 4, 4), dtype=np.float64)
    v_inh = np.zeros((4, 2, 4, 4), dtype=np.float64)
    lam_inh = np.zeros((4, 2, 3), dtype=np.float64)
    with pytest.raises(ValueError, match="conventions_version"):
        QSLInhomogeneousPortState(
            U_v2=u_inh, V_v2=v_inh, lambdas_v2=lam_inh, conventions_version="legacy-v0"
        )

    # Reject NaN / inf in U_v2, V_v2, lambdas_v2 (previously NaN > 1e-12 evaluated to False)
    for bad_val in (float("nan"), float("inf"), float("-inf")):
        u_bad = u_ok.copy()
        u_bad[0, 0, 0] = bad_val
        with pytest.raises(ValueError, match="finite"):
            QSLMeanFieldPortState(U_v2=u_bad, V_v2=v_ok, lambdas_v2=lam_ok)

        v_bad = v_ok.copy()
        v_bad[0, 0, 1] = bad_val
        v_bad[0, 1, 0] = -bad_val
        with pytest.raises(ValueError, match="finite"):
            QSLMeanFieldPortState(U_v2=u_ok, V_v2=v_bad, lambdas_v2=lam_ok)

        lam_bad = lam_ok.copy()
        lam_bad[0, 0] = bad_val
        with pytest.raises(ValueError, match="finite"):
            QSLMeanFieldPortState(U_v2=u_ok, V_v2=v_ok, lambdas_v2=lam_bad)

        ui_bad = u_inh.copy()
        ui_bad[0, 0, 0, 0] = bad_val
        with pytest.raises(ValueError, match="finite"):
            QSLInhomogeneousPortState(U_v2=ui_bad, V_v2=v_inh, lambdas_v2=lam_inh)

        vi_bad = v_inh.copy()
        vi_bad[0, 0, 0, 1] = bad_val
        vi_bad[0, 0, 1, 0] = -bad_val
        with pytest.raises(ValueError, match="finite"):
            QSLInhomogeneousPortState(U_v2=u_inh, V_v2=vi_bad, lambdas_v2=lam_inh)

        li_bad = lam_inh.copy()
        li_bad[0, 0, 0] = bad_val
        with pytest.raises(ValueError, match="finite"):
            QSLInhomogeneousPortState(U_v2=u_inh, V_v2=v_inh, lambdas_v2=li_bad)

        v_site_bad = np.zeros((4, 4), dtype=np.float64)
        v_site_bad[0, 1] = bad_val
        v_site_bad[1, 0] = -bad_val
        with pytest.raises(ValueError, match="finite"):
            extract_site_su2_expectations(v_site_bad)

    # Reject non-finite u_gauge and reversed B->A bond orientation (UOW-02 R4 Lakatosian probe)
    cell4 = make_honeycomb_4site_rectangular()
    u_nan = np.ones(len(cell4.bonds), dtype=np.float64)
    u_nan[0] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        compute_plaquette_wilson_loops_from_u(u_nan, cell4)

    b0 = cell4.bonds[0]
    reversed_bonds = (
        Bond(
            sublattice_i=b0.sublattice_j,
            sublattice_j=b0.sublattice_i,
            cell_offset=(-b0.cell_offset[0], -b0.cell_offset[1]),
            bond_type=b0.bond_type,
        ),
        *cell4.bonds[1:],
    )
    reversed_cell = UnitCell(
        name="reversed_bond0_4site",
        num_sublattices=cell4.num_sublattices,
        a1=cell4.a1,
        a2=cell4.a2,
        basis_positions=cell4.basis_positions,
        bonds=reversed_bonds,
        plaquettes=cell4.plaquettes,
    )
    with pytest.raises(ValueError, match="A -> B"):
        compute_plaquette_wilson_loops_from_u(
            np.ones(len(cell4.bonds), dtype=np.float64), reversed_cell
        )

    # Reject unbalanced plaquette oriented_bonds
    bad_plaq = Plaquette(
        name="p_unbalanced",
        bond_indices=cell4.plaquettes[0].bond_indices,
        oriented_bonds=((0, 1), (3, 1), (1, 1), (2, -1), (4, 1), (5, -1)),
    )
    unbalanced_cell = UnitCell(
        name="unbalanced_plaq_4site",
        num_sublattices=cell4.num_sublattices,
        a1=cell4.a1,
        a2=cell4.a2,
        basis_positions=cell4.basis_positions,
        bonds=cell4.bonds,
        plaquettes=(bad_plaq,),
    )
    with pytest.raises(ValueError, match="oriented_bonds"):
        compute_plaquette_wilson_loops_from_u(
            np.ones(len(cell4.bonds), dtype=np.float64), unbalanced_cell
        )
