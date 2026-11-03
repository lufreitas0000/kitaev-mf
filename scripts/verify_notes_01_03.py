"""Independent numerical verification of the statements in Notes 01-03 (revised).

Self-contained (NumPy only). Run with `pytest -q scripts/verify_notes_01_03.py` or
`python scripts/verify_notes_01_03.py`.
Every test corresponds to a numbered statement in the revised notes.
Conventions: flavor index mu in {0,1,2,3}; bond type gamma in {1,2,3} <-> 0-based g in {0,1,2}.
"""

import sys

import numpy as np

TOL = 1e-12

I2 = np.eye(2)
sx = np.array([[0.0, 1.0], [1.0, 0.0]])
sz = np.diag([1.0, -1.0])
ep = np.array([[0.0, 1.0], [-1.0, 0.0]])
e11 = np.diag([1.0, 0.0])

G = [np.kron(I2, ep), np.kron(ep, sz), np.kron(ep, sx)]
M = [np.kron(sz, ep), np.kron(ep, I2), np.kron(sx, ep)]
N = [np.kron(e11, ep), np.kron(ep, e11), 0.5 * (np.kron(sx, ep) + np.kron(ep, sx))]

LC = {
    (0, 1, 2): 1,
    (1, 2, 0): 1,
    (2, 0, 1): 1,
    (0, 2, 1): -1,
    (2, 1, 0): -1,
    (1, 0, 2): -1,
}


def com(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Matrix commutator [A, B] = A @ B - B @ A."""
    return A @ B - B @ A


# ---------- Jordan-Wigner Majoranas
_X = np.array([[0, 1], [1, 0j]])
_Y = np.array([[0, -1j], [1j, 0]])
_Z = np.diag([1, -1]).astype(complex)
_I = np.eye(2, dtype=complex)


def majoranas(nq: int) -> list[np.ndarray]:
    out = []
    for m in range(nq):
        for P in (_X, _Y):
            op = np.array([[1]], dtype=complex)
            for idx_l in range(nq):
                op = np.kron(op, _Z if idx_l < m else (P if idx_l == m else _I))
            out.append(op)
    return out


def op_quad(c: list[np.ndarray], A: np.ndarray, norm: float) -> np.ndarray:
    """Quadratic Majorana form (i / norm) * c^T A c."""
    n = A.shape[0]
    return 1j / norm * sum(A[a, b] * c[a] @ c[b] for a in range(n) for b in range(n))


# ====================================================== Note 01
def test_matrices_and_index_formulas() -> None:
    """Def. 2.2 / 4.1: explicit Kronecker forms equal the index formulas."""
    for g in (1, 2, 3):
        Nf = np.zeros((4, 4))
        Nf[0, g] = 1
        Nf[g, 0] = -1
        a, b = [(2, 3), (3, 1), (1, 2)][g - 1]
        Gf = Nf.copy()
        Gf[a, b] += 1
        Gf[b, a] -= 1
        assert np.allclose(N[g - 1], Nf)
        assert np.allclose(G[g - 1], Gf)
        assert np.allclose(N[g - 1], (M[g - 1] + G[g - 1]) / 2)


def test_so4_algebra_and_traces() -> None:
    """Prop. 2.3 (matrix level)."""
    for (a, b, c_idx), s in LC.items():
        assert np.allclose(com(G[a], G[b]), -2 * s * G[c_idx])
        assert np.allclose(com(M[a], M[b]), 2 * s * M[c_idx])
    assert all(np.allclose(com(M[a], G[b]), 0) for a in range(3) for b in range(3))
    assert all(
        np.isclose(np.trace(G[a].T @ G[b]), 4 * (a == b))
        and np.isclose(np.trace(M[a].T @ M[b]), 4 * (a == b))
        and np.isclose(np.trace(N[a].T @ N[b]), 2 * (a == b))
        and np.isclose(np.trace(M[a].T @ G[b]), 0)
        for a in range(3)
        for b in range(3)
    )


def test_lemma_2_1_corrected() -> None:
    """Lemma 2.1: [O_A, O_B] = i O_[A,B]."""
    c = majoranas(2)
    A = G[0]
    B = M[0] + G[1]
    lhs = op_quad(c, A, 4) @ op_quad(c, B, 4) - op_quad(c, B, 4) @ op_quad(c, A, 4)
    assert np.allclose(lhs, 1j * op_quad(c, com(A, B), 4))
    assert not np.allclose(lhs, op_quad(c, com(A, B), 4))


def test_casimir_parity_and_pauli() -> None:
    """Prop. 2.3 (operator level), Thm 3.1, Prop. 3.4."""
    c = majoranas(2)
    D = c[0] @ c[1] @ c[2] @ c[3]
    Id = np.eye(4)
    S = [op_quad(c, m, 8) for m in M]
    T = [op_quad(c, g, 8) for g in G]

    for (a, b, cc), s in LC.items():
        assert np.allclose(com(S[a], S[b]), 1j * s * S[cc])
        assert np.allclose(com(T[a], T[b]), -1j * s * T[cc])
    assert all(np.allclose(com(S[a], T[b]), 0) for a in range(3) for b in range(3))
    assert np.allclose(sum(t @ t for t in T), 3 / 8 * (Id - D))
    assert np.allclose(sum(s @ s for s in S), 3 / 8 * (Id + D))

    sg = [1j * c[0] @ c[g] for g in (1, 2, 3)]
    assert all(
        np.allclose(sg[g], 2 * (S[g] + T[g])) and np.allclose(sg[g], op_quad(c, N[g], 2))
        for g in range(3)
    )
    for (a, b, cc), s in LC.items():
        assert np.allclose(sg[a] @ sg[b], 1j * s * sg[cc] - 4j * s * T[cc])

    w, v = np.linalg.eigh(D)
    P = v[:, w > 0]
    s1 = [P.conj().T @ s @ P for s in sg]
    assert np.allclose(s1[0] @ s1[1], 1j * s1[2])


def _quad(c: list[np.ndarray], A: np.ndarray, i: int, j: int) -> np.ndarray:
    return sum(A[a, b] * c[4 * i + a] @ c[4 * j + b] for a in range(4) for b in range(4))


def test_wick_and_hfb_two_site() -> None:
    """Lemma 5.3, Prop. 5.5: Wick formula, <H_MF>=E, and tangent property."""
    rng = np.random.default_rng(1)
    c = majoranas(4)
    Xr = rng.normal(size=(8, 8))
    Xr = Xr - Xr.T
    Hr = 1j / 4 * sum(Xr[a, b] * c[a] @ c[b] for a in range(8) for b in range(8))
    psi = np.linalg.eigh(Hr)[1][:, 0]

    def ev(op: np.ndarray) -> complex:
        return psi.conj() @ op @ psi

    Cov = np.array([[ev(1j * c[a] @ c[b]) for b in range(8)] for a in range(8)])
    U = Cov[:4, 4:].real
    V = [(Cov[:4, :4] - 1j * np.eye(4)).real, (Cov[4:, 4:] - 1j * np.eye(4)).real]
    assert np.allclose(V[0], -V[0].T)

    J = rng.normal(size=(3, 3))
    h = rng.normal(size=3)
    lam = [rng.normal(size=3), rng.normal(size=3)]

    def sig(s: int, g: int) -> np.ndarray:
        return 0.5j * sum(
            N[g][a, b] * c[4 * s + a] @ c[4 * s + b] for a in range(4) for b in range(4)
        )

    tr = np.trace

    def Q(s: int, g: int) -> np.ndarray:
        return 2 * sig(s, g)

    for a in range(3):
        for b in range(3):
            val = ev(Q(0, a) @ Q(1, b))
            expected = tr(V[0].T @ N[a]) * tr(V[1].T @ N[b]) + 2 * tr(U.T @ N[a] @ U @ N[b])
            assert abs(val - expected) < 1e-12
    assert all(abs(ev(sig(0, g)) - V[0][0, g + 1]) < 1e-12 for g in range(3))

    def E(U_mat: np.ndarray, V_list: list[np.ndarray]) -> float:
        e = 0.0
        for Jo, Uo, Vi, Vj in ((J, U_mat, V_list[0], V_list[1]), (J.T, -U_mat.T, V_list[1], V_list[0])):
            e += sum(
                Jo[a, b] * (tr(Vi.T @ N[a]) * tr(Vj.T @ N[b]) + 2 * tr(Uo.T @ N[a] @ Uo @ N[b]))
                for a in range(3)
                for b in range(3)
            ) / 8
        for i in range(2):
            e += sum(
                -0.5 * h[g] * tr(V_list[i].T @ N[g]) + 0.25 * lam[i][g] * tr(G[g].T @ V_list[i])
                for g in range(3)
            )
        return float(e)

    def blocks(U_mat: np.ndarray, V_list: list[np.ndarray]) -> tuple[np.ndarray, list[np.ndarray]]:
        Aij = 2 * sum(J[a, b] * N[a] @ U_mat @ N[b] for a in range(3) for b in range(3))
        Bs = [
            sum(J[a, b] * N[a] * tr(V_list[1].T @ N[b]) for a in range(3) for b in range(3)),
            sum(J[b, a] * N[a] * tr(V_list[0].T @ N[b]) for a in range(3) for b in range(3)),
        ]
        for i in range(2):
            Bs[i] = Bs[i] + sum(lam[i][g] * G[g] - 2 * h[g] * N[g] for g in range(3))
        return Aij, Bs

    def Lin(Up: np.ndarray, Vp: list[np.ndarray], U_mat: np.ndarray, V_list: list[np.ndarray]) -> float:
        Aij, Bs = blocks(U_mat, V_list)
        return float(0.25 * (2 * tr(Aij.T @ Up) + sum(tr(Bs[i].T @ Vp[i]) for i in range(2))))

    Cc = Lin(U, V, U, V) - E(U, V)
    C_closed = (
        sum(
            Jo[a, b] * (tr(Vi.T @ N[a]) * tr(Vj.T @ N[b]) + 2 * tr(Uo.T @ N[a] @ Uo @ N[b]))
            for (Jo, Uo, Vi, Vj) in ((J, U, V[0], V[1]), (J.T, -U.T, V[1], V[0]))
            for a in range(3)
            for b in range(3)
        )
        / 8
    )
    Aij, Bs = blocks(U, V)
    HMF = (1j / 4) * (
        _quad(c, Aij, 0, 1)
        + _quad(c, -Aij.T, 1, 0)
        + _quad(c, Bs[0], 0, 0)
        + _quad(c, Bs[1], 1, 1)
    )
    assert abs(ev(HMF).real - C_closed - E_exact(ev, c, J, h, lam, sig)) < 1e-10
    assert abs(Cc - C_closed) < 1e-10

    Up = rng.normal(size=(4, 4))
    Vp = [x - x.T for x in rng.normal(size=(2, 4, 4))]
    t = 1e-3
    dE = (
        E(U + t * (Up - U), [V[i] + t * (Vp[i] - V[i]) for i in range(2)])
        - E(U - t * (Up - U), [V[i] - t * (Vp[i] - V[i]) for i in range(2)])
    ) / (2 * t)
    assert abs(Lin(Up, Vp, U, V) - Cc - (E(U, V) + dE)) < 1e-8


def E_exact(
    ev: callable,  # type: ignore[valid-type]
    c: list[np.ndarray],
    J: np.ndarray,
    h: np.ndarray,
    lam: list[np.ndarray],
    sig: callable,  # type: ignore[valid-type]
) -> float:
    Hb = sum(J[a, b] * sig(0, a) @ sig(1, b) for a in range(3) for b in range(3))
    Hz = -sum(h[g] * (sig(0, g) + sig(1, g)) for g in range(3))
    Hl = sum(
        0.25 * lam[s][g] * 1j * sum(
            G[g][a, b] * c[4 * s + a] @ c[4 * s + b] for a in range(4) for b in range(4)
        )
        for s in range(2)
        for g in range(3)
    )
    return float(ev(Hb + Hz + Hl).real)


# ====================================================== Notes 02/03: Bloch machinery
def Jmat(
    g: int,
    Jv: float = 0.0,
    K: float = 0.0,
    Gm: float = 0.0,
    Gp: float = 0.0,
    D: float = 0.0,
) -> np.ndarray:
    a, b = (g + 1) % 3, (g + 2) % 3
    Jm = Jv * np.eye(3)
    Jm[g, g] += K
    Jm[a, b] += Gm + D
    Jm[b, a] += Gm - D
    for x in (a, b):
        Jm[x, g] += Gp
        Jm[g, x] += Gp
    return Jm


BONDS = [
    (0, 1, 0, -1, 0),
    (0, 1, -1, 0, 1),
    (0, 1, 0, 0, 2),
]  # (s_i, s_j, n1, n2, gamma-1)


def blocks_k(
    bonds: list[tuple[int, int, int, int, int]],
    nsub: int,
    Ub: list[np.ndarray],
    Vs: list[np.ndarray],
    Jb: list[np.ndarray],
    lam: list[np.ndarray],
    h: np.ndarray,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    A = [
        2 * sum(Jb[b][a, bb] * N[a] @ Ub[b] @ N[bb] for a in range(3) for bb in range(3))
        for b in range(len(bonds))
    ]
    B = [np.zeros((4, 4)) for _ in range(nsub)]
    for ib, (si, sj, n1, n2, g) in enumerate(bonds):
        B[si] += sum(
            Jb[ib][a, bb] * N[a] * np.trace(Vs[sj].T @ N[bb]) for a in range(3) for bb in range(3)
        )
        B[sj] += sum(
            Jb[ib][bb, a] * N[a] * np.trace(Vs[si].T @ N[bb]) for a in range(3) for bb in range(3)
        )
    for s in range(nsub):
        B[s] += sum(lam[s][g] * G[g] - 2 * h[g] * N[g] for g in range(3))
    return A, B


def Hk(
    bonds: list[tuple[int, int, int, int, int]],
    nsub: int,
    A: list[np.ndarray],
    B: list[np.ndarray],
    k1: float,
    k2: float,
) -> np.ndarray:
    Mk = np.zeros((4 * nsub, 4 * nsub), complex)
    for s in range(nsub):
        Mk[4 * s : 4 * s + 4, 4 * s : 4 * s + 4] += B[s]
    for ib, (si, sj, n1, n2, g) in enumerate(bonds):
        ph = np.exp(2j * np.pi * (k1 * n1 + k2 * n2))
        Mk[4 * si : 4 * si + 4, 4 * sj : 4 * sj + 4] += A[ib] * ph
        Mk[4 * sj : 4 * sj + 4, 4 * si : 4 * si + 4] -= A[ib].T * np.conj(ph)
    return 0.25j * Mk


def Pk(Hh: np.ndarray) -> np.ndarray:
    """P = 2 Theta(-H) with Theta(0)=1/2 (Note 03, Thm 1.2)."""
    w, v = np.linalg.eigh(Hh)
    wt = np.where(w < -1e-12, 1.0, np.where(w > 1e-12, 0.0, 0.5))
    return 2 * (v * wt) @ v.conj().T


def correlators(
    bonds: list[tuple[int, int, int, int, int]],
    nsub: int,
    A: list[np.ndarray],
    B: list[np.ndarray],
    L1: int,
    L2: int,
) -> tuple[list[np.ndarray], list[np.ndarray]]:
    Ub = [np.zeros((4, 4)) for _ in bonds]
    Vn = [np.zeros((4, 4)) for _ in range(nsub)]
    n = 0
    for m1 in range(L1):
        for m2 in range(L2):
            k1, k2 = m1 / L1, m2 / L2
            P = Pk(Hk(bonds, nsub, A, B, k1, k2))
            n += 1
            for ib, (si, sj, n1, n2, g) in enumerate(bonds):
                ph = 2 * np.pi * (k1 * n1 + k2 * n2)
                blk = P[4 * si : 4 * si + 4, 4 * sj : 4 * sj + 4]
                Ub[ib] += blk.imag * np.cos(ph) - blk.real * np.sin(ph)
            for s in range(nsub):
                Vn[s] += P[4 * s : 4 * s + 4, 4 * s : 4 * s + 4].imag
    return [u / n for u in Ub], [v / n for v in Vn]


def realspace_M(
    bonds: list[tuple[int, int, int, int, int]],
    nsub: int,
    A: list[np.ndarray],
    B: list[np.ndarray],
    L: int,
) -> tuple[np.ndarray, callable]:  # type: ignore[valid-type]
    def idx(s: int, r1: int, r2: int) -> int:
        return (((r1 % L) * L + (r2 % L)) * nsub + s) * 4

    Mr = np.zeros((4 * nsub * L * L, 4 * nsub * L * L))
    for r1 in range(L):
        for r2 in range(L):
            for s in range(nsub):
                i = idx(s, r1, r2)
                Mr[i : i + 4, i : i + 4] += B[s]
            for ib, (si, sj, n1, n2, g) in enumerate(bonds):
                i = idx(si, r1, r2)
                j = idx(sj, r1 + n1, r2 + n2)
                Mr[i : i + 4, j : j + 4] += A[ib]
                Mr[j : j + 4, i : i + 4] -= A[ib].T
    return Mr, idx


def generic_model(
    seed: int = 7, field: bool = True
) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray], list[np.ndarray], np.ndarray]:
    rng = np.random.default_rng(seed)
    Jb = [Jmat(g, Jv=-0.05, K=-1.0, Gm=0.2, Gp=0.03, D=0.04) for (*_, g) in BONDS]
    Ub = [np.diag([0.5, 0.9, 0.9, 0.9]) + 0.05 * rng.normal(size=(4, 4)) for _ in BONDS]
    Vr = rng.normal(size=(4, 4))
    Vr = 0.2 * (Vr - Vr.T)
    h = np.array([0.2, 0.1, 0.3]) if field else np.zeros(3)
    return Jb, Ub, [Vr, Vr], [rng.normal(size=3) * 0.1] * 2, h


def test_hermiticity_and_ph_symmetry() -> None:
    """Props. 3.1, 3.2 (Note 02): H(k) Hermitian, H(-k) = -H(k)^*."""
    Jb, Ub, Vs, lam, h = generic_model()
    A, B = blocks_k(BONDS, 2, Ub, Vs, Jb, lam, h)
    H = Hk(BONDS, 2, A, B, 0.17, 0.36)
    assert np.allclose(H, H.conj().T)
    assert np.allclose(Hk(BONDS, 2, A, B, -0.17, -0.36), -H.conj())

    e = np.linalg.eigvalsh(H)
    assert np.allclose(np.linalg.eigvalsh(Hk(BONDS, 2, A, B, -0.17, -0.36)), np.sort(-e))
    assert not np.allclose(e, np.sort(-e), atol=1e-3)

    Jb2 = [Jmat(g, K=-1.0, Gm=0.2) for (*_, g) in BONDS]
    A2, B2 = blocks_k(
        BONDS, 2, Ub, [np.zeros((4, 4))] * 2, Jb2, [np.zeros(3)] * 2, np.zeros(3)
    )
    e2 = np.linalg.eigvalsh(Hk(BONDS, 2, A2, B2, 0.17, 0.36))
    assert np.allclose(e2, np.sort(-e2))


def test_factor_four_and_k_vs_real_space() -> None:
    """Prop. 3.5, Thm 1.2/2.1 of Note 03 vs real space."""
    Jb, Ub, Vs, lam, h = generic_model()
    A, B = blocks_k(BONDS, 2, Ub, Vs, Jb, lam, h)
    L = 8
    Mr, idx = realspace_M(BONDS, 2, A, B, L)
    lam_, W = np.linalg.eigh(1j * Mr)
    ek = np.sort(
        np.concatenate(
            [
                4 * np.linalg.eigvalsh(Hk(BONDS, 2, A, B, m1 / L, m2 / L))
                for m1 in range(L)
                for m2 in range(L)
            ]
        )
    )
    assert np.allclose(np.sort(lam_), ek, atol=1e-10)

    cov = 2 * W[:, lam_ > 0] @ W[:, lam_ > 0].conj().T
    Uk, Vk = correlators(BONDS, 2, A, B, L, L)
    i0 = idx(0, 0, 0)
    assert np.allclose(Vk[0], (1j * cov[i0 : i0 + 4, i0 : i0 + 4] - 1j * np.eye(4)).real, atol=1e-10)
    for ib, (si, sj, n1, n2, g) in enumerate(BONDS):
        j0 = idx(sj, n1, n2)
        assert np.allclose(
            Uk[ib],
            (1j * cov[idx(si, 0, 0) : idx(si, 0, 0) + 4, j0 : j0 + 4]).real,
            atol=1e-10,
        )


def test_ground_state_correlator_is_2Pi_plus_in_real_space() -> None:
    """Thm 1.2 (Note 03): <c c^T> = 2 Theta(+h) for H = (1/4) c^T (i M) c."""
    rng = np.random.default_rng(3)
    n = 8
    Mr = rng.normal(size=(n, n))
    Mr = Mr - Mr.T
    c = majoranas(n // 2)
    H = 1j / 4 * sum(Mr[a, b] * c[a] @ c[b] for a in range(n) for b in range(n))
    w, v = np.linalg.eigh(H)
    ps = v[:, 0]
    cov = np.array([[ps.conj() @ c[a] @ c[b] @ ps for b in range(n)] for a in range(n)])
    lam_, W = np.linalg.eigh(1j * Mr)
    assert np.allclose(cov, 2 * W[:, lam_ > 0] @ W[:, lam_ > 0].conj().T, atol=1e-10)
    assert abs(w[0] + 0.5 * lam_[lam_ > 0].sum()) < 1e-10


def test_pure_kitaev_dispersion_signs_and_energy() -> None:
    """Thm 4.1 / Cor 4.2 (Note 02), Thm 3.1 / Lemma 3.2 / Prop 3.3 (Note 03)."""
    K = -1.0
    Jb = [Jmat(g, K=K) for (*_, g) in BONDS]
    z = [np.zeros((4, 4))] * 2
    zl = [np.zeros(3)] * 2
    h0 = np.zeros(3)

    def Ub(um: float, ug: float) -> list[np.ndarray]:
        out = []
        for *_, g in BONDS:
            U_mat = np.zeros((4, 4))
            U_mat[0, 0] = um
            U_mat[g + 1, g + 1] = ug
            out.append(U_mat)
        return out

    um = 0.5248657
    A, B = blocks_k(BONDS, 2, Ub(um, 1.0), z, Jb, zl, h0)
    for k1, k2 in [(0.0, 0.0), (0.13, 0.41), (0.3, 0.9)]:
        f = 1 + np.exp(-2j * np.pi * k1) + np.exp(-2j * np.pi * k2)
        e = np.sort(np.abs(np.linalg.eigvalsh(Hk(BONDS, 2, A, B, k1, k2))))
        assert np.allclose(e, np.sort([abs(K) * um / 2] * 6 + [abs(K) * abs(f) / 2] * 2))
    assert abs(abs(1 + 1 + 1) * abs(K) / 2 - 1.5) < 1e-15

    for ug, um_in in ((1.0, 0.5248657), (1.0, -0.5248657)):
        A, B = blocks_k(BONDS, 2, Ub(um_in, ug), z, Jb, zl, h0)
        Uc, _ = correlators(BONDS, 2, A, B, 45, 45)
        assert np.isclose(Uc[2][3, 3], np.sign(K * um_in), atol=1e-12)
        assert np.isclose(Uc[2][0, 0], np.sign(K * ug) * 0.52486575, atol=3e-3)

    kk = (np.arange(1500) + 0.5) / 1500
    K1, K2 = np.meshgrid(kk, kk)
    f = 1 + np.exp(-2j * np.pi * K1) + np.exp(-2j * np.pi * K2)
    fbar = abs(f).mean()
    assert abs(fbar - 1.5745972) < 2e-7
    assert abs(fbar / 3 - 0.5248657) < 1e-7

    um_ = -fbar / 3
    C = -3 * K * um_ * 1.0
    e_cell = -2 * (abs(K) * fbar / 2 + 3 * abs(K) * abs(um_) / 2) - C
    assert abs(e_cell + abs(K) * fbar) < 1e-12


def test_unimodular_relabelling_gives_equal_bond_correlators() -> None:
    """Lemma 3.2: f(k) e^{i2pi k2} = f(k1-k2, -k2)."""
    def f(a: float, b: float) -> complex:
        return 1 + np.exp(-2j * np.pi * a) + np.exp(-2j * np.pi * b)

    r = np.random.default_rng(0).random((50, 2))
    assert np.allclose([f(a, b) * np.exp(2j * np.pi * b) for a, b in r], [f(a - b, -b) for a, b in r])
    assert np.allclose([f(a, b) * np.exp(2j * np.pi * a) for a, b in r], [f(-a, b - a) for a, b in r])
    assert abs(f(1 / 3, 2 / 3)) < 1e-14 and abs(f(2 / 3, 1 / 3)) < 1e-14


def test_bz_folding_is_exact_unitary_equivalence() -> None:
    """Thm 5.1 (Note 02): H_sup(q) = W [H(k1) + H(k2)] W^dag."""
    Jb, Ub, Vs, lam, h = generic_model()
    A, B = blocks_k(BONDS, 2, Ub, Vs, Jb, lam, h)
    reps = [(0, 0), (1, 0)]

    def locate(n1: int, n2: int) -> tuple[int, int, int]:
        for ri, (r1, r2) in enumerate(reps):
            d1, d2 = n1 - r1, n2 - r2
            if (d1 + d2) % 2 == 0:
                return ri, (d1 - d2) // 2, (d1 + d2) // 2
        raise RuntimeError("Cell representation not found")

    sb, sA = [], []
    for ri, (r1, r2) in enumerate(reps):
        for ib, (si, sj, n1, n2, g) in enumerate(BONDS):
            rj, m1, m2 = locate(r1 + n1, r2 + n2)
            sb.append((2 * ri + si, 2 * rj + sj, m1, m2, g))
            sA.append(A[ib])
    sB = [B[s] for _ in range(2) for s in range(2)]

    def Hsup(q1: float, q2: float) -> np.ndarray:
        return Hk(sb, 4, sA, sB, q1, q2)

    for q1, q2 in [(0.21, 0.67), (0.0, 0.0), (0.5, 0.25)]:
        k1 = (((q1 + q2) / 2) % 1, ((q2 - q1) / 2) % 1)
        k2 = (((q1 + q2 + 1) / 2) % 1, ((q2 - q1 + 1) / 2) % 1)
        Hd = np.zeros((16, 16), complex)
        Hd[:8, :8] = Hk(BONDS, 2, A, B, *k1)
        Hd[8:, 8:] = Hk(BONDS, 2, A, B, *k2)
        W = np.zeros((16, 16), complex)
        for ri, r in enumerate(reps):
            for ka, kk_ in enumerate((k1, k2)):
                ph = np.exp(2j * np.pi * (kk_[0] * r[0] + kk_[1] * r[1])) / np.sqrt(2)
                for s in range(2):
                    for mu in range(4):
                        W[4 * (2 * ri + s) + mu, 8 * ka + 4 * s + mu] = ph
        assert np.allclose(W @ W.conj().T, np.eye(16), atol=1e-13)
        assert np.allclose(Hsup(q1, q2), W @ Hd @ W.conj().T, atol=1e-13)


def test_dirac_points_and_geometry() -> None:
    """Def. 1.1 (Note 02): geometry and rectangular supercell."""
    s3 = np.sqrt(3)
    dx = np.array([0.5, 1 / (2 * s3)])
    dy = np.array([-0.5, 1 / (2 * s3)])
    dz = np.array([0.0, -1 / s3])
    a1 = dz - dy
    a2 = dz - dx
    assert np.allclose(dx - dz, -a2)
    assert np.allclose(dy - dz, -a1)
    assert np.isclose(a1 @ a2, 0.5)
    assert abs((a1 - a2) @ (a1 + a2)) < 1e-15
    assert np.isclose(np.linalg.norm(a1 - a2), 1)
    assert np.isclose(np.linalg.norm(a1 + a2), s3)
    b1 = 2 * np.pi * np.array([1.0, -1 / s3])
    b2 = 2 * np.pi * np.array([-1.0, -1 / s3])
    assert np.allclose([[a1 @ b1, a1 @ b2], [a2 @ b1, a2 @ b2]], 2 * np.pi * np.eye(2))


# ====================================================== Note 03, Sec. 4
def test_gauge_residual_sign_and_lambda_update() -> None:
    """Lemma 4.2: <i c^T G c> = -g, <T> = -g/8; dE0/dlambda = -g/4."""
    rng = np.random.default_rng(5)
    c = majoranas(2)
    X = rng.normal(size=(4, 4))
    X = X - X.T
    H = 1j / 4 * sum(X[a, b] * c[a] @ c[b] for a in range(4) for b in range(4))
    psi = np.linalg.eigh(H)[1][:, 0]

    def ev(op: np.ndarray) -> float:
        return float((psi.conj() @ op @ psi).real)

    V = np.array([[ev(1j * c[a] @ c[b]) for b in range(4)] for a in range(4)])
    for g in range(3):
        gg = np.trace(G[g] @ V)
        quad_val = ev(1j * sum(G[g][a, b] * c[a] @ c[b] for a in range(4) for b in range(4)))
        assert abs(quad_val + gg) < 1e-12
        assert abs(ev(op_quad(c, G[g], 8)) + gg / 8) < 1e-12


if __name__ == "__main__":
    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("PASS", name)
            except Exception as exc:
                fails += 1
                print("FAIL", name, repr(exc))
    sys.exit(fails)
