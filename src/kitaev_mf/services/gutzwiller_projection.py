"""Wave 1 HC-02: Exact and 50-seed VMC Gutzwiller projection via Thouless Pfaffians."""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt

import numpy as np
import torch
from numpy.typing import NDArray
from scipy.special import ndtr, ndtri

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.hamiltonian import build_bond_hopping_matrices, build_on_site_matrices


@dataclass(frozen=True, slots=True)
class ExactGutzwillerResult:
    """Exact 2^N Gutzwiller-projected observables on a finite periodic honeycomb torus."""

    nx: int
    ny: int
    n_sites: int
    h: float
    bc_twist: tuple[int, int]
    ref_condition_number: float
    e_proj_per_site: float
    magnetization_111: float
    plaquette_flux: float
    nn_bond_correlator: float


@dataclass(frozen=True, slots=True)
class VMCGutzwillerResult:
    """50-chain Metropolis-Hastings Gutzwiller VMC result with BCa CIs and n_eff."""

    nx: int
    ny: int
    n_sites: int
    h: float
    bc_twist: tuple[int, int]
    n_chains: int
    samples_per_chain: int
    e_proj_mean: float
    e_proj_ci_low: float
    e_proj_ci_high: float
    e_loc_std: float
    tau_int_energy: float
    n_eff_total: float
    magnetization_111_mean: float
    magnetization_111_ci_low: float
    magnetization_111_ci_high: float


def bca_confidence_interval(
    values: NDArray[np.float64],
    *,
    confidence: float = 0.95,
    n_boot: int = 2000,
    seed: int = 12345,
) -> tuple[float, float]:
    """Compute non-parametric BCa confidence interval with degenerate-array guard."""
    arr = np.asarray(values, dtype=np.float64)
    n = arr.size
    if n < 2 or float(np.ptp(arr)) <= 1e-13 * max(1.0, float(np.abs(np.mean(arr)))):
        m = float(np.mean(arr))
        return m, m

    alpha = 0.5 * (1.0 - confidence)
    theta_hat = float(np.mean(arr))
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, n, size=(n_boot, n))
    boot_means = np.mean(arr[idx], axis=1)

    prop_less = float(np.mean(boot_means < theta_hat))
    prop_less = min(max(prop_less, 0.5 / n_boot), 1.0 - 0.5 / n_boot)
    z0 = float(ndtri(prop_less))

    total = float(np.sum(arr))
    jack = (total - arr) / float(n - 1)
    jack_mean = float(np.mean(jack))
    diff = jack_mean - jack
    denom = float(np.sum(diff * diff)) ** 1.5
    acc = float(np.sum(diff**3)) / (6.0 * denom) if denom > 1e-30 else 0.0

    def _adj_q(z_a: float) -> float:
        numer = z0 + z_a
        denom_q = 1.0 - acc * numer
        if abs(denom_q) < 1e-12:
            return alpha if z_a < 0.0 else 1.0 - alpha
        return float(ndtr(z0 + numer / denom_q))

    q_low = float(np.clip(_adj_q(float(ndtri(alpha))), 0.0, 1.0))
    q_high = float(np.clip(_adj_q(float(ndtri(1.0 - alpha))), 0.0, 1.0))
    lo, hi = np.quantile(boot_means, [q_low, q_high])
    return float(lo), float(hi)


def complex_pfaffian(matrix: NDArray[np.complex128], tol: float = 1e-14) -> complex:
    """Compute the complex Pfaffian Pf(M) via skew-symmetric Parlett-Reid pivoting."""
    n = matrix.shape[0]
    if n == 0:
        return 1.0 + 0.0j
    if n % 2 != 0:
        return 0.0 + 0.0j
    m = matrix.astype(np.complex128, copy=True)
    pf = 1.0 + 0.0j
    for k in range(0, n, 2):
        col = int(np.argmax(np.abs(m[k + 1 :, k]))) + k + 1
        pivot = m[k, col]
        if abs(pivot) < tol:
            return 0.0 + 0.0j
        if col != k + 1:
            m[[k + 1, col], k:] = m[[col, k + 1], k:]
            m[k:, [k + 1, col]] = m[k:, [col, k + 1]]
            pf = -pf
        pf *= pivot
        if k + 2 < n:
            b0 = m[k, k + 2 :]
            b1 = m[k + 1, k + 2 :]
            m[k + 2 :, k + 2 :] += (np.outer(b1, b0) - np.outer(b0, b1)) / pivot
    return complex(pf)


def torus_site(nx: int, ny: int, x: int, y: int, sublattice: int) -> int:
    """Site index matching kitaev-numerics HoneycombTorus(nx, ny)."""
    return 2 * ((x % nx) * ny + (y % ny)) + sublattice


def torus_bonds(nx: int, ny: int) -> list[tuple[int, int, str]]:
    """Return ordered list of `(site_A, site_B, axis)` bonds on `HoneycombTorus(nx, ny)`."""
    bonds: list[tuple[int, int, str]] = []
    for x in range(nx):
        for y in range(ny):
            a = torus_site(nx, ny, x, y, 0)
            bonds.append((a, torus_site(nx, ny, x, y, 1), "z"))
            bonds.append((a, torus_site(nx, ny, x + 1, y, 1), "x"))
            bonds.append((a, torus_site(nx, ny, x, y + 1, 1), "y"))
    return bonds


def torus_plaquettes(nx: int, ny: int) -> list[tuple[tuple[int, str], ...]]:
    """Return six `(site, axis)` Pauli factors for each plaquette matching `kitaev-numerics`."""
    plaquettes: list[tuple[tuple[int, str], ...]] = []
    for x in range(nx):
        for y in range(ny):
            plaquettes.append(
                (
                    (torus_site(nx, ny, x, y, 0), "x"),
                    (torus_site(nx, ny, x, y, 1), "y"),
                    (torus_site(nx, ny, x - 1, y, 0), "z"),
                    (torus_site(nx, ny, x - 1, y + 1, 1), "x"),
                    (torus_site(nx, ny, x - 1, y + 1, 0), "y"),
                    (torus_site(nx, ny, x, y + 1, 1), "z"),
                )
            )
    return plaquettes


def assemble_torus_majorana_matrix(
    nx: int,
    ny: int,
    params: ModelParameters,
    state: MeanFieldState,
    bc_twist: tuple[int, int] = (-1, -1),
) -> NDArray[np.float64]:
    """Assemble real skew-symmetric 4N x 4N matrix `A` (H = (i/4) c^T A c) on `(nx, ny)`."""
    cell = make_honeycomb_2site()
    b_mat = build_on_site_matrices(cell, params, state, torch.device("cpu")).numpy()
    a_mat = build_bond_hopping_matrices(cell, params, state, torch.device("cpu")).numpy()

    n_sites = 2 * nx * ny
    mat = np.zeros((4 * n_sites, 4 * n_sites), dtype=np.float64)
    for x in range(nx):
        for y in range(ny):
            for sub in (0, 1):
                s = torus_site(nx, ny, x, y, sub)
                mat[4 * s : 4 * s + 4, 4 * s : 4 * s + 4] = b_mat[sub]
            sa = torus_site(nx, ny, x, y, 0)
            sz = torus_site(nx, ny, x, y, 1)
            mat[4 * sa : 4 * sa + 4, 4 * sz : 4 * sz + 4] += a_mat[0]
            mat[4 * sz : 4 * sz + 4, 4 * sa : 4 * sa + 4] -= a_mat[0].T

            sgn_x = float(bc_twist[0]) if x == nx - 1 else 1.0
            sx = torus_site(nx, ny, x + 1, y, 1)
            mat[4 * sa : 4 * sa + 4, 4 * sx : 4 * sx + 4] += sgn_x * a_mat[1]
            mat[4 * sx : 4 * sx + 4, 4 * sa : 4 * sa + 4] -= sgn_x * a_mat[1].T

            sgn_y = float(bc_twist[1]) if y == ny - 1 else 1.0
            sy = torus_site(nx, ny, x, y + 1, 1)
            mat[4 * sa : 4 * sa + 4, 4 * sy : 4 * sy + 4] += sgn_y * a_mat[2]
            mat[4 * sy : 4 * sy + 4, 4 * sa : 4 * sa + 4] -= sgn_y * a_mat[2].T

    return np.asarray(0.5 * (mat - mat.T), dtype=np.float64)


def build_thouless_data(
    a_mat: NDArray[np.float64],
    n_sites: int,
    *,
    seed: int = 42,
) -> tuple[NDArray[np.int8], NDArray[np.complex128], float] | None:
    """Build skew-symmetric Thouless matrix Z_eff and canonical reference configuration."""
    n_jw = 2 * n_sites
    rng = np.random.default_rng(seed)
    best_cond = float("inf")
    best_ref: NDArray[np.int8] | None = None
    best_z: NDArray[np.complex128] | None = None

    for trial in range(16):
        ref = (
            np.zeros(n_sites, dtype=np.int8)
            if trial == 0
            else rng.integers(0, 2, size=n_sites, dtype=np.int8)
        )
        t_cf = np.zeros((2 * n_jw, 2 * n_jw), dtype=np.complex128)
        for s in range(n_sites):
            m0, m1 = 2 * s, 2 * s + 1
            if ref[s] == 0:
                t_cf[2 * m0, m0] = 1.0
                t_cf[2 * m0, n_jw + m0] = 1.0
                t_cf[2 * m0 + 1, m0] = -1.0j
                t_cf[2 * m0 + 1, n_jw + m0] = 1.0j
                t_cf[2 * m1, m1] = 1.0
                t_cf[2 * m1, n_jw + m1] = 1.0
                t_cf[2 * m1 + 1, m1] = 1.0j
                t_cf[2 * m1 + 1, n_jw + m1] = -1.0j
            else:
                t_cf[2 * m0, m0] = 1.0
                t_cf[2 * m0, n_jw + m0] = 1.0
                t_cf[2 * m0 + 1, m0] = 1.0j
                t_cf[2 * m0 + 1, n_jw + m0] = -1.0j
                t_cf[2 * m1, m1] = 1.0
                t_cf[2 * m1, n_jw + m1] = 1.0
                t_cf[2 * m1 + 1, m1] = -1.0j
                t_cf[2 * m1 + 1, n_jw + m1] = 1.0j
        h_bdg = 0.5j * (t_cf.conj().T @ a_mat @ t_cf)
        evals, evecs = np.linalg.eigh(h_bdg)
        pos_idx = np.where(evals > 0.0)[0]
        if pos_idx.size != n_jw:
            continue
        row_mat = evecs[:, pos_idx].conj().T
        a_blk = row_mat[:, :n_jw]
        cond = float(np.linalg.cond(a_blk))
        if cond < best_cond:
            best_cond = cond
            best_ref = ref
            best_z = -np.linalg.solve(a_blk, row_mat[:, n_jw:])
            if cond < 50.0:
                break

    if best_ref is None or best_z is None or best_cond > 1e8:
        return None

    # Symmetrize-skew and permute (2s <-> 2s+1) whenever ref[s] == 1 so canonical pairs (2s, 2s+1) work uniformly
    z_skew = 0.5 * (best_z - best_z.T)
    perm = np.arange(n_jw, dtype=np.intp)
    for s in range(n_sites):
        if best_ref[s] == 1:
            perm[2 * s], perm[2 * s + 1] = perm[2 * s + 1], perm[2 * s]
    z_eff = z_skew[np.ix_(perm, perm)]
    return best_ref, np.asarray(0.5 * (z_eff - z_eff.T), dtype=np.complex128), best_cond


def evaluate_exact_fock_equivalence(a_mat: NDArray[np.float64], n_sites: int = 2) -> float:
    """Verify 2^N Thouless Pfaffian amplitudes equal full 2^{2N} Fock ground state."""
    if n_sites > 3:
        raise ValueError("Full 2^{2N} Fock check is intended for N <= 3.")
    n_jw = 2 * n_sites
    i2 = np.eye(2, dtype=np.complex128)
    px = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    py = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
    pz = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
    majs: list[NDArray[np.complex128]] = []
    for mode in range(n_jw):
        for pmat in (px, py):
            op = np.array([[1.0 + 0.0j]], dtype=np.complex128)
            for s in range(n_jw):
                op = np.asarray(
                    np.kron(op, pmat if s == mode else (pz if s < mode else i2)),
                    dtype=np.complex128,
                )
            majs.append(op)

    dim_f = 1 << n_jw
    h_fock = np.zeros((dim_f, dim_f), dtype=np.complex128)
    for p in range(2 * n_jw):
        for q in range(2 * n_jw):
            if a_mat[p, q] != 0.0:
                h_fock += 0.25j * a_mat[p, q] * (majs[p] @ majs[q])
    _, evecs = np.linalg.eigh(h_fock)
    psi0 = evecs[:, 0]

    phys_idx: list[int] = []
    for tau in range(1 << n_sites):
        bits: list[int] = []
        for s in range(n_sites):
            t_s = (tau >> s) & 1
            bits.extend([1, 0] if t_s == 1 else [0, 1])
        st_int = sum(b << (n_jw - 1 - m) for m, b in enumerate(bits))
        phys_idx.append(st_int)
    vec_exact = np.asarray(psi0[phys_idx], dtype=np.complex128)
    norm_exact = float(np.linalg.norm(vec_exact))
    if norm_exact < 1e-12:
        raise RuntimeError("Fock ground state has zero physical weight in parity check.")
    vec_exact /= norm_exact

    tdata = build_thouless_data(a_mat, n_sites)
    if tdata is None:
        raise RuntimeError("Thouless construction failed on even-parity test state.")
    ref, z_eff, _ = tdata

    vec_pf = np.empty(1 << n_sites, dtype=np.complex128)
    for tau in range(1 << n_sites):
        j_idx = [
            m for s in range(n_sites) if ((tau >> s) & 1) != int(ref[s]) for m in (2 * s, 2 * s + 1)
        ]
        vec_pf[tau] = complex_pfaffian(z_eff[np.ix_(j_idx, j_idx)])
    vec_pf /= np.linalg.norm(vec_pf)
    ov = np.vdot(vec_pf, vec_exact)
    phase = ov / abs(ov) if abs(ov) > 1e-12 else 1.0 + 0.0j
    return float(np.max(np.abs(vec_exact - phase * vec_pf)))


def _exact_spin_observables(
    nx: int,
    ny: int,
    h: float,
    psi: NDArray[np.complex128],
) -> tuple[float, float, float, float]:
    """Evaluate exact `(e_per_site, m111, wp_mean, nn_bond)` on normalized 2^N vector `psi`."""
    n_sites = 2 * nx * ny
    dim = 1 << n_sites
    idx = np.arange(dim, dtype=np.int32)
    spins_x = [np.where((idx >> s) & 1, 1.0, -1.0) for s in range(n_sites)]
    phases_y = [np.where((idx >> s) & 1, 1.0j, -1.0j) for s in range(n_sites)]

    bond_vals: list[float] = []
    e_k = 0.0
    for u, v, ax in torus_bonds(nx, ny):
        if ax == "x":
            val = float(np.vdot(psi, spins_x[u] * spins_x[v] * psi).real)
        elif ax == "z":
            val = float(np.vdot(psi, psi[idx ^ ((1 << u) | (1 << v))]).real)
        else:
            val = float(
                np.vdot(psi, phases_y[u] * phases_y[v] * psi[idx ^ ((1 << u) | (1 << v))]).real
            )
        bond_vals.append(val)
        e_k -= val

    mag_sum = 0.0
    for s in range(n_sites):
        ex_x = float(np.vdot(psi, spins_x[s] * psi).real)
        flipped = psi[idx ^ (1 << s)]
        ex_z = float(np.vdot(psi, flipped).real)
        ex_y = float(np.vdot(psi, phases_y[s] * flipped).real)
        mag_sum += ex_x + ex_y + ex_z
    comp = h / sqrt(3.0)
    e_f = -comp * mag_sum
    m_111 = mag_sum / (sqrt(3.0) * n_sites)

    wp_vals: list[float] = []
    for factors in torus_plaquettes(nx, ny):
        flip_mask = 0
        ph: NDArray[np.complex128] = np.ones(dim, dtype=np.complex128)
        for site_f, ax in factors:
            if ax == "x":
                ph *= spins_x[site_f]
            elif ax == "z":
                flip_mask |= 1 << site_f
            else:
                flip_mask |= 1 << site_f
                ph *= phases_y[site_f]
        wp_vals.append(float(np.vdot(psi, ph * psi[idx ^ flip_mask]).real))

    return (
        (e_k + e_f) / float(n_sites),
        float(m_111),
        float(np.mean(wp_vals)),
        float(np.mean(bond_vals)),
    )


def evaluate_exact_gutzwiller(
    nx: int,
    ny: int,
    h: float,
    state: MeanFieldState,
    *,
    twists: tuple[tuple[int, int], ...] = ((-1, -1), (-1, 1), (1, -1), (1, 1)),
) -> ExactGutzwillerResult:
    """Evaluate exact 2^N Gutzwiller-projected observables selecting lowest-energy twist."""
    n_sites = 2 * nx * ny
    if n_sites > 14:
        raise ValueError(
            "Exact 2^N Gutzwiller contraction bounded to N <= 14; use VMC for N >= 16."
        )
    comp = h / sqrt(3.0)
    params = ModelParameters(K=-1.0, h=(comp, comp, comp))

    best_res: ExactGutzwillerResult | None = None
    for bc in twists:
        a_mat = assemble_torus_majorana_matrix(nx, ny, params, state, bc_twist=bc)
        tdata = build_thouless_data(a_mat, n_sites)
        if tdata is None:
            continue
        ref, z_eff, cond = tdata
        psi = np.empty(1 << n_sites, dtype=np.complex128)
        for tau in range(1 << n_sites):
            j_idx = [
                m
                for s in range(n_sites)
                if ((tau >> s) & 1) != int(ref[s])
                for m in (2 * s, 2 * s + 1)
            ]
            psi[tau] = complex_pfaffian(z_eff[np.ix_(j_idx, j_idx)])
        norm_val = float(np.linalg.norm(psi))
        if norm_val < 1e-12:
            continue
        psi /= norm_val
        e_site, m111, wp, nn_bond = _exact_spin_observables(nx, ny, h, psi)
        cand = ExactGutzwillerResult(
            nx=nx,
            ny=ny,
            n_sites=n_sites,
            h=float(h),
            bc_twist=bc,
            ref_condition_number=cond,
            e_proj_per_site=e_site,
            magnetization_111=m111,
            plaquette_flux=wp,
            nn_bond_correlator=nn_bond,
        )
        if best_res is None or cand.e_proj_per_site < best_res.e_proj_per_site:
            best_res = cand

    if best_res is None:
        raise RuntimeError(f"No parity-admissible topological twist found on ({nx},{ny}) at h={h}.")
    return best_res


class _ThoulessStateEvaluator:
    """Fast O(N) Pfaffian ratio & local-observable evaluator for spin bitstrings."""

    def __init__(self, ref: NDArray[np.int8], z_eff: NDArray[np.complex128]) -> None:
        self.ref = ref
        self.n_sites = ref.size
        self.z_eff = z_eff

    def amplitude(self, tau: NDArray[np.int8]) -> complex:
        active = [s for s in range(self.n_sites) if int(tau[s]) != int(self.ref[s])]
        if not active:
            return 1.0 + 0.0j
        j_idx = [m for s in active for m in (2 * s, 2 * s + 1)]
        return complex_pfaffian(self.z_eff[np.ix_(j_idx, j_idx)])

    def flip_ratio_with_inv(
        self,
        tau: NDArray[np.int8],
        active: list[int],
        pos_map: dict[int, int],
        j_idx: list[int],
        z_inv: NDArray[np.complex128] | None,
        pf_val: complex,
        flip_site: int,
    ) -> complex:
        p, q = 2 * flip_site, 2 * flip_site + 1
        if flip_site in pos_map:
            pos = pos_map[flip_site]
            assert z_inv is not None
            return -complex(z_inv[2 * pos, 2 * pos + 1])
        if not active:
            return complex(self.z_eff[p, q])
        assert z_inv is not None
        z_aj = self.z_eff[np.ix_([p, q], j_idx)]
        schur = self.z_eff[np.ix_([p, q], [p, q])] + z_aj @ z_inv @ z_aj.T
        return complex(schur[0, 1])

    def two_flip_ratio_with_inv(
        self,
        active: list[int],
        pos_map: dict[int, int],
        j_idx: list[int],
        z_inv: NDArray[np.complex128] | None,
        pf_val: complex,
        u: int,
        v: int,
    ) -> complex:
        u_act = u in pos_map
        v_act = v in pos_map
        if u_act and v_act:
            assert z_inv is not None
            pu, pv = pos_map[u], pos_map[v]
            k4 = [2 * pu, 2 * pu + 1, 2 * pv, 2 * pv + 1]
            return complex_pfaffian(-z_inv[np.ix_(k4, k4)])
        if (not u_act) and (not v_act):
            k4 = [2 * u, 2 * u + 1, 2 * v, 2 * v + 1]
            if not active:
                return complex_pfaffian(self.z_eff[np.ix_(k4, k4)])
            assert z_inv is not None
            z_aj = self.z_eff[np.ix_(k4, j_idx)]
            schur4 = self.z_eff[np.ix_(k4, k4)] + z_aj @ z_inv @ z_aj.T
            return complex_pfaffian(schur4)
        # One active spin removed, one inactive spin added
        assert z_inv is not None
        rem_site, add_site = (u, v) if u_act else (v, u)
        pos = pos_map[rem_site]
        k2 = [2 * pos, 2 * pos + 1]
        w = z_inv[k2[0], k2[1]]
        if abs(w) < 1e-12:
            # Exact direct Pfaffian fallback near accidental node
            rem_set = {rem_site}
            new_act = sorted([s for s in active if s not in rem_set] + [add_site])
            j_new = [m for s in new_act for m in (2 * s, 2 * s + 1)]
            return complex_pfaffian(self.z_eff[np.ix_(j_new, j_new)]) / pf_val
        rem_c = [i for i in range(len(j_idx)) if i not in (k2[0], k2[1])]
        j_1r = [j_idx[i] for i in rem_c]
        pq = [2 * add_site, 2 * add_site + 1]
        if not j_1r:
            return (-complex(w)) * complex(self.z_eff[pq[0], pq[1]])
        u_r = self.z_eff[np.ix_(pq, j_1r)]
        # 2x2 skew inverse of [[0, w], [-w, 0]] is (1/w) * [[0, -1], [1, 0]]
        s_inv = np.array([[0.0, -1.0 / w], [1.0 / w, 0.0]], dtype=np.complex128)
        t1 = u_r @ z_inv[np.ix_(rem_c, rem_c)] @ u_r.T
        c_mat = u_r @ z_inv[np.ix_(rem_c, k2)]
        z_kr_u = z_inv[np.ix_(k2, rem_c)] @ u_r.T
        t2 = c_mat @ s_inv @ z_kr_u
        schur2 = self.z_eff[np.ix_(pq, pq)] + (t1 - t2)
        return (-complex(w)) * complex(schur2[0, 1])


def _integrated_autocorr(series: NDArray[np.float64]) -> float:
    """Positive initial-sequence estimate of integrated autocorrelation time tau_int."""
    n = series.size
    if n < 8:
        return 1.0
    var = float(np.var(series))
    if var <= 1e-24:
        return 1.0
    centered = series - float(np.mean(series))
    max_lag = min(n // 4, 25)
    tau_int = 1.0
    for lag in range(1, max_lag + 1):
        rho = float(np.dot(centered[:-lag], centered[lag:]) / ((n - lag) * var))
        if rho <= 0.0:
            break
        tau_int += 2.0 * rho
    return max(1.0, tau_int)


def select_admissible_twist(
    nx: int,
    ny: int,
    h: float,
    state: MeanFieldState,
    *,
    preferred_bc: tuple[int, int] | None = None,
) -> tuple[tuple[int, int], NDArray[np.int8], NDArray[np.complex128], float]:
    """Pick an even-parity topological twist `(bc_x, bc_y)` minimizing condition number/energy."""
    n_sites = 2 * nx * ny
    comp = h / sqrt(3.0)
    params = ModelParameters(K=-1.0, h=(comp, comp, comp))
    candidates: tuple[tuple[int, int], ...] = ((-1, -1), (-1, 1), (1, -1), (1, 1))
    if preferred_bc is not None and preferred_bc in candidates:
        candidates = (preferred_bc,) + tuple(c for c in candidates if c != preferred_bc)

    if n_sites <= 12 and preferred_bc is None:
        exact = evaluate_exact_gutzwiller(nx, ny, h, state, twists=candidates)
        a_mat = assemble_torus_majorana_matrix(nx, ny, params, state, bc_twist=exact.bc_twist)
        tdata = build_thouless_data(a_mat, n_sites)
        assert tdata is not None
        return exact.bc_twist, tdata[0], tdata[1], tdata[2]

    for bc in candidates:
        a_mat = assemble_torus_majorana_matrix(nx, ny, params, state, bc_twist=bc)
        tdata = build_thouless_data(a_mat, n_sites)
        if tdata is not None:
            return bc, tdata[0], tdata[1], tdata[2]
    raise RuntimeError(f"No parity-admissible topological twist found on ({nx},{ny}) at h={h}.")


def evaluate_vmc_gutzwiller(
    nx: int,
    ny: int,
    h: float,
    state: MeanFieldState,
    *,
    n_chains: int = 50,
    burn_in_sweeps: int = 12,
    samples_per_chain: int = 40,
    seed_base: int = 900,
    preferred_bc: tuple[int, int] | None = None,
) -> VMCGutzwillerResult:
    """Run `n_chains` independent Metropolis-Hastings Gutzwiller VMC chains with BCa 95% CIs."""
    n_sites = 2 * nx * ny
    comp = h / sqrt(3.0)
    bc, ref, z_eff, _ = select_admissible_twist(nx, ny, h, state, preferred_bc=preferred_bc)
    evaluator = _ThoulessStateEvaluator(ref, z_eff)
    bonds = torus_bonds(nx, ny)

    chain_e_means = np.empty(n_chains, dtype=np.float64)
    chain_m_means = np.empty(n_chains, dtype=np.float64)
    all_e_samples = np.empty((n_chains, samples_per_chain), dtype=np.float64)
    taus: list[float] = []

    for cid in range(n_chains):
        rng = np.random.default_rng(seed_base + cid)
        # Initialize from a well-conditioned bitstring with non-zero weight
        tau = ref.copy()
        for _ in range(n_sites // 2):
            s_try = int(rng.integers(0, n_sites))
            cand_tau = tau.copy()
            cand_tau[s_try] = 1 - cand_tau[s_try]
            if abs(evaluator.amplitude(cand_tau)) > 1e-6:
                tau = cand_tau
        pf_curr = evaluator.amplitude(tau)

        def _sync_state(
            curr_tau: NDArray[np.int8],
        ) -> tuple[list[int], dict[int, int], list[int], NDArray[np.complex128] | None, complex]:
            act = [s for s in range(n_sites) if int(curr_tau[s]) != int(ref[s])]
            pm = {s: idx_a for idx_a, s in enumerate(act)}
            j_lst = [m for s in act for m in (2 * s, 2 * s + 1)]
            if not act:
                return act, pm, j_lst, None, 1.0 + 0.0j
            sub = evaluator.z_eff[np.ix_(j_lst, j_lst)]
            return act, pm, j_lst, np.linalg.inv(sub), complex_pfaffian(sub)

        active, pos_map, j_idx, z_inv, pf_curr = _sync_state(tau)
        total_sweeps = burn_in_sweeps + samples_per_chain
        e_series = np.empty(samples_per_chain, dtype=np.float64)
        m_series = np.empty(samples_per_chain, dtype=np.float64)

        for sw in range(total_sweeps):
            for _step in range(n_sites):
                if rng.random() < 0.4:
                    # Two-spin bond move keeps VMC ergodic inside fixed spin-parity sectors at h=0
                    bu, bv, _ = bonds[int(rng.integers(0, len(bonds)))]
                    r2 = evaluator.two_flip_ratio_with_inv(
                        active, pos_map, j_idx, z_inv, pf_curr, bu, bv
                    )
                    prob2 = float(abs(r2) ** 2)
                    if prob2 >= 1.0 or rng.random() < prob2:
                        tau[bu] = 1 - tau[bu]
                        tau[bv] = 1 - tau[bv]
                        active, pos_map, j_idx, z_inv, pf_curr = _sync_state(tau)
                else:
                    s_flip = int(rng.integers(0, n_sites))
                    ratio = evaluator.flip_ratio_with_inv(
                        tau, active, pos_map, j_idx, z_inv, pf_curr, s_flip
                    )
                    prob = float(abs(ratio) ** 2)
                    if prob >= 1.0 or rng.random() < prob:
                        tau[s_flip] = 1 - tau[s_flip]
                        active, pos_map, j_idx, z_inv, pf_curr = _sync_state(tau)

            if sw >= burn_in_sweeps:
                m_idx = sw - burn_in_sweeps
                # Local observables in |01>, |10> basis: sigma_x=-Z, sigma_z=X, sigma_y=Y (<1|Y|0>=+i, <0|Y|1>=-i)
                sx_vals = np.where(tau == 1, 1.0, -1.0)
                ph_y_dst = np.where(tau == 1, 1.0j, -1.0j)  # <tau| sigma_y |tau^flip>

                flip1 = [
                    evaluator.flip_ratio_with_inv(tau, active, pos_map, j_idx, z_inv, pf_curr, s)
                    for s in range(n_sites)
                ]
                m_tot = 0.0
                for s in range(n_sites):
                    # sigma_x is diagonal; Re(<tau|sigma_z|psi>/<tau|psi>) + Re(ph_y*r)
                    r1 = flip1[s]
                    m_tot += float(sx_vals[s] + r1.real + (ph_y_dst[s] * r1).real)

                e_k = 0.0
                for u, v, ax in bonds:
                    if ax == "x":
                        e_k -= float(sx_vals[u] * sx_vals[v])
                    else:
                        r2 = evaluator.two_flip_ratio_with_inv(
                            active, pos_map, j_idx, z_inv, pf_curr, u, v
                        )
                        if ax == "z":
                            e_k -= float(r2.real)
                        else:
                            e_k -= float((ph_y_dst[u] * ph_y_dst[v] * r2).real)

                e_loc = (e_k - comp * m_tot) / float(n_sites)
                m_loc = m_tot / (sqrt(3.0) * n_sites)
                e_series[m_idx] = e_loc
                m_series[m_idx] = m_loc

        chain_e_means[cid] = float(np.mean(e_series))
        chain_m_means[cid] = float(np.mean(m_series))
        all_e_samples[cid, :] = e_series
        taus.append(_integrated_autocorr(e_series))

    e_lo, e_hi = bca_confidence_interval(chain_e_means, seed=seed_base)
    m_lo, m_hi = bca_confidence_interval(chain_m_means, seed=seed_base + 1)
    tau_mean = float(np.mean(taus))
    n_eff = float(n_chains * samples_per_chain) / tau_mean

    return VMCGutzwillerResult(
        nx=nx,
        ny=ny,
        n_sites=n_sites,
        h=float(h),
        bc_twist=bc,
        n_chains=n_chains,
        samples_per_chain=samples_per_chain,
        e_proj_mean=float(np.mean(chain_e_means)),
        e_proj_ci_low=e_lo,
        e_proj_ci_high=e_hi,
        e_loc_std=float(np.std(all_e_samples)),
        tau_int_energy=tau_mean,
        n_eff_total=n_eff,
        magnetization_111_mean=float(np.mean(chain_m_means)),
        magnetization_111_ci_low=m_lo,
        magnetization_111_ci_high=m_hi,
    )
