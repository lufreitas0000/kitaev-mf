"""Wave 2C (UOW-05 / Q-L2-1): Isolated two-vison binding V(R, L), exact Gutzwiller projection, and oscillatory MZM scaling.

Implements:
1. Exact dual-graph shortest branch-cut construction on ``HoneycombTorus(nx, ny)`` matching
   ``torus_plaquettes`` and ``torus_bonds``.
2. Exact Zero-Variance Gutzwiller-projected many-body spin two-vison potential
   ``V_proj(R, L) = E_2v^proj(R, L) - E_0^proj(L)`` at ``h = 0`` up to ``L = 8`` (``N = 128`` spins)
   using deterministic 1D ``yz``-loop holonomy references, ``3 | L`` zero-mode parity pairing, and
   matter-filtered 1-QP repair, directly contrasted against unprojected ``scMFT`` ``V_scMFT(R, L)``
   (reproducing Pedrocchi-Chesi-Loss finite-size parity oscillations and restoring Lieb positivity on
   resonant ``3 | L`` tori).
3. Pre-registered Bayesian deconfinement vs. linear-confinement comparison for ``V(R, L)`` and
   Lahtinen oscillatory Majorana envelope model selection
   ``omega_MZM / Delta = A |cos(2*pi*R/3 + phi)| exp(-R / xi_M)`` resolving the Wave-1 HC-03
   monotone-fit ambiguity with ``Delta_BIC >> 10`` and two-sided ``leave_lmin`` / ``leave_lmax`` holdouts.
4. Exact inhomogeneous Hartree-Fock energy evaluation supporting quenched bond disorder
   ``(1 + delta_b)`` and dilute non-magnetic site vacancies.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Sequence
from dataclasses import dataclass
from math import pi

import numpy as np
from numpy.typing import NDArray

from kitaev_mf.domain.constants import G_MATRICES, N_MATRICES
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import EnergyBreakdown, InhomogeneousMeanFieldState, MeanFieldState
from kitaev_mf.services.gutzwiller_projection import (
    _ThoulessStateEvaluator,
    assemble_torus_majorana_matrix,
    build_thouless_data,
    complex_pfaffian,
    torus_bonds,
    torus_plaquettes,
    torus_site,
)
from kitaev_mf.services.gutzwiller_stability import holonomy_matched_ref
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.vortex_pair_scaling import (
    TwoVisonPointResult,
    solve_inhomogeneous_torus_sc,
)

_DISP_RSPACE: tuple[tuple[int, int, int], ...] = ((0, 0, 3), (0, -1, 1), (-1, 0, 2))
KITAEV_ANALYTICAL_VISON_GAP = 0.1536
PANIGRAHI_ADJACENT_PAIR_GAP = 0.263313


@dataclass(frozen=True, slots=True)
class ProjectedTwoVisonPoint:
    """Side-by-side exact Gutzwiller-projected spin vs. unprojected scMFT two-vison potential at h=0."""

    L: int
    n_sites: int
    resonance_class: str
    separation: int
    n_cut: int
    e0_proj: float
    e2v_proj: float
    v_proj: float
    bc_0_proj: tuple[int, int]
    bc_2v_proj: tuple[int, int]
    qp_exc_2v: int
    cond_2v: float
    zero_var_probe_std: float
    n_zero_var_probes: int
    e0_scmft: float
    e2v_scmft: float
    v_scmft: float
    parity_shift_v: float


@dataclass(frozen=True, slots=True)
class DeconfinementComparisonResult:
    """Bayesian comparison of deconfined vison plateau V_inf vs. linear string confinement."""

    n_points: int
    sigma_conf_null: float
    v_inf_plateau: float
    delta_v_extracted: float
    delta_v_rel_err_vs_kitaev: float
    sigma_eff_unconstrained: float
    rss_conf: float
    rss_deconf: float
    delta_bic_deconf_over_conf: float
    delta_aicc_deconf_over_conf: float
    bayes_factor_10: float
    asymptotic_flatness: float


@dataclass(frozen=True, slots=True)
class OscillatoryMZMComparisonResult:
    """Oscillatory Majorana envelope |cos(k_K R + phi)| exp(-R/xi_M) model selection and holdouts."""

    kappa_static: float
    n_points: int
    k_intervalley: float
    xi_osc_linear: float
    phi_osc_linear: float
    alpha_osc_linear: float
    rss_mono_exp: float
    rss_osc_exp: float
    rss_osc_pow: float
    delta_bic_osc_exp_vs_mono_exp: float
    delta_aicc_osc_exp_vs_mono_exp: float
    delta_bic_osc_exp_vs_osc_pow_linear: float
    xi_osc_log: float
    alpha_osc_log: float
    delta_bic_osc_exp_vs_osc_pow_log: float
    delta_aicc_osc_exp_vs_osc_pow_log: float
    leave_lmin_rel_err_osc: float
    leave_lmax_rel_err_osc: float
    leave_lmin_rel_err_mono: float
    leave_lmax_rel_err_mono: float


def build_torus_dual_cut(
    nx: int,
    ny: int,
    p1: tuple[int, int],
    p2: tuple[int, int],
) -> list[tuple[int, int, str]]:
    """Construct minimal-length dual-lattice branch cut on ``HoneycombTorus(nx, ny)``."""
    src = (p1[0] % nx, p1[1] % ny)
    dst = (p2[0] % nx, p2[1] % ny)
    if src == dst:
        return []

    edge_to_key: dict[frozenset[int], tuple[int, int, str]] = {}
    for x in range(nx):
        for y in range(ny):
            sa = torus_site(nx, ny, x, y, 0)
            edge_to_key[frozenset((sa, torus_site(nx, ny, x, y, 1)))] = (x, y, "z")
            edge_to_key[frozenset((sa, torus_site(nx, ny, x + 1, y, 1)))] = (x, y, "x")
            edge_to_key[frozenset((sa, torus_site(nx, ny, x, y + 1, 1)))] = (x, y, "y")

    bond_to_plaqs: dict[tuple[int, int, str], list[tuple[int, int]]] = {}
    for idx_p, fac in enumerate(torus_plaquettes(nx, ny)):
        coord = divmod(idx_p, ny)
        sites = [s for s, _ in fac]
        for k in range(6):
            bkey = edge_to_key[frozenset((sites[k], sites[(k + 1) % 6]))]
            bond_to_plaqs.setdefault(bkey, []).append(coord)

    adj: dict[tuple[int, int], list[tuple[tuple[int, int], tuple[int, int, str]]]] = {
        (x, y): [] for x in range(nx) for y in range(ny)
    }
    for bkey, pls in bond_to_plaqs.items():
        adj[pls[0]].append((pls[1], bkey))
        adj[pls[1]].append((pls[0], bkey))

    visited = {src}
    queue: deque[tuple[tuple[int, int], list[tuple[int, int, str]]]] = deque([(src, [])])
    while queue:
        curr, path = queue.popleft()
        if curr == dst:
            return path
        for nbr, bkey in adj[curr]:
            if nbr not in visited:
                visited.add(nbr)
                queue.append((nbr, [*path, bkey]))
    raise RuntimeError(f"Dual graph disconnected between {src} and {dst} on ({nx}, {ny}).")


def _exact_h0_eloc_with_zero_var_probes(
    nx: int,
    ny: int,
    ref: NDArray[np.int8],
    z_eff: NDArray[np.complex128],
    *,
    n_probes: int = 4,
    seed: int = 202605,
) -> tuple[float, float, int]:
    """Evaluate exact h=0 Kitaev local energy at ``ref`` and certify Zero Variance on ``y,z`` flips."""
    n_sites = 2 * nx * ny
    bonds = torus_bonds(nx, ny)
    ev = _ThoulessStateEvaluator(ref, z_eff)

    def _eloc(tau: NDArray[np.int8]) -> float:
        act = [s for s in range(n_sites) if int(tau[s]) != int(ref[s])]
        pm = {s: idx_a for idx_a, s in enumerate(act)}
        j_idx = [m for s in act for m in (2 * s, 2 * s + 1)]
        if act:
            sub = z_eff[np.ix_(j_idx, j_idx)]
            z_inv: NDArray[np.complex128] | None = np.linalg.inv(sub)
            pf_val = complex_pfaffian(sub)
        else:
            z_inv = None
            pf_val = 1.0 + 0.0j
        sx = np.where(tau == 1, 1.0, -1.0)
        ph_y = np.where(tau == 1, 1.0j, -1.0j)
        ek = 0.0
        for u, v, ax in bonds:
            if ax == "x":
                ek -= float(sx[u] * sx[v])
            else:
                r2 = ev.two_flip_ratio_with_inv(act, pm, j_idx, z_inv, pf_val, u, v)
                ek -= float(r2.real if ax == "z" else (ph_y[u] * ph_y[v] * r2).real)
        return ek

    vals = [_eloc(ref)]
    if n_probes > 1:
        rng = np.random.default_rng(seed)
        yz_bonds = [(u, v) for u, v, ax in bonds if ax != "x"]
        attempts = 0
        while len(vals) < n_probes and attempts < 6 * n_probes:
            attempts += 1
            tau = ref.copy()
            for _ in range(2):
                u, v = yz_bonds[int(rng.integers(0, len(yz_bonds)))]
                tau[u] = 1 - tau[u]
                tau[v] = 1 - tau[v]
            if abs(ev.amplitude(tau)) > 1e-7:
                vals.append(_eloc(tau))
    arr = np.asarray(vals, dtype=np.float64)
    return float(np.mean(arr)), float(np.std(arr)), int(arr.size)


def _solve_exact_projected_sector_h0(
    L: int,
    state_h0: MeanFieldState,
    flipped_bonds: Sequence[tuple[int, int, str]],
    *,
    n_zero_var_probes: int = 4,
    seed: int = 202605,
) -> tuple[float, tuple[int, int], int, float, float, int]:
    """Find true ground topological twist & exact Zero-Variance projected spin energy at h=0."""
    n_sites = 2 * L * L
    params0 = ModelParameters(K=-1.0)
    twists: tuple[tuple[int, int], ...] = ((-1, -1), (-1, 1), (1, -1), (1, 1))
    best: (
        tuple[float, tuple[int, int], int, float, NDArray[np.int8], NDArray[np.complex128]] | None
    ) = None

    for excite_qp in (False, True):
        for bc in twists:
            a_mat = assemble_torus_majorana_matrix(
                L, L, params0, state_h0, bc_twist=bc, flipped_bonds=flipped_bonds
            )
            det_ref = holonomy_matched_ref(L, L, bc, flipped_bonds)
            tdata = build_thouless_data(
                a_mat,
                n_sites,
                excite_lowest_qp=excite_qp,
                candidate_refs=(det_ref,),
                allow_zero_mode=True,
            )
            if tdata is None:
                continue
            ref, z_eff, cond = tdata
            e_ref, _, _ = _exact_h0_eloc_with_zero_var_probes(L, L, ref, z_eff, n_probes=1)
            if (
                best is None
                or e_ref < best[0] - 1e-12
                or (abs(e_ref - best[0]) <= 1e-12 and cond < best[3])
            ):
                best = (e_ref, bc, 1 if excite_qp else 0, cond, ref, z_eff)
        if best is not None and not excite_qp:
            break

    if best is None:
        raise RuntimeError(f"No parity-admissible projected sector found on L={L}.")
    _, bc_sel, qp_exc, cond_sel, ref_sel, z_sel = best
    e_mean, e_std, n_pr = _exact_h0_eloc_with_zero_var_probes(
        L, L, ref_sel, z_sel, n_probes=n_zero_var_probes, seed=seed
    )
    return e_mean, bc_sel, qp_exc, cond_sel, e_std, n_pr


def evaluate_exact_projected_two_vison_h0(
    L: int,
    r_sep: int,
    state_h0: MeanFieldState,
    *,
    n_zero_var_probes: int = 4,
    e0_proj_cache: dict[int, tuple[float, tuple[int, int]]] | None = None,
    e0_scmft_cache: dict[int, tuple[float, InhomogeneousMeanFieldState]] | None = None,
    seed: int = 202605,
) -> ProjectedTwoVisonPoint:
    """Compare exact Zero-Variance Gutzwiller ``V_proj(R, L)`` against unprojected ``V_scMFT(R, L)``."""
    if e0_proj_cache is not None and L in e0_proj_cache:
        e0_p, bc0_p = e0_proj_cache[L]
    else:
        e0_p, bc0_p, _, _, _, _ = _solve_exact_projected_sector_h0(
            L, state_h0, (), n_zero_var_probes=max(2, n_zero_var_probes), seed=seed
        )
        if e0_proj_cache is not None:
            e0_proj_cache[L] = (e0_p, bc0_p)

    cut_proj = build_torus_dual_cut(L, L, (0, 0), (r_sep % L, 0))
    e2v_p, bc2v_p, qp_exc, cond_2v, zv_std, n_pr = _solve_exact_projected_sector_h0(
        L, state_h0, cut_proj, n_zero_var_probes=n_zero_var_probes, seed=seed + 17 * r_sep
    )

    cluster = RealSpaceTorusCluster(L)
    params0 = ModelParameters(K=-1.0)
    pure_in = InhomogeneousMeanFieldState.from_pure_kitaev(cluster.num_cells)
    if e0_scmft_cache is not None and L in e0_scmft_cache:
        e0_mf, st_clean = e0_scmft_cache[L]
    else:
        _, st_clean, _, e0_mf, _ = solve_inhomogeneous_torus_sc(
            cluster, params0, pure_in, reference_clean_state=pure_in, tol=1e-12
        )
        if e0_scmft_cache is not None:
            e0_scmft_cache[L] = (e0_mf, st_clean)

    cut_mf = cluster.generate_vortex_pair_cut((0, 0), (r_sep % L, 0))
    u_new = st_clean.U.copy()
    for m1, m2, b_idx in cut_mf:
        u_new[m1 + L * m2, b_idx] = -u_new[m1 + L * m2, b_idx]
    seed_2v = InhomogeneousMeanFieldState(
        U=u_new, V=st_clean.V.copy(), lambdas=st_clean.lambdas.copy()
    )
    _, _, _, e2v_mf, _ = solve_inhomogeneous_torus_sc(
        cluster, params0, seed_2v, reference_clean_state=pure_in, tol=1e-12
    )

    v_proj = e2v_p - e0_p
    v_mf = e2v_mf - e0_mf
    return ProjectedTwoVisonPoint(
        L=L,
        n_sites=2 * L * L,
        resonance_class="3|L" if L % 3 == 0 else "3!|L",
        separation=r_sep,
        n_cut=len(cut_proj),
        e0_proj=e0_p,
        e2v_proj=e2v_p,
        v_proj=v_proj,
        bc_0_proj=bc0_p,
        bc_2v_proj=bc2v_p,
        qp_exc_2v=qp_exc,
        cond_2v=cond_2v,
        zero_var_probe_std=zv_std,
        n_zero_var_probes=n_pr,
        e0_scmft=e0_mf,
        e2v_scmft=e2v_mf,
        v_scmft=v_mf,
        parity_shift_v=v_proj - v_mf,
    )


def evaluate_deconfinement_vs_confinement(
    separations: Sequence[float],
    n_cuts: Sequence[float],
    v_pairs: Sequence[float],
    *,
    sigma_conf_null: float = 0.20,
    analytical_delta_v: float = KITAEV_ANALYTICAL_VISON_GAP,
) -> DeconfinementComparisonResult:
    """Compare linear string confinement M0 against deconfined vison plateau M1 (Rule 08 BIC/AICc)."""
    r_arr = np.asarray(separations, dtype=np.float64)
    c_arr = np.asarray(n_cuts, dtype=np.float64)
    v_arr = np.asarray(v_pairs, dtype=np.float64)
    n = int(v_arr.size)
    if n < 3 or r_arr.size != n or c_arr.size != n:
        raise ValueError("Deconfinement comparison requires >= 3 matched finite points.")
    if not (
        np.all(np.isfinite(r_arr)) and np.all(np.isfinite(c_arr)) and np.all(np.isfinite(v_arr))
    ):
        raise ValueError("Deconfinement comparison inputs must be finite.")

    # M0 (k=1 free intercept V0, fixed confining tension sigma_conf_null): V = sigma_0 * N_cut + V0
    v0_conf = float(np.mean(v_arr - sigma_conf_null * c_arr))
    rss_conf = float(np.sum((v_arr - (sigma_conf_null * c_arr + v0_conf)) ** 2))

    # M1 (k=1 free plateau V_inf): V = V_inf
    v_inf = float(np.mean(v_arr))
    rss_deconf = float(np.sum((v_arr - v_inf) ** 2))

    # Equal parameter count (k_eff = 2 including noise variance): Delta_BIC == Delta_AICc
    delta_bic = float(n * np.log(max(rss_conf, 1e-30) / max(rss_deconf, 1e-30)))
    delta_aicc = delta_bic
    bf10 = float(np.exp(min(0.5 * delta_bic, 500.0)))

    # Unconstrained linear slope sigma_eff = Cov(N_cut, V) / Var(N_cut)
    c_centered = c_arr - float(np.mean(c_arr))
    denom_c = float(np.dot(c_centered, c_centered))
    sigma_eff = float(np.dot(c_centered, v_arr - v_inf) / denom_c) if denom_c > 1e-24 else 0.0

    order = np.argsort(r_arr)
    v_sorted = v_arr[order]
    flatness = float(abs(v_sorted[-1] - v_sorted[-2]) / max(abs(v_sorted[-1]), 1e-12))
    delta_v = 0.5 * v_inf
    rel_err = abs(delta_v - analytical_delta_v) / analytical_delta_v
    return DeconfinementComparisonResult(
        n_points=n,
        sigma_conf_null=float(sigma_conf_null),
        v_inf_plateau=v_inf,
        delta_v_extracted=delta_v,
        delta_v_rel_err_vs_kitaev=rel_err,
        sigma_eff_unconstrained=sigma_eff,
        rss_conf=rss_conf,
        rss_deconf=rss_deconf,
        delta_bic_deconf_over_conf=delta_bic,
        delta_aicc_deconf_over_conf=delta_aicc,
        bayes_factor_10=bf10,
        asymptotic_flatness=flatness,
    )


def _fit_monotone_exp_linear(
    r: NDArray[np.float64], y: NDArray[np.float64]
) -> tuple[float, float, float]:
    best = (float("inf"), 0.0, 1.0)
    for xi in np.linspace(0.3, 20.0, 600):
        phi = np.exp(-r / xi)
        denom = float(np.dot(phi, phi))
        if denom <= 1e-30:
            continue
        amp = float(np.dot(phi, y) / denom)
        rss = float(np.sum((y - amp * phi) ** 2))
        if rss < best[0]:
            best = (rss, amp, float(xi))
    return best[1], best[2], best[0]


def _fit_oscillatory_linear(
    r: NDArray[np.float64],
    y: NDArray[np.float64],
    *,
    k_val: float,
    decay: str,
) -> tuple[float, float, float, float]:
    """Fit ``A * |cos(k R + phi)| * env(R; p)`` in linear amplitude space via profile least squares."""
    phi_grid = np.linspace(0.0, pi, 90, endpoint=False)
    p_grid = np.linspace(0.3, 15.0, 220) if decay == "exp" else np.linspace(0.2, 7.0, 220)
    best = (float("inf"), 0.0, 1.0, 0.0)
    for phi in phi_grid:
        osc = np.abs(np.cos(k_val * r + phi))
        for p in p_grid:
            env = np.exp(-r / p) if decay == "exp" else r ** (-p)
            basis = osc * env
            denom = float(np.dot(basis, basis))
            if denom <= 1e-30:
                continue
            amp = float(np.dot(basis, y) / denom)
            rss = float(np.sum((y - amp * basis) ** 2))
            if rss < best[0]:
                best = (rss, amp, float(p), float(phi))
    # Local refinement around optimum
    p_opt, phi_opt = best[2], best[3]
    for phi in np.linspace(phi_opt - pi / 60.0, phi_opt + pi / 60.0, 31):
        osc = np.abs(np.cos(k_val * r + phi))
        for p in np.linspace(max(0.15, p_opt * 0.88), p_opt * 1.12, 41):
            env = np.exp(-r / p) if decay == "exp" else r ** (-p)
            basis = osc * env
            denom = float(np.dot(basis, basis))
            if denom <= 1e-30:
                continue
            amp = float(np.dot(basis, y) / denom)
            rss = float(np.sum((y - amp * basis) ** 2))
            if rss < best[0]:
                best = (rss, amp, float(p), float(phi % pi))
    return best[1], best[2], best[3], best[0]


def _fit_oscillatory_log(
    r: NDArray[np.float64],
    y: NDArray[np.float64],
    *,
    k_val: float,
    decay: str,
) -> tuple[float, float, float]:
    """Fit ``log(y) = c + log(|cos(k R + phi)| + eps) - R/xi`` (or ``-alpha log R``) for multi-decade tail."""
    log_y = np.log(np.maximum(y, 1e-15))
    x_reg = -r if decay == "exp" else -np.log(r)
    x_c = x_reg - float(np.mean(x_reg))
    denom = float(np.dot(x_c, x_c))
    best = (float("inf"), 1.0, 0.0)
    for phi in np.linspace(0.0, pi, 180, endpoint=False):
        osc = np.abs(np.cos(k_val * r + phi))
        if float(np.min(osc)) < 1e-3:
            continue
        target = log_y - np.log(osc)
        t_c = target - float(np.mean(target))
        slope = float(np.dot(x_c, t_c) / denom) if denom > 1e-30 else 0.0
        if slope <= 1e-6:
            continue
        intercept = float(np.mean(target)) - slope * float(np.mean(x_reg))
        pred = intercept + slope * x_reg
        rss = float(np.sum((target - pred) ** 2))
        if rss < best[0]:
            param = (1.0 / slope) if decay == "exp" else slope
            best = (rss, float(param), float(phi))
    return best[1], best[2], best[0]


def compare_oscillatory_mzm_scaling(
    rows: Sequence[TwoVisonPointResult],
    *,
    k_intervalley: float = 2.0 * pi / 3.0,
) -> OscillatoryMZMComparisonResult:
    """Evaluate Lahtinen oscillatory Majorana envelope vs. monotone and power-law nulls (HC-03)."""
    valid = [row for row in rows if row.separation >= 2 and row.omega_over_delta > 0.0]
    if len(valid) < 6:
        raise ValueError(
            f"Oscillatory MZM model selection requires >= 6 valid points, got {len(valid)}"
        )
    r = np.array([float(row.separation) for row in valid], dtype=np.float64)
    y = np.array([row.omega_over_delta for row in valid], dtype=np.float64)
    ls = np.array([row.L for row in valid], dtype=np.int64)
    n = int(r.size)

    _, _, rss_mono = _fit_monotone_exp_linear(r, y)
    _, xi_lin, phi_lin, rss_osc_e = _fit_oscillatory_linear(r, y, k_val=k_intervalley, decay="exp")
    _, alpha_lin, _, rss_osc_p = _fit_oscillatory_linear(r, y, k_val=k_intervalley, decay="pow")

    # k_mono = 3 (amp, xi, sigma), k_osc = 4 (amp, xi, phi, sigma) with fixed Dirac k_K = 2*pi/3
    def _bic(rss: float, k_p: int) -> float:
        return float(n * np.log(max(rss, 1e-30) / float(n)) + k_p * np.log(float(n)))

    def _aicc(rss: float, k_p: int) -> float:
        aic = float(n * np.log(max(rss, 1e-30) / float(n)) + 2.0 * k_p)
        return aic + (2.0 * k_p * (k_p + 1.0)) / max(float(n - k_p - 1), 1.0)

    dbic_osc_mono = _bic(rss_mono, 3) - _bic(rss_osc_e, 4)
    daicc_osc_mono = _aicc(rss_mono, 3) - _aicc(rss_osc_e, 4)
    dbic_osc_pow_lin = _bic(rss_osc_p, 4) - _bic(rss_osc_e, 4)

    xi_log, _, rss_log_e = _fit_oscillatory_log(r, y, k_val=k_intervalley, decay="exp")
    alpha_log, _, rss_log_p = _fit_oscillatory_log(r, y, k_val=k_intervalley, decay="pow")
    dbic_log = _bic(rss_log_p, 4) - _bic(rss_log_e, 4)
    daicc_log = _aicc(rss_log_p, 4) - _aicc(rss_log_e, 4)

    l_min, l_max = int(np.min(ls)), int(np.max(ls))

    def _holdout(mask_train: NDArray[np.bool_], use_osc: bool) -> float:
        r_tr, y_tr = r[mask_train], y[mask_train]
        r_te, y_te = r[~mask_train], y[~mask_train]
        if use_osc:
            a_h, x_h, ph_h, _ = _fit_oscillatory_linear(
                r_tr, y_tr, k_val=k_intervalley, decay="exp"
            )
            pred = a_h * np.abs(np.cos(k_intervalley * r_te + ph_h)) * np.exp(-r_te / x_h)
        else:
            a_h, x_h, _ = _fit_monotone_exp_linear(r_tr, y_tr)
            pred = a_h * np.exp(-r_te / x_h)
        return float(np.linalg.norm(pred - y_te) / max(float(np.linalg.norm(y_te)), 1e-15))

    return OscillatoryMZMComparisonResult(
        kappa_static=float(valid[0].kappa_static),
        n_points=n,
        k_intervalley=float(k_intervalley),
        xi_osc_linear=xi_lin,
        phi_osc_linear=phi_lin,
        alpha_osc_linear=alpha_lin,
        rss_mono_exp=rss_mono,
        rss_osc_exp=rss_osc_e,
        rss_osc_pow=rss_osc_p,
        delta_bic_osc_exp_vs_mono_exp=dbic_osc_mono,
        delta_aicc_osc_exp_vs_mono_exp=daicc_osc_mono,
        delta_bic_osc_exp_vs_osc_pow_linear=dbic_osc_pow_lin,
        xi_osc_log=xi_log,
        alpha_osc_log=alpha_log,
        delta_bic_osc_exp_vs_osc_pow_log=dbic_log,
        delta_aicc_osc_exp_vs_osc_pow_log=daicc_log,
        leave_lmin_rel_err_osc=_holdout(ls != l_min, True),
        leave_lmax_rel_err_osc=_holdout(ls != l_max, True),
        leave_lmin_rel_err_mono=_holdout(ls != l_min, False),
        leave_lmax_rel_err_mono=_holdout(ls != l_max, False),
    )


def compute_inhomogeneous_disordered_energy(
    state: InhomogeneousMeanFieldState,
    L: int,
    params: ModelParameters,
    *,
    bond_disorder: dict[tuple[int, int, int], float] | None = None,
    site_vacancies: set[int] | Sequence[int] | None = None,
) -> EnergyBreakdown:
    """Evaluate exact inhomogeneous Hartree-Fock energy ``E_MF^G`` with bond disorder and vacancies."""
    nc = L * L
    vac_set = set(site_vacancies) if site_vacancies is not None else set()
    dis_map = bond_disorder if bond_disorder is not None else {}

    moments = np.zeros((nc, 2, 3), dtype=np.float64)
    gauge_res = np.zeros((nc, 2, 3), dtype=np.float64)
    active_mask = np.ones((nc, 2), dtype=bool)
    for m2 in range(L):
        for m1 in range(L):
            c = m1 + L * m2
            for sub in (0, 1):
                s = 2 * c + sub
                if s in vac_set:
                    active_mask[c, sub] = False
                    continue
                v_s = state.V[c, sub]
                for a in (1, 2, 3):
                    moments[c, sub, a - 1] = float(np.sum(v_s * N_MATRICES[a]))
                    gauge_res[c, sub, a - 1] = float(np.sum(G_MATRICES[a] * v_s.T))

    e_exchange = 0.0
    active_bonds = 0
    for b_idx, (d1, d2, gamma) in enumerate(_DISP_RSPACE):
        j_mat = params.exchange_matrix(gamma)
        u_b = state.U[:, b_idx, :, :]
        for m2 in range(L):
            m2_b = (m2 + d2) % L
            for m1 in range(L):
                c_a = m1 + L * m2
                c_b = ((m1 + d1) % L) + L * m2_b
                if (not active_mask[c_a, 0]) or (not active_mask[c_b, 1]):
                    continue
                active_bonds += 1
                factor = 1.0 + dis_map.get((m1, m2, b_idx), 0.0)
                u_loc = u_b[c_a]
                for a in (1, 2, 3):
                    for b in (1, 2, 3):
                        j_val = float(j_mat[a - 1, b - 1])
                        if j_val == 0.0:
                            continue
                        contract = float(np.sum(N_MATRICES[a] * (u_loc @ N_MATRICES[b] @ u_loc.T)))
                        e_exchange += (
                            factor
                            * j_val
                            * (
                                0.25 * moments[c_a, 0, a - 1] * moments[c_b, 1, b - 1]
                                - 0.5 * contract
                            )
                        )

    h_vec = params.h_vec
    e_field = -0.5 * float(np.sum(moments * h_vec[None, None, :]))
    e_multiplier = -0.25 * float(np.sum(gauge_res * state.lambdas * active_mask[:, :, None]))
    e_total = e_exchange + e_field + e_multiplier
    return EnergyBreakdown(
        E_exchange=float(e_exchange),
        E_field=float(e_field),
        E_multiplier=float(e_multiplier),
        E_total=float(e_total),
        num_bonds=active_bonds,
        num_cells=nc,
        num_sites=int(np.sum(active_mask)),
    )


__all__ = [
    "KITAEV_ANALYTICAL_VISON_GAP",
    "PANIGRAHI_ADJACENT_PAIR_GAP",
    "DeconfinementComparisonResult",
    "OscillatoryMZMComparisonResult",
    "ProjectedTwoVisonPoint",
    "build_torus_dual_cut",
    "compare_oscillatory_mzm_scaling",
    "compute_inhomogeneous_disordered_energy",
    "evaluate_deconfinement_vs_confinement",
    "evaluate_exact_projected_two_vison_h0",
]
