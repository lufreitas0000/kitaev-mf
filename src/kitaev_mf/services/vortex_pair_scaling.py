"""Wave 1 HC-03: Self-consistent multi-L two-vison energetics, chiral NNN coupling, and MZM scaling."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from kitaev_mf.domain.constants import G_MATRICES, N_MATRICES
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import EnergyBreakdown, InhomogeneousMeanFieldState
from kitaev_mf.services.fixed_point import JointAndersonSolver
from kitaev_mf.services.purification import purify_newton_schulz
from kitaev_mf.services.real_space import RealSpaceTorusCluster

_DISP: tuple[tuple[int, int, int], ...] = ((0, 0, 3), (0, -1, 1), (-1, 0, 2))


@dataclass(frozen=True, slots=True)
class TorusVortexSolveMetrics:
    """Typed solver convergence and spectral diagnostics for an inhomogeneous torus solve."""

    converged: bool
    iterations: int
    final_residual: float
    max_gauge_residual: float
    w_pos: np.ndarray
    evecs: np.ndarray


@dataclass(frozen=True, slots=True)
class TwoVisonPointResult:
    """Converged self-consistent two-vison torus observable summary."""

    L: int
    resonance_class: str
    separation: int
    kappa_static: float
    h_field: float
    converged: bool
    iterations: int
    max_gauge_residual: float
    defect_count: int
    min_abs_pfaffian_wp: float
    e_total: float
    v_pair: float
    omega_mzm_can: float
    delta_bulk_can: float
    e0_over_delta: float
    omega_over_delta: float
    mzm_core_weight: float
    bulk_core_weight: float


@dataclass(frozen=True, slots=True)
class ScalingComparisonResult:
    """Multi-L model comparison (exp(-R/xi_M) vs R^-alpha) and two-sided holdout errors."""

    kappa_static: float
    n_points: int
    xi_fit: float
    alpha_fit: float
    rss_exp: float
    rss_pow: float
    delta_bic_exp_over_pow: float
    delta_aicc_exp_over_pow: float
    leave_lmin_rel_err: float
    leave_lmax_rel_err: float


def compute_inhomogeneous_gaussian_energy(
    state: InhomogeneousMeanFieldState,
    L: int,
    params: ModelParameters,
) -> EnergyBreakdown:
    """Evaluate exact inhomogeneous Hartree-Fock energy E_MF^G[U, V, lambda] on an LxL torus."""
    nc = L * L
    moments = np.zeros((nc, 2, 3), dtype=np.float64)
    gauge_res = np.zeros((nc, 2, 3), dtype=np.float64)
    for a in (1, 2, 3):
        n_a = N_MATRICES[a]
        g_a = G_MATRICES[a]
        # tr(V^T N^a) = sum(V * N^a); tr(G^gamma V) = sum_{ij} G_{ij} V_{ji}
        moments[:, :, a - 1] = np.einsum("csij,ij->cs", state.V, n_a)
        gauge_res[:, :, a - 1] = np.einsum("ij,csji->cs", g_a, state.V)

    e_exchange = 0.0
    for b_idx, (d1, d2, gamma) in enumerate(_DISP):
        j_mat = params.exchange_matrix(gamma)
        u_b = state.U[:, b_idx, :, :]
        for m2 in range(L):
            m2_b = (m2 + d2) % L
            for m1 in range(L):
                c_a = m1 + L * m2
                c_b = ((m1 + d1) % L) + L * m2_b
                u_loc = u_b[c_a]
                for a in (1, 2, 3):
                    for b in (1, 2, 3):
                        j_val = float(j_mat[a - 1, b - 1])
                        if j_val == 0.0:
                            continue
                        contract = float(np.sum(N_MATRICES[a] * (u_loc @ N_MATRICES[b] @ u_loc.T)))
                        e_exchange += j_val * (
                            0.25 * moments[c_a, 0, a - 1] * moments[c_b, 1, b - 1] - 0.5 * contract
                        )

    h_vec = params.h_vec
    e_field = -0.5 * float(np.sum(moments * h_vec[None, None, :]))
    e_multiplier = -0.25 * float(np.sum(gauge_res * state.lambdas))
    e_total = e_exchange + e_field + e_multiplier
    return EnergyBreakdown(
        E_exchange=float(e_exchange),
        E_field=float(e_field),
        E_multiplier=float(e_multiplier),
        E_total=float(e_total),
        num_bonds=3 * nc,
        num_cells=nc,
        num_sites=2 * nc,
    )


def bond_signs_from_state_gauge(
    state: InhomogeneousMeanFieldState,
    L: int,
    reference_state: InhomogeneousMeanFieldState,
    flipped_bonds: Sequence[tuple[int, int, int]] | None = None,
) -> np.ndarray:
    """Compute relative bond signs u_{c,b} / u_{c,b}^{(ref)} for gauge-covariant NNN coupling."""
    if flipped_bonds is not None:
        fl = set(flipped_bonds)
        signs = np.ones((L, L, 3), dtype=np.float64)
        for m1, m2, b_idx in fl:
            signs[m1 % L, m2 % L, b_idx] = -1.0
        return signs

    signs = np.ones((L, L, 3), dtype=np.float64)
    disp_gamma = (3, 1, 2)
    for m2 in range(L):
        for m1 in range(L):
            c = m1 + L * m2
            for b_idx, gamma in enumerate(disp_gamma):
                u_val = float(state.U[c, b_idx, gamma, gamma])
                u_ref = float(reference_state.U[c, b_idx, gamma, gamma])
                if u_val * u_ref < 0.0:
                    signs[m1, m2, b_idx] = -1.0
    return signs


def build_covariant_nnn_matrix(
    cluster: RealSpaceTorusCluster,
    kappa_static: float,
    bond_signs: np.ndarray,
) -> np.ndarray:
    """Build gauge-covariant c_0 Haldane NNN skew-matrix A_NNN (canonical 4*H_B units).

    Normalization: hopping element is `c_nnn = 2.0 * kappa_static`, matching
    `Delta_can^MF = 12 sqrt(3) |kappa_static| = 6 sqrt(3) |kappa_KWG|` when
    `kappa_static = kappa_KWG / 2`.
    """
    M = cluster.num_majoranas
    L = cluster.L
    A = np.zeros((M, M), dtype=np.float64)
    if kappa_static == 0.0:
        return A

    c_nnn = 2.0 * float(kappa_static)
    nn_sign: dict[tuple[int, int], float] = {}
    for m2 in range(L):
        for m1 in range(L):
            s_a = cluster.site_index(m1, m2, 0)
            for b_idx, (d1, d2, _) in enumerate(_DISP):
                s_b = cluster.site_index(m1 + d1, m2 + d2, 1)
                sgn = float(bond_signs[m1, m2, b_idx])
                nn_sign[(s_a, s_b)] = sgn
                nn_sign[(s_b, s_a)] = sgn

    for m2 in range(L):
        for m1 in range(L):
            s_b = cluster.site_index(m1, m2, 1)
            a_ring = (
                cluster.site_index(m1, m2, 0),
                cluster.site_index(m1, m2 + 1, 0),
                cluster.site_index(m1 + 1, m2, 0),
            )
            for k in range(3):
                s1 = a_ring[k]
                s2 = a_ring[(k + 1) % 3]
                cov = nn_sign[(s1, s_b)] * nn_sign[(s_b, s2)]
                A[4 * s1, 4 * s2] += c_nnn * cov
                A[4 * s2, 4 * s1] -= c_nnn * cov

            s_a = cluster.site_index(m1, m2, 0)
            b_ring = (
                cluster.site_index(m1, m2, 1),
                cluster.site_index(m1, m2 - 1, 1),
                cluster.site_index(m1 - 1, m2, 1),
            )
            for k in range(3):
                s1 = b_ring[k]
                s2 = b_ring[(k + 1) % 3]
                cov = nn_sign[(s1, s_a)] * nn_sign[(s_a, s2)]
                A[4 * s1, 4 * s2] += c_nnn * cov
                A[4 * s2, 4 * s1] -= c_nnn * cov

    return np.asarray(0.5 * (A - A.T), dtype=np.float64)


def _compute_gauge_residual_tensor(v_tensor: np.ndarray) -> np.ndarray:
    g_res = np.zeros((v_tensor.shape[0], 2, 3), dtype=np.float64)
    for gamma in (1, 2, 3):
        g_res[:, :, gamma - 1] = np.einsum("ij,csji->cs", G_MATRICES[gamma], v_tensor)
    return g_res


def solve_inhomogeneous_torus_sc(
    cluster: RealSpaceTorusCluster,
    params: ModelParameters,
    initial_state: InhomogeneousMeanFieldState,
    *,
    reference_clean_state: InhomogeneousMeanFieldState,
    kappa_static: float = 0.0,
    use_eigh: bool = True,
    max_iter: int = 60,
    tol: float = 1e-12,
    mixing_alpha: float = 0.5,
    lambda_step: float = 0.6,
    diis_history_depth: int = 8,
    flipped_bonds: Sequence[tuple[int, int, int]] | None = None,
) -> tuple[np.ndarray, InhomogeneousMeanFieldState, np.ndarray, float, TorusVortexSolveMetrics]:
    """Run DIIS-accelerated self-consistent inhomogeneous Hartree-Fock on the torus."""
    curr_state = initial_state
    accel = JointAndersonSolver(history_depth=diis_history_depth, beta=mixing_alpha)
    prev_omega: np.ndarray | None = None
    converged = False
    iters = 0
    max_g = 1.0
    final_res = 1.0
    last_a_nnn = np.zeros((cluster.num_majoranas, cluster.num_majoranas), dtype=np.float64)
    cached_signs: np.ndarray | None = None
    last_w_pos: np.ndarray | None = None
    last_vecs: np.ndarray | None = None

    for it in range(max_iter):
        iters = it + 1
        a_mf = cluster.assemble_skew_hamiltonian(params, curr_state, flipped_bonds=flipped_bonds)
        if kappa_static != 0.0:
            signs = bond_signs_from_state_gauge(
                curr_state, cluster.L, reference_clean_state, flipped_bonds=flipped_bonds
            )
            if cached_signs is None or not np.array_equal(signs, cached_signs):
                last_a_nnn = build_covariant_nnn_matrix(cluster, kappa_static, signs)
                cached_signs = signs
            a_tot = a_mf + last_a_nnn
        else:
            a_tot = a_mf

        if use_eigh:
            s_mat = 0.5 * ((-a_tot @ a_tot) + (-a_tot @ a_tot).T)
            w2, vecs = np.linalg.eigh(s_mat)
            w_pos = np.sqrt(np.maximum(w2, 0.0))
            inv_w = np.where(w_pos > 1e-11, 1.0 / np.maximum(w_pos, 1e-30), 0.0)
            omega_next = -a_tot @ (vecs @ (inv_w[:, None] * vecs.T))
            omega_next = 0.5 * (omega_next - omega_next.T)
            last_w_pos = w_pos
            last_vecs = vecs
        else:
            omega_next, _ = purify_newton_schulz(
                a_tot, warm_start=prev_omega, tol=1e-13, max_iter=35
            )
        prev_omega = omega_next

        raw_state = cluster.extract_inhomogeneous_state(omega_next, curr_state.lambdas)
        g_res = _compute_gauge_residual_tensor(raw_state.V)
        lam_target = curr_state.lambdas - lambda_step * g_res
        diff_u = float(np.max(np.abs(raw_state.U - curr_state.U)))
        diff_v = float(np.max(np.abs(raw_state.V - curr_state.V)))
        max_g = float(np.max(np.abs(g_res)))
        final_res = max(diff_u, diff_v, max_g)

        if final_res < tol:
            converged = True
            curr_state = InhomogeneousMeanFieldState(
                U=raw_state.U, V=raw_state.V, lambdas=lam_target
            )
            break

        phi_state = InhomogeneousMeanFieldState(U=raw_state.U, V=raw_state.V, lambdas=lam_target)
        x_next, _ = accel.step(curr_state.flatten(), phi_state.flatten())
        curr_state = InhomogeneousMeanFieldState.from_flattened(x_next, cluster.num_cells)

    assert prev_omega is not None
    if last_w_pos is None or last_vecs is None or not use_eigh:
        a_final = cluster.assemble_skew_hamiltonian(params, curr_state, flipped_bonds=flipped_bonds)
        if kappa_static != 0.0:
            signs = bond_signs_from_state_gauge(
                curr_state, cluster.L, reference_clean_state, flipped_bonds=flipped_bonds
            )
            if cached_signs is None or not np.array_equal(signs, cached_signs):
                last_a_nnn = build_covariant_nnn_matrix(cluster, kappa_static, signs)
            a_final = a_final + last_a_nnn
        s_mat = 0.5 * ((-a_final @ a_final) + (-a_final @ a_final).T)
        w2, vecs = np.linalg.eigh(s_mat)
        last_w_pos = np.sqrt(np.maximum(w2, 0.0))
        last_vecs = vecs
        inv_w = np.where(last_w_pos > 1e-11, 1.0 / np.maximum(last_w_pos, 1e-30), 0.0)
        omega_final = -a_final @ (last_vecs @ (inv_w[:, None] * last_vecs.T))
        omega_final = 0.5 * (omega_final - omega_final.T)
    else:
        omega_final = prev_omega

    final_raw = cluster.extract_inhomogeneous_state(omega_final, curr_state.lambdas)
    max_g_final = float(np.max(np.abs(_compute_gauge_residual_tensor(final_raw.V))))
    e_hf = compute_inhomogeneous_gaussian_energy(curr_state, cluster.L, params).E_total
    e_nnn = 0.25 * float(np.sum(last_a_nnn * omega_final)) if kappa_static != 0.0 else 0.0
    e_tot = e_hf + e_nnn
    metrics = TorusVortexSolveMetrics(
        converged=converged,
        iterations=iters,
        final_residual=final_res,
        max_gauge_residual=max_g_final,
        w_pos=last_w_pos,
        evecs=last_vecs,
    )
    return omega_final, curr_state, last_w_pos, e_tot, metrics


def measure_two_vison_point(
    L: int,
    r_sep: int,
    *,
    kappa_static: float = 0.0,
    h_field: float = 0.0,
    e0_cache: dict[tuple[int, float, float], tuple[float, float, InhomogeneousMeanFieldState]]
    | None = None,
    use_eigh: bool = True,
    tol: float = 1e-12,
) -> TwoVisonPointResult:
    """Solve clean reference + two-vison separation R on an LxL torus."""
    cluster = RealSpaceTorusCluster(L)
    comp = h_field / float(np.sqrt(3.0))
    params = ModelParameters(K=-1.0, h=(comp, comp, comp))
    pure_in = InhomogeneousMeanFieldState.from_pure_kitaev(cluster.num_cells)

    key = (L, kappa_static, h_field)
    if e0_cache is not None and key in e0_cache:
        e0_clean, delta_clean_torus, st_clean = e0_cache[key]
    else:
        _, st_clean, w_pos_clean, e0_clean, _ = solve_inhomogeneous_torus_sc(
            cluster,
            params,
            pure_in,
            reference_clean_state=pure_in,
            kappa_static=kappa_static,
            use_eigh=use_eigh,
            tol=tol,
        )
        delta_clean_torus = 2.0 * float(np.min(w_pos_clean))
        if e0_cache is not None:
            e0_cache[key] = (e0_clean, delta_clean_torus, st_clean)

    p1 = (2, 2)
    p2 = ((2 + r_sep) % L, 2)
    if h_field == 0.0:
        cut = cluster.generate_vortex_pair_cut(p1, p2)
        u_new = st_clean.U.copy()
        for m1, m2, b_idx in cut:
            u_new[m1 + L * m2, b_idx] = -u_new[m1 + L * m2, b_idx]
        seed_2v = InhomogeneousMeanFieldState(
            U=u_new, V=st_clean.V.copy(), lambdas=st_clean.lambdas.copy()
        )
    else:
        seed_2v, _ = cluster.create_unpinned_vortex_pair_seed(pure_in, p1, p2)
    omega, _, w_pos, e_2v, metrics = solve_inhomogeneous_torus_sc(
        cluster,
        params,
        seed_2v,
        reference_clean_state=pure_in,
        kappa_static=kappa_static,
        use_eigh=use_eigh,
        tol=tol,
    )

    fluxes = cluster.compute_plaquette_fluxes(omega)
    defs = sum(1 for wp in fluxes.values() if wp < 0.0)
    min_abs_pf = float(
        min(abs(cluster.compute_plaquette_pfaffian_flux(omega, p)) for p in (p1, p2))
    )

    w_sorted = np.sort(w_pos)
    eps_min = float(w_sorted[0])
    omega_mzm_can = 2.0 * eps_min
    delta_bulk = (
        12.0 * float(np.sqrt(3.0)) * abs(kappa_static) if kappa_static != 0.0 else delta_clean_torus
    )
    e0_over_delta = (eps_min / delta_bulk) if delta_bulk > 0.0 else 0.0
    omega_over_delta = (omega_mzm_can / delta_bulk) if delta_bulk > 0.0 else 0.0

    c_mzm, c_bulk, _ = cluster.compute_vortex_mode_core_weight(metrics.evecs[:, 0], p1, p2)

    return TwoVisonPointResult(
        L=L,
        resonance_class="3|L" if L % 3 == 0 else "3!|L",
        separation=r_sep,
        kappa_static=kappa_static,
        h_field=h_field,
        converged=metrics.converged,
        iterations=metrics.iterations,
        max_gauge_residual=metrics.max_gauge_residual,
        defect_count=defs,
        min_abs_pfaffian_wp=min_abs_pf,
        e_total=e_2v,
        v_pair=e_2v - e0_clean,
        omega_mzm_can=omega_mzm_can,
        delta_bulk_can=delta_bulk,
        e0_over_delta=e0_over_delta,
        omega_over_delta=omega_over_delta,
        mzm_core_weight=c_mzm,
        bulk_core_weight=c_bulk,
    )


def _fit_exp_profile(r: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Fit y(R) ~ A * exp(-R / xi) minimizing sum of squared residuals via exact xi profile."""
    xi_grid = np.linspace(0.4, 25.0, 2000)
    best_rss = float("inf")
    best_xi = 1.0
    best_amp = 0.0
    for xi in xi_grid:
        phi = np.exp(-r / xi)
        denom = float(np.dot(phi, phi))
        if denom <= 0.0:
            continue
        amp = float(np.dot(phi, y) / denom)
        rss = float(np.sum((y - amp * phi) ** 2))
        if rss < best_rss:
            best_rss = rss
            best_xi = float(xi)
            best_amp = amp
    # Refine locally around best_xi
    lo = max(0.2, best_xi * 0.9)
    hi = best_xi * 1.1
    for xi in np.linspace(lo, hi, 400):
        phi = np.exp(-r / xi)
        amp = float(np.dot(phi, y) / float(np.dot(phi, phi)))
        rss = float(np.sum((y - amp * phi) ** 2))
        if rss < best_rss:
            best_rss = rss
            best_xi = float(xi)
            best_amp = amp
    return best_amp, best_xi, best_rss


def _fit_pow_profile(r: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """Fit y(R) ~ C * R^(-alpha) minimizing sum of squared residuals via exact alpha profile."""
    alpha_grid = np.linspace(0.05, 8.0, 2000)
    best_rss = float("inf")
    best_alpha = 1.0
    best_amp = 0.0
    for alpha in alpha_grid:
        phi = r ** (-alpha)
        denom = float(np.dot(phi, phi))
        if denom <= 0.0:
            continue
        amp = float(np.dot(phi, y) / denom)
        rss = float(np.sum((y - amp * phi) ** 2))
        if rss < best_rss:
            best_rss = rss
            best_alpha = float(alpha)
            best_amp = amp
    lo = max(0.01, best_alpha * 0.9)
    hi = best_alpha * 1.1
    for alpha in np.linspace(lo, hi, 400):
        phi = r ** (-alpha)
        amp = float(np.dot(phi, y) / float(np.dot(phi, phi)))
        rss = float(np.sum((y - amp * phi) ** 2))
        if rss < best_rss:
            best_rss = rss
            best_alpha = float(alpha)
            best_amp = amp
    return best_amp, best_alpha, best_rss


def compare_exp_vs_powerlaw_scaling(rows: Sequence[TwoVisonPointResult]) -> ScalingComparisonResult:
    """Compare exp(-R/xi_M) vs R^(-alpha) and compute leave-L_min / leave-L_max holdout errors."""
    valid = [row for row in rows if row.separation >= 2 and row.omega_over_delta > 0.0]
    r = np.array([float(row.separation) for row in valid], dtype=np.float64)
    y = np.array([row.omega_over_delta for row in valid], dtype=np.float64)
    ls = np.array([row.L for row in valid], dtype=np.int64)

    amp_e, xi_fit, rss_exp = _fit_exp_profile(r, y)
    _, alpha_fit, rss_pow = _fit_pow_profile(r, y)
    n = len(valid)
    # Equal k parameters (amp, decay + noise variance => k_eff = 3): Delta_BIC == Delta_AICc
    delta_bic = float(n * np.log(max(rss_pow, 1e-30) / max(rss_exp, 1e-30)))
    delta_aicc = delta_bic

    l_min = int(np.min(ls))
    l_max = int(np.max(ls))

    mask_no_min = ls != l_min
    amp_no_min, xi_no_min, _ = _fit_exp_profile(r[mask_no_min], y[mask_no_min])
    pred_min = amp_no_min * np.exp(-r[~mask_no_min] / xi_no_min)
    err_lmin = float(np.linalg.norm(pred_min - y[~mask_no_min]) / np.linalg.norm(y[~mask_no_min]))

    mask_no_max = ls != l_max
    amp_no_max, xi_no_max, _ = _fit_exp_profile(r[mask_no_max], y[mask_no_max])
    pred_max = amp_no_max * np.exp(-r[~mask_no_max] / xi_no_max)
    err_lmax = float(np.linalg.norm(pred_max - y[~mask_no_max]) / np.linalg.norm(y[~mask_no_max]))
    _ = amp_e

    return ScalingComparisonResult(
        kappa_static=valid[0].kappa_static,
        n_points=n,
        xi_fit=xi_fit,
        alpha_fit=alpha_fit,
        rss_exp=rss_exp,
        rss_pow=rss_pow,
        delta_bic_exp_over_pow=delta_bic,
        delta_aicc_exp_over_pow=delta_aicc,
        leave_lmin_rel_err=err_lmin,
        leave_lmax_rel_err=err_lmax,
    )
