"""Wave 2C (UOW-05 / Q-L2-1): Self-consistent inhomogeneous disorder, vacancy stability & P_0 parity.

Implements:
1. Full inhomogeneous ``D = 80 L^2`` DIIS ``scMFT`` solver supporting Kenney-Laub Pade ``[3/2]``
   rational sign purification, cubic Newton-Schulz purification, and exact skew-symmetric ``eigh``
   under simultaneous quenched bond disorder ``delta_b``, non-magnetic site vacancies, covariant
   chiral ``kappa_static``, and isolated two-vison branch cuts.
2. Multi-seed quenched bond-disorder two-vison stability ensembles with BCa 95% CIs, exact
   12-Majorana Pfaffian flux tracking, KKT gauge residuals, and Pade ``[3/2]`` vs. Newton-Schulz
   iteration benchmarks.
3. Dilute vacancy stability benchmarks across both odd (``n_vac = 1``, ``N_act = 2 L^2 - 1``) and
   even (``n_vac = 0, 2``) active-site clusters, directly verifying the UOW-02 ``F2`` improper
   ``O(4)`` Pfaffian parity law ``Pf(P_0^{oplus N_act} Omega P_0^{oplus N_act}) = (-1)^{N_act} Pf(Omega)``,
   the exact single-site positivity diamond ``||V_s||_2 = 2(||<S_s>|| + ||<G_s^{v2}>||) <= 1``,
   and Willans-Chalker-Moessner / Kao-Knolle Eq. (8) vacancy-pair flux binding ``E_bind < 0``.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Sequence
from dataclasses import dataclass
from math import sqrt
from typing import Literal

import numpy as np

from kitaev_mf.domain.constants import G_MATRICES
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.qsl_bridge import P0_PARITY_BRIDGE, extract_site_su2_expectations
from kitaev_mf.domain.state import InhomogeneousMeanFieldState
from kitaev_mf.services.fixed_point import JointAndersonSolver
from kitaev_mf.services.gutzwiller_projection import bca_confidence_interval
from kitaev_mf.services.purification import purify_newton_schulz, purify_pade_sign
from kitaev_mf.services.real_space import RealSpaceTorusCluster, compute_skew_pfaffian
from kitaev_mf.services.vison_binding_wave2c import compute_inhomogeneous_disordered_energy
from kitaev_mf.services.vortex_pair_scaling import (
    bond_signs_from_state_gauge,
    build_covariant_nnn_matrix,
)

SolverMode = Literal["pade32", "newton_schulz", "eigh"]
_DISP_RSPACE: tuple[tuple[int, int, int], ...] = ((0, 0, 3), (0, -1, 1), (-1, 0, 2))
KAO_KNOLLE_VACANCY_BINDING_THERMO = -0.0268


@dataclass(frozen=True, slots=True)
class InhomogeneousDisorderSolveResult:
    """Full inhomogeneous scMFT solve diagnostics under bond disorder, vacancies, and flux defects."""

    L: int
    n_active_sites: int
    solver_mode: SolverMode
    converged: bool
    iterations: int
    final_residual: float
    max_gauge_residual: float
    kappa_active: float
    n_active_zero_modes: int
    e_total: float
    omega: np.ndarray
    state: InhomogeneousMeanFieldState
    w_pos: np.ndarray
    evecs: np.ndarray
    last_purify_iters: int


@dataclass(frozen=True, slots=True)
class DisorderedVisonEnsembleRow:
    """Multi-seed quenched bond-disorder summary for self-consistent two-vison trapping."""

    L: int
    separation: int
    kappa_static: float
    h_field: float
    delta_disorder: float
    n_seeds: int
    n_converged: int
    n_trapped_2v: int
    max_gauge_residual: float
    min_abs_pfaffian_wp: float
    v_pair_mean: float
    v_pair_ci_low: float
    v_pair_ci_high: float
    omega_mzm_mean: float
    omega_mzm_ci_low: float
    omega_mzm_ci_high: float
    core_enrichment_mean: float
    kappa_active_mean: float
    pade32_iters_mean: float
    ns_iters_mean: float
    max_pade_vs_eigh_diff: float


@dataclass(frozen=True, slots=True)
class VacancyStabilityRow:
    """Inhomogeneous scMFT vacancy stability, Pade vs. NS conditioning, and P_0 Pfaffian parity check."""

    L: int
    h_field: float
    delta_disorder: float
    n_vacancies: int
    n_active_sites: int
    converged_pade: bool
    sc_iters_pade: int
    max_gauge_residual: float
    kappa_active: float
    n_active_zero_modes: int
    pade32_purify_iters: int
    ns_purify_iters: int
    pade_vs_eigh_max_diff: float
    ns_vs_eigh_max_diff: float
    e_total: float
    e_double_counting_err: float
    pfaffian_kmf: float
    pfaffian_v2: float
    expected_p0_sign: int
    p0_pfaffian_parity_ok: bool
    max_positivity_diamond_err: float
    max_site_spectral_norm: float
    max_neighbor_spin_norm: float


@dataclass(frozen=True, slots=True)
class VacancyFluxBindingRow:
    """Willans-Chalker-Moessner / Kao-Knolle Eq. (8) A-B vacancy-pair bound-flux energy check."""

    L: int
    resonance_class: str
    vac_a_site: int
    vac_b_site: int
    cut_length: int
    e_zero_flux: float
    e_bound_flux: float
    e_bind_per_vacancy: float
    converged_zero: bool
    converged_bound: bool
    max_gauge_residual: float


@dataclass(frozen=True, slots=True)
class SingleVacancyCrossingDiagnostic:
    """Lakatosian diagnostic of pure-state T=0 KKT jump discontinuity on an unpaired vacancy."""

    L: int
    h_field: float
    lambda_critical_bare: float
    g_z_below: float
    g_z_above: float
    gauge_jump: float
    min_active_sv: float


def _active_gauge_residual_tensor(v_tensor: np.ndarray, active_mask: np.ndarray) -> np.ndarray:
    g_res = np.zeros((v_tensor.shape[0], 2, 3), dtype=np.float64)
    for gamma in (1, 2, 3):
        g_res[:, :, gamma - 1] = np.einsum("ij,csji->cs", G_MATRICES[gamma], v_tensor)
    return np.asarray(g_res * active_mask[:, :, None], dtype=np.float64)


def _active_spectral_diagnostics(
    a_mat: np.ndarray, active_sites: Sequence[int]
) -> tuple[float, int]:
    act_idx = [4 * s + mu for s in active_sites for mu in range(4)]
    a_act = a_mat[np.ix_(act_idx, act_idx)]
    sv = np.linalg.svd(a_act, compute_uv=False)
    n_zero = int(np.sum(sv < 1e-10))
    pos_sv = sv[sv >= 1e-10]
    kappa = float(pos_sv[0] / pos_sv[-1]) if pos_sv.size > 0 else float("inf")
    return kappa, n_zero


def _build_vacant_covariant_nnn(
    cluster: RealSpaceTorusCluster,
    kappa_static: float,
    bond_signs: np.ndarray,
    vac_set: set[int],
) -> np.ndarray:
    a_nnn = build_covariant_nnn_matrix(cluster, kappa_static, bond_signs)
    if vac_set:
        for s in vac_set:
            sl = cluster.majorana_slice(s)
            a_nnn[sl, :] = 0.0
            a_nnn[:, sl] = 0.0
    return np.asarray(0.5 * (a_nnn - a_nnn.T), dtype=np.float64)


def solve_inhomogeneous_sc_general(
    cluster: RealSpaceTorusCluster,
    params: ModelParameters,
    initial_state: InhomogeneousMeanFieldState,
    *,
    reference_clean_state: InhomogeneousMeanFieldState | None = None,
    solver_mode: SolverMode = "pade32",
    kappa_static: float = 0.0,
    bond_disorder: dict[tuple[int, int, int], float] | None = None,
    site_vacancies: set[int] | Sequence[int] | None = None,
    flipped_bonds: Sequence[tuple[int, int, int]] | None = None,
    max_iter: int = 80,
    tol: float = 1e-10,
    mixing_alpha: float = 0.5,
    lambda_step: float = 0.5,
    diis_history_depth: int = 8,
    purify_tol: float = 1e-12,
    purify_max_iter: int = 40,
) -> InhomogeneousDisorderSolveResult:
    """Solve full inhomogeneous D = 80 L^2 scMFT via Pade [3/2], Newton-Schulz, or eigh + DIIS."""
    vac_set = set(site_vacancies) if site_vacancies is not None else set()
    ref_state = (
        reference_clean_state
        if reference_clean_state is not None
        else InhomogeneousMeanFieldState.from_pure_kitaev(cluster.num_cells)
    )
    active_mask = np.ones((cluster.num_cells, 2), dtype=np.float64)
    active_sites: list[int] = []
    for s in range(cluster.num_sites):
        if s in vac_set:
            active_mask[s // 2, s % 2] = 0.0
        else:
            active_sites.append(s)

    curr_state = initial_state
    accel = JointAndersonSolver(history_depth=diis_history_depth, beta=mixing_alpha)
    converged = False
    iters = 0
    final_res = 1.0
    last_purify_iters = 0
    last_a_nnn = np.zeros((cluster.num_majoranas, cluster.num_majoranas), dtype=np.float64)
    cached_signs: np.ndarray | None = None
    prev_omega: np.ndarray | None = None

    for it in range(max_iter):
        iters = it + 1
        a_mf = cluster.assemble_skew_hamiltonian(
            params,
            curr_state,
            flipped_bonds=flipped_bonds,
            bond_disorder=bond_disorder,
            site_vacancies=vac_set,
        )
        if kappa_static != 0.0:
            signs = bond_signs_from_state_gauge(
                curr_state, cluster.L, ref_state, flipped_bonds=flipped_bonds
            )
            if cached_signs is None or not np.array_equal(signs, cached_signs):
                last_a_nnn = _build_vacant_covariant_nnn(cluster, kappa_static, signs, vac_set)
                cached_signs = signs
            a_tot = a_mf + last_a_nnn
        else:
            a_tot = a_mf

        if solver_mode == "pade32":
            omega_next, pmet = purify_pade_sign(
                a_tot, tol=purify_tol, max_iter=purify_max_iter, order="32"
            )
            last_purify_iters = pmet.iterations
        elif solver_mode == "newton_schulz":
            omega_next, pmet = purify_newton_schulz(a_tot, tol=purify_tol, max_iter=purify_max_iter)
            last_purify_iters = pmet.iterations
        else:
            omega_next, _ = cluster.solve_covariance_complex_eigh(a_tot)
            last_purify_iters = 0
        prev_omega = omega_next

        raw_state = cluster.extract_inhomogeneous_state(omega_next, curr_state.lambdas)
        g_res = _active_gauge_residual_tensor(raw_state.V, active_mask)
        lam_target = (curr_state.lambdas - lambda_step * g_res) * active_mask[:, :, None]
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
        next_st = InhomogeneousMeanFieldState.from_flattened(x_next, cluster.num_cells)
        curr_state = InhomogeneousMeanFieldState(
            U=next_st.U,
            V=next_st.V,
            lambdas=next_st.lambdas * active_mask[:, :, None],
        )

    assert prev_omega is not None
    a_final = cluster.assemble_skew_hamiltonian(
        params,
        curr_state,
        flipped_bonds=flipped_bonds,
        bond_disorder=bond_disorder,
        site_vacancies=vac_set,
    )
    if kappa_static != 0.0:
        signs = bond_signs_from_state_gauge(
            curr_state, cluster.L, ref_state, flipped_bonds=flipped_bonds
        )
        last_a_nnn = _build_vacant_covariant_nnn(cluster, kappa_static, signs, vac_set)
        a_final = a_final + last_a_nnn

    s_mat = 0.5 * ((-a_final @ a_final) + (-a_final @ a_final).T)
    w2, vecs = np.linalg.eigh(s_mat)
    w_pos = np.sqrt(np.maximum(w2, 0.0))
    kappa_act, n_zero_act = _active_spectral_diagnostics(a_final, active_sites)
    final_raw = cluster.extract_inhomogeneous_state(prev_omega, curr_state.lambdas)
    max_g_final = float(np.max(np.abs(_active_gauge_residual_tensor(final_raw.V, active_mask))))

    e_hf = compute_inhomogeneous_disordered_energy(
        curr_state,
        cluster.L,
        params,
        bond_disorder=bond_disorder,
        site_vacancies=vac_set,
    ).E_total
    e_nnn = 0.25 * float(np.sum(last_a_nnn * prev_omega)) if kappa_static != 0.0 else 0.0
    return InhomogeneousDisorderSolveResult(
        L=cluster.L,
        n_active_sites=len(active_sites),
        solver_mode=solver_mode,
        converged=converged,
        iterations=iters,
        final_residual=final_res,
        max_gauge_residual=max_g_final,
        kappa_active=kappa_act,
        n_active_zero_modes=n_zero_act,
        e_total=e_hf + e_nnn,
        omega=prev_omega,
        state=curr_state,
        w_pos=w_pos,
        evecs=vecs,
        last_purify_iters=last_purify_iters,
    )


def evaluate_bond_disordered_two_vison_ensemble(
    L: int,
    r_sep: int,
    delta_levels: Sequence[float],
    *,
    kappa_static: float = 0.06,
    h_field: float = 0.0,
    n_seeds: int = 8,
    seed_base: int = 202611,
    tol: float = 1e-10,
) -> list[DisorderedVisonEnsembleRow]:
    """Audit self-consistent two-vison trapping & MZM stability under quenched bond disorder."""
    cluster = RealSpaceTorusCluster(L)
    comp = h_field / sqrt(3.0)
    params = ModelParameters(K=-1.0, h=(comp, comp, comp))
    pure_in = InhomogeneousMeanFieldState.from_pure_kitaev(cluster.num_cells)
    p1 = (1, 1)
    p2 = ((1 + r_sep) % L, 1)
    cut = cluster.generate_vortex_pair_cut(p1, p2)

    rows: list[DisorderedVisonEnsembleRow] = []
    for d_idx, delta in enumerate(delta_levels):
        actual_seeds = 1 if delta == 0.0 else n_seeds
        v_vals: list[float] = []
        om_vals: list[float] = []
        core_ratios: list[float] = []
        kappas: list[float] = []
        pade_its: list[int] = []
        ns_its: list[int] = []
        n_conv = 0
        n_trapped = 0
        max_g_all = 0.0
        min_pf_all = 1.0
        max_pd_err = 0.0

        for s_idx in range(actual_seeds):
            rng = np.random.default_rng(seed_base + 1000 * d_idx + s_idx)
            b_dis = (
                {
                    (m1, m2, b_i): float(rng.uniform(-delta, delta))
                    for m2 in range(L)
                    for m1 in range(L)
                    for b_i in range(3)
                }
                if delta > 0.0
                else {}
            )
            res0 = solve_inhomogeneous_sc_general(
                cluster,
                params,
                pure_in,
                reference_clean_state=pure_in,
                solver_mode="pade32",
                kappa_static=kappa_static,
                bond_disorder=b_dis,
                tol=tol,
            )
            if h_field == 0.0:
                u_2v = res0.state.U.copy()
                for m1, m2, b_i in cut:
                    u_2v[m1 + L * m2, b_i] = -u_2v[m1 + L * m2, b_i]
                seed_2v = InhomogeneousMeanFieldState(
                    U=u_2v, V=res0.state.V.copy(), lambdas=res0.state.lambdas.copy()
                )
            else:
                seed_2v, _ = cluster.create_unpinned_vortex_pair_seed(res0.state, p1, p2)

            res2v = solve_inhomogeneous_sc_general(
                cluster,
                params,
                seed_2v,
                reference_clean_state=pure_in,
                solver_mode="pade32",
                kappa_static=kappa_static,
                bond_disorder=b_dis,
                tol=tol,
            )

            # Direct Pade [3/2] vs. Newton-Schulz vs. eigh check on converged disordered two-vison A
            a_conv = cluster.assemble_skew_hamiltonian(params, res2v.state, bond_disorder=b_dis)
            if kappa_static != 0.0:
                sgn = bond_signs_from_state_gauge(res2v.state, L, pure_in)
                a_conv = a_conv + build_covariant_nnn_matrix(cluster, kappa_static, sgn)
            om_pd, m_pd = purify_pade_sign(a_conv, tol=1e-12, max_iter=40, order="32")
            _, m_ns = purify_newton_schulz(a_conv, tol=1e-12, max_iter=50)
            om_eg, _ = cluster.solve_covariance_complex_eigh(a_conv)
            max_pd_err = max(max_pd_err, float(np.max(np.abs(om_pd - om_eg))))
            pade_its.append(m_pd.iterations)
            ns_its.append(m_ns.iterations)

            fluxes = cluster.compute_plaquette_fluxes(res2v.omega)
            defs = [p for p, wp in fluxes.items() if wp < 0.0]
            min_pf = float(
                min(abs(cluster.compute_plaquette_pfaffian_flux(res2v.omega, p)) for p in (p1, p2))
            )
            if res0.converged and res2v.converged:
                n_conv += 1
            if set(defs) == {p1, p2}:
                n_trapped += 1
            max_g_all = max(max_g_all, res0.max_gauge_residual, res2v.max_gauge_residual)
            min_pf_all = min(min_pf_all, min_pf)

            v_vals.append(res2v.e_total - res0.e_total)
            om_vals.append(2.0 * float(np.min(res2v.w_pos)))
            _, _, enrich = cluster.compute_vortex_mode_core_weight(res2v.evecs[:, 0], p1, p2)
            core_ratios.append(enrich)
            kappas.append(res2v.kappa_active)

        v_arr = np.asarray(v_vals, dtype=np.float64)
        om_arr = np.asarray(om_vals, dtype=np.float64)
        v_lo, v_hi = bca_confidence_interval(v_arr, seed=seed_base + d_idx)
        om_lo, om_hi = bca_confidence_interval(om_arr, seed=seed_base + 100 + d_idx)
        rows.append(
            DisorderedVisonEnsembleRow(
                L=L,
                separation=r_sep,
                kappa_static=float(kappa_static),
                h_field=float(h_field),
                delta_disorder=float(delta),
                n_seeds=actual_seeds,
                n_converged=n_conv,
                n_trapped_2v=n_trapped,
                max_gauge_residual=max_g_all,
                min_abs_pfaffian_wp=min_pf_all,
                v_pair_mean=float(np.mean(v_arr)),
                v_pair_ci_low=v_lo,
                v_pair_ci_high=v_hi,
                omega_mzm_mean=float(np.mean(om_arr)),
                omega_mzm_ci_low=om_lo,
                omega_mzm_ci_high=om_hi,
                core_enrichment_mean=float(np.mean(core_ratios)),
                kappa_active_mean=float(np.mean(kappas)),
                pade32_iters_mean=float(np.mean(pade_its)),
                ns_iters_mean=float(np.mean(ns_its)),
                max_pade_vs_eigh_diff=max_pd_err,
            )
        )
    return rows


def evaluate_vacancy_stability_and_p0_parity(
    L: int,
    h_values: Sequence[float],
    vacancy_counts: Sequence[int] = (0, 1, 2),
    *,
    delta_disorder: float = 0.0,
    seed: int = 202611,
    tol: float = 1e-10,
) -> list[VacancyStabilityRow]:
    """Audit inhomogeneous vacancy scMFT, Pade vs. NS conditioning, P_0 Pfaffian parity & positivity."""
    cluster = RealSpaceTorusCluster(L)
    pure_in = InhomogeneousMeanFieldState.from_pure_kitaev(cluster.num_cells)
    rng = np.random.default_rng(seed)
    b_dis = (
        {
            (m1, m2, b_i): float(rng.uniform(-delta_disorder, delta_disorder))
            for m2 in range(L)
            for m1 in range(L)
            for b_i in range(3)
        }
        if delta_disorder > 0.0
        else {}
    )

    rows: list[VacancyStabilityRow] = []
    for h_val in h_values:
        comp = h_val / sqrt(3.0)
        params = ModelParameters(K=-1.0, h=(comp, comp, comp))
        for n_vac in vacancy_counts:
            if n_vac == 0:
                vacs: set[int] = set()
            elif n_vac == 1:
                vacs = {0}
            elif n_vac == 2:
                vacs = {0, cluster.site_index(L // 2, L // 2, 1)}
            else:
                vacs = set(range(n_vac))

            res = solve_inhomogeneous_sc_general(
                cluster,
                params,
                pure_in,
                reference_clean_state=pure_in,
                solver_mode="pade32",
                bond_disorder=b_dis,
                site_vacancies=vacs,
                max_iter=120,
                tol=tol,
                mixing_alpha=0.5,
                lambda_step=0.55,
                purify_tol=1e-12,
                purify_max_iter=45,
            )
            a_conv = cluster.assemble_skew_hamiltonian(
                params, res.state, bond_disorder=b_dis, site_vacancies=vacs
            )
            om_pd, m_pd = purify_pade_sign(a_conv, tol=1e-12, max_iter=45, order="32")
            om_ns, m_ns = purify_newton_schulz(a_conv, tol=1e-12, max_iter=45)
            om_eg, _ = cluster.solve_covariance_complex_eigh(a_conv)

            # Exact K-007 double-counting check: E_HF == <H_quad> - E_exchange
            eb = compute_inhomogeneous_disordered_energy(
                res.state, L, params, bond_disorder=b_dis, site_vacancies=vacs
            )
            e_quad = 0.25 * float(np.sum(a_conv * om_eg))
            dc_err = abs((e_quad - eb.E_exchange) - eb.E_total)

            # Active-site O(4) improper Pfaffian parity law: Pf(P0^N Omega P0^N) = (-1)^N Pf(Omega)
            act_sites = [s for s in range(cluster.num_sites) if s not in vacs]
            act_idx = [4 * s + mu for s in act_sites for mu in range(4)]
            om_act = om_eg[np.ix_(act_idx, act_idx)]
            p0_big = np.kron(np.eye(len(act_sites), dtype=np.float64), P0_PARITY_BRIDGE)
            om_act_v2 = p0_big @ om_act @ p0_big
            pf_kmf = compute_skew_pfaffian(om_act)
            pf_v2 = compute_skew_pfaffian(om_act_v2)
            exp_sign = -1 if (len(act_sites) % 2 == 1) else 1
            if abs(pf_kmf) > 1e-10:
                parity_ok = abs(pf_v2 - exp_sign * pf_kmf) < 1e-10 * abs(pf_kmf)
            else:
                parity_ok = abs(pf_v2) < 1e-10

            # Site positivity diamond ||V_s||_2 = 2*(||<S_s>|| + ||<G_s^v2>||) <= 1
            max_diam_err = 0.0
            max_spec_norm = 0.0
            max_spin_norm = 0.0
            for s in act_sites:
                v_s = res.state.V[s // 2, s % 2]
                dec = extract_site_su2_expectations(v_s)
                s_norm = float(np.linalg.norm(dec.spin_s))
                g_norm = float(np.linalg.norm(dec.gauge_g))
                v_norm = float(np.linalg.norm(v_s, ord=2))
                max_diam_err = max(max_diam_err, abs(v_norm - 2.0 * (s_norm + g_norm)))
                max_spec_norm = max(max_spec_norm, v_norm)
                max_spin_norm = max(max_spin_norm, s_norm)

            rows.append(
                VacancyStabilityRow(
                    L=L,
                    h_field=float(h_val),
                    delta_disorder=float(delta_disorder),
                    n_vacancies=len(vacs),
                    n_active_sites=len(act_sites),
                    converged_pade=res.converged,
                    sc_iters_pade=res.iterations,
                    max_gauge_residual=res.max_gauge_residual,
                    kappa_active=res.kappa_active,
                    n_active_zero_modes=res.n_active_zero_modes,
                    pade32_purify_iters=m_pd.iterations,
                    ns_purify_iters=m_ns.iterations,
                    pade_vs_eigh_max_diff=float(np.max(np.abs(om_pd - om_eg))),
                    ns_vs_eigh_max_diff=float(np.max(np.abs(om_ns - om_eg))),
                    e_total=eb.E_total,
                    e_double_counting_err=dc_err,
                    pfaffian_kmf=pf_kmf,
                    pfaffian_v2=pf_v2,
                    expected_p0_sign=exp_sign,
                    p0_pfaffian_parity_ok=parity_ok,
                    max_positivity_diamond_err=max_diam_err,
                    max_site_spectral_norm=max_spec_norm,
                    max_neighbor_spin_norm=max_spin_norm,
                )
            )
    return rows


def evaluate_vacancy_pair_flux_binding(
    L: int,
    *,
    tol: float = 1e-11,
) -> VacancyFluxBindingRow:
    """Measure Willans-Chalker-Moessner / Kao-Knolle Eq. (8) A-B vacancy-pair bound-flux energy."""
    cluster = RealSpaceTorusCluster(L)
    params0 = ModelParameters(K=-1.0)
    pure_in = InhomogeneousMeanFieldState.from_pure_kitaev(cluster.num_cells)

    vac_a = cluster.site_index(0, 0, 0)
    vac_b = cluster.site_index(L // 2, L // 2, 1)
    vacs = {vac_a, vac_b}

    def _plaqs_touching(site: int) -> list[tuple[int, int]]:
        return [
            (m1, m2)
            for m2 in range(L)
            for m1 in range(L)
            if site in cluster.compute_plaquette_boundary_sites((m1, m2))
        ]

    def _bond_touches_vac(bk: tuple[int, int, int]) -> bool:
        m1, m2, b_i = bk
        d1, d2, _ = _DISP_RSPACE[b_i]
        return (
            cluster.site_index(m1, m2, 0) in vacs or cluster.site_index(m1 + d1, m2 + d2, 1) in vacs
        )

    p_start = _plaqs_touching(vac_a)[0]
    p_targets = set(_plaqs_touching(vac_b))
    adj = cluster.build_dual_graph()
    visited = {p_start}
    queue: deque[tuple[tuple[int, int], list[tuple[int, int, int]]]] = deque([(p_start, [])])
    cut_bound: list[tuple[int, int, int]] = []
    while queue:
        curr, path = queue.popleft()
        if curr in p_targets:
            cut_bound = path
            break
        for nbr, bk in adj[curr]:
            if nbr not in visited and not _bond_touches_vac(bk):
                visited.add(nbr)
                queue.append((nbr, [*path, bk]))
    if not cut_bound:
        raise RuntimeError(f"Failed to construct active dual cut between vacancies on L={L}.")

    res_zero = solve_inhomogeneous_sc_general(
        cluster,
        params0,
        pure_in,
        reference_clean_state=pure_in,
        solver_mode="pade32",
        site_vacancies=vacs,
        tol=tol,
    )
    u_bound = res_zero.state.U.copy()
    for m1, m2, b_i in cut_bound:
        u_bound[m1 + L * m2, b_i] = -u_bound[m1 + L * m2, b_i]
    seed_bound = InhomogeneousMeanFieldState(
        U=u_bound, V=res_zero.state.V.copy(), lambdas=res_zero.state.lambdas.copy()
    )
    res_bound = solve_inhomogeneous_sc_general(
        cluster,
        params0,
        seed_bound,
        reference_clean_state=pure_in,
        solver_mode="pade32",
        site_vacancies=vacs,
        tol=tol,
    )
    return VacancyFluxBindingRow(
        L=L,
        resonance_class="3|L" if L % 3 == 0 else "3!|L",
        vac_a_site=vac_a,
        vac_b_site=vac_b,
        cut_length=len(cut_bound),
        e_zero_flux=res_zero.e_total,
        e_bound_flux=res_bound.e_total,
        e_bind_per_vacancy=0.5 * (res_bound.e_total - res_zero.e_total),
        converged_zero=res_zero.converged,
        converged_bound=res_bound.converged,
        max_gauge_residual=max(res_zero.max_gauge_residual, res_bound.max_gauge_residual),
    )


def probe_single_vacancy_kkt_level_crossing(
    L: int = 4,
    h_field: float = 0.10,
    delta_lambda: float = 0.01,
) -> SingleVacancyCrossingDiagnostic:
    """Demonstrate pure-state T=0 KKT discontinuity across ``lambda_z ~ 2 h / sqrt(3)`` on a single vacancy.

    Removing a single A-site ``s = 0`` leaves a dangling ``b^z`` Majorana on neighbor site ``1``
    coupled to ``c^0`` with coefficient ``lambda_z - 2 h_z + Hartree``. On gapped non-resonant
    (``3 !| L``) or C_3-broken disordered clusters at ``h > 0``, ``<G_1^z>(lambda_z)`` exhibits a
    discontinuous sign jump (> 1.0) bracketing zero, so no pure ``T = 0`` Gaussian state
    (``Omega^2 = -I``) satisfies ``<G_1^z> = 0``.
    """
    cluster = RealSpaceTorusCluster(L)
    pure_in = InhomogeneousMeanFieldState.from_pure_kitaev(cluster.num_cells)
    comp = h_field / sqrt(3.0)
    params = ModelParameters(K=-1.0, h=(comp, comp, comp))
    vacs = {0}
    res = solve_inhomogeneous_sc_general(
        cluster,
        params,
        pure_in,
        reference_clean_state=pure_in,
        solver_mode="eigh",
        site_vacancies=vacs,
        max_iter=50,
        tol=1e-10,
        mixing_alpha=0.5,
        lambda_step=0.2,
    )
    lam_crit = 2.0 * comp
    act_idx = [4 * s + mu for s in range(1, cluster.num_sites) for mu in range(4)]

    def _eval_at(lam_z: float) -> tuple[float, float]:
        lam_try = res.state.lambdas.copy()
        lam_try[0, 1, 2] = lam_z
        st_try = InhomogeneousMeanFieldState(U=res.state.U, V=res.state.V, lambdas=lam_try)
        a_mat = cluster.assemble_skew_hamiltonian(params, st_try, site_vacancies=vacs)
        om, _ = cluster.solve_covariance_complex_eigh(a_mat)
        g_z = float(np.sum(G_MATRICES[3] * om[4:8, 4:8].T))
        sv_min = float(np.min(np.linalg.svd(a_mat[np.ix_(act_idx, act_idx)], compute_uv=False)))
        return g_z, sv_min

    g_below, sv_below = _eval_at(lam_crit - delta_lambda)
    g_above, sv_above = _eval_at(lam_crit + delta_lambda)
    return SingleVacancyCrossingDiagnostic(
        L=L,
        h_field=float(h_field),
        lambda_critical_bare=lam_crit,
        g_z_below=g_below,
        g_z_above=g_above,
        gauge_jump=g_above - g_below,
        min_active_sv=min(sv_below, sv_above),
    )


__all__ = [
    "KAO_KNOLLE_VACANCY_BINDING_THERMO",
    "DisorderedVisonEnsembleRow",
    "InhomogeneousDisorderSolveResult",
    "SingleVacancyCrossingDiagnostic",
    "SolverMode",
    "VacancyFluxBindingRow",
    "VacancyStabilityRow",
    "evaluate_bond_disordered_two_vison_ensemble",
    "evaluate_vacancy_pair_flux_binding",
    "evaluate_vacancy_stability_and_p0_parity",
    "probe_single_vacancy_kkt_level_crossing",
    "solve_inhomogeneous_sc_general",
]
