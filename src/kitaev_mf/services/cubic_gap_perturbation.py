"""Wave 1 HC-01: Cubic Dirac gap Feshbach-Schur decomposition and continuum extrapolation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field

K_DIRAC = np.array([[1.0 / 3.0, 2.0 / 3.0]], dtype=np.float64)


@dataclass(frozen=True, slots=True)
class CubicGapDecomposition:
    """Single-(h, N_k) Feshbach-Schur breakdown of the canonical Dirac gap 4*Delta_K^MF."""

    h: float
    nk: int
    converged: bool
    norm_lambda: float
    norm_lambda_over_h: float
    gauge_dimer_gap_can: float
    delta_can_direct: float
    delta_schur_full: float
    delta_schur_bare_zeeman: float
    delta_3rd_dressed: float
    delta_3rd_qq_lambda_only: float
    delta_3rd_zeeman_uv_no_lambda: float
    delta_3rd_lambda_only: float
    delta_direct_over_h3: float
    norm_dv_pp: float


@dataclass(frozen=True, slots=True)
class ContinuumCubicSummary:
    """Continuum N_k -> infinity extrapolated prefactor and empirical exponent."""

    nk_grid: tuple[int, ...]
    c_mf_by_nk: dict[int, float]
    c_mf_continuum: float
    lambda_over_h_continuum: float
    exponent_fit: float


def _bare_zeeman_state_with_lambda_only(st0: MeanFieldState, lambdas: np.ndarray) -> MeanFieldState:
    return MeanFieldState(U=st0.U.copy(), V=st0.V.copy(), lambdas=lambdas.copy())


def analyze_cubic_gap_at_point(
    h: float,
    nk: int,
    *,
    initial_state: MeanFieldState | None = None,
    st0: MeanFieldState | None = None,
    tol: float = 1e-12,
) -> tuple[CubicGapDecomposition, MeanFieldState]:
    """Solve scMFT at `(h, nk x nk)` and decompose Delta_can(K) via Feshbach-Schur."""
    cell = make_honeycomb_2site()
    cfg = SolverConfig(
        k_grid_shape=(nk, nk),
        max_iterations=80,
        tol_correlator=tol,
        tol_gauge=tol,
        method="diis",
        mixing_alpha=0.6,
        skew_mode="retract_step",
    )
    if st0 is None:
        st0_sol, _ = solve_self_consistent_mean_field(
            cell, ModelParameters(K=-1.0), make_qsl_seed(cell), cfg
        )
    else:
        st0_sol = st0

    h0 = 4.0 * build_bloch_hamiltonian(K_DIRAC, cell, ModelParameters(K=-1.0), st0_sol)[0].numpy()
    evals0, evecs0 = np.linalg.eigh(h0)
    zero_mask = np.abs(evals0) < 1e-9
    p_idx = np.where(zero_mask)[0]
    q_idx = np.where(~zero_mask)[0]
    if p_idx.size != 2 or q_idx.size != 6:
        raise RuntimeError(f"Unexpected Dirac zero-space partitioning: P={p_idx}, Q={q_idx}")
    gauge_gap = float(2.0 * np.min(np.abs(evals0[q_idx])))

    comp = h / float(np.sqrt(3.0))
    params = ModelParameters(K=-1.0, h=(comp, comp, comp))
    seed = st0_sol if initial_state is None else initial_state
    st_h, met = solve_self_consistent_mean_field(cell, params, seed, cfg)

    hh = 4.0 * build_bloch_hamiltonian(K_DIRAC, cell, params, st_h)[0].numpy()
    evals_h = np.linalg.eigvalsh(hh)
    delta_can = float(evals_h[4] - evals_h[3])

    hh_rot = evecs0.conj().T @ hh @ evecs0
    h_pp = hh_rot[np.ix_(p_idx, p_idx)]
    h_pq = hh_rot[np.ix_(p_idx, q_idx)]
    h_qq = hh_rot[np.ix_(q_idx, q_idx)]
    s_full = h_pp - h_pq @ np.linalg.solve(h_qq, h_pq.conj().T)
    delta_schur_full = float(np.ptp(np.linalg.eigvalsh(s_full)))

    # Bare Zeeman perturbation with unperturbed U^(0), V^(0)=0, lambda=0
    hb = 4.0 * build_bloch_hamiltonian(K_DIRAC, cell, params, st0_sol)[0].numpy()
    hb_rot = evecs0.conj().T @ hb @ evecs0
    s_bare = hb_rot[np.ix_(p_idx, p_idx)] - hb_rot[np.ix_(p_idx, q_idx)] @ np.linalg.solve(
        hb_rot[np.ix_(q_idx, q_idx)], hb_rot[np.ix_(q_idx, p_idx)]
    )
    delta_schur_bare = float(np.ptp(np.linalg.eigvalsh(s_bare)))

    # Leading 3rd-order dressed Feshbach-Schur:
    # H_FS^(3) = P dV_PP P + P dV_PQ D0^-1 dV_QQ D0^-1 dV_QP P
    dv = hh_rot - np.diag(evals0)
    dv_pp = dv[np.ix_(p_idx, p_idx)]
    dv_pq = dv[np.ix_(p_idx, q_idx)]
    dv_qq = dv[np.ix_(q_idx, q_idx)]
    d0_inv = np.diag(1.0 / evals0[q_idx])
    h3_dressed = dv_pp + dv_pq @ d0_inv @ dv_qq @ d0_inv @ dv_pq.conj().T
    delta_3rd_dressed = float(np.ptp(np.linalg.eigvalsh(h3_dressed)))

    # Pure KKT lambda* operator inside Q-space paired with full PQ vertex
    st_only_lam = MeanFieldState(
        U=st0_sol.U.copy(), V=st0_sol.V.copy(), lambdas=st_h.lambdas.copy()
    )
    h_only_lam = (
        4.0
        * build_bloch_hamiltonian(K_DIRAC, cell, ModelParameters(K=-1.0), st_only_lam)[0].numpy()
    )
    dv_only_lam = evecs0.conj().T @ (h_only_lam - h0) @ evecs0
    h3_qq_lam = dv_pq @ d0_inv @ dv_only_lam[np.ix_(q_idx, q_idx)] @ d0_inv @ dv_pq.conj().T
    delta_3rd_qq_lam = float(np.ptp(np.linalg.eigvalsh(h3_qq_lam)))

    # Full Zeeman + Hartree (U, V) without KKT lambda* (verifies Hartree alone generates zero cubic gap)
    st_no_lam = MeanFieldState(
        U=st_h.U.copy(), V=st_h.V.copy(), lambdas=np.zeros_like(st_h.lambdas)
    )
    h_no_lam = 4.0 * build_bloch_hamiltonian(K_DIRAC, cell, params, st_no_lam)[0].numpy()
    dv_no_lam = evecs0.conj().T @ (h_no_lam - h0) @ evecs0
    h3_no_lam = (
        dv_no_lam[np.ix_(p_idx, q_idx)]
        @ d0_inv
        @ dv_no_lam[np.ix_(q_idx, q_idx)]
        @ d0_inv
        @ dv_no_lam[np.ix_(q_idx, p_idx)]
    )
    delta_3rd_no_lam = float(np.ptp(np.linalg.eigvalsh(h3_no_lam)))

    # Isolate strictly Zeeman + KKT lambda* G^gamma (unscreened by Hartree dU/dV)
    st_lam_only = _bare_zeeman_state_with_lambda_only(st0_sol, st_h.lambdas)
    hl = 4.0 * build_bloch_hamiltonian(K_DIRAC, cell, params, st_lam_only)[0].numpy()
    dvl = evecs0.conj().T @ hl @ evecs0 - np.diag(evals0)
    h3_lam = (
        dvl[np.ix_(p_idx, q_idx)]
        @ d0_inv
        @ dvl[np.ix_(q_idx, q_idx)]
        @ d0_inv
        @ dvl[np.ix_(q_idx, p_idx)]
    )
    delta_3rd_lam = float(np.ptp(np.linalg.eigvalsh(h3_lam)))

    norm_lam = float(np.linalg.norm(st_h.lambdas))
    return (
        CubicGapDecomposition(
            h=float(h),
            nk=nk,
            converged=bool(met.converged),
            norm_lambda=norm_lam,
            norm_lambda_over_h=norm_lam / h if h > 0.0 else 0.0,
            gauge_dimer_gap_can=gauge_gap,
            delta_can_direct=delta_can,
            delta_schur_full=delta_schur_full,
            delta_schur_bare_zeeman=delta_schur_bare,
            delta_3rd_dressed=delta_3rd_dressed,
            delta_3rd_qq_lambda_only=delta_3rd_qq_lam,
            delta_3rd_zeeman_uv_no_lambda=delta_3rd_no_lam,
            delta_3rd_lambda_only=delta_3rd_lam,
            delta_direct_over_h3=delta_can / (h**3) if h > 0.0 else 0.0,
            norm_dv_pp=float(np.max(np.abs(dv_pp))),
        ),
        st_h,
    )


def extrapolate_cubic_prefactor_continuum(
    h_ref: float = 0.01,
    nk_grid: tuple[int, ...] = (24, 36, 48, 72),
    h_exponent_grid: tuple[float, ...] = (0.002, 0.005, 0.008, 0.012, 0.016, 0.020, 0.030),
) -> tuple[list[CubicGapDecomposition], ContinuumCubicSummary]:
    """Run deterministic k-mesh refinement and small-h power-law exponent check."""
    rows: list[CubicGapDecomposition] = []
    c_by_nk: dict[int, float] = {}
    lam_by_nk: list[float] = []

    for nk in nk_grid:
        decomp, _ = analyze_cubic_gap_at_point(h_ref, nk)
        rows.append(decomp)
        c_by_nk[nk] = decomp.delta_direct_over_h3
        lam_by_nk.append(decomp.norm_lambda_over_h)

    inv_nk = np.array([1.0 / float(nk) for nk in nk_grid], dtype=np.float64)
    # Use quadratic-in-(1/Nk) Richardson/Vandermonde extrapolation to 1/Nk -> 0 on finest grids
    x_ext = inv_nk[-3:]
    y_ext = np.array([c_by_nk[nk] for nk in nk_grid[-3:]], dtype=np.float64)
    c_cont = float(np.polyval(np.polyfit(x_ext, y_ext, 2), 0.0))
    lam_cont = float(np.polyval(np.polyfit(x_ext, np.array(lam_by_nk[-3:]), 2), 0.0))

    # Measure log-log slope on finest nk (nk_grid[-2]) across small-h grid
    nk_slope = nk_grid[-2]
    hs = []
    gaps = []
    curr_st = None
    for h in h_exponent_grid:
        decomp, curr_st = analyze_cubic_gap_at_point(h, nk_slope, initial_state=curr_st)
        rows.append(decomp)
        hs.append(h)
        gaps.append(decomp.delta_can_direct)
    slope = float(np.polyfit(np.log(np.array(hs)), np.log(np.array(gaps)), 1)[0])

    return rows, ContinuumCubicSummary(
        nk_grid=nk_grid,
        c_mf_by_nk=c_by_nk,
        c_mf_continuum=c_cont,
        lambda_over_h_continuum=lam_cont,
        exponent_fit=slope,
    )
