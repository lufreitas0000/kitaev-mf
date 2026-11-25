"""Wave 2B (UOW-03 / HC-02): Gutzwiller D_j=+1 VMC stability across staggered SU(2)_R & flux sectors.

Implements:
1. Exact site-local SU(2)_R gauge rotation engine ``O = bigoplus_s exp(theta_s . SO4_L)`` in the
   ``(-c, b)`` basis and direct verification of the Thouless block-diagonality theorem
   ``U_cf^dagger O U_cf = diag(U_O, U_O^*)`` (``U_O in SU(2)^N``, ``det(u_s) = +1``), proving exact
   reference condition-number invariance and bitstring-by-bitstring Pfaffian invariance.
2. Explicit Lakatosian contrast between direct quadratic conjugation ``O A O^T`` / pure-spin
   ``SO4_R`` hopping (continuous ``SU(2)_R`` covariant) and ``kitaev-mf``'s ``N^alpha = R^alpha - L^alpha``
   decoupling (breaks continuous ``SU(2)_R -> Z_2``, preserving only ``exp(2 pi n.L) = -I_4``).
3. Static Z_2 flux / vison sector Gutzwiller projection (``zero_flux``, ``two_vison``, ``pi_flux``)
   with automatic 1-quasiparticle particle-hole-paired parity repair resolving Pedrocchi-Chesi-Loss
   odd-parity obstructions (e.g. 2x3 pi-flux) and exact Zero-Variance verification.
4. Multi-chain MCMC autocorrelation & convergence metrology (Sokal inefficiency ``g_sokal``,
   Wolff ``tau_int_wolff = g_sokal / 2``, grand-mean Geyer Initial Positive Sequence ``g_ips_geyer``,
   split-Rhat, acceptance rates, and chain-vs-autocorr MCSE cross-check).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, sin, sqrt
from typing import Literal

import numpy as np
import torch
from numpy.typing import NDArray

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.qsl_bridge import SO4_L, SO4_R
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.gutzwiller_projection import (
    _build_physical_ref_bdg_transform,
    _exact_spin_observables,
    _integrated_autocorr,
    _run_vmc_from_thouless,
    _ThoulessStateEvaluator,
    assemble_torus_majorana_matrix,
    bca_confidence_interval,
    build_thouless_data,
    exact_wavefunction_from_thouless,
    select_admissible_twist,
    torus_bonds,
    torus_site,
)
from kitaev_mf.services.hamiltonian import build_on_site_matrices

FluxSectorName = Literal["zero_flux", "two_vison", "pi_flux"]


@dataclass(frozen=True, slots=True)
class AutocorrelationDiagnostics:
    """Multi-chain MCMC autocorrelation, equilibration, and acceptance diagnostics.

    Convention reconciliation with ``contracts/qsl_metrology_guidelines.md``:
    - ``g_sokal = 1 + 2 sum_{k>=1} rho(k) = 2 * tau_int_wolff`` (statistical inefficiency >= 1).
    - ``tau_int_wolff = g_sokal / 2`` so ``n_eff_sokal = N_tot / (2 * tau_int_wolff) = N_tot / g_sokal``.
    - ``g_ips_geyer`` uses Geyer's Initial Positive Sequence on the grand-mean PSD autocovariance.
    """

    g_sokal: float
    tau_int_wolff: float
    g_ips_geyer: float
    n_eff_sokal: float
    n_eff_ips: float
    split_rhat: float
    mcse_chains: float
    mcse_autocorr: float
    acc_rate_1flip: float
    acc_rate_2flip: float


@dataclass(frozen=True, slots=True)
class ThoulessSU2BlockCheck:
    """Direct algebraic check of ``U_cf^dagger O U_cf = diag(U_O, U_O^*)`` with ``U_O in SU(2)^N``."""

    n_sites: int
    max_off_block: float
    max_conj_err: float
    max_unitarity_err: float
    max_det_err: float


@dataclass(frozen=True, slots=True)
class SU2GaugeStabilityResult:
    """Exact and 50-chain VMC Gutzwiller stability record under a staggered SU(2)_R gauge rotation."""

    nx: int
    ny: int
    n_sites: int
    h: float
    bc_twist: tuple[int, int]
    pattern_name: str
    matrix_rotation_norm: float
    z_eff_diff_norm: float
    cond_base: float
    cond_rot: float
    cond_rel_diff: float
    max_probe_amp_diff: float
    max_probe_eloc_diff: float
    exact_fidelity: float | None
    delta_e_exact: float | None
    delta_m111_exact: float | None
    delta_wp_exact: float | None
    exact_energy_variance: float | None
    vmc_e_mean: float
    vmc_e_ci_low: float
    vmc_e_ci_high: float
    vmc_e_loc_std: float
    vmc_delta_e_vs_base: float
    vmc_m111_mean: float
    vmc_m111_ci_low: float
    vmc_m111_ci_high: float
    diagnostics: AutocorrelationDiagnostics


@dataclass(frozen=True, slots=True)
class NDecouplingBreakingResult:
    """Lakatosian contrast: ``O A O^T`` & pure-spin ``R^a`` covariance vs ``N^a`` breaking ``SU(2)_R -> Z_2``."""

    nx: int
    ny: int
    n_sites: int
    h: float
    direct_conjugation_fidelity: float
    pure_spin_r_decoupling_fidelity: float
    n_decoupling_staggered_fidelity: float
    n_decoupling_z2_fidelity: float


@dataclass(frozen=True, slots=True)
class FluxSectorStabilityResult:
    """Exact and VMC Gutzwiller projection record across static Z_2 flux / vison sectors."""

    nx: int
    ny: int
    n_sites: int
    h: float
    sector_name: str
    bc_twist: tuple[int, int]
    n_flipped_bonds: int
    parity_repaired: bool
    n_excited_qp: int
    ref_condition_number: float
    e_proj_exact: float | None
    exact_energy_variance: float | None
    plaquette_flux_exact: float | None
    magnetization_111_exact: float | None
    vmc_e_mean: float
    vmc_e_ci_low: float
    vmc_e_ci_high: float
    vmc_e_loc_std: float
    diagnostics: AutocorrelationDiagnostics


def compute_chain_diagnostics(
    samples_2d: NDArray[np.float64],
    *,
    acc_rate_1flip: float = 0.0,
    acc_rate_2flip: float = 0.0,
) -> AutocorrelationDiagnostics:
    """Compute Sokal/Wolff/Geyer-IPS autocorrelation, split-Rhat, and MCSE diagnostics."""
    arr = np.asarray(samples_2d, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[0] < 1 or arr.shape[1] < 2:
        raise ValueError(f"Expected 2D array (n_chains >= 1, n_samples >= 2), got {arr.shape}")
    if not np.all(np.isfinite(arr)):
        raise ValueError("Chain sample array must contain finite values.")

    n_chains, n_samples = arr.shape
    n_tot = float(n_chains * n_samples)
    grand_mean = float(np.mean(arr))
    overall_std = float(np.std(arr))

    if overall_std <= 1e-13 * max(1.0, abs(grand_mean)):
        return AutocorrelationDiagnostics(
            g_sokal=1.0,
            tau_int_wolff=0.5,
            g_ips_geyer=1.0,
            n_eff_sokal=n_tot,
            n_eff_ips=n_tot,
            split_rhat=1.0,
            mcse_chains=0.0,
            mcse_autocorr=0.0,
            acc_rate_1flip=float(acc_rate_1flip),
            acc_rate_2flip=float(acc_rate_2flip),
        )

    g_sokal = float(np.mean([_integrated_autocorr(arr[c]) for c in range(n_chains)]))
    tau_wolff = 0.5 * g_sokal

    # Grand-mean multi-chain 1/L PSD autocovariance + Geyer Initial Positive Sequence (IPS)
    centered = arr - grand_mean
    gamma0 = float(np.mean(centered * centered))
    max_lag = min(n_samples - 1, max(2, n_samples // 2))
    rho_psd = np.empty(max_lag + 1, dtype=np.float64)
    rho_psd[0] = 1.0
    for lag in range(1, max_lag + 1):
        cov_lag = float(np.mean(np.sum(centered[:, :-lag] * centered[:, lag:], axis=1) / n_samples))
        rho_psd[lag] = cov_lag / gamma0

    gamma_sum = 0.0
    m_pair = 0
    while 2 * m_pair + 1 <= max_lag:
        pair_val = float(rho_psd[2 * m_pair] + rho_psd[2 * m_pair + 1])
        if pair_val <= 0.0:
            break
        gamma_sum += pair_val
        m_pair += 1
    g_ips = max(1.0, -1.0 + 2.0 * gamma_sum)

    half = n_samples // 2
    if half >= 2:
        split_chains = np.vstack([arr[:, :half], arr[:, n_samples - half :]])
        m_split = split_chains.shape[0]
        w_var = float(np.mean(np.var(split_chains, axis=1, ddof=1)))
        b_var = float(half * np.var(np.mean(split_chains, axis=1), ddof=1)) if m_split > 1 else 0.0
        if w_var <= 1e-28:
            split_rhat = 1.0 if b_var <= 1e-28 else float("inf")
        else:
            var_plus = ((half - 1.0) / float(half)) * w_var + b_var / float(half)
            split_rhat = float(sqrt(max(1.0, var_plus / w_var)))
    else:
        split_rhat = 1.0

    chain_means = np.mean(arr, axis=1)
    mcse_chains = (
        float(np.std(chain_means, ddof=1) / sqrt(float(n_chains))) if n_chains > 1 else 0.0
    )
    n_eff_sokal = n_tot / g_sokal
    n_eff_ips = n_tot / g_ips
    return AutocorrelationDiagnostics(
        g_sokal=g_sokal,
        tau_int_wolff=tau_wolff,
        g_ips_geyer=g_ips,
        n_eff_sokal=n_eff_sokal,
        n_eff_ips=n_eff_ips,
        split_rhat=split_rhat,
        mcse_chains=mcse_chains,
        mcse_autocorr=overall_std / sqrt(max(1.0, n_eff_ips)),
        acc_rate_1flip=float(acc_rate_1flip),
        acc_rate_2flip=float(acc_rate_2flip),
    )


def build_site_su2_gauge_orthogonal(thetas: NDArray[np.float64]) -> NDArray[np.float64]:
    """Build block-diagonal ``O = bigoplus_s exp(theta_s . SO4_L) in SO(4)^N`` in ``(-c, b)`` basis."""
    th = np.asarray(thetas, dtype=np.float64)
    if th.ndim != 2 or th.shape[1] != 3 or th.shape[0] < 1:
        raise ValueError(f"Expected thetas shape (n_sites, 3), got {th.shape}")
    if not np.all(np.isfinite(th)):
        raise ValueError("Rotation angles must be finite.")

    n_sites = th.shape[0]
    o_mat = np.zeros((4 * n_sites, 4 * n_sites), dtype=np.float64)
    eye4 = np.eye(4, dtype=np.float64)
    lx, ly, lz = SO4_L["x"], SO4_L["y"], SO4_L["z"]

    for s in range(n_sites):
        tx, ty, tz = float(th[s, 0]), float(th[s, 1]), float(th[s, 2])
        half = 0.5 * sqrt(tx * tx + ty * ty + tz * tz)
        coeff = sin(half) / half if half > 1e-8 else 1.0 - (half * half) / 6.0
        gen = tx * lx + ty * ly + tz * lz
        o_mat[4 * s : 4 * s + 4, 4 * s : 4 * s + 4] = cos(half) * eye4 + coeff * gen
    return o_mat


def verify_thouless_su2_block_diagonality(
    thetas: NDArray[np.float64],
    ref: NDArray[np.int8],
) -> ThoulessSU2BlockCheck:
    """Verify ``U_cf^dagger O U_cf = diag(U_O, U_O^*)`` with site blocks ``u_s in SU(2)`` (``det = +1``)."""
    n_sites = ref.size
    if thetas.shape != (n_sites, 3):
        raise ValueError(f"Expected thetas shape ({n_sites}, 3), got {thetas.shape}")
    n_jw = 2 * n_sites
    u_cf = _build_physical_ref_bdg_transform(ref) / sqrt(2.0)
    m_bdg = u_cf.conj().T @ build_site_su2_gauge_orthogonal(thetas) @ u_cf

    u_blk = m_bdg[:n_jw, :n_jw]
    eye2 = np.eye(2, dtype=np.complex128)
    unit_errs: list[float] = []
    det_errs: list[float] = []
    for s in range(n_sites):
        u_s = u_blk[2 * s : 2 * s + 2, 2 * s : 2 * s + 2]
        unit_errs.append(float(np.max(np.abs(u_s.conj().T @ u_s - eye2))))
        det_errs.append(float(abs(np.linalg.det(u_s) - 1.0)))

    return ThoulessSU2BlockCheck(
        n_sites=n_sites,
        max_off_block=float(np.max(np.abs(m_bdg[:n_jw, n_jw:]))),
        max_conj_err=float(np.max(np.abs(m_bdg[n_jw:, n_jw:] - u_blk.conj()))),
        max_unitarity_err=max(unit_errs),
        max_det_err=max(det_errs),
    )


def make_staggered_su2_patterns(
    n_sites: int,
    *,
    seed: int = 2026,
) -> dict[str, NDArray[np.float64]]:
    """Construct canonical uniform, sublattice-staggered, pi-flip, Z_2, and random SU(2)_R patterns."""
    rng = np.random.default_rng(seed)
    th_a = np.array([0.7, -1.1, 0.5], dtype=np.float64)
    th_b = np.array([-0.9, 0.4, 1.3], dtype=np.float64)
    axes = np.eye(3, dtype=np.float64)
    return {
        "identity": np.zeros((n_sites, 3), dtype=np.float64),
        "sublattice_staggered": np.array(
            [th_a if s % 2 == 0 else th_b for s in range(n_sites)], dtype=np.float64
        ),
        "pi_flips_xyz": np.array([np.pi * axes[s % 3] for s in range(n_sites)], dtype=np.float64),
        "z2_2pi_staggered": np.array(
            [
                np.array([2.0 * np.pi, 0.0, 0.0]) if s % 2 == 0 else np.zeros(3)
                for s in range(n_sites)
            ],
            dtype=np.float64,
        ),
        "random_inhomogeneous": rng.uniform(-np.pi, np.pi, size=(n_sites, 3)),
    }


def build_flux_sector_bonds(nx: int, ny: int, sector: FluxSectorName) -> list[tuple[int, int, str]]:
    """Return ``flipped_bonds`` coordinates ``(x, y, axis)`` generating the requested Z_2 flux sector."""
    if sector == "zero_flux":
        return []
    if sector == "two_vison":
        return [(0, 0, "z")]
    if sector == "pi_flux":
        if nx % 2 != 0:
            raise ValueError(
                f"Uniform pi-flux (all W_p = -1) is topologically forbidden on odd-nx torus ({nx}, {ny})."
            )
        return [(x, y, "z") for x in range(1, nx, 2) for y in range(ny)]
    raise ValueError(f"Unsupported flux sector {sector!r}.")


def exact_energy_variance(nx: int, ny: int, h: float, psi: NDArray[np.complex128]) -> float:
    """Evaluate exact energy-per-site variance ``Var(H / N) = <(H/N)^2> - <H/N>^2`` on normalized 2^N ``psi``."""
    n_sites = 2 * nx * ny
    dim = 1 << n_sites
    idx = np.arange(dim, dtype=np.int32)
    spins_x = [np.where((idx >> s) & 1, 1.0, -1.0) for s in range(n_sites)]
    phases_y = [np.where((idx >> s) & 1, 1.0j, -1.0j) for s in range(n_sites)]

    h_psi = np.zeros(dim, dtype=np.complex128)
    for u, v, ax in torus_bonds(nx, ny):
        if ax == "x":
            h_psi -= spins_x[u] * spins_x[v] * psi
        elif ax == "z":
            h_psi -= psi[idx ^ ((1 << u) | (1 << v))]
        else:
            h_psi -= phases_y[u] * phases_y[v] * psi[idx ^ ((1 << u) | (1 << v))]

    if h != 0.0:
        comp = h / sqrt(3.0)
        for s in range(n_sites):
            flipped = psi[idx ^ (1 << s)]
            h_psi -= comp * (spins_x[s] * psi + flipped + phases_y[s] * flipped)

    e_mean = float(np.vdot(psi, h_psi).real) / float(n_sites)
    e2_mean = float(np.vdot(h_psi, h_psi).real) / float(n_sites * n_sites)
    return max(0.0, e2_mean - e_mean * e_mean)


def _probe_bitstring_invariance(
    nx: int,
    ny: int,
    h: float,
    ref_base: NDArray[np.int8],
    z_base: NDArray[np.complex128],
    ref_rot: NDArray[np.int8],
    z_rot: NDArray[np.complex128],
    *,
    n_probes: int = 24,
    seed: int = 777,
) -> tuple[float, float]:
    """Check bitstring-by-bitstring normalized amplitude and local-energy invariance for any N (incl. N=18)."""
    n_sites = 2 * nx * ny
    comp = h / sqrt(3.0)
    bonds = torus_bonds(nx, ny)
    ev0 = _ThoulessStateEvaluator(ref_base, z_base)
    ev1 = _ThoulessStateEvaluator(ref_rot, z_rot)
    rng = np.random.default_rng(seed)

    def _local_e(ev: _ThoulessStateEvaluator, tau: NDArray[np.int8]) -> tuple[complex, float]:
        pf = ev.amplitude(tau)
        sx = np.where(tau == 1, 1.0, -1.0)
        ph_y = np.where(tau == 1, 1.0j, -1.0j)
        m_tot = 0.0
        for s in range(n_sites):
            t1 = tau.copy()
            t1[s] = 1 - t1[s]
            r1 = ev.amplitude(t1) / pf
            m_tot += float(sx[s] + r1.real + (ph_y[s] * r1).real)
        e_k = 0.0
        for u, v, ax in bonds:
            if ax == "x":
                e_k -= float(sx[u] * sx[v])
            else:
                t2 = tau.copy()
                t2[u] = 1 - t2[u]
                t2[v] = 1 - t2[v]
                r2 = ev.amplitude(t2) / pf
                e_k -= float(r2.real if ax == "z" else (ph_y[u] * ph_y[v] * r2).real)
        return pf, (e_k - comp * m_tot) / float(n_sites)

    ratios: list[complex] = []
    e_diffs: list[float] = []
    for _ in range(n_probes * 4):
        if len(ratios) >= n_probes:
            break
        tau = rng.integers(0, 2, size=n_sites, dtype=np.int8)
        pf0 = ev0.amplitude(tau)
        if abs(pf0) < 1e-6:
            continue
        pf0_v, e0 = _local_e(ev0, tau)
        pf1_v, e1 = _local_e(ev1, tau)
        ratios.append(pf1_v / pf0_v)
        e_diffs.append(abs(e1 - e0))

    if not ratios:
        return 0.0, 0.0
    r_arr = np.asarray(ratios, dtype=np.complex128)
    rel_amp = float(np.max(np.abs(r_arr - r_arr[0])) / max(1e-12, float(abs(r_arr[0]))))
    return rel_amp, float(max(e_diffs))


def evaluate_su2_gauge_stability(
    nx: int,
    ny: int,
    h: float,
    state: MeanFieldState,
    *,
    patterns: dict[str, NDArray[np.float64]] | None = None,
    n_chains: int = 50,
    burn_in_sweeps: int = 15,
    samples_per_chain: int = 50,
    seed_base: int = 3100,
    preferred_bc: tuple[int, int] | None = None,
    independent_vmc_seeds: bool = False,
) -> list[SU2GaugeStabilityResult]:
    """Audit exact & VMC Gutzwiller stability under staggered SU(2)_R gauge conjugations ``A -> O A O^T``."""
    n_sites = 2 * nx * ny
    comp = h / sqrt(3.0)
    params = ModelParameters(K=-1.0, h=(comp, comp, comp))
    bc, ref_base, z_base, cond_base = select_admissible_twist(
        nx, ny, h, state, preferred_bc=preferred_bc
    )
    a_base = assemble_torus_majorana_matrix(nx, ny, params, state, bc_twist=bc)
    psi_base = exact_wavefunction_from_thouless(ref_base, z_base) if n_sites <= 12 else None
    obs_base = _exact_spin_observables(nx, ny, h, psi_base) if psi_base is not None else None

    gauge_map = patterns if patterns is not None else make_staggered_su2_patterns(n_sites)
    base_vmc_e: float | None = None
    results: list[SU2GaugeStabilityResult] = []

    for p_idx, (name, thetas) in enumerate(gauge_map.items()):
        o_mat = build_site_su2_gauge_orthogonal(thetas)
        a_rot_raw = o_mat @ a_base @ o_mat.T
        a_rot = np.asarray(0.5 * (a_rot_raw - a_rot_raw.T), dtype=np.float64)
        rot_norm = float(np.linalg.norm(a_rot - a_base))

        tdata = build_thouless_data(a_rot, n_sites)
        if tdata is None:
            raise RuntimeError(
                f"Thouless construction failed under SU(2)_R gauge pattern {name!r}."
            )
        ref_rot, z_rot, cond_rot = tdata
        amp_diff, eloc_diff = _probe_bitstring_invariance(
            nx, ny, h, ref_base, z_base, ref_rot, z_rot, n_probes=16, seed=777 + p_idx
        )

        exact_fid: float | None = None
        d_e_ex: float | None = None
        d_m_ex: float | None = None
        d_wp_ex: float | None = None
        ex_var: float | None = None
        if psi_base is not None and obs_base is not None:
            psi_rot = exact_wavefunction_from_thouless(ref_rot, z_rot)
            if psi_rot is None:
                raise RuntimeError(f"Projected state vanished under SU(2)_R pattern {name!r}.")
            exact_fid = float(abs(np.vdot(psi_base, psi_rot)))
            obs_rot = _exact_spin_observables(nx, ny, h, psi_rot)
            d_e_ex = abs(obs_rot[0] - obs_base[0])
            d_m_ex = abs(obs_rot[1] - obs_base[1])
            d_wp_ex = abs(obs_rot[2] - obs_base[2])
            ex_var = exact_energy_variance(nx, ny, h, psi_rot)

        run_seed = seed_base + (1000 * p_idx if independent_vmc_seeds else 0)
        batch = _run_vmc_from_thouless(
            nx,
            ny,
            h,
            ref_rot,
            z_rot,
            n_chains=n_chains,
            burn_in_sweeps=burn_in_sweeps,
            samples_per_chain=samples_per_chain,
            seed_base=run_seed,
        )
        e_mean = float(np.mean(batch.chain_e_means))
        if base_vmc_e is None:
            base_vmc_e = e_mean
        e_lo, e_hi = bca_confidence_interval(batch.chain_e_means, seed=run_seed)
        m_lo, m_hi = bca_confidence_interval(batch.chain_m_means, seed=run_seed + 1)
        diag = compute_chain_diagnostics(
            batch.all_e_samples,
            acc_rate_1flip=batch.acc_rate_1flip,
            acc_rate_2flip=batch.acc_rate_2flip,
        )

        results.append(
            SU2GaugeStabilityResult(
                nx=nx,
                ny=ny,
                n_sites=n_sites,
                h=float(h),
                bc_twist=bc,
                pattern_name=name,
                matrix_rotation_norm=rot_norm,
                z_eff_diff_norm=float(np.linalg.norm(z_rot - z_base)),
                cond_base=cond_base,
                cond_rot=cond_rot,
                cond_rel_diff=abs(cond_rot - cond_base) / max(1.0, cond_base),
                max_probe_amp_diff=amp_diff,
                max_probe_eloc_diff=eloc_diff,
                exact_fidelity=exact_fid,
                delta_e_exact=d_e_ex,
                delta_m111_exact=d_m_ex,
                delta_wp_exact=d_wp_ex,
                exact_energy_variance=ex_var,
                vmc_e_mean=e_mean,
                vmc_e_ci_low=e_lo,
                vmc_e_ci_high=e_hi,
                vmc_e_loc_std=float(np.std(batch.all_e_samples)),
                vmc_delta_e_vs_base=abs(e_mean - base_vmc_e),
                vmc_m111_mean=float(np.mean(batch.chain_m_means)),
                vmc_m111_ci_low=m_lo,
                vmc_m111_ci_high=m_hi,
                diagnostics=diag,
            )
        )
    return results


def evaluate_n_decoupling_su2_breaking(
    nx: int,
    ny: int,
    h: float,
    state: MeanFieldState,
    *,
    bc_twist: tuple[int, int] = (-1, -1),
) -> NDecouplingBreakingResult:
    """Contrast direct ``O A O^T`` & pure-spin ``R^a`` covariance against ``N^a`` breaking ``SU(2)_R -> Z_2``."""
    n_sites = 2 * nx * ny
    comp = h / sqrt(3.0)
    params = ModelParameters(K=-1.0, h=(comp, comp, comp))
    a_base = assemble_torus_majorana_matrix(nx, ny, params, state, bc_twist=bc_twist)
    td_base = build_thouless_data(a_base, n_sites)
    if td_base is None:
        raise RuntimeError("Base Thouless state inadmissible for N-decoupling contrast.")
    psi_base = exact_wavefunction_from_thouless(td_base[0], td_base[1])
    assert psi_base is not None

    th_a = np.array([0.7, -1.1, 0.5], dtype=np.float64)
    th_b = np.array([-0.9, 0.4, 1.3], dtype=np.float64)
    o_a = build_site_su2_gauge_orthogonal(th_a[None, :])
    o_b = build_site_su2_gauge_orthogonal(th_b[None, :])
    o_full = build_site_su2_gauge_orthogonal(
        np.array([th_a if s % 2 == 0 else th_b for s in range(n_sites)], dtype=np.float64)
    )

    def _fidelity_from_a(a_cand: NDArray[np.float64]) -> float:
        td = build_thouless_data(a_cand, n_sites)
        if td is None:
            raise RuntimeError("Candidate Thouless state inadmissible.")
        psi = exact_wavefunction_from_thouless(td[0], td[1])
        assert psi is not None
        return float(abs(np.vdot(psi_base, psi)))

    a_conj_raw = o_full @ a_base @ o_full.T
    fid_conj = _fidelity_from_a(np.asarray(0.5 * (a_conj_raw - a_conj_raw.T), dtype=np.float64))

    rot_state = MeanFieldState(
        U=np.array([o_a @ state.U[b] @ o_b.T for b in range(state.U.shape[0])], dtype=np.float64),
        V=np.array([o_a @ state.V[0] @ o_a.T, o_b @ state.V[1] @ o_b.T], dtype=np.float64),
        lambdas=state.lambdas.copy(),
    )
    fid_n_stag = _fidelity_from_a(
        assemble_torus_majorana_matrix(nx, ny, params, rot_state, bc_twist=bc_twist)
    )
    z2_state = MeanFieldState(U=-state.U.copy(), V=state.V.copy(), lambdas=state.lambdas.copy())
    fid_n_z2 = _fidelity_from_a(
        assemble_torus_majorana_matrix(nx, ny, params, z2_state, bc_twist=bc_twist)
    )

    # Pure-spin SO4_R hopping decoupling: 8 * K_b * R^g U R^g (since [L^a, R^g] = 0!)
    cell = make_honeycomb_2site()
    b_base = build_on_site_matrices(cell, params, state, torch.device("cpu")).numpy()
    k_bonds = [float(params.exchange_matrix(bt)[bt - 1, bt - 1]) for bt in (3, 1, 2)]

    def _assemble_pure_spin_r(
        u_arr: NDArray[np.float64], b_arr: NDArray[np.float64]
    ) -> NDArray[np.float64]:
        mat = np.zeros((4 * n_sites, 4 * n_sites), dtype=np.float64)
        a_r = [
            8.0 * k_bonds[b_i] * (SO4_R[ax] @ u_arr[b_i] @ SO4_R[ax])
            for b_i, ax in enumerate(("z", "x", "y"))
        ]
        for x in range(nx):
            for y in range(ny):
                for sub in (0, 1):
                    s = torus_site(nx, ny, x, y, sub)
                    mat[4 * s : 4 * s + 4, 4 * s : 4 * s + 4] = b_arr[sub]
                sa = torus_site(nx, ny, x, y, 0)
                offsets = (
                    (0, 0, 1.0),
                    (1, 0, float(bc_twist[0]) if x == nx - 1 else 1.0),
                    (0, 1, float(bc_twist[1]) if y == ny - 1 else 1.0),
                )
                for b_i, (dx, dy, sgn) in enumerate(offsets):
                    sb = torus_site(nx, ny, x + dx, y + dy, 1)
                    mat[4 * sa : 4 * sa + 4, 4 * sb : 4 * sb + 4] += sgn * a_r[b_i]
                    mat[4 * sb : 4 * sb + 4, 4 * sa : 4 * sa + 4] -= sgn * a_r[b_i].T
        return np.asarray(0.5 * (mat - mat.T), dtype=np.float64)

    a_r_base = _assemble_pure_spin_r(state.U, b_base)
    b_rot_exact = np.array([o_a @ b_base[0] @ o_a.T, o_b @ b_base[1] @ o_b.T], dtype=np.float64)
    a_r_rot = _assemble_pure_spin_r(rot_state.U, b_rot_exact)
    td_r0 = build_thouless_data(a_r_base, n_sites)
    td_r1 = build_thouless_data(a_r_rot, n_sites)
    assert td_r0 is not None and td_r1 is not None
    psi_r0 = exact_wavefunction_from_thouless(td_r0[0], td_r0[1])
    psi_r1 = exact_wavefunction_from_thouless(td_r1[0], td_r1[1])
    assert psi_r0 is not None and psi_r1 is not None

    return NDecouplingBreakingResult(
        nx=nx,
        ny=ny,
        n_sites=n_sites,
        h=float(h),
        direct_conjugation_fidelity=fid_conj,
        pure_spin_r_decoupling_fidelity=float(abs(np.vdot(psi_r0, psi_r1))),
        n_decoupling_staggered_fidelity=fid_n_stag,
        n_decoupling_z2_fidelity=fid_n_z2,
    )


def evaluate_flux_sector_gutzwiller(
    nx: int,
    ny: int,
    h: float,
    state: MeanFieldState,
    sector: FluxSectorName,
    *,
    n_chains: int = 50,
    burn_in_sweeps: int = 12,
    samples_per_chain: int = 40,
    seed_base: int = 4200,
    twists: tuple[tuple[int, int], ...] = ((-1, -1), (-1, 1), (1, -1), (1, 1)),
) -> FluxSectorStabilityResult:
    """Evaluate exact & VMC Gutzwiller projection in a static Z_2 flux sector with 1-QP parity repair."""
    n_sites = 2 * nx * ny
    comp = h / sqrt(3.0)
    params = ModelParameters(K=-1.0, h=(comp, comp, comp))
    flipped = build_flux_sector_bonds(nx, ny, sector)

    best_exact: (
        tuple[
            float,
            tuple[int, int],
            NDArray[np.int8],
            NDArray[np.complex128],
            float,
            bool,
            float | None,
            float | None,
            float | None,
        ]
        | None
    ) = None

    for excite_qp in (False, True):
        for bc in twists:
            a_mat = assemble_torus_majorana_matrix(
                nx, ny, params, state, bc_twist=bc, flipped_bonds=flipped
            )
            tdata = build_thouless_data(a_mat, n_sites, excite_lowest_qp=excite_qp)
            if tdata is None:
                continue
            ref, z_eff, cond = tdata
            if n_sites <= 12:
                psi = exact_wavefunction_from_thouless(ref, z_eff)
                if psi is None:
                    continue
                e_site, m111, wp, _ = _exact_spin_observables(nx, ny, h, psi)
                e_var = exact_energy_variance(nx, ny, h, psi)
                if best_exact is None or e_site < best_exact[0]:
                    best_exact = (e_site, bc, ref, z_eff, cond, excite_qp, e_var, wp, m111)
            else:
                if best_exact is None or cond < best_exact[4]:
                    best_exact = (0.0, bc, ref, z_eff, cond, excite_qp, None, None, None)
        if best_exact is not None:
            break

    if best_exact is None:
        raise RuntimeError(
            f"No admissible 0-QP or 1-QP Gutzwiller twist found on ({nx},{ny}) sector={sector!r}."
        )

    e_ex_val, bc_sel, ref_sel, z_sel, cond_sel, repaired, e_var_ex, wp_ex, m_ex = best_exact
    batch = _run_vmc_from_thouless(
        nx,
        ny,
        h,
        ref_sel,
        z_sel,
        n_chains=n_chains,
        burn_in_sweeps=burn_in_sweeps,
        samples_per_chain=samples_per_chain,
        seed_base=seed_base,
    )
    e_lo, e_hi = bca_confidence_interval(batch.chain_e_means, seed=seed_base)
    diag = compute_chain_diagnostics(
        batch.all_e_samples,
        acc_rate_1flip=batch.acc_rate_1flip,
        acc_rate_2flip=batch.acc_rate_2flip,
    )

    return FluxSectorStabilityResult(
        nx=nx,
        ny=ny,
        n_sites=n_sites,
        h=float(h),
        sector_name=sector,
        bc_twist=bc_sel,
        n_flipped_bonds=len(flipped),
        parity_repaired=repaired,
        n_excited_qp=1 if repaired else 0,
        ref_condition_number=cond_sel,
        e_proj_exact=e_ex_val if n_sites <= 12 else None,
        exact_energy_variance=e_var_ex,
        plaquette_flux_exact=wp_ex,
        magnetization_111_exact=m_ex,
        vmc_e_mean=float(np.mean(batch.chain_e_means)),
        vmc_e_ci_low=e_lo,
        vmc_e_ci_high=e_hi,
        vmc_e_loc_std=float(np.std(batch.all_e_samples)),
        diagnostics=diag,
    )


__all__ = [
    "AutocorrelationDiagnostics",
    "FluxSectorName",
    "FluxSectorStabilityResult",
    "NDecouplingBreakingResult",
    "SU2GaugeStabilityResult",
    "ThoulessSU2BlockCheck",
    "build_flux_sector_bonds",
    "build_site_su2_gauge_orthogonal",
    "compute_chain_diagnostics",
    "evaluate_flux_sector_gutzwiller",
    "evaluate_n_decoupling_su2_breaking",
    "evaluate_su2_gauge_stability",
    "exact_energy_variance",
    "make_staggered_su2_patterns",
    "verify_thouless_su2_block_diagonality",
]
