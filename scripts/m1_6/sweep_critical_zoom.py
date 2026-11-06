"""High-resolution critical region zoom sweep (M1.6 Phase 3).

Maps the first-order thermodynamic crossing h_c^* and forward spinodal h_c2^+
across h in [0.280, 0.520] with fine resolution Delta h = 0.002.
Tracks dual forward (QSL) and backward (polarized) adiabatic continuations,
energy gap Delta E(h), Wilson loop derivative dW_p/dh, and contraction factor mu_k.
"""

from __future__ import annotations

import time
from typing import Any

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.observables import (
    compute_plaquette_fluxes,
    compute_sublattice_magnetizations,
)
from kitaev_mf.services.seeds import make_polarized_seed, make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field
from scripts.m1_6._common import (
    append_jsonl_record,
    write_dat_benchmark,
)


def run_critical_zoom_sweep() -> None:
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()

    h_min = 0.280
    h_max = 0.520
    delta_h = 0.002
    h_grid = np.arange(h_min, h_max + 1e-6, delta_h)

    config = SolverConfig(
        k_grid_shape=(24, 24),
        max_iterations=60,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        method="diis",
        mixing_alpha=0.6,
        skew_mode="retract_step",
    )

    headers = [
        "Description: High-resolution zoom sweep h in [0.280, 0.520] with Delta h = 0.002 (M1.6 Phase 3)",
        "Columns: h  branch  conv  n_phi  wall_s  E_site  m_111  W_p  norm_lam  mu_final  eps_comm  eps_gauge",
    ]
    data_lines = []

    print(
        f"[Critical Zoom] Starting fine sweep across {len(h_grid)} field points (Delta h = {delta_h:.3f})..."
    )

    # 1. Forward continuation: start from converged QSL state at h_min
    print("[Critical Zoom] Initializing forward QSL branch at h =", h_min)
    h_vec_init = (h_min / np.sqrt(3.0), h_min / np.sqrt(3.0), h_min / np.sqrt(3.0))
    params_init = ModelParameters(K=-1.0, h=h_vec_init)
    qsl_seed = make_qsl_seed(cell)
    curr_qsl, m_init = solve_self_consistent_mean_field(
        cell, params_init, qsl_seed, config, backend
    )
    assert m_init.converged, "Failed to initialize forward QSL branch"

    fwd_records: list[dict[str, Any]] = []

    for h_val in h_grid:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)

        t0 = time.perf_counter()
        state_final, metrics = solve_self_consistent_mean_field(
            cell, params, curr_qsl, config, backend
        )
        wall_s = time.perf_counter() - t0

        eb = compute_gaussian_energy(state_final, cell, params)
        sig_mag, _ = compute_sublattice_magnetizations(state_final, cell)
        m_vec = 0.5 * (sig_mag[0] + sig_mag[1])
        m_111 = float(np.dot(m_vec, np.ones(3) / np.sqrt(3.0)))
        w_p = float(compute_plaquette_fluxes(state_final, cell)[0])
        norm_lam = float(np.linalg.norm(state_final.lambdas))
        eps_comm = metrics.comm_residuals[-1] if metrics.comm_residuals else float("nan")

        mu_final = metrics.contraction_factors[-1] if metrics.contraction_factors else float("nan")

        row = (
            f"{h_val:6.4f}  fwd_qsl    {str(metrics.converged):5s}  {metrics.map_evaluations:5d}  {wall_s:8.4f}  "
            f"{eb.per_site:12.6f}  {m_111:9.5f}  {w_p:8.4f}  {norm_lam:10.6e}  {mu_final:9.4f}  "
            f"{eps_comm:12.4e}  {metrics.final_gauge_residual:12.4e}"
        )
        data_lines.append(row)

        rec = {
            "h": float(h_val),
            "branch": "fwd_qsl",
            "conv": metrics.converged,
            "n_phi": metrics.map_evaluations,
            "wall_s": wall_s,
            "e_site": eb.per_site,
            "m_111": m_111,
            "w_p": w_p,
            "norm_lam": norm_lam,
            "mu_final": mu_final,
            "eps_comm": eps_comm,
            "eps_gauge": metrics.final_gauge_residual,
        }
        fwd_records.append(rec)
        append_jsonl_record("runs_critical_zoom.jsonl", rec)

        if metrics.converged:
            curr_qsl = state_final

    # 2. Backward continuation: start from converged polarized state at h_max
    print("[Critical Zoom] Initializing backward Polarized branch at h =", h_max)
    h_vec_max = (h_max / np.sqrt(3.0), h_max / np.sqrt(3.0), h_max / np.sqrt(3.0))
    params_max = ModelParameters(K=-1.0, h=h_vec_max)
    pol_seed = make_polarized_seed(cell, axis="111")
    curr_pol, m_pol = solve_self_consistent_mean_field(cell, params_max, pol_seed, config, backend)
    assert m_pol.converged, "Failed to initialize backward polarized branch"

    bwd_records: list[dict[str, Any]] = []

    for h_val in reversed(h_grid):
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)

        t0 = time.perf_counter()
        state_final, metrics = solve_self_consistent_mean_field(
            cell, params, curr_pol, config, backend
        )
        wall_s = time.perf_counter() - t0

        eb = compute_gaussian_energy(state_final, cell, params)
        sig_mag, _ = compute_sublattice_magnetizations(state_final, cell)
        m_vec = 0.5 * (sig_mag[0] + sig_mag[1])
        m_111 = float(np.dot(m_vec, np.ones(3) / np.sqrt(3.0)))
        w_p = float(compute_plaquette_fluxes(state_final, cell)[0])
        norm_lam = float(np.linalg.norm(state_final.lambdas))
        eps_comm = metrics.comm_residuals[-1] if metrics.comm_residuals else float("nan")
        mu_final = metrics.contraction_factors[-1] if metrics.contraction_factors else float("nan")

        row = (
            f"{h_val:6.4f}  bwd_pol    {str(metrics.converged):5s}  {metrics.map_evaluations:5d}  {wall_s:8.4f}  "
            f"{eb.per_site:12.6f}  {m_111:9.5f}  {w_p:8.4f}  {norm_lam:10.6e}  {mu_final:9.4f}  "
            f"{eps_comm:12.4e}  {metrics.final_gauge_residual:12.4e}"
        )
        data_lines.append(row)

        rec = {
            "h": float(h_val),
            "branch": "bwd_pol",
            "conv": metrics.converged,
            "n_phi": metrics.map_evaluations,
            "wall_s": wall_s,
            "e_site": eb.per_site,
            "m_111": m_111,
            "w_p": w_p,
            "norm_lam": norm_lam,
            "mu_final": mu_final,
            "eps_comm": eps_comm,
            "eps_gauge": metrics.final_gauge_residual,
        }
        bwd_records.append(rec)
        append_jsonl_record("runs_critical_zoom.jsonl", rec)

        if metrics.converged:
            curr_pol = state_final

    write_dat_benchmark("critical_zoom_111.dat", headers, data_lines)
    print("[Critical Zoom] Completed fine zoom sweep and exported critical_zoom_111.dat")

    # Locate thermodynamic crossing h_c* where Delta E(h) = E_QSL(h) - E_pol(h) = 0
    fwd_dict = {r["h"]: r for r in fwd_records}
    bwd_dict = {r["h"]: r for r in bwd_records}
    common_h = sorted(set(fwd_dict.keys()) & set(bwd_dict.keys()))

    crossings = []
    for i in range(len(common_h) - 1):
        h1, h2 = common_h[i], common_h[i + 1]
        de1 = fwd_dict[h1]["e_site"] - bwd_dict[h1]["e_site"]
        de2 = fwd_dict[h2]["e_site"] - bwd_dict[h2]["e_site"]
        if de1 * de2 <= 0.0 and abs(de1 - de2) > 1e-12:
            # Linear interpolation for zero crossing
            h_cross = h1 - de1 * (h2 - h1) / (de2 - de1)
            crossings.append(h_cross)

    if crossings:
        print(f"[Critical Zoom] >> Precision First-Order Crossing: h_c* = {crossings[0]:.4f} <<")
    else:
        print("[Critical Zoom] >> No crossing detected in the tested range <<")


if __name__ == "__main__":
    run_critical_zoom_sweep()
