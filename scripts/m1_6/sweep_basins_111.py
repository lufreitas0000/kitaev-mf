#!/usr/bin/env python3
"""Wide-range field multistability and adiabatic continuation sweep (h in [0.0, 5.0]).

Usage:
    python3 scripts/m1_6/sweep_basins_111.py
"""

from __future__ import annotations

import time

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.observables import (
    compute_plaquette_fluxes,
    compute_sublattice_magnetizations,
)
from kitaev_mf.services.seeds import (
    add_controlled_perturbation,
    make_polarized_seed,
    make_qsl_seed,
)
from kitaev_mf.services.solver import solve_self_consistent_mean_field
from scripts.m1_6._common import (
    append_jsonl_record,
    setup_cpu_environment,
    write_dat_benchmark,
)


def run_wide_range_basin_sweep() -> None:
    """Execute wide-range sweep mapping multistability and hysteresis across h in [0.0, 5.0]."""
    setup_cpu_environment(1)
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()
    rng = np.random.default_rng(2026)

    # 31 field points: dense in [0.0, 1.0], coarse up to 5.0
    h_dense = np.linspace(0.0, 1.0, 21)
    h_coarse = np.array([1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.5, 5.0])
    field_grid = np.concatenate([h_dense, h_coarse])

    headers = [
        "Description: Wide-range field multistability and hysteresis benchmark h in [0, 5] (M1.6)",
        "Columns: h  sweep_type  seed_name  delta  conv  n_phi  wall_s  E_site  m_x  m_y  m_z  m_111  W_p  norm_lam  eps_comm  eps_gauge",
    ]
    data_lines = []

    config = SolverConfig(
        k_grid_shape=(24, 24),
        max_iterations=50,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        method="diis",
        mixing_alpha=0.6,
        skew_mode="retract_step",
    )

    print("[Basin Sweep] 1. Static seed matrix sweep across all field points...")

    static_seeds_def = [
        ("qsl_exact", 0.0),
        ("qsl_pert_med", 0.05),
        ("qsl_pert_large", 0.20),
        ("pol_exact", 0.0),
        ("pol_pert_med", 0.05),
        ("pol_pert_large", 0.20),
    ]

    for h_val in field_grid:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)

        for seed_tag, delta in static_seeds_def:
            if "qsl" in seed_tag:
                base = make_qsl_seed(cell)
            else:
                base = make_polarized_seed(cell, axis="111")

            init_state = add_controlled_perturbation(base, delta=delta, rng=rng, skew_only=True)

            t0 = time.perf_counter()
            state_final, metrics = solve_self_consistent_mean_field(
                cell, params, init_state, config, backend
            )
            wall_s = time.perf_counter() - t0

            eb = compute_gaussian_energy(state_final, cell, params)
            e_site = eb.per_site
            sig_mag, _ = compute_sublattice_magnetizations(state_final, cell)
            m_vec = 0.5 * (sig_mag[0] + sig_mag[1])
            m_111 = float(np.dot(m_vec, np.ones(3) / np.sqrt(3.0)))
            w_p = float(compute_plaquette_fluxes(state_final, cell)[0])
            norm_lam = float(np.linalg.norm(state_final.lambdas))
            eps_comm = metrics.comm_residuals[-1] if metrics.comm_residuals else float("nan")

            row = (
                f"{h_val:6.3f}  static        {seed_tag:16s}  {delta:6.3f}  "
                f"{str(metrics.converged):5s}  {metrics.map_evaluations:6d}  {wall_s:8.4f}  "
                f"{e_site:12.6f}  {m_vec[0]:9.5f}  {m_vec[1]:9.5f}  {m_vec[2]:9.5f}  "
                f"{m_111:9.5f}  {w_p:8.4f}  {norm_lam:10.6e}  {eps_comm:12.4e}  {metrics.final_gauge_residual:12.4e}"
            )
            data_lines.append(row)

            append_jsonl_record(
                "runs_basin_sweep.jsonl",
                {
                    "h": h_val,
                    "sweep": "static",
                    "seed": seed_tag,
                    "delta": delta,
                    "conv": metrics.converged,
                    "n_phi": metrics.map_evaluations,
                    "wall_s": wall_s,
                    "e_site": e_site,
                    "m_vec": m_vec.tolist(),
                    "m_111": m_111,
                    "w_p": w_p,
                    "norm_lam": norm_lam,
                    "eps_comm": eps_comm,
                    "eps_gauge": metrics.final_gauge_residual,
                },
            )

    print("[Basin Sweep] 2. Adiabatic Forward Continuation (0.0 -> 5.0)...")
    curr_state_fwd: MeanFieldState = make_qsl_seed(cell)

    for h_val in field_grid:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)

        t0 = time.perf_counter()
        state_final, metrics = solve_self_consistent_mean_field(
            cell, params, curr_state_fwd, config, backend
        )
        wall_s = time.perf_counter() - t0
        curr_state_fwd = state_final

        eb = compute_gaussian_energy(state_final, cell, params)
        e_site = eb.per_site
        sig_mag, _ = compute_sublattice_magnetizations(state_final, cell)
        m_vec = 0.5 * (sig_mag[0] + sig_mag[1])
        m_111 = float(np.dot(m_vec, np.ones(3) / np.sqrt(3.0)))
        w_p = float(compute_plaquette_fluxes(state_final, cell)[0])
        norm_lam = float(np.linalg.norm(state_final.lambdas))
        eps_comm = metrics.comm_residuals[-1] if metrics.comm_residuals else float("nan")

        row = (
            f"{h_val:6.3f}  fwd_contin    {'qsl_continued':16s}  {0.0:6.3f}  "
            f"{str(metrics.converged):5s}  {metrics.map_evaluations:6d}  {wall_s:8.4f}  "
            f"{e_site:12.6f}  {m_vec[0]:9.5f}  {m_vec[1]:9.5f}  {m_vec[2]:9.5f}  "
            f"{m_111:9.5f}  {w_p:8.4f}  {norm_lam:10.6e}  {eps_comm:12.4e}  {metrics.final_gauge_residual:12.4e}"
        )
        data_lines.append(row)

        append_jsonl_record(
            "runs_basin_sweep.jsonl",
            {
                "h": h_val,
                "sweep": "fwd_contin",
                "seed": "qsl_continued",
                "delta": 0.0,
                "conv": metrics.converged,
                "n_phi": metrics.map_evaluations,
                "wall_s": wall_s,
                "e_site": e_site,
                "m_vec": m_vec.tolist(),
                "m_111": m_111,
                "w_p": w_p,
                "norm_lam": norm_lam,
                "eps_comm": eps_comm,
                "eps_gauge": metrics.final_gauge_residual,
            },
        )

    print("[Basin Sweep] 3. Adiabatic Backward Continuation (5.0 -> 0.0)...")
    curr_state_bwd: MeanFieldState = make_polarized_seed(cell, axis="111")

    for h_val in reversed(field_grid):
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)

        t0 = time.perf_counter()
        state_final, metrics = solve_self_consistent_mean_field(
            cell, params, curr_state_bwd, config, backend
        )
        wall_s = time.perf_counter() - t0
        curr_state_bwd = state_final

        eb = compute_gaussian_energy(state_final, cell, params)
        e_site = eb.per_site
        sig_mag, _ = compute_sublattice_magnetizations(state_final, cell)
        m_vec = 0.5 * (sig_mag[0] + sig_mag[1])
        m_111 = float(np.dot(m_vec, np.ones(3) / np.sqrt(3.0)))
        w_p = float(compute_plaquette_fluxes(state_final, cell)[0])
        norm_lam = float(np.linalg.norm(state_final.lambdas))
        eps_comm = metrics.comm_residuals[-1] if metrics.comm_residuals else float("nan")

        row = (
            f"{h_val:6.3f}  bwd_contin    {'pol_continued':16s}  {0.0:6.3f}  "
            f"{str(metrics.converged):5s}  {metrics.map_evaluations:6d}  {wall_s:8.4f}  "
            f"{e_site:12.6f}  {m_vec[0]:9.5f}  {m_vec[1]:9.5f}  {m_vec[2]:9.5f}  "
            f"{m_111:9.5f}  {w_p:8.4f}  {norm_lam:10.6e}  {eps_comm:12.4e}  {metrics.final_gauge_residual:12.4e}"
        )
        data_lines.append(row)

        append_jsonl_record(
            "runs_basin_sweep.jsonl",
            {
                "h": h_val,
                "sweep": "bwd_contin",
                "seed": "pol_continued",
                "delta": 0.0,
                "conv": metrics.converged,
                "n_phi": metrics.map_evaluations,
                "wall_s": wall_s,
                "e_site": e_site,
                "m_vec": m_vec.tolist(),
                "m_111": m_111,
                "w_p": w_p,
                "norm_lam": norm_lam,
                "eps_comm": eps_comm,
                "eps_gauge": metrics.final_gauge_residual,
            },
        )

    write_dat_benchmark("basin_sweep_111.dat", headers, data_lines)
    print("[Basin Sweep] Sweep completed and exported.")


if __name__ == "__main__":
    run_wide_range_basin_sweep()
