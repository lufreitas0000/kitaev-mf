#!/usr/bin/env python3
"""Weak-field exploratory solver and seed multistability runner for h in {0.01, 0.05}.

Usage:
    python3 scripts/m1_6/run_weak_field_batch.py
"""

from __future__ import annotations

import time

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
from kitaev_mf.services.seeds import (
    add_controlled_perturbation,
    apply_gauge_transformation_sublattice_b,
    make_polarized_seed,
    make_qsl_seed,
)
from kitaev_mf.services.solver import solve_self_consistent_mean_field
from scripts.m1_6._common import (
    append_jsonl_record,
    setup_cpu_environment,
    write_dat_benchmark,
)


def run_weak_field_batch() -> None:
    """Execute exploratory batch comparing solvers and seeds under h || [111]."""
    setup_cpu_environment(1)
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()

    fields = [0.01, 0.05]
    solvers = ["fp", "picard", "diis", "broyden", "jfnk"]
    deltas = [0.0, 0.001, 0.01, 0.05]

    rng = np.random.default_rng(123)

    headers = [
        "Description: Weak-field solver convergence and seed multistability benchmark (M1.6)",
        "Columns: h  solver  seed_name  delta  converged  n_phi  wall_s  E_site  m_x  m_y  m_z  m_111  W_p  norm_lam  eps_comm  eps_gauge",
    ]
    data_lines = []

    print("[Weak Field Batch] Running solver x seed exploratory sweep...")

    for h_val in fields:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)

        for seed_base_name in ["qsl", "polarized", "gauge_flip_b"]:
            if seed_base_name == "qsl":
                base_state = make_qsl_seed(cell)
            elif seed_base_name == "polarized":
                base_state = make_polarized_seed(cell, axis="111")
            elif seed_base_name == "gauge_flip_b":
                base_state = apply_gauge_transformation_sublattice_b(make_qsl_seed(cell), cell)
            else:
                continue

            for delta in deltas:
                if delta == 0.0:
                    init_state = base_state
                else:
                    init_state = add_controlled_perturbation(
                        base_state, delta=delta, rng=rng, skew_only=True
                    )

                for solver in solvers:
                    config = SolverConfig(
                        k_grid_shape=(16, 16),
                        max_iterations=60,
                        tol_correlator=1e-8,
                        tol_gauge=1e-8,
                        method=solver,  # type: ignore[arg-type]
                        mixing_alpha=0.5,
                        broyden_alpha=0.5,
                        skew_mode="retract_step",
                    )

                    t0 = time.perf_counter()
                    try:
                        state_final, metrics = solve_self_consistent_mean_field(
                            cell, params, init_state, config, backend
                        )
                        wall_s = time.perf_counter() - t0
                        conv = metrics.converged
                        n_phi = metrics.map_evaluations

                        # Evaluate physical observables
                        eb = compute_gaussian_energy(state_final, cell, params)
                        e_site = eb.per_site
                        sig_mag, _ = compute_sublattice_magnetizations(state_final, cell)
                        # Average over sublattices
                        m_vec = 0.5 * (sig_mag[0] + sig_mag[1])
                        m_111 = float(np.dot(m_vec, np.ones(3) / np.sqrt(3.0)))
                        fluxes = compute_plaquette_fluxes(state_final, cell)
                        w_p = float(fluxes[0])
                        norm_lam = float(np.linalg.norm(state_final.lambdas))
                        eps_comm = (
                            metrics.comm_residuals[-1] if metrics.comm_residuals else float("nan")
                        )
                        eps_gauge = metrics.final_gauge_residual
                    except Exception as exc:
                        wall_s = time.perf_counter() - t0
                        conv = False
                        n_phi = 0
                        e_site = float("nan")
                        m_vec = np.zeros(3)
                        m_111 = float("nan")
                        w_p = float("nan")
                        norm_lam = float("nan")
                        eps_comm = float("nan")
                        eps_gauge = float("nan")
                        print(f"Failed {h_val} {seed_base_name} {delta} {solver}: {exc}")

                    row = (
                        f"{h_val:6.3f}  {solver:8s}  {seed_base_name:12s}  {delta:6.3f}  "
                        f"{str(conv):5s}  {n_phi:6d}  {wall_s:8.4f}  {e_site:12.6f}  "
                        f"{m_vec[0]:9.5f}  {m_vec[1]:9.5f}  {m_vec[2]:9.5f}  {m_111:9.5f}  "
                        f"{w_p:8.4f}  {norm_lam:10.6e}  {eps_comm:12.4e}  {eps_gauge:12.4e}"
                    )
                    data_lines.append(row)
                    print(f"  -> {row}")

                    append_jsonl_record(
                        "runs_weak_field.jsonl",
                        {
                            "h": h_val,
                            "solver": solver,
                            "seed": seed_base_name,
                            "delta": delta,
                            "converged": conv,
                            "n_phi": n_phi,
                            "wall_s": wall_s,
                            "e_site": e_site,
                            "m_vec": m_vec.tolist(),
                            "m_111": m_111,
                            "w_p": w_p,
                            "norm_lam": norm_lam,
                            "eps_comm": eps_comm,
                            "eps_gauge": eps_gauge,
                            "residual_history": list(metrics.residual_history),
                            "contraction_factors": list(metrics.contraction_factors),
                        },
                    )

    write_dat_benchmark("weak_field_solvers.dat", headers, data_lines)
    print("[Weak Field Batch] Finished successfully.")


if __name__ == "__main__":
    run_weak_field_batch()
