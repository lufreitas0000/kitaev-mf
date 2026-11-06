#!/usr/bin/env python3
"""A/B Benchmark runner evaluating skew-symmetry handling arms (AH-04).

Usage:
    python3 scripts/m1_6/ab_skew_retraction.py
"""

from __future__ import annotations

import time

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.infra.linalg_backend import TorchBackend
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


def run_skew_ab_benchmark() -> None:
    """Execute A/B testing matrix across skew-symmetry arms."""
    setup_cpu_environment(1)
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()

    fields = [0.01, 0.05]
    solvers = ["fp", "diis", "broyden", "jfnk"]
    arms = ["raw", "symmetrize_output", "retract_step"]

    rng = np.random.default_rng(77)

    headers = [
        "Description: Skew-symmetrization A/B benchmark (AH-04) on weak field h || [111]",
        "Columns: arm  solver  h  seed_type  converged  n_phi  wall_s  skew_err  eps_comm  eps_gauge",
    ]
    data_lines = []

    print("[Skew A/B] Starting benchmark sweep...")

    for h_val in fields:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)

        # Base seeds
        qsl_base = make_qsl_seed(cell)
        qsl_pert = add_controlled_perturbation(qsl_base, delta=0.01, rng=rng)
        pol_base = make_polarized_seed(cell, axis="111")
        pol_pert = add_controlled_perturbation(pol_base, delta=0.01, rng=rng)

        seed_dict = {
            "qsl_exact": qsl_base,
            "qsl_pert": qsl_pert,
            "pol_exact": pol_base,
            "pol_pert": pol_pert,
        }

        for seed_name, init_state in seed_dict.items():
            for solver in solvers:
                for arm in arms:
                    config = SolverConfig(
                        k_grid_shape=(16, 16),
                        max_iterations=40,
                        tol_correlator=1e-8,
                        tol_gauge=1e-8,
                        method=solver,  # type: ignore[arg-type]
                        skew_mode=arm,  # type: ignore[arg-type]
                    )

                    t0 = time.perf_counter()
                    try:
                        state_final, metrics = solve_self_consistent_mean_field(
                            cell, params, init_state, config, backend
                        )
                        wall_s = time.perf_counter() - t0
                        conv = metrics.converged
                        n_phi = metrics.map_evaluations
                        skew_err = float(
                            np.max(
                                [
                                    np.max(np.abs(state_final.V[s] + state_final.V[s].T))
                                    for s in range(2)
                                ]
                            )
                        )
                        eps_comm = (
                            metrics.comm_residuals[-1] if metrics.comm_residuals else float("nan")
                        )
                        eps_gauge = metrics.final_gauge_residual
                    except Exception as exc:
                        wall_s = time.perf_counter() - t0
                        conv = False
                        n_phi = 0
                        skew_err = float("nan")
                        eps_comm = float("nan")
                        eps_gauge = float("nan")
                        print(f"Failed {arm} {solver} {h_val} {seed_name}: {exc}")

                    row = (
                        f"{arm:18s}  {solver:8s}  {h_val:6.3f}  {seed_name:12s}  "
                        f"{str(conv):5s}  {n_phi:6d}  {wall_s:8.4f}  "
                        f"{skew_err:12.4e}  {eps_comm:12.4e}  {eps_gauge:12.4e}"
                    )
                    data_lines.append(row)
                    print(f"  -> {row}")

                    append_jsonl_record(
                        "runs_skew_ab.jsonl",
                        {
                            "arm": arm,
                            "solver": solver,
                            "h": h_val,
                            "seed": seed_name,
                            "converged": conv,
                            "n_phi": n_phi,
                            "wall_s": wall_s,
                            "skew_err": skew_err,
                            "eps_comm": eps_comm,
                            "eps_gauge": eps_gauge,
                        },
                    )

    write_dat_benchmark("ab_skew_retraction.dat", headers, data_lines)
    print("[Skew A/B] Benchmark finished successfully.")


if __name__ == "__main__":
    run_skew_ab_benchmark()
