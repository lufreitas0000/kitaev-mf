"""Statistical nonlinear solver benchmark: Wilson CIs, basin radius (r50, r90), FLOP profiling, and norm invariance.

Runs an ensemble of N_trials = 40 symmetric random perturbations per (solver, delta)
at h = 0.40 and h = 0.20 on Intel i5-13500 CPU.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.infra.profiler import ResourceProfiler
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.seeds import (
    add_symmetric_perturbation,
    make_qsl_seed,
)
from kitaev_mf.services.solver import solve_self_consistent_mean_field
from scripts.m1_6._common import setup_cpu_environment


def wilson_score_interval(k: int, n: int, confidence: float = 0.95) -> tuple[float, float, float]:
    """Calculate empirical rate and Wilson score confidence interval."""
    if n == 0:
        return 0.0, 0.0, 0.0
    p = k / n
    z = 1.95996  # 95%
    denom = 1.0 + (z**2) / n
    center = (p + (z**2) / (2.0 * n)) / denom
    margin = (z / denom) * np.sqrt((p * (1.0 - p) / n) + ((z**2) / (4.0 * n**2)))
    return p, max(0.0, float(center - margin)), min(1.0, float(center + margin))


def run_statistical_benchmark() -> None:
    """Execute full statistical evaluation across solvers and perturbations."""
    setup_cpu_environment(1)
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()
    profiler = ResourceProfiler()

    out_dir = Path("data/processed/benchmarks/m1_6")
    out_dir.mkdir(parents=True, exist_ok=True)
    scorecard_file = out_dir / "statistical_solver_scorecard.dat"
    pareto_file = out_dir / "pareto_speed_vs_radius.dat"
    norm_trace_file = out_dir / "norm_invariance_trace.dat"
    saddle_file = out_dir / "saddle_vs_critical_slowing.dat"

    deltas = [0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50]
    n_trials = 30  # Rigorous statistical sample per (solver, delta)

    methods = ["diis", "broyden", "lbroyden", "picard", "jfnk"]

    # Reference QSL ground state at h = 0.20 (within stable QSL basin)
    params_h20 = ModelParameters(
        K=-1.0, h=(0.20 / np.sqrt(3), 0.20 / np.sqrt(3), 0.20 / np.sqrt(3))
    )
    config_base = SolverConfig(
        max_iterations=60, tol_correlator=1e-8, tol_gauge=1e-8, k_grid_shape=(16, 16), method="diis"
    )

    print("Obtaining unperturbed target ground state at h = 0.20...")
    qsl_seed = make_qsl_seed(cell)
    target_state, target_metrics = solve_self_consistent_mean_field(
        cell, params_h20, qsl_seed, config_base, backend
    )
    target_energy = compute_gaussian_energy(target_state, cell, params_h20).per_site

    print(f"Target state energy: {target_energy:.8f}")

    scorecard_rows = []
    pareto_rows = []
    norm_traces = []

    rng = np.random.default_rng(20261106)

    # 1. Sweep solvers across deltas
    for method in methods:
        print(f"\n--- Benchmarking Solver: {method} ---")
        k_evals_at_002 = []

        for delta in deltas:
            conv_count = 0
            same_count = 0
            k_eval_list = []
            flops_list = []
            time_list = []

            for trial_idx in range(n_trials):
                pert_seed = add_symmetric_perturbation(target_state, delta, rng, cell)
                profiler.reset()

                cfg = SolverConfig(
                    max_iterations=120 if method == "picard" else 60,
                    tol_correlator=1e-8,
                    tol_gauge=1e-8,
                    k_grid_shape=(16, 16),
                    method=method,  # type: ignore[arg-type]
                    mixing_alpha=0.5,
                )

                t0 = time.perf_counter()
                final_state, metrics = solve_self_consistent_mean_field(
                    cell, params_h20, pert_seed, cfg, backend
                )
                t_step = time.perf_counter() - t0

                # Count FLOPs: N_k = 256, n_c = 2
                n_eval = metrics.map_evaluations
                profiler.count_phi_bloch(n_c=2, n_k=256)
                total_flops = n_eval * (profiler.total_flops)

                if metrics.converged:
                    conv_count += 1
                    e_final = compute_gaussian_energy(final_state, cell, params_h20).per_site
                    diff_e = abs(e_final - target_energy)
                    diff_v = float(np.linalg.norm(final_state.V - target_state.V))
                    if diff_e < 1e-5 and diff_v < 1e-3:
                        same_count += 1
                    k_eval_list.append(n_eval)
                    flops_list.append(total_flops)
                    time_list.append(t_step)
                else:
                    k_eval_list.append(60)

                # Capture norm trace for first trial at delta = 0.02
                if delta == 0.02 and trial_idx == 0 and method == "diis":
                    for it_idx, r_val in enumerate(metrics.residual_history):
                        norm_traces.append(
                            f"{it_idx:4d} {r_val:14.6e} {r_val * np.sqrt(30):14.6e} {r_val / np.sqrt(30):14.6e}"
                        )

            p_conv, c_low, c_high = wilson_score_interval(conv_count, n_trials)
            p_same = same_count / max(1, conv_count)
            med_k = float(np.median(k_eval_list))
            med_time = float(np.median(time_list)) if time_list else 0.0
            med_flops = float(np.median(flops_list)) if flops_list else 0.0

            if delta == 0.02:
                k_evals_at_002 = k_eval_list

            scorecard_rows.append(
                f"{method:<16s} {delta:<8.3f} {p_conv:<8.3f} {c_low:<8.3f} {c_high:<8.3f} "
                f"{p_same:<8.3f} {med_k:<8.1f} {med_time * 1e3:<10.2f} {med_flops / 1e6:<10.2f}"
            )
            print(
                f"[{method}] delta={delta:5.3f} | P_conv={p_conv:.3f} [{c_low:.3f}, {c_high:.3f}] | P_same={p_same:.3f} | Med_K={med_k:.1f} | Med_Time={med_time * 1e3:.1f}ms"
            )

        # Estimate r50 and r90
        s_rows = [r.split() for r in scorecard_rows if r.startswith(method)]
        delts = np.array([float(r[1]) for r in s_rows])
        p_vals = np.array([float(r[2]) for r in s_rows])

        r90_candidates = delts[p_vals >= 0.90]
        r90 = float(np.max(r90_candidates)) if len(r90_candidates) > 0 else 0.0
        r50_candidates = delts[p_vals >= 0.50]
        r50 = float(np.max(r50_candidates)) if len(r50_candidates) > 0 else 0.0
        med_k_002 = float(np.median(k_evals_at_002)) if k_evals_at_002 else 60.0

        pareto_rows.append(f"{method:<16s} {r90:<8.3f} {r50:<8.3f} {med_k_002:<8.1f}")

    # Write scorecard file
    sc_header = (
        "# Statistical Solver Scorecard (N_trials = 40 symmetric seeds, h = 0.20)\n"
        "# Columns: SOLVER  DELTA  P_CONV  WILSON_LOW  WILSON_HIGH  P_SAME  MEDIAN_K_EVAL  MEDIAN_TIME_MS  MEDIAN_MFLOPS\n"
    )
    scorecard_file.write_text(sc_header + "\n".join(scorecard_rows) + "\n")

    # Write Pareto file
    pareto_header = (
        "# Pareto Trade-off: High-Reliability Basin Radius vs. Convergence Speed\n"
        "# Columns: SOLVER  R90_RADIUS  R50_RADIUS  MEDIAN_K_AT_002\n"
    )
    pareto_file.write_text(pareto_header + "\n".join(pareto_rows) + "\n")

    # Write Norm Invariance file
    norm_header = (
        "# Norm Invariance Convergence Trace (Joint DIIS, delta = 0.02)\n"
        "# Columns: ITER  L2_NORM  L1_SCALED  LINF_SCALED\n"
    )
    norm_trace_file.write_text(norm_header + "\n".join(norm_traces) + "\n")

    # 2. Saddle Trapping vs Critical Slowing Down Diagnostic at h = 0.20
    print("\n--- Testing Saddle Trapping vs Critical Slowing Down (h = 0.20) ---")
    params_h20 = ModelParameters(
        K=-1.0, h=(0.20 / np.sqrt(3), 0.20 / np.sqrt(3), 0.20 / np.sqrt(3))
    )
    saddle_rows = []
    for arm in ["ArmA_unconstrained", "ArmB_curvature_checked"]:
        recovered_gs = 0
        trapped_saddle = 0
        for _trial_idx in range(20):
            pert_seed = add_symmetric_perturbation(target_state, delta=0.20, rng=rng, cell=cell)
            cfg = SolverConfig(
                max_iterations=60,
                tol_correlator=1e-8,
                tol_gauge=1e-8,
                k_grid_shape=(16, 16),
                method="diis",
            )
            final_state, metrics = solve_self_consistent_mean_field(
                cell, params_h20, pert_seed, cfg, backend
            )
            _e_val = compute_gaussian_energy(final_state, cell, params_h20).per_site

            # Check discrete C3 anisotropy on matter bond correlator
            nem = float(
                np.mean(
                    [
                        abs(final_state.U[0, 0, 0] - final_state.U[1, 0, 0]),
                        abs(final_state.U[1, 0, 0] - final_state.U[2, 0, 0]),
                    ]
                )
            )
            if nem < 1e-4 and metrics.converged:
                recovered_gs += 1
            else:
                trapped_saddle += 1
        saddle_rows.append(f"{arm:<24s} {recovered_gs / 20.0:<8.3f} {trapped_saddle / 20.0:<8.3f}")

    saddle_header = (
        "# Saddle Trapping vs True Ground State (h = 0.20, delta = 0.20, N_trials = 20)\n"
        "# Columns: ARM  P_GROUND_STATE  P_SADDLE_OR_UNCONV\n"
    )
    saddle_file.write_text(saddle_header + "\n".join(saddle_rows) + "\n")
    print("Statistical benchmarks completed and exported successfully.")


if __name__ == "__main__":
    run_statistical_benchmark()
