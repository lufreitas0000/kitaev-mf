#!/usr/bin/env python3
"""Benchmark script for Task T3.W1.1: Scalable Canonical Skew-Schur Solver and SVD Overlap Engine.

Executes systematic runtime and memory profiling across cluster sizes L in {4, 6, 8, 12, 16, 20},
exporting certified ACID-provenance data tables.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

# Enforce single-threaded BLAS to avoid OpenBLAS thread contention
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import numpy as np

# Ensure src is on path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from kitaev_mf.services.phase_3.gauge_theory import (
    build_honeycomb_torus_cluster,
    build_honeycomb_cylinder_cluster,
    build_uniform_gauge_configuration,
    flip_bond_gauge_link,
    assemble_matter_skew_matrix,
)
from kitaev_mf.services.phase_3.bogoliubov import (
    extract_mode_wavefunctions,
    apply_particle_hole_interchange,
    find_optimal_particle_hole_mode_svd,
)
from kitaev_mf.services.phase_3.scalable_skew import (
    canonical_skew_schur,
    solve_scalable_quench_pair,
)


def run_benchmark(output_path: Path) -> None:
    """Run scaling benchmark across cluster sizes L."""
    print("=" * 88)
    print("TASK T3.W1.1: SCALABLE SKEW-SCHUR SOLVER & SVD OVERLAP BENCHMARK")
    print("=" * 88)
    print(f"{'L':>4} | {'N_sites':>7} | {'t_Schur (s)':>12} | {'Ortho Err':>12} | {'Recon Err':>12} | {'t_SVD (ms)':>10} | {'t_Loop (ms)':>11} | {'Speedup':>8}")
    print("-" * 88)

    records: list[dict[str, float]] = []
    L_values = [4, 6, 8, 10, 12]

    for L in L_values:
        num_sites, bonds, triplets = build_honeycomb_torus_cluster(L)
        g0 = build_uniform_gauge_configuration(num_sites, bonds)
        A0 = assemble_matter_skew_matrix(g0, K=1.0, kappa=0.1, triplets=triplets)

        # 1. Benchmark Schur
        schur_res = canonical_skew_schur(A0)
        t_schur = schur_res.elapsed_seconds
        ortho_err = schur_res.orthogonality_error
        recon_err = schur_res.reconstruction_error

        # 2. Benchmark Quench and SVD mode selection vs loop on singular cylinder boundary quench
        # Cylinder boundary quenches guarantee genuine det(X) ~ 0 singular Bogoliubov transitions
        ns_cyl, bonds_cyl, trips_cyl = build_honeycomb_cylinder_cluster(Lx=L, Ly=L)
        gc0 = build_uniform_gauge_configuration(ns_cyl, bonds_cyl)
        bc0 = bonds_cyl[1]
        gc1 = flip_bond_gauge_link(gc0, bc0.site_a, bc0.site_b)
        Ac0 = assemble_matter_skew_matrix(gc0, K=1.0, kappa=0.1, triplets=trips_cyl)
        Ac1 = assemble_matter_skew_matrix(gc1, K=1.0, kappa=0.1, triplets=trips_cyl)
        schur_c0 = canonical_skew_schur(Ac0)
        schur_c1 = canonical_skew_schur(Ac1)

        phic0 = extract_mode_wavefunctions(schur_c0.Q)
        phic1 = extract_mode_wavefunctions(schur_c1.Q)

        X = 2.0 * phic1 @ phic0.conj().T
        Y = 2.0 * phic1 @ phic0.T

        # Time SVD mode selection
        t0 = time.perf_counter()
        svd_mode, svd_det = find_optimal_particle_hole_mode_svd(X, Y)
        t_svd = time.perf_counter() - t0

        # Time loop (for L <= 8 to avoid exhaustive O(N^4) scaling bottleneck)
        if L <= 8:
            t0 = time.perf_counter()
            best_det = 0.0
            for m in range(X.shape[0]):
                X_m, _ = apply_particle_hole_interchange(X, Y, mode_idx=m)
                d = float(np.abs(np.linalg.det(X_m)))
                if d > best_det:
                    best_det = d
            t_loop = time.perf_counter() - t0
            speedup = t_loop / max(t_svd, 1e-9)
        else:
            t_loop = float("nan")
            speedup = float("nan")

        records.append({
            "L": float(L),
            "num_sites": float(num_sites),
            "t_schur": t_schur,
            "ortho_err": ortho_err,
            "recon_err": recon_err,
            "t_svd": t_svd * 1000.0,
            "t_loop": t_loop * 1000.0,
            "speedup": speedup,
        })

        t_loop_str = f"{t_loop * 1000.0:11.2f}" if not np.isnan(t_loop) else f"{'N/A':>11}"
        speedup_str = f"{speedup:7.1f}x" if not np.isnan(speedup) else f"{'N/A':>8}"
        print(f"{L:4d} | {num_sites:7d} | {t_schur:12.4f} | {ortho_err:12.2e} | {recon_err:12.2e} | {t_svd * 1000.0:10.2f} | {t_loop_str} | {speedup_str}")

    print("-" * 88)

    # Export to .dat file with ACID header
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("# ACID Provenance Header: Task T3.W1.1 Scalable Skew Solver Benchmark\n")
        f.write("# Author: Lu Freitas\n")
        f.write(f"# Timestamp: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")
        f.write("# Model: Pure FM Kitaev (K=1.0, kappa=0.1) on Honeycomb Torus\n")
        f.write("# Columns: L num_sites t_schur_s ortho_err recon_err t_svd_ms t_loop_ms speedup\n")
        for r in records:
            f.write(f"{int(r['L']):4d}  {int(r['num_sites']):7d}  {r['t_schur']:12.6f}  {r['ortho_err']:12.4e}  {r['recon_err']:12.4e}  {r['t_svd']:10.4f}  {r['t_loop']:10.4f}  {r['speedup']:8.2f}\n")

    print(f"Exported benchmark table: {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark scalable skew solver and SVD overlap engine.")
    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "data/processed/benchmarks/phase_3/skew_solver_scaling.dat",
        help="Target .dat output path.",
    )
    args = parser.parse_args()
    run_benchmark(args.output)


if __name__ == "__main__":
    main()
