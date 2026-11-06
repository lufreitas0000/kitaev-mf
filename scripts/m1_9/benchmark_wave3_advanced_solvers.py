"""Benchmark script for Wave 3: Advanced Solvers, Purification & Disordered Scaling.

Evaluates:
1. Task 3.1: Self-Consistent Matrix Sign Purification across accelerators (DIIS, L-Broyden, Picard)
   on real-space torus clusters (L in {4, 6, 8}).
2. Task 3.2: Multi-Secant L-Broyden (alpha=0.5) vs DIIS across Zeeman fields h in [0.01, 0.30].
3. Task 3.3: Disordered and defect cluster scaling under bond disorder delta in [0, 0.20] and vacancies.
"""

from __future__ import annotations

import os
import time

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


def run_wave3_benchmarks() -> None:
    output_dir = "data/processed/benchmarks/m1_9"
    os.makedirs(output_dir, exist_ok=True)
    out_file1 = os.path.join(output_dir, "wave3_purification_and_solvers.dat")
    out_file2 = os.path.join(output_dir, "wave3_disorder_scaling.dat")

    print("[Wave 3 Benchmark] Starting advanced solver & disorder benchmarks...")

    cell = make_honeycomb_2site()
    seed = make_qsl_seed(cell)
    params_base = ModelParameters(K=-1.0, h=(0.02, 0.02, 0.02))

    # =========================================================================
    # 1. Real-Space Cluster Scaling across Accelerators (Task 3.1 & 3.2)
    # =========================================================================
    print("\n--- 1. Real-Space Cluster Scaling (Purification + Accelerators) ---")
    cluster_results: list[dict[str, float | int | str]] = []

    for L in (4, 6, 8):
        cluster = RealSpaceTorusCluster(L=L)
        M = cluster.num_majoranas
        print(f"  Cluster L = {L} (N_sites = {cluster.num_sites}, M = {M}):")

        for acc in ("diis", "lbroyden", "picard"):
            t0 = time.perf_counter()
            _, _, met = cluster.solve_real_space_purification(
                params_base,
                seed,
                max_iter=30,
                tol=1e-5,
                accelerator=acc,  # type: ignore[arg-type]
                mixing_alpha=0.5,
            )
            wall_ms = (time.perf_counter() - t0) * 1000.0
            conv = met["converged"]
            it = met["iterations"]
            gemm = met["total_gemm_calls"]
            res = met["final_residual"]
            print(
                f"    {acc:8s} | Conv: {conv!s:5s} | Iter: {it:2d} | "
                f"GEMMs: {gemm:4d} | Wall: {wall_ms:7.1f} ms | Res: {res:.2e}"
            )
            cluster_results.append({
                "L": L,
                "M": M,
                "acc": acc,
                "conv": 1 if conv else 0,
                "it": it,
                "gemm": gemm,
                "wall_ms": wall_ms,
                "res": res,
            })

    # =========================================================================
    # 2. Damped L-Broyden (alpha=0.5) vs DIIS Field Robustness (Task 3.2)
    # =========================================================================
    print("\n--- 2. L-Broyden (alpha=0.5) vs DIIS Field Robustness (L = 4) ---")
    cluster4 = RealSpaceTorusCluster(L=4)
    h_values = [0.01, 0.05, 0.10, 0.20, 0.30]
    field_results: list[tuple[float, int, int, float, float]] = []

    for h_val in h_values:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params_h = ModelParameters(K=-1.0, h=h_vec)

        # DIIS
        _, _, met_diis = cluster4.solve_real_space_purification(
            params_h, seed, max_iter=35, tol=1e-5, accelerator="diis", mixing_alpha=0.5
        )
        # L-Broyden
        _, _, met_lb = cluster4.solve_real_space_purification(
            params_h, seed, max_iter=35, tol=1e-5, accelerator="lbroyden", mixing_alpha=0.5
        )

        it_diis = met_diis["iterations"]
        it_lb = met_lb["iterations"]
        res_diis = met_diis["final_residual"]
        res_lb = met_lb["final_residual"]
        print(
            f"  h = {h_val:4.2f} | DIIS: {it_diis:2d} iters (res {res_diis:.1e}) | "
            f"L-Broyden: {it_lb:2d} iters (res {res_lb:.1e})"
        )
        field_results.append((h_val, it_diis, it_lb, res_diis, res_lb))

    # =========================================================================
    # 3. Disordered and Defect Cluster Solver Scaling (Task 3.3)
    # =========================================================================
    print("\n--- 3. Disordered and Defect Cluster Scaling (L = 4) ---")
    disorder_levels = [0.00, 0.05, 0.10, 0.15, 0.20]
    rng = np.random.default_rng(20261106)
    disorder_results: list[tuple[float, int, float, int, float, float]] = []

    for d_level in disorder_levels:
        bond_disorder = {
            (m1, m2, b_idx): float(rng.uniform(-d_level, d_level))
            for m2 in range(cluster4.L)
            for m1 in range(cluster4.L)
            for b_idx in range(3)
        }
        for n_vac in (0, 2):
            vacancies = set(range(n_vac)) if n_vac > 0 else set()
            A_init = cluster4.assemble_skew_hamiltonian(
                params_base, seed, bond_disorder=bond_disorder, site_vacancies=vacancies
            )
            # Evaluate condition number kappa(A) on non-zero subspace
            sv = np.linalg.svd(A_init, compute_uv=False)
            active_sv = sv[sv > 1e-10]
            kappa = float(active_sv[0] / active_sv[-1]) if len(active_sv) > 0 else 1.0

            _, _, met_dis = cluster4.solve_real_space_purification(
                params_base,
                seed,
                max_iter=45,
                tol=1e-5,
                accelerator="diis",
                bond_disorder=bond_disorder,
                site_vacancies=vacancies,
            )
            it_dis = met_dis["iterations"]
            res_dis = met_dis["final_residual"]
            print(
                f"  Disorder delta = {d_level:4.2f}, Vacancies = {n_vac} | "
                f"kappa(A) = {kappa:6.2f} | DIIS iters: {it_dis:2d} | Res: {res_dis:.1e}"
            )
            disorder_results.append((d_level, n_vac, kappa, it_dis, res_dis, float(met_dis["converged"])))

    # =========================================================================
    # Export Benchmark Data
    # =========================================================================
    with open(out_file1, "w") as f:
        f.write("# Description: Wave 3 real-space cluster scaling across accelerators and fields\n")
        f.write("# Generated: " + time.strftime("%Y-%m-%d %H:%M:%S") + "\n")
        f.write("# Section 1: Cluster size and accelerator comparison\n")
        f.write("# L  M  accelerator  converged  iterations  gemm_calls  wall_time_ms  residual\n")
        for r in cluster_results:
            f.write(
                f"{r['L']:2d}  {r['M']:4d}  {str(r['acc']):8s}  {r['conv']:1d}  "
                f"{r['it']:2d}  {r['gemm']:5d}  {r['wall_ms']:8.2f}  {r['res']:12.4e}\n"
            )
        f.write("\n\n# Section 2: L-Broyden (alpha=0.5) vs DIIS field robustness\n")
        f.write("# h  diis_iter  lbroyden_iter  diis_res  lbroyden_res\n")
        for h_val, it_d, it_lb, r_d, r_lb in field_results:
            f.write(f"{h_val:5.2f}  {it_d:2d}  {it_lb:2d}  {r_d:12.4e}  {r_lb:12.4e}\n")

    with open(out_file2, "w") as f:
        f.write("# Description: Wave 3 disordered and defect cluster solver scaling\n")
        f.write("# Generated: " + time.strftime("%Y-%m-%d %H:%M:%S") + "\n")
        f.write("# delta_disorder  n_vacancies  kappa_A  iterations  residual  converged\n")
        for d_l, n_v, kap, it_d, r_d, c_d in disorder_results:
            f.write(f"{d_l:5.2f}  {n_v:2d}  {kap:10.4f}  {it_d:2d}  {r_d:12.4e}  {int(c_d):1d}\n")

    print(f"\n[Wave 3 Benchmark] Exported successfully to:\n  - {out_file1}\n  - {out_file2}")


if __name__ == "__main__":
    run_wave3_benchmarks()
