"""Benchmark real-space torus clusters: Track A (Dense Real vs Complex) and Track B (Eigenvector-Free Purification).

Clusters:
  L in {6, 12, 18, 24}
  N_sites = 2 * L^2 in {72, 288, 648, 1152}
  M = 8 * L^2 in {288, 1152, 2592, 4608}
"""

import time
from pathlib import Path

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.purification import purify_newton_schulz
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


def run_real_space_benchmarks() -> None:
    out_dir = Path("data/processed/benchmarks/m1_8")
    out_dir.mkdir(parents=True, exist_ok=True)

    params = ModelParameters(K=-1.0, h=(0.20, 0.20, 0.20))
    cell2 = make_honeycomb_2site()
    state = make_qsl_seed(cell2)

    sizes = [4, 6, 8, 10, 12]
    # Include L=24 if feasible, or test L=6, 12, 18 first
    print("=== Track A: Real Skew-Symmetric (dsyevd) vs. Complex Hermitian (zheevd) ===")
    track_a_rows = []

    for L in sizes:
        cluster = RealSpaceTorusCluster(L=L)
        M = cluster.num_majoranas
        n_sites = cluster.num_sites
        print(f"\n--- Cluster L = {L} (N_sites = {n_sites}, M = {M} Majoranas) ---")

        # Assemble Hamiltonian
        A = cluster.assemble_skew_hamiltonian(params, state)

        # Memory footprints:
        # Complex matrix: M x M x 16 bytes
        mem_complex_mb = (M * M * 16) / (1024 * 1024)
        # Real matrix: M x M x 8 bytes
        mem_real_mb = (M * M * 8) / (1024 * 1024)

        # Complex Hermitian eigensolve (zheevd)
        t_c0 = time.perf_counter()
        omega_c, e_c = cluster.solve_covariance_complex_eigh(A)
        time_c_ms = (time.perf_counter() - t_c0) * 1000.0

        # Real symmetric eigensolve (dsyevd on A^2)
        t_r0 = time.perf_counter()
        omega_r, e_r = cluster.solve_covariance_real_eigh(A)
        time_r_ms = (time.perf_counter() - t_r0) * 1000.0

        # Flop counts:
        # Complex: 16/3 M^3
        flops_c = (16.0 / 3.0) * (M**3)
        # Real: 4/3 M^3
        flops_r = (4.0 / 3.0) * (M**3)

        # Throughputs:
        gflops_c = (flops_c / (time_c_ms * 1e-3)) / 1e9
        gflops_r = (flops_r / (time_r_ms * 1e-3)) / 1e9

        diff_frob = float(np.linalg.norm(omega_c - omega_r)) / np.sqrt(M)
        diff_energy = abs(e_c - e_r)

        track_a_rows.append(
            (
                L,
                n_sites,
                M,
                mem_complex_mb,
                mem_real_mb,
                flops_c / 1e9,
                flops_r / 1e9,
                time_c_ms,
                time_r_ms,
                gflops_c,
                gflops_r,
                diff_frob,
                diff_energy,
            )
        )
        print(
            f"Complex zheevd: {time_c_ms:7.2f} ms ({gflops_c:5.2f} GFLOP/s) | Mem: {mem_complex_mb:5.2f} MB"
        )
        print(
            f"Real dsyevd:    {time_r_ms:7.2f} ms ({gflops_r:5.2f} GFLOP/s) | Mem: {mem_real_mb:5.2f} MB"
        )
        print(f"Parity: diff_frob = {diff_frob:.2e}, diff_energy = {diff_energy:.2e}")

    # Export Track A file
    track_a_file = out_dir / "real_space_dense_vs_complex.dat"
    with open(track_a_file, "w") as f:
        f.write("# Real-Space Track A: Real Skew (dsyevd) vs Complex (zheevd)\n")
        f.write(
            "# Columns: L  N_SITES  M  MEM_C_MB  MEM_R_MB  GFLOP_C  GFLOP_R  TIME_C_MS  TIME_R_MS  GFLOPS_C  GFLOPS_R  DIFF_FROB  DIFF_ENERGY\n"
        )
        for r in track_a_rows:
            f.write(
                f"{r[0]:<3d} {r[1]:<7d} {r[2]:<6d} {r[3]:<10.3f} {r[4]:<10.3f} {r[5]:<9.3f} {r[6]:<9.3f} {r[7]:<11.2f} {r[8]:<11.2f} {r[9]:<10.2f} {r[10]:<10.2f} {r[11]:<12.3e} {r[12]:<12.3e}\n"
            )
    print(f"\nExported {track_a_file}")

    print("\n=== Track B: Eigenvector-Free Matrix Sign Purification vs. Eigensolver ===")
    track_b_rows: list[tuple[int, int, float, float, float, int, int, float, float, float]] = []

    for L in sizes:
        cluster = RealSpaceTorusCluster(L=L)
        M = cluster.num_majoranas
        A = cluster.assemble_skew_hamiltonian(params, state)

        # Baseline eigensolver time
        t_eig0 = time.perf_counter()
        omega_eig, e_eig = cluster.solve_covariance_complex_eigh(A)
        time_eig_ms = (time.perf_counter() - t_eig0) * 1000.0

        # Track B: Newton-Schulz purification (pure real dgemm)
        omega_purif, met_purif = purify_newton_schulz(A, tol=1e-11, max_iter=30)
        time_purif_ms = met_purif.wall_time_ms
        speedup = time_eig_ms / max(time_purif_ms, 1e-4)

        # Accuracy
        diff_purif = float(np.linalg.norm(omega_purif - omega_eig)) / np.sqrt(M)
        e_diff = abs(met_purif.ground_state_energy - e_eig)

        # Total dgemm FLOPs: 2 * iter * (2 * M^3) = 4 * iter * M^3
        total_purif_gflops_ops = (4.0 * met_purif.iterations * (M**3)) / 1e9
        purif_throughput_gflops = total_purif_gflops_ops / (time_purif_ms * 1e-3)

        track_b_rows.append(
            (
                L,
                M,
                time_eig_ms,
                time_purif_ms,
                speedup,
                met_purif.iterations,
                met_purif.gemm_calls,
                purif_throughput_gflops,
                diff_purif,
                e_diff,
            )
        )
        print(
            f"[L={L:2d}, M={M:4d}] Eigensolve: {time_eig_ms:7.2f} ms | Purification: {time_purif_ms:7.2f} ms ({met_purif.iterations} iters, {purif_throughput_gflops:5.1f} GFLOP/s) | Speedup: {speedup:4.2f}x | diff: {diff_purif:.2e}"
        )

    # Export Track B file
    track_b_file = out_dir / "real_space_purification_vs_eig.dat"
    with open(track_b_file, "w") as f:
        f.write("# Real-Space Track B: Eigenvector-Free Matrix Sign Purification vs Eigensolvers\n")
        f.write(
            "# Columns: L  M  TIME_EIG_MS  TIME_PURIF_MS  SPEEDUP  ITERS  GEMM_CALLS  THROUGHPUT_GFLOPS  DIFF_NORM  DIFF_ENERGY\n"
        )
        for rb in track_b_rows:
            f.write(
                f"{rb[0]:<3d} {rb[1]:<6d} {rb[2]:<12.2f} {rb[3]:<14.2f} {rb[4]:<9.2f} {rb[5]:<7d} {rb[6]:<11d} {rb[7]:<18.2f} {rb[8]:<12.3e} {rb[9]:<12.3e}\n"
            )
    print(f"Exported {track_b_file}")


if __name__ == "__main__":
    run_real_space_benchmarks()
