"""Micro-benchmark to calibrate host CPU sustained throughput (GFLOP/s).

Measures empirical performance of real/complex GEMM and eigensolvers across matrix sizes
M in {8, 16, 32, 64, 128, 288, 1152, 2592, 4608} on Intel i5-13500 CPU.
"""

from __future__ import annotations

import time
from pathlib import Path

import torch


def calibrate_cpu_flops() -> None:
    """Run GEMM and eigensolver benchmarks and export calibration table."""
    sizes = [8, 16, 32, 64, 128, 288, 1152]

    out_dir = Path("data/processed/benchmarks/m1_6")
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "hardware_flops_calibration.dat"

    print(f"Calibrating CPU GFLOP/s on host across M = {sizes}...")

    rows = []
    for m in sizes:
        # 1. Real GEMM
        n_repeats = max(2, int(5e7 / (2.0 * m**3))) if m <= 288 else 3
        a_real = torch.randn(m, m, dtype=torch.float64, device="cpu")
        b_real = torch.randn(m, m, dtype=torch.float64, device="cpu")

        # Warmup
        _ = a_real @ b_real

        t0 = time.perf_counter()
        for _ in range(n_repeats):
            _ = a_real @ b_real
        t_gemm_real = (time.perf_counter() - t0) / n_repeats
        gflops_gemm_real = (2.0 * m**3) / (t_gemm_real * 1e9)

        # 2. Complex GEMM
        a_c = torch.randn(m, m, dtype=torch.complex128, device="cpu")
        b_c = torch.randn(m, m, dtype=torch.complex128, device="cpu")
        _ = a_c @ b_c
        t0 = time.perf_counter()
        for _ in range(n_repeats):
            _ = a_c @ b_c
        t_gemm_c = (time.perf_counter() - t0) / n_repeats
        gflops_gemm_c = (8.0 * m**3) / (t_gemm_c * 1e9)

        # 3. Real Eigh (Symmetric)
        n_rep_eig = max(2, int(1e7 / (1.33 * m**3))) if m <= 288 else 2
        sym_real = a_real + a_real.T
        _ = torch.linalg.eigh(sym_real)
        t0 = time.perf_counter()
        for _ in range(n_rep_eig):
            _ = torch.linalg.eigh(sym_real)
        t_eig_real = (time.perf_counter() - t0) / n_rep_eig
        gflops_eig_real = ((4.0 / 3.0) * m**3) / (t_eig_real * 1e9)

        # 4. Complex Eigh (Hermitian)
        herm_c = a_c + a_c.mH
        _ = torch.linalg.eigh(herm_c)
        t0 = time.perf_counter()
        for _ in range(n_rep_eig):
            _ = torch.linalg.eigh(herm_c)
        t_eig_c = (time.perf_counter() - t0) / n_rep_eig
        gflops_eig_c = ((16.0 / 3.0) * m**3) / (t_eig_c * 1e9)

        rows.append(
            f"{m:<8d} {gflops_gemm_real:<14.2f} {gflops_gemm_c:<14.2f} "
            f"{gflops_eig_real:<14.2f} {gflops_eig_c:<14.2f} "
            f"{t_gemm_real * 1e3:<12.3f} {t_eig_c * 1e3:<12.3f}"
        )
        print(
            f"M={m:4d} | GEMM real: {gflops_gemm_real:6.1f} GFLOP/s, complex: {gflops_gemm_c:6.1f} GFLOP/s | Eigh real: {gflops_eig_real:6.1f} GFLOP/s, complex: {gflops_eig_c:6.1f} GFLOP/s"
        )

    header = (
        "# Host Hardware FLOP/s Calibration (Intel i5-13500 CPU)\n"
        "# Columns: M  GFLOPS_GEMM_REAL  GFLOPS_GEMM_COMPLEX  GFLOPS_EIGH_REAL  GFLOPS_EIGH_COMPLEX  TIME_GEMM_REAL_MS  TIME_EIGH_C_MS\n"
    )
    content = header + "\n".join(rows) + "\n"
    out_file.write_text(content)
    print(f"Calibration data saved to {out_file}")


if __name__ == "__main__":
    calibrate_cpu_flops()
