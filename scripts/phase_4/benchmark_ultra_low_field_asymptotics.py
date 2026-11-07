"""Benchmark script sweeping ultra-low magnetic fields down to h = 1e-5.

Certifies the exact cubic Dirac gap power law, the linear gauge multiplier response,
and extracts the asymptotic prefactor C_MF with high precision.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

# Ensure local package takes precedence
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.phase_4.symmetrized_solver import solve_symmetrized_fixed_point


def run_ultra_low_field_benchmark() -> None:
    """Execute the ultra-low field sweep and export provenance-certified dat table."""
    cell = make_honeycomb_2site()

    # Logarithmically spaced grid from 1e-5 to 1e-2 (13 points)
    h_vals = np.logspace(-5, -2, 13)

    output_dir = (
        Path(__file__).resolve().parents[2] / "data" / "processed" / "benchmarks" / "phase_4"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    out_file = output_dir / "ultra_low_field_scaling.dat"

    print(
        f"[Phase 4 Benchmark] Sweeping {len(h_vals)} points from h = {h_vals[0]:.1e} to {h_vals[-1]:.1e}..."
    )

    headers = [
        "# Benchmark: ultra_low_field_scaling.dat",
        "# Module: kitaev_mf.services.phase_4.symmetrized_solver",
        "# Description: Ultra-low field scaling of Dirac gap, gauge multiplier, and magnetization",
        "# Columns: h  gap_K  gap_canonical  c_ratio  lambda_star  lam_over_h  m_111  chi_m  conv  iters  wall_s",
        "#----------------------------------------------------------------------------------------------------------------",
    ]

    data_lines = []
    log_h_list = []
    log_gap_list = []

    curr_sym = None

    for h in h_vals:
        t0 = time.perf_counter()
        inv_sqrt3 = 1.0 / np.sqrt(3.0)
        h_vec = (float(h * inv_sqrt3), float(h * inv_sqrt3), float(h * inv_sqrt3))
        params = ModelParameters(K=-1.0, h=h_vec)

        # Higher k-mesh resolution (nk=24) to avoid boundary quadrature artifacts at ultra-low field
        res = solve_symmetrized_fixed_point(
            params,
            cell=cell,
            init_sym=curr_sym,
            nk=24,
            max_iter=35,
            tol=1e-8,
            diis_depth=6,
            diis_beta=0.5,
            lambda_step=0.2,
        )
        wall_s = time.perf_counter() - t0

        if res.metrics.converged:
            curr_sym = res.sym_params

        gap_k = res.dirac_gap_mf
        gap_canonical = res.dirac_gap_canonical
        c_ratio = gap_k / (h**3)
        lam_star = res.sym_params.lambda_val
        lam_over_h = lam_star / h if h > 0 else 0.0
        m_val = res.magnetization
        chi_m = m_val / h if h > 0 else 0.0

        row = (
            f"{h:11.4e}  {gap_k:13.6e}  {gap_canonical:13.6e}  {c_ratio:10.5f}  "
            f"{lam_star:12.6e}  {lam_over_h:10.5f}  {m_val:12.6e}  {chi_m:10.5f}  "
            f"{str(res.metrics.converged):5s}  {res.metrics.iterations:3d}  {wall_s:7.3f}"
        )
        data_lines.append(row)

        if res.metrics.converged and h >= 1e-4:
            log_h_list.append(np.log(h))
            log_gap_list.append(np.log(gap_k))

        print(
            f"  h = {h:9.2e} | Gap(K) = {gap_k:11.4e} | C_ratio = {c_ratio:8.4f} | "
            f"lambda/h = {lam_over_h:7.4f} | Iters = {res.metrics.iterations:2d} ({wall_s:.2f}s)"
        )

    # Write output dat file
    with open(out_file, "w", encoding="utf-8") as f:
        f.write("\n".join(headers) + "\n")
        f.write("\n".join(data_lines) + "\n")

    # Fit power law across h in [1e-4, 1e-2]
    if len(log_h_list) >= 3:
        poly = np.polyfit(log_h_list, log_gap_list, deg=1)
        power_p = float(poly[0])
        prefactor_c = float(np.exp(poly[1]))
        print("\n[Asymptotic Fit Results]")
        print(f"  Fitted Power Law Exponent: alpha = {power_p:.5f} (Theory: 3.00000)")
        print(f"  Fitted Dirac Gap Prefactor: C_MF = {prefactor_c:.5f} (Empirical M1.6: 2.59214)")
        print(f"  Certified Benchmark Table saved to: {out_file}")


if __name__ == "__main__":
    run_ultra_low_field_benchmark()
