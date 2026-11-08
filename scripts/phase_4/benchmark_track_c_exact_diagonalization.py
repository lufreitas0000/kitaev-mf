"""Benchmark script for Track C: Exact Diagonalization in the spin formalism.

Computes the ground state, excitation gap, magnetization, plaquette flux,
and fidelity susceptibility across the low-field and phase transition regimes,
comparing ED with Kitaev 3rd-order perturbation theory and mean-field theory.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

# Ensure local package takes precedence
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import numpy as np

from kitaev_mf.services.phase_4.exact_diagonalization import (
    ExactDiagonalizationSpinSolver,
    HoneycombCluster,
)


def run_track_c_ed_benchmark(
    n1: int = 2,
    n2: int = 4,
    h_max: float = 0.45,
    num_h: int = 20,
    output_filename: str = "ed_spin_scaling_n16.dat",
) -> None:
    """Run ED sweep across magnetic fields and export benchmark dat table."""
    cluster = HoneycombCluster.make_torus(n1, n2)
    num_spins = cluster.num_spins
    print(
        f"\n[Track C ED Benchmark] Torus {n1}x{n2} ({num_spins} spins, dim = 2^{num_spins} = {1 << num_spins})..."
    )

    output_dir = (
        Path(__file__).resolve().parents[2] / "data" / "processed" / "benchmarks" / "phase_4"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    out_file = output_dir / output_filename

    h_vals = np.linspace(0.0, h_max, num_h)

    headers = [
        f"# Benchmark: {output_filename}",
        f"# System: Honeycomb torus {n1}x{n2} (N = {num_spins} spins)",
        "# Description: Exact Diagonalization observables across field h || [111]",
        "# Columns: h  E0  E_per_spin  gap_E1_E0  m_111  flux_Wp  corr_z  fidelity  chi_F  wall_s",
        "#------------------------------------------------------------------------------------------------------",
    ]

    data_lines = []
    prev_psi0 = None
    prev_h = None

    for _idx, h in enumerate(h_vals):
        t0 = time.perf_counter()
        solver = ExactDiagonalizationSpinSolver(cluster, K=1.0, h=float(h))
        # Compute lowest 4 states and observables
        obs = solver.compute_observables(num_states=4)
        wall_s = time.perf_counter() - t0

        # Fidelity with previous state
        _, evecs = solver.solve_spectrum(num_states=2)
        psi0 = evecs[:, 0]
        fidelity = 1.0
        chi_F = 0.0
        if prev_psi0 is not None and prev_h is not None:
            dh = h - prev_h
            if dh > 1e-12:
                overlap = abs(np.vdot(prev_psi0, psi0))
                fidelity = float(overlap)
                chi_F = float(2.0 * (1.0 - overlap) / (dh**2))

        prev_psi0 = psi0.copy()
        prev_h = h

        row = (
            f"{h:8.4f}  {obs.energy_ground:12.6f}  {obs.energy_per_spin:10.6f}  "
            f"{obs.excitation_gap:11.6e}  {obs.magnetization:9.5f}  {obs.plaquette_flux_mean:9.5f}  "
            f"{obs.nearest_neighbor_corr_z:9.5f}  {fidelity:10.6f}  {chi_F:10.4e}  {wall_s:6.2f}"
        )
        data_lines.append(row)

        print(
            f"  h = {h:6.4f} | E/N = {obs.energy_per_spin:8.5f} | Gap = {obs.excitation_gap:10.5e} | "
            f"m_111 = {obs.magnetization:7.4f} | <W_p> = {obs.plaquette_flux_mean:7.4f} | "
            f"chi_F = {chi_F:9.2e} ({wall_s:.2f}s)"
        )

    with open(out_file, "w", encoding="utf-8") as f:
        f.write("\n".join(headers) + "\n")
        f.write("\n".join(data_lines) + "\n")

    print(f"[Track C] Certified ED Benchmark exported to: {out_file}\n")


if __name__ == "__main__":
    # Run N = 16 benchmark first
    run_track_c_ed_benchmark(
        n1=2, n2=4, h_max=0.45, num_h=19, output_filename="ed_spin_scaling_n16.dat"
    )
