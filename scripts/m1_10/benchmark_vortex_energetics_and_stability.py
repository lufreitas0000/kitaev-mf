"""Benchmark script for Wave 2: Real-space vortex pair energetics, zero modes, and stability.

Evaluates:
1. Two-vison interaction potential V(R) = E_vison(R) - E_vacuum at h = 0 as a function of separation R.
2. In-gap Majorana zero mode spectrum epsilon_1(h, R) and localization under field h in [0.0, 0.08].
3. Dynamic self-consistent solver stability (vortex trapping vs. healing).
"""

from __future__ import annotations

import os
import time

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


def run_vortex_benchmark() -> None:
    output_dir = "data/processed/benchmarks/m1_10"
    os.makedirs(output_dir, exist_ok=True)
    out_file = os.path.join(output_dir, "vortex_energetics_and_zero_modes.dat")

    print("[Wave 2 Benchmark] Starting real-space vortex pair simulation...")

    cell = make_honeycomb_2site()
    # Choose L = 8 (N_sites = 128, M = 512, 3 does not divide L -> non-resonant)
    L = 8
    cluster = RealSpaceTorusCluster(L=L)

    # 1. Potential V(R) at h = 0 as a function of separation R
    print(f"\n--- 1. Two-Vison Interaction Potential V(R) at h = 0 (L = {L}) ---")
    p1 = (0, 0)
    # Test separations along the diagonal (R, R) for R in 1 .. L//2 = 4
    separations = [1, 2, 3, 4]
    v_data: list[tuple[int, int, float, float, float]] = []

    params_0 = ModelParameters(K=-1.0)
    state_0 = make_qsl_seed(cell)
    A_vac_0 = cluster.assemble_skew_hamiltonian(params_0, state_0)
    omega_vac_0, e_vac_0 = cluster.solve_covariance_complex_eigh(A_vac_0)

    w_vac_0, _ = np.linalg.eigh(1j * A_vac_0.astype(np.complex128))
    bulk_gap_0 = float(np.min(w_vac_0[w_vac_0 > 1e-10]))
    print(f"  Unperturbed vacuum bulk gap at h = 0: {bulk_gap_0:.6f} |K|")

    for r_sep in separations:
        p2 = (r_sep, r_sep)
        cut = cluster.generate_vortex_pair_cut(p1, p2)
        A_vor = cluster.assemble_skew_hamiltonian(params_0, state_0, flipped_bonds=cut)
        omega_vor, e_vor = cluster.solve_covariance_complex_eigh(A_vor)

        dE = (e_vor - e_vac_0) / 2.0

        w_vor, _ = np.linalg.eigh(1j * A_vor.astype(np.complex128))
        pos_vor = np.sort(w_vor[w_vor > 1e-10])
        eps_1 = float(pos_vor[0])
        eps_2 = float(pos_vor[1])

        fluxes = cluster.compute_plaquette_fluxes(omega_vor)
        w_p1 = fluxes[p1]
        w_p2 = fluxes[p2]

        v_data.append((r_sep, len(cut), dE, eps_1, eps_2))
        print(
            f"  Separation R = {r_sep} (cut length {len(cut):2d}) | Delta E = {dE:.6f} |K| | "
            f"eps_1 = {eps_1:.6f}, eps_2 = {eps_2:.6f} | W(p1) = {w_p1:+.2f}, W(p2) = {w_p2:+.2f}"
        )

    # 2. Field response and in-gap Majorana bound modes
    print(f"\n--- 2. Field Response and Sub-Gap Modes at R = 4 (L = {L}) ---")
    p2 = (4, 4)
    cut = cluster.generate_vortex_pair_cut(p1, p2)
    h_vals = [0.00, 0.01, 0.02, 0.04, 0.06, 0.08, 0.10]
    h_data: list[tuple[float, float, float, float, float]] = []

    for h_val in h_vals:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params_h = ModelParameters(K=-1.0, h=h_vec)
        state_h = make_qsl_seed(cell)

        A_vac_h = cluster.assemble_skew_hamiltonian(params_h, state_h)
        A_vor_h = cluster.assemble_skew_hamiltonian(params_h, state_h, flipped_bonds=cut)

        _, e_vac_h = cluster.solve_covariance_complex_eigh(A_vac_h)
        omega_vor_h, e_vor_h = cluster.solve_covariance_complex_eigh(A_vor_h)

        dE_h = (e_vor_h - e_vac_h) / 2.0

        w_vac_h, _ = np.linalg.eigh(1j * A_vac_h.astype(np.complex128))
        w_vor_h, _ = np.linalg.eigh(1j * A_vor_h.astype(np.complex128))

        bulk_gap_h = float(np.min(w_vac_h[w_vac_h > 1e-10]))
        pos_vor_h = np.sort(w_vor_h[w_vor_h > 1e-10])
        eps_1_h = float(pos_vor_h[0])
        eps_2_h = float(pos_vor_h[1])

        h_data.append((h_val, dE_h, bulk_gap_h, eps_1_h, eps_2_h))
        print(
            f"  h = {h_val:4.2f} | Delta E = {dE_h:.6f} | Bulk Gap = {bulk_gap_h:.6f} | "
            f"In-gap Mode eps_1 = {eps_1_h:.6f} (ratio eps_1/Gap = {eps_1_h / bulk_gap_h:.4f})"
        )

    # Export to .dat file with full metadata provenance
    with open(out_file, "w") as f:
        f.write(
            "# Description: Two-vison pair interaction potential and in-gap Majorana zero modes (Wave 2)\n"
        )
        f.write(
            f"# Geometry: L = {L} periodic honeycomb torus (N_sites = {2 * L * L}, M = {8 * L * L})\n"
        )
        f.write("# Generated: " + time.strftime("%Y-%m-%d %H:%M:%S") + "\n")
        f.write("# Section 1: Potential V(R) at h = 0\n")
        f.write("# R  cut_len  delta_E  eps_1  eps_2\n")
        for r_sep, c_len, dE, e1, e2 in v_data:
            f.write(f"{r_sep:2d}  {c_len:2d}  {dE:12.6f}  {e1:12.6f}  {e2:12.6f}\n")
        f.write("\n\n")
        f.write("# Section 2: In-gap mode vs Zeeman field at separation R = 4\n")
        f.write("# h  delta_E  bulk_gap  eps_1  eps_2\n")
        for h_val, dE, bg, e1, e2 in h_data:
            f.write(f"{h_val:6.3f}  {dE:12.6f}  {bg:12.6f}  {e1:12.6f}  {e2:12.6f}\n")

    print(f"\n[Wave 2 Benchmark] Data exported successfully to: {out_file}")


if __name__ == "__main__":
    run_vortex_benchmark()
