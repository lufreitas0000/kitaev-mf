"""Multi-site supercell benchmark comparing 2-site vs 4-site cell (M1.6 Phase 3).

Evaluates energy per physical site, discrete C_3 nematicity, plaquette flux uniformity,
and O(N_site^3) Hamiltonian diagonalization scaling.
"""

from __future__ import annotations

import time

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site, make_honeycomb_4site_rectangular
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.observables import (
    compute_plaquette_fluxes,
    compute_sublattice_magnetizations,
)
from kitaev_mf.services.seeds import make_polarized_seed, make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field
from scripts.m1_6._common import (
    append_jsonl_record,
    write_dat_benchmark,
)


def run_supercell_benchmark() -> None:
    backend = TorchBackend(device="cpu")
    cell_2 = make_honeycomb_2site()
    cell_4 = make_honeycomb_4site_rectangular()

    h_points = [0.00, 0.05, 0.10, 0.20, 0.30, 0.34, 0.40, 0.50, 1.00, 2.00]

    # Use comparable total BZ sampling: 24x24 for 2-site (area 1), 16x16 for 4-site (area 2)
    config_2 = SolverConfig(
        k_grid_shape=(24, 24),
        max_iterations=40,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        method="diis",
        mixing_alpha=0.6,
        skew_mode="retract_step",
    )
    config_4 = SolverConfig(
        k_grid_shape=(16, 16),
        max_iterations=40,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        method="diis",
        mixing_alpha=0.6,
        skew_mode="retract_step",
    )

    headers = [
        "Description: Supercell comparison (2-site vs 4-site) across magnetic field h in [0, 2] (M1.6 Phase 3)",
        "Columns: h  cell_type  seed_type  conv  n_phi  wall_s  ms_per_phi  E_site  m_111  nematicity  W_p_avg  delta_W_p",
    ]
    data_lines = []

    print("[Supercell Benchmark] Comparing 2-site vs 4-site cell...")

    for h_val in h_points:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)

        for cell_name, cell, config in [("2site", cell_2, config_2), ("4site", cell_4, config_4)]:
            for seed_type in ["qsl", "pol"]:
                if seed_type == "qsl":
                    seed = make_qsl_seed(cell)
                else:
                    seed = make_polarized_seed(cell, axis="111")

                t0 = time.perf_counter()
                state_final, metrics = solve_self_consistent_mean_field(
                    cell, params, seed, config, backend
                )
                wall_s = time.perf_counter() - t0

                ms_per_phi = (wall_s / max(1, metrics.map_evaluations)) * 1000.0
                eb = compute_gaussian_energy(state_final, cell, params)
                e_site = eb.per_site

                sig_mag, _ = compute_sublattice_magnetizations(state_final, cell)
                m_avg = np.mean(sig_mag, axis=0)
                m_111 = float(np.dot(m_avg, np.ones(3) / np.sqrt(3.0)))

                # Nematicity: deviation from C3 isotropy
                nematicity = float(
                    np.mean([np.abs(s[0] - s[1]) + np.abs(s[1] - s[2]) for s in sig_mag])
                )

                fluxes = compute_plaquette_fluxes(state_final, cell)
                w_p_avg = float(np.mean(fluxes))
                delta_w_p = float(np.max(fluxes) - np.min(fluxes)) if len(fluxes) > 1 else 0.0

                row = (
                    f"{h_val:6.3f}  {cell_name:8s}  {seed_type:6s}  {str(metrics.converged):5s}  "
                    f"{metrics.map_evaluations:5d}  {wall_s:8.4f}  {ms_per_phi:10.2f}  "
                    f"{e_site:12.6f}  {m_111:9.5f}  {nematicity:12.4e}  {w_p_avg:8.4f}  {delta_w_p:10.4e}"
                )
                data_lines.append(row)

                rec = {
                    "h": float(h_val),
                    "cell": cell_name,
                    "seed": seed_type,
                    "conv": metrics.converged,
                    "n_phi": metrics.map_evaluations,
                    "wall_s": wall_s,
                    "ms_per_phi": ms_per_phi,
                    "e_site": e_site,
                    "m_111": m_111,
                    "nematicity": nematicity,
                    "w_p_avg": w_p_avg,
                    "delta_w_p": delta_w_p,
                }
                append_jsonl_record("runs_supercell.jsonl", rec)

                print(
                    f"  [{cell_name:5s}|{seed_type:3s}] h={h_val:4.2f} -> E={e_site:9.5f}, W_p={w_p_avg:6.4f}, nem={nematicity:.1e}, ms/phi={ms_per_phi:.1f}"
                )

    write_dat_benchmark("supercell_comparison.dat", headers, data_lines)
    print("[Supercell Benchmark] Finished and exported supercell_comparison.dat")


if __name__ == "__main__":
    run_supercell_benchmark()
