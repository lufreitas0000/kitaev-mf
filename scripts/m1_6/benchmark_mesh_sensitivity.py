#!/usr/bin/env python3
"""Benchmark Brillouin-zone mesh sensitivity and Dirac-point sampling (AH-03, AH-04).

Usage:
    python3 scripts/m1_6/benchmark_mesh_sensitivity.py
"""

from __future__ import annotations

import time

import numpy as np
import torch

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.observables import (
    compute_plaquette_fluxes,
    compute_sublattice_magnetizations,
)
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field
from scripts.m1_6._common import (
    append_jsonl_record,
    setup_cpu_environment,
    write_dat_benchmark,
)


def compute_spectral_gap(h_batch: torch.Tensor, backend: TorchBackend) -> float:
    """Compute the single-particle excitation gap across the Brillouin zone."""
    eigvals, _ = backend.batched_eigh(h_batch)
    # Dimension is 8: occupied bands are 0, 1, 2, 3; unoccupied are 4, 5, 6, 7
    gaps_per_k = eigvals[:, 4] - eigvals[:, 3]
    return float(torch.min(gaps_per_k).item())


def run_mesh_sensitivity_benchmark() -> None:
    """Evaluate solver contraction, gap, and energy density across Nk in {12, 16, 24, 32}."""
    setup_cpu_environment(1)
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()

    mesh_sizes = [12, 16, 24, 32]
    fields = [0.0, 0.05, 0.5, 2.0]
    solvers = ["diis", "broyden"]

    headers = [
        "Description: Brillouin zone mesh sensitivity and Dirac sampling benchmark (M1.6)",
        "Columns: Nk  div3  h  solver  conv  n_phi  wall_s  ms_per_phi  E_site  m_111  W_p  norm_lam  gap  eps_comm",
    ]
    data_lines = []

    print("[Mesh Benchmark] Evaluating grid sensitivity across Nk and fields...")

    for nk in mesh_sizes:
        is_div3 = nk % 3 == 0
        k_grid = cell.make_k_grid(nk, nk)

        for h_val in fields:
            if h_val == 0.0:
                params = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))
            else:
                h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
                params = ModelParameters(K=-1.0, h=h_vec)

            init_state = make_qsl_seed(cell)

            for solver in solvers:
                config = SolverConfig(
                    k_grid_shape=(nk, nk),
                    max_iterations=40,
                    tol_correlator=1e-8,
                    tol_gauge=1e-8,
                    method=solver,  # type: ignore[arg-type]
                    mixing_alpha=0.6,
                    broyden_alpha=0.6,
                    skew_mode="retract_step",
                )

                t0 = time.perf_counter()
                try:
                    state_final, metrics = solve_self_consistent_mean_field(
                        cell, params, init_state, config, backend
                    )
                    wall_s = time.perf_counter() - t0
                    conv = metrics.converged
                    n_phi = metrics.map_evaluations
                    ms_per_phi = (wall_s / max(n_phi, 1)) * 1000.0

                    eb = compute_gaussian_energy(state_final, cell, params)
                    e_site = eb.per_site
                    sig_mag, _ = compute_sublattice_magnetizations(state_final, cell)
                    m_vec = 0.5 * (sig_mag[0] + sig_mag[1])
                    m_111 = float(np.dot(m_vec, np.ones(3) / np.sqrt(3.0)))
                    fluxes = compute_plaquette_fluxes(state_final, cell)
                    w_p = float(fluxes[0])
                    norm_lam = float(np.linalg.norm(state_final.lambdas))

                    # Compute final spectral gap
                    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state_final, backend)
                    gap_val = compute_spectral_gap(h_batch, backend)
                    eps_comm = (
                        metrics.comm_residuals[-1] if metrics.comm_residuals else float("nan")
                    )

                except Exception as exc:
                    wall_s = time.perf_counter() - t0
                    conv = False
                    n_phi = 0
                    ms_per_phi = float("nan")
                    e_site = float("nan")
                    m_111 = float("nan")
                    w_p = float("nan")
                    norm_lam = float("nan")
                    gap_val = float("nan")
                    eps_comm = float("nan")
                    print(f"Failed Nk={nk} h={h_val} solver={solver}: {exc}")

                row = (
                    f"{nk:4d}  {str(is_div3):5s}  {h_val:6.3f}  {solver:8s}  {str(conv):5s}  "
                    f"{n_phi:6d}  {wall_s:8.4f}  {ms_per_phi:10.2f}  {e_site:12.6f}  "
                    f"{m_111:9.5f}  {w_p:8.4f}  {norm_lam:10.6e}  {gap_val:10.6f}  {eps_comm:12.4e}"
                )
                data_lines.append(row)
                print(f"  -> {row}")

                append_jsonl_record(
                    "runs_mesh_sensitivity.jsonl",
                    {
                        "nk": nk,
                        "div3": is_div3,
                        "h": h_val,
                        "solver": solver,
                        "converged": conv,
                        "n_phi": n_phi,
                        "wall_s": wall_s,
                        "ms_per_phi": ms_per_phi,
                        "e_site": e_site,
                        "m_111": m_111,
                        "w_p": w_p,
                        "norm_lam": norm_lam,
                        "gap": gap_val,
                        "eps_comm": eps_comm,
                    },
                )

    write_dat_benchmark("mesh_sensitivity.dat", headers, data_lines)
    print("[Mesh Benchmark] Finished successfully.")


if __name__ == "__main__":
    run_mesh_sensitivity_benchmark()
