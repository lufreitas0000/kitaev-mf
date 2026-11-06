"""Benchmark measuring the Majorana excitation gap and cubic scaling (M1.6 Phase 3).

Computes the single-particle Dirac excitation gap Delta_K(h) = 2 epsilon_5(K)
and global bandgap Delta_min(h) across h in [0.0, 0.35].
Fits the power law Delta_K(h) = C * h^p and compares to Kitaev chiral perturbation theory (kappa ~ h^3 / K^2).
"""

from __future__ import annotations

import time

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field
from scripts.m1_6._common import (
    append_jsonl_record,
    write_dat_benchmark,
)


def run_majorana_gap_benchmark() -> None:
    backend = TorchBackend(device="cpu")
    cell = make_honeycomb_2site()

    h_vals = np.array(
        [
            0.000,
            0.005,
            0.010,
            0.015,
            0.020,
            0.030,
            0.040,
            0.050,
            0.075,
            0.100,
            0.125,
            0.150,
            0.175,
            0.200,
            0.250,
            0.300,
            0.350,
        ],
        dtype=np.float64,
    )

    config = SolverConfig(
        k_grid_shape=(24, 24),
        max_iterations=40,
        tol_correlator=1e-8,
        tol_gauge=1e-8,
        method="diis",
        mixing_alpha=0.6,
        skew_mode="retract_step",
    )

    k_dirac = np.array([[1.0 / 3.0, 2.0 / 3.0]], dtype=np.float64)

    # Dense path Gamma -> K -> M -> Gamma to find global minimum bandgap
    path_fracs = []
    # Gamma -> K
    for t in np.linspace(0.0, 1.0, 25, endpoint=False):
        path_fracs.append(t * k_dirac[0])
    # K -> M (M = (1/2, 1/2))
    k_m = np.array([0.5, 0.5])
    for t in np.linspace(0.0, 1.0, 25, endpoint=False):
        path_fracs.append((1.0 - t) * k_dirac[0] + t * k_m)
    # M -> Gamma
    for t in np.linspace(0.0, 1.0, 25, endpoint=True):
        path_fracs.append((1.0 - t) * k_m)
    dense_k_path = np.array(path_fracs, dtype=np.float64)

    headers = [
        "Description: Majorana single-particle bandgap and cubic power-law scaling (M1.6 Phase 3)",
        "Columns: h  gap_K  gap_global  e5_K  e4_K  conv  n_phi  wall_s",
    ]
    data_lines = []

    print(f"[Majorana Gap] Evaluating excitation gaps across {len(h_vals)} field points...")

    # Forward continuation from zero field
    curr_state = make_qsl_seed(cell)

    fit_h: list[float] = []
    fit_gap: list[float] = []

    for h_val in h_vals:
        h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)

        t0 = time.perf_counter()
        state_final, metrics = solve_self_consistent_mean_field(
            cell, params, curr_state, config, backend
        )
        wall_s = time.perf_counter() - t0

        if metrics.converged:
            curr_state = state_final

        # 1. Direct gap at Dirac point K
        h_k = build_bloch_hamiltonian(k_dirac, cell, params, state_final, backend)
        evals_k, _ = backend.eigh(h_k)
        evals_k_np = evals_k.cpu().numpy()[0]
        e4_k = float(evals_k_np[3])
        e5_k = float(evals_k_np[4])
        gap_k = float(e5_k - e4_k)

        # 2. Global gap along high-symmetry path
        h_path = build_bloch_hamiltonian(dense_k_path, cell, params, state_final, backend)
        evals_path, _ = backend.eigh(h_path)
        evals_path_np = evals_path.cpu().numpy()  # (N_pts, 8)
        direct_gaps = evals_path_np[:, 4] - evals_path_np[:, 3]
        gap_global = float(np.min(direct_gaps))

        row = (
            f"{h_val:6.3f}  {gap_k:12.6e}  {gap_global:12.6e}  {e5_k:12.6e}  {e4_k:12.6e}  "
            f"{str(metrics.converged):5s}  {metrics.map_evaluations:5d}  {wall_s:8.4f}"
        )
        data_lines.append(row)

        rec = {
            "h": float(h_val),
            "gap_k": gap_k,
            "gap_global": gap_global,
            "e5_k": e5_k,
            "e4_k": e4_k,
            "conv": metrics.converged,
            "n_phi": metrics.map_evaluations,
            "wall_s": wall_s,
        }
        append_jsonl_record("runs_majorana_gap.jsonl", rec)

        if 0.005 <= h_val <= 0.150:
            fit_h.append(float(h_val))
            fit_gap.append(gap_k)

        print(
            f"  -> h = {h_val:5.3f} | Gap(K) = {gap_k:10.4e} | Gap(global) = {gap_global:10.4e} | Iter = {metrics.map_evaluations:2d}"
        )

    write_dat_benchmark("majorana_gap_scaling.dat", headers, data_lines)

    # Power law fit: log(gap) = p * log(h) + log(C)
    log_h = np.log(fit_h)
    log_g = np.log(fit_gap)
    p_fit, log_c = np.polyfit(log_h, log_g, deg=1)
    c_fit = float(np.exp(log_c))

    print("\n[Majorana Gap] ===============================================")
    print(f"[Majorana Gap] Fitted Power Law: Delta_K(h) = {c_fit:.4f} * h^{p_fit:.4f}")
    print(f"[Majorana Gap] Exponent p = {p_fit:.4f} (Expected cubic: p = 3.000)")
    print("[Majorana Gap] ===============================================\n")


if __name__ == "__main__":
    run_majorana_gap_benchmark()
