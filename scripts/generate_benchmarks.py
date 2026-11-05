#!/usr/bin/env python3
"""Export verified numerical benchmark data to structured plain text .dat files.

Usage:
    python3 scripts/generate_benchmarks.py --all
    python3 scripts/generate_benchmarks.py --benchmark k006
    python3 scripts/generate_benchmarks.py --benchmark hellmann-feynman
    python3 scripts/generate_benchmarks.py --benchmark susceptibility
    python3 scripts/generate_benchmarks.py --benchmark klein-duality
"""

from __future__ import annotations

import argparse
import datetime
import subprocess
from pathlib import Path

import numpy as np

from kitaev_mf.domain.lattice import (
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.observables import compute_sublattice_magnetizations
from kitaev_mf.services.quadrature import midpoint_convergence
from kitaev_mf.services.response import (
    evaluate_finite_difference_magnetization,
    evaluate_klein_dual_pair,
    evaluate_thermodynamic_susceptibility,
)
from kitaev_mf.services.solver import solve_self_consistent_mean_field

BENCHMARK_DIR = Path("data/processed/benchmarks")


def get_git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def write_dat_file(filename: str, header_lines: list[str], data_lines: list[str]) -> Path:
    BENCHMARK_DIR.mkdir(parents=True, exist_ok=True)
    out_path = BENCHMARK_DIR / filename
    timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    commit = get_git_commit()

    content = [
        f"# Benchmark: {filename}",
        f"# Generated: {timestamp}",
        f"# Git Commit: {commit}",
    ]
    content.extend([f"# {h}" for h in header_lines])
    content.append("#" + "-" * 78)
    content.extend(data_lines)
    content.append("")

    out_path.write_text("\n".join(content), encoding="utf-8")
    print(f"Exported benchmark data -> {out_path}")
    return out_path


def export_k006_quadrature() -> None:
    meshes = (96, 192, 384, 768, 1536, 3072, 6144, 12288)
    estimates = midpoint_convergence(meshes, row_chunk_size=128)

    headers = [
        "Description: K-006 pure-Kitaev matter correlator midpoint BZ quadrature",
        "Columns: mesh_size  U_bar  observed_order  estimated_error  richardson_extrapolated",
    ]
    rows = []
    for est in estimates:
        order_val = est.observed_order if est.observed_order is not None else float("nan")
        err_val = est.estimated_error if est.estimated_error is not None else float("nan")
        rich_val = est.richardson_value if est.richardson_value is not None else float("nan")
        rows.append(
            f"{est.mesh_size:6d}  {est.value:18.13f}  {order_val:10.4f}  {err_val:14.6e}  {rich_val:18.13f}"
        )

    write_dat_file("k006_quadrature.dat", headers, rows)


def export_hellmann_feynman_scaling(backend: TorchBackend) -> None:
    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, J=0.04, h=(0.02, 0.02, 0.03))
    config = SolverConfig(
        k_grid_shape=(16, 16), tol_correlator=1e-9, tol_gauge=1e-9, max_iterations=60
    )
    deltas = [0.008, 0.004, 0.002, 0.001]

    headers = [
        "Description: Variational Hellmann-Feynman derivative error scaling vs central difference step delta",
        "Columns: delta  de_numerical  neg_sigma_tot  discrepancy  observed_order",
    ]
    results = [
        evaluate_finite_difference_magnetization(
            cell=cell, params=params, h_dir=2, delta=d, config=config, backend=backend
        )
        for d in deltas
    ]

    rows = []
    for i, (d, res) in enumerate(zip(deltas, results, strict=True)):
        if i > 0 and results[i - 1].discrepancy > 0.0 and res.discrepancy > 0.0:
            order = np.log2(results[i - 1].discrepancy / res.discrepancy)
        else:
            order = float("nan")
        rows.append(
            f"{d:8.5f}  {res.de_numerical:14.8f}  {res.neg_sigma_tot:14.8f}  {res.discrepancy:14.6e}  {order:10.4f}"
        )

    write_dat_file("hellmann_feynman_scaling.dat", headers, rows)


def export_susceptibility_curve(backend: TorchBackend) -> None:
    cell = make_honeycomb_2site()
    h_z_vals = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10]
    config = SolverConfig(
        k_grid_shape=(16, 16), tol_correlator=1e-8, tol_gauge=1e-8, max_iterations=50
    )

    headers = [
        "Description: Thermodynamic susceptibility chi^zz along external Zeeman field [0, 0, hz]",
        "Columns: hz  sigma_tot_z  susceptibility_chi_zz",
    ]
    rows = []
    for hz in h_z_vals:
        params = ModelParameters(K=-1.0, J=0.03, h=(0.0, 0.0, hz))
        chi = evaluate_thermodynamic_susceptibility(
            cell=cell, params=params, h_dir=2, delta=1e-3, config=config, backend=backend
        )
        init_s = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.4)
        state, _ = solve_self_consistent_mean_field(cell, params, init_s, config, backend)
        sigma_mag, _ = compute_sublattice_magnetizations(state, cell)
        tot_mag = float(np.sum(sigma_mag[:, 2]))
        rows.append(f"{hz:8.4f}  {tot_mag:14.8f}  {chi:14.8f}")

    write_dat_file("susceptibility_curve.dat", headers, rows)


def export_klein_duality_splitting(backend: TorchBackend) -> None:
    cell = make_honeycomb_4site_rectangular()
    config = SolverConfig(
        k_grid_shape=(12, 12), tol_correlator=1e-8, tol_gauge=1e-8, max_iterations=50
    )
    epsilons = [0.04, 0.02, 0.01, 0.005, 0.002]

    headers = [
        "Description: Symmetric Klein duality finite-cell boundary splitting on 4-sublattice cell",
        "Columns: eps_J  eps_Gamma  E_site  E_site_dual  energy_diff",
    ]
    rows = []
    for eps in epsilons:
        res = evaluate_klein_dual_pair(
            cell=cell,
            K0=-1.0,
            eps_J=eps,
            eps_Gamma=eps,
            eps_Gamma_prime=-0.01,
            config=config,
            backend=backend,
        )
        rows.append(
            f"{eps:8.4f}  {eps:8.4f}  {res.energy_site:14.8f}  {res.energy_dual_site:14.8f}  {res.energy_diff:14.6e}"
        )

    write_dat_file("klein_duality_splitting.dat", headers, rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Export MMFT benchmark data files")
    parser.add_argument(
        "--benchmark",
        choices=["all", "k006", "hellmann-feynman", "susceptibility", "klein-duality"],
        default="all",
        help="Benchmark to generate",
    )
    args = parser.parse_args()
    backend = TorchBackend(device="cpu")

    if args.benchmark in ("all", "k006"):
        export_k006_quadrature()
    if args.benchmark in ("all", "hellmann-feynman"):
        export_hellmann_feynman_scaling(backend)
    if args.benchmark in ("all", "susceptibility"):
        export_susceptibility_curve(backend)
    if args.benchmark in ("all", "klein-duality"):
        export_klein_duality_splitting(backend)


if __name__ == "__main__":
    main()
