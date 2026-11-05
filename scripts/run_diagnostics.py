#!/usr/bin/env python3
"""Unified numerical diagnostics for variational response, KKT stationarity, and dualities.

Usage:
    python3 scripts/run_diagnostics.py --suite all
    python3 scripts/run_diagnostics.py --suite hellmann-feynman
    python3 scripts/run_diagnostics.py --suite klein-duality
"""

from __future__ import annotations

import argparse
from math import log2

import numpy as np

from kitaev_mf.domain.lattice import (
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.correlators import compute_spectral_projector
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.observables import compute_sublattice_magnetizations
from kitaev_mf.services.response import (
    compute_riemannian_gradient_norm,
    evaluate_finite_difference_magnetization,
    evaluate_klein_dual_pair,
    evaluate_thermodynamic_susceptibility,
)
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def run_kkt_stationarity_diagnostic(backend: TorchBackend) -> None:
    print("\n### 1. KKT Riemannian Gradient Norm at Self-Consistency")
    print(
        "| State / Regimes | Iterations | Correlator Res. | Gauge Res. | Max ||[H(k), P(k)]||_F |"
    )
    print("|---|---:|---:|---:|---:|")

    cell = make_honeycomb_2site()
    cases = [
        ("Pure Kitaev (K=-1)", ModelParameters(K=-1.0)),
        ("Heisenberg-Kitaev (J=0.05)", ModelParameters(K=-1.0, J=0.05)),
        ("Field-Polarized (h=0.04)", ModelParameters(K=-1.0, J=0.03, h=(0.04, 0.02, 0.03))),
    ]
    config = SolverConfig(
        k_grid_shape=(16, 16), tol_correlator=1e-9, tol_gauge=1e-9, max_iterations=60
    )
    k_grid = cell.make_k_grid(16, 16)

    for label, params in cases:
        init_s = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.4)
        state, metrics = solve_self_consistent_mean_field(cell, params, init_s, config, backend)

        h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, backend)
        p_batch = compute_spectral_projector(h_batch, backend)
        grad_norm = compute_riemannian_gradient_norm(h_batch, p_batch)

        print(
            f"| {label} | {metrics.iterations} | {metrics.final_correlator_residual:.2e} | "
            f"{metrics.final_gauge_residual:.2e} | {grad_norm:.2e} |"
        )


def run_hellmann_feynman_diagnostic(backend: TorchBackend) -> None:
    print("\n### 2. Hellmann-Feynman Error Scaling vs. Step Size delta")
    print("| Step delta | dE_num | -<sigma_tot^z> | Discrepancy Delta | Observed Order p |")
    print("|---:|---:|---:|---:|---:|")

    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, J=0.04, h=(0.02, 0.02, 0.03))
    config = SolverConfig(
        k_grid_shape=(16, 16), tol_correlator=1e-9, tol_gauge=1e-9, max_iterations=60
    )

    deltas = [0.008, 0.004, 0.002, 0.001]
    results = []

    for delta in deltas:
        res = evaluate_finite_difference_magnetization(
            cell=cell,
            params=params,
            h_dir=2,
            delta=delta,
            config=config,
            backend=backend,
        )
        results.append(res)

    for i, (delta, res) in enumerate(zip(deltas, results, strict=True)):
        if i > 0 and results[i - 1].discrepancy > 0.0 and res.discrepancy > 0.0:
            order = log2(results[i - 1].discrepancy / res.discrepancy)
            order_str = f"{order:.3f}"
        else:
            order_str = "--"

        print(
            f"| {delta:.4f} | {res.de_numerical:.6f} | {res.neg_sigma_tot:.6f} | "
            f"{res.discrepancy:.2e} | {order_str} |"
        )


def run_susceptibility_diagnostic(backend: TorchBackend) -> None:
    print("\n### 3. Thermodynamic Susceptibility along [0, 0, 1]")
    print("| Field h_z | Total Magnetization <sigma^z> | Susceptibility chi^zz | Convexity Check |")
    print("|---:|---:|---:|:---:|")

    cell = make_honeycomb_2site()
    h_z_vals = [0.01, 0.02, 0.03, 0.04, 0.05]
    config = SolverConfig(
        k_grid_shape=(16, 16), tol_correlator=1e-8, tol_gauge=1e-8, max_iterations=50
    )

    for hz in h_z_vals:
        params = ModelParameters(K=-1.0, J=0.03, h=(0.0, 0.0, hz))
        chi = evaluate_thermodynamic_susceptibility(
            cell=cell,
            params=params,
            h_dir=2,
            delta=1e-3,
            config=config,
            backend=backend,
        )
        init_s = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.4)
        state, _ = solve_self_consistent_mean_field(cell, params, init_s, config, backend)
        sigma_mag, _ = compute_sublattice_magnetizations(state, cell)
        tot_mag = float(np.sum(sigma_mag[:, 2]))
        check = "PASS (chi >= 0)" if chi >= -1e-6 else "FAIL"

        print(f"| {hz:.2f} | {tot_mag:.6f} | {chi:.4f} | {check} |")


def run_klein_duality_diagnostic(backend: TorchBackend) -> None:
    print("\n### 4. Symmetrically Parameterized Klein Duality on 4-Site Cell")
    print("| eps_J | eps_Gamma | E_site(p) | E_site(p_dual) | Splitting |delta E| |")
    print("|---:|---:|---:|---:|---:|")

    cell = make_honeycomb_4site_rectangular()
    config = SolverConfig(
        k_grid_shape=(12, 12), tol_correlator=1e-8, tol_gauge=1e-8, max_iterations=50
    )
    epsilons = [0.04, 0.02, 0.01, 0.005]

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
        print(
            f"| {eps:.3f} | {eps:.3f} | {res.energy_site:.6f} | {res.energy_dual_site:.6f} | "
            f"{res.energy_diff:.2e} |"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description="MMFT response and duality diagnostics")
    parser.add_argument(
        "--suite",
        choices=["all", "kkt", "hellmann-feynman", "susceptibility", "klein-duality"],
        default="all",
        help="Diagnostic suite to execute",
    )
    args = parser.parse_args()
    backend = TorchBackend(device="cpu")

    if args.suite in ("all", "kkt"):
        run_kkt_stationarity_diagnostic(backend)
    if args.suite in ("all", "hellmann-feynman"):
        run_hellmann_feynman_diagnostic(backend)
    if args.suite in ("all", "susceptibility"):
        run_susceptibility_diagnostic(backend)
    if args.suite in ("all", "klein-duality"):
        run_klein_duality_diagnostic(backend)


if __name__ == "__main__":
    main()
