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
import torch

from kitaev_mf.domain.lattice import (
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.hamiltonian import (
    build_bloch_hamiltonian,
    build_bond_hopping_matrices,
    build_on_site_matrices,
)
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


def export_solver_convergence(backend: TorchBackend) -> None:
    cell = make_honeycomb_2site()
    params_pure = ModelParameters(K=(-1.0, -1.0, -1.0))
    params_field = ModelParameters(K=(-1.0, -1.0, -1.0), h=(0.1, 0.1, 0.1))

    # Run Picard (depth=1), DIIS m=3, DIIS m=6
    conf_picard = SolverConfig(
        k_grid_shape=(12, 12),
        tol_correlator=1e-12,
        tol_gauge=1e-12,
        max_iterations=30,
        diis_history_depth=1,
    )
    conf_diis3 = SolverConfig(
        k_grid_shape=(12, 12),
        tol_correlator=1e-12,
        tol_gauge=1e-12,
        max_iterations=30,
        diis_history_depth=3,
    )
    conf_diis6 = SolverConfig(
        k_grid_shape=(12, 12),
        tol_correlator=1e-12,
        tol_gauge=1e-12,
        max_iterations=30,
        diis_history_depth=6,
    )

    init_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0)
    # add small seed perturbation to test convergence
    rng = np.random.default_rng(42)
    init_state_pert = MeanFieldState(
        U=init_state.U + 0.05 * rng.standard_normal(init_state.U.shape),
        V=init_state.V,
        lambdas=init_state.lambdas,
    )

    _, res_picard = solve_self_consistent_mean_field(
        cell, params_pure, init_state_pert, conf_picard, backend
    )
    _, res_diis3 = solve_self_consistent_mean_field(
        cell, params_pure, init_state_pert, conf_diis3, backend
    )
    _, res_diis6 = solve_self_consistent_mean_field(
        cell, params_pure, init_state_pert, conf_diis6, backend
    )
    _, res_field = solve_self_consistent_mean_field(
        cell, params_field, init_state_pert, conf_diis6, backend
    )

    hist_p = list(res_picard.residual_history)
    hist_d3 = list(res_diis3.residual_history)
    hist_d6 = list(res_diis6.residual_history)
    hist_f = list(res_field.residual_history)

    max_len = max(len(hist_p), len(hist_d3), len(hist_d6), len(hist_f), 15)
    headers = [
        "Description: Fixed-point residual history for Picard vs DIIS acceleration",
        "Columns: iter  res_picard  res_diis3  res_diis6  res_field",
    ]
    rows = []
    for it in range(max_len):
        vp = hist_p[it] if it < len(hist_p) else hist_p[-1]
        vd3 = hist_d3[it] if it < len(hist_d3) else hist_d3[-1]
        vd6 = hist_d6[it] if it < len(hist_d6) else hist_d6[-1]
        vf = hist_f[it] if it < len(hist_f) else hist_f[-1]
        rows.append(f"{it + 1:4d}  {vp:14.6e}  {vd3:14.6e}  {vd6:14.6e}  {vf:14.6e}")

    write_dat_file("solver_convergence.dat", headers, rows)


def export_dispersion_path(backend: TorchBackend) -> None:
    cell = make_honeycomb_2site()
    params = ModelParameters(K=(-1.0, -1.0, -1.0))
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0)

    # Path: Gamma (0,0) -> K (1/3, -1/3) -> M (1/2, 0) -> Gamma (0,0)
    # High-symmetry nodes in fractional reciprocal coordinates
    n_seg = 40
    # Seg 1: Gamma to K
    k_seg1 = [np.array([s / (3 * n_seg), -s / (3 * n_seg)]) for s in range(n_seg)]
    # Seg 2: K to M
    k_seg2 = [
        np.array([1 / 3 + s * (1 / 2 - 1 / 3) / n_seg, -1 / 3 + s * (0 - (-1 / 3)) / n_seg])
        for s in range(n_seg)
    ]
    # Seg 3: M to Gamma
    k_seg3 = [np.array([1 / 2 * (1 - s / n_seg), 0.0]) for s in range(n_seg + 1)]
    k_path = k_seg1 + k_seg2 + k_seg3

    headers = [
        "Description: Pure-Kitaev quasiparticle band dispersion along Gamma-K-M-Gamma",
        "Columns: s_path  k1  k2  e1  e2  e3  e4  e5  e6  e7  e8",
    ]
    rows = []
    for s_idx, k_pt in enumerate(k_path):
        k_tensor = torch.tensor([[k_pt[0], k_pt[1]]], dtype=torch.float64, device=backend.device)
        h_k = build_bloch_hamiltonian(k_tensor, cell, params, state, backend)
        evals, _ = backend.eigh(h_k)
        ev = evals[0].detach().cpu().numpy()
        row_str = f"{s_idx:4d}  {k_pt[0]:8.4f}  {k_pt[1]:8.4f}  " + "  ".join(
            f"{val:10.5f}" for val in ev
        )
        rows.append(row_str)

    write_dat_file("dispersion_path.dat", headers, rows)


def export_parseval_sum_rule(backend: TorchBackend) -> None:
    cell = make_honeycomb_2site()
    params = ModelParameters(
        J=0.08, K=(-1.0, -1.05, -0.95), Gamma=0.12, Gamma_prime=-0.04, D=0.03, h=(0.04, 0.02, 0.05)
    )
    rng = np.random.default_rng(101)
    v_raw = 0.15 * rng.standard_normal((2, 4, 4))
    v_skew = 0.5 * (v_raw - np.swapaxes(v_raw, -1, -2))
    state = MeanFieldState(
        U=0.25 * rng.standard_normal((3, 4, 4)),
        V=v_skew,
        lambdas=0.1 * rng.standard_normal((2, 3)),
    )

    b_tensors = build_on_site_matrices(cell, params, state, backend.device).detach().cpu().numpy()
    a_tensors = (
        build_bond_hopping_matrices(cell, params, state, backend.device).detach().cpu().numpy()
    )
    sum_b_sq = float(np.sum([np.trace(b @ b.T) for b in b_tensors]))
    sum_a_sq = float(np.sum([np.trace(a @ a.T) for a in a_tensors]))
    rhs_frobenius = (1.0 / 16.0) * sum_b_sq + (1.0 / 8.0) * sum_a_sq

    meshes = [4, 6, 8, 12, 16, 24, 32]
    headers = [
        "Description: Parseval-Frobenius spectral sum rule verification across BZ mesh sizes",
        "Columns: mesh_size  lhs_spectral  rhs_frobenius  abs_residual",
    ]
    rows = []
    for m in meshes:
        k_grid = cell.make_k_grid(m, m)
        h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, backend)
        evals, _ = backend.eigh(h_batch)
        ev_np = evals.detach().cpu().numpy()
        lhs_spectral = float(np.mean(np.sum(ev_np**2, axis=-1)))
        diff = abs(lhs_spectral - rhs_frobenius)
        rows.append(f"{m:4d}  {lhs_spectral:14.8f}  {rhs_frobenius:14.8f}  {diff:14.6e}")

    write_dat_file("parseval_sum_rule.dat", headers, rows)


def export_wick_cumulant_check() -> None:
    # 4-Majorana Jordan-Wigner evaluation
    dim = 4
    sigma_x = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
    sigma_y = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
    sigma_z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
    eye = np.eye(2, dtype=np.complex128)

    d1 = np.kron(0.5 * (sigma_x - 1.0j * sigma_y), eye)
    d2 = np.kron(sigma_z, 0.5 * (sigma_x - 1.0j * sigma_y))
    c_ops = [
        d1 + d1.conj().T,
        -1.0j * (d1 - d1.conj().T),
        d2 + d2.conj().T,
        -1.0j * (d2 - d2.conj().T),
    ]

    def measure_state(rho: np.ndarray) -> tuple[float, float, float]:
        omega = np.zeros((4, 4), dtype=np.float64)
        for a in range(4):
            for b in range(4):
                if a != b:
                    omega[a, b] = float(np.imag(np.trace(rho @ c_ops[a] @ c_ops[b])))
        op4 = c_ops[0] @ c_ops[1] @ c_ops[2] @ c_ops[3]
        exp_4 = float(np.real(np.trace(rho @ op4)))
        pf = float(
            omega[0, 1] * omega[2, 3] - omega[0, 2] * omega[1, 3] + omega[0, 3] * omega[1, 2]
        )
        cumulant = abs(exp_4 - (-pf))
        return exp_4, pf, cumulant

    # State 1: Pure Gaussian ground state (vacuum)
    rho_gauss = np.zeros((dim, dim), dtype=np.complex128)
    rho_gauss[0, 0] = 1.0  # vacuum is Gaussian

    # State 2: Singlet state (entangled, non-Gaussian)
    psi_singlet = np.zeros(dim, dtype=np.complex128)
    psi_singlet[1] = 1.0 / np.sqrt(2.0)
    psi_singlet[2] = -1.0 / np.sqrt(2.0)
    rho_singlet = np.outer(psi_singlet, psi_singlet.conj())

    # State 3: Statistical mixture of two Gaussians
    rho_mix = 0.5 * np.diag([1.0, 0.0, 0.0, 0.0]) + 0.5 * np.diag([0.0, 0.0, 0.0, 1.0])

    headers = [
        "Description: Fourth cumulant evaluation kappa_4 = |<c0 c1 c2 c3> - (-Pf(Omega))|",
        "Columns: state_idx  state_label  exp_c0123  pfaffian  kappa_4",
    ]
    rows = []
    cases = [
        (1, "Pure_Gaussian", rho_gauss),
        (2, "Singlet_Entangled", rho_singlet),
        (3, "Gaussian_Mixture", rho_mix),
    ]
    for idx, label, r_mat in cases:
        e4, pf, k4 = measure_state(r_mat)
        rows.append(f"{idx:4d}  {label:20s}  {e4:12.6f}  {pf:12.6f}  {k4:14.6e}")

    write_dat_file("wick_cumulant_check.dat", headers, rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Export MMFT benchmark data files")
    parser.add_argument(
        "--benchmark",
        choices=[
            "all",
            "k006",
            "hellmann-feynman",
            "susceptibility",
            "klein-duality",
            "solver-convergence",
            "dispersion-path",
            "parseval-sum-rule",
            "wick-cumulant",
        ],
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
    if args.benchmark in ("all", "solver-convergence"):
        export_solver_convergence(backend)
    if args.benchmark in ("all", "dispersion-path"):
        export_dispersion_path(backend)
    if args.benchmark in ("all", "parseval-sum-rule"):
        export_parseval_sum_rule(backend)
    if args.benchmark in ("all", "wick-cumulant"):
        export_wick_cumulant_check()


if __name__ == "__main__":
    main()
