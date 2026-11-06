"""Benchmark generalized Bloch supercells n_c in {2, 4, 8, 16, 32}.

Measures:
1. Scaling of evaluation FLOPs and wall time (O(N_k * n_c^3))
2. Ground-state intensive energy parity across folded Brillouin zones
3. Unconstrained random seed stability and absence of spontaneous symmetry breaking
"""

import time
from pathlib import Path

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site, make_honeycomb_supercell
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.observables import compute_sublattice_magnetizations
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def run_supercell_benchmarks() -> None:
    backend = TorchBackend(device="cpu")
    out_dir = Path("data/processed/benchmarks/m1_7")
    out_dir.mkdir(parents=True, exist_ok=True)

    h_val = 0.20
    params = ModelParameters(K=-1.0, h=(h_val, h_val, h_val))

    supercell_configs = [
        (1, 1, "2site_1x1"),
        (2, 1, "4site_2x1"),
        (2, 2, "8site_2x2"),
        (4, 2, "16site_4x2"),
        (4, 4, "32site_4x4"),
    ]

    print("=== 1. Benchmarking Bloch Supercell Scaling and BZ Folding Parity ===")
    scaling_data = []

    # Reference 2-site energy on a 24x24 mesh
    cell2 = make_honeycomb_2site()
    seed2 = make_qsl_seed(cell2)
    cfg2 = SolverConfig(
        max_iterations=60, tol_correlator=1e-8, tol_gauge=1e-8, k_grid_shape=(24, 24), method="diis"
    )
    st2, met2 = solve_self_consistent_mean_field(cell2, params, seed2, cfg2, backend)
    e_ref = compute_gaussian_energy(st2, cell2, params).per_site
    print(f"Reference 2-site energy (24x24): {e_ref:.10f}")

    for n1, n2, name in supercell_configs:
        cell = make_honeycomb_supercell(n1, n2)
        nc = cell.num_sublattices
        M = 4 * nc
        dim_x = 16 * len(cell.bonds) + 16 * nc + 3 * nc

        # Folded k-mesh matching 24x24 primitive sampling: (24/n1, 24/n2)
        k1 = 24 // n1
        k2 = 24 // n2
        N_k = k1 * k2

        seed = make_qsl_seed(cell)
        cfg = SolverConfig(
            max_iterations=60,
            tol_correlator=1e-8,
            tol_gauge=1e-8,
            k_grid_shape=(k1, k2),
            method="diis",
        )

        t0 = time.perf_counter()
        st, met = solve_self_consistent_mean_field(cell, params, seed, cfg, backend)
        wall_time_ms = (time.perf_counter() - t0) * 1000.0

        e_site = compute_gaussian_energy(st, cell, params).per_site
        e_err = abs(e_site - e_ref)

        # Flops per map evaluation: N_k * (16/3 * M^3) for complex Hermitian eigh
        flops_per_eval = (16.0 / 3.0) * (M**3) * N_k
        mflops_per_eval = flops_per_eval / 1e6

        time_per_eval_ms = wall_time_ms / max(met.iterations, 1)

        # Polarization nematicity
        sigma, _ = compute_sublattice_magnetizations(st, cell)
        nematicity = float(
            np.mean(np.abs(sigma[:, 0] - sigma[:, 1]) + np.abs(sigma[:, 1] - sigma[:, 2]))
        )

        scaling_data.append(
            (
                nc,
                M,
                dim_x,
                N_k,
                met.iterations,
                wall_time_ms,
                time_per_eval_ms,
                mflops_per_eval,
                e_site,
                e_err,
                nematicity,
            )
        )
        print(
            f"[{name}] nc={nc:2d} | M={M:3d} | D={dim_x:4d} | N_k={N_k:3d} | K={met.iterations:2d} | Time={wall_time_ms:6.1f}ms | per_eval={time_per_eval_ms:5.2f}ms | MFLOP/eval={mflops_per_eval:6.2f} | diff={e_err:.2e} | nem={nematicity:.2e}"
        )

    # Export scaling table
    scaling_file = out_dir / "bloch_supercell_scaling.dat"
    with open(scaling_file, "w") as f:
        f.write("# Bloch Supercell Scaling and BZ Folding Parity (h = 0.20)\n")
        f.write(
            "# Columns: N_C  M_DIM  X_DIM  N_K  ITERS  WALL_MS  MS_PER_EVAL  MFLOPS_PER_EVAL  E_SITE  E_DIFF  NEMATICITY\n"
        )
        for row in scaling_data:
            f.write(
                f"{row[0]:<4d} {row[1]:<6d} {row[2]:<6d} {row[3]:<5d} {row[4]:<6d} {row[5]:<9.2f} {row[6]:<12.3f} {row[7]:<16.3f} {row[8]:<14.10f} {row[9]:<12.4e} {row[10]:<12.4e}\n"
            )
    print(f"Exported {scaling_file}")

    print("\n=== 2. Testing Unconstrained Random Seed Symmetry Stability ===")
    # For nc in {2, 4, 8, 16}, run N_trials = 10 unconstrained random seeds
    stability_data: list[tuple[int, int, float, float, float, float]] = []
    rng = np.random.default_rng(2026)

    for n1, n2, name in [(1, 1, "2site"), (2, 1, "4site"), (2, 2, "8site"), (4, 2, "16site")]:
        cell = make_honeycomb_supercell(n1, n2)
        nc = cell.num_sublattices
        k1, k2 = 12 // n1, 12 // n2
        N_trials = 10
        converged_count = 0
        isotropic_count = 0
        max_nem = 0.0
        max_e_diff = 0.0

        # Unperturbed reference energy on this mesh
        seed_ref = make_qsl_seed(cell)
        cfg_ref = SolverConfig(
            max_iterations=60,
            tol_correlator=1e-7,
            tol_gauge=1e-7,
            k_grid_shape=(k1, k2),
            method="diis",
        )
        st_ref, _ = solve_self_consistent_mean_field(cell, params, seed_ref, cfg_ref, backend)
        e0 = compute_gaussian_energy(st_ref, cell, params).per_site

        for _ in range(N_trials):
            # Unconstrained random seed: perturb all bond correlators independently
            pert_seed = make_qsl_seed(cell)
            noise_u = rng.standard_normal(pert_seed.U.shape) * 0.05
            u_pert = pert_seed.U + noise_u
            # Antisymmetrize off-diagonal if needed or clip
            state_pert = MeanFieldState(u_pert, pert_seed.V, pert_seed.lambdas)

            cfg = SolverConfig(
                max_iterations=70,
                tol_correlator=1e-7,
                tol_gauge=1e-7,
                k_grid_shape=(k1, k2),
                method="diis",
            )
            st_out, met_out = solve_self_consistent_mean_field(
                cell, params, state_pert, cfg, backend
            )

            if met_out.converged:
                converged_count += 1
                e_run = compute_gaussian_energy(st_out, cell, params).per_site
                diff_e = abs(e_run - e0)
                sigma, _ = compute_sublattice_magnetizations(st_out, cell)
                nem = float(
                    np.mean(np.abs(sigma[:, 0] - sigma[:, 1]) + np.abs(sigma[:, 1] - sigma[:, 2]))
                )

                max_nem = max(max_nem, nem)
                max_e_diff = max(max_e_diff, diff_e)
                if nem < 1e-3 and diff_e < 1e-4:
                    isotropic_count += 1

        p_conv = converged_count / N_trials
        p_iso = isotropic_count / max(converged_count, 1)
        stability_data.append((nc, N_trials, p_conv, p_iso, max_nem, max_e_diff))
        print(
            f"[{name}] nc={nc:2d} | P_conv={p_conv:.2f} | P_iso={p_iso:.2f} | max_nem={max_nem:.2e} | max_dE={max_e_diff:.2e}"
        )

    stability_file = out_dir / "supercell_symmetry_stability.dat"
    with open(stability_file, "w") as f:
        f.write("# Supercell Symmetry Stability Under Unconstrained Random Seeds (h = 0.20)\n")
        f.write("# Columns: N_C  N_TRIALS  P_CONV  P_ISOTROPIC  MAX_NEMATICITY  MAX_ENERGY_DIFF\n")
        for stab_row in stability_data:
            f.write(
                f"{stab_row[0]:<4d} {stab_row[1]:<8d} {stab_row[2]:<8.3f} {stab_row[3]:<12.3f} {stab_row[4]:<16.4e} {stab_row[5]:<16.4e}\n"
            )
    print(f"Exported {stability_file}")


if __name__ == "__main__":
    run_supercell_benchmarks()
