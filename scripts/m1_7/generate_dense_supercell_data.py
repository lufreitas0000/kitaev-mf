"""Generate dense supercell scaling benchmarks with intermediate sizes.

Supercell sizes: n_c in {2, 4, 6, 8, 12, 16, 24, 32}.
Measures intensive ground-state energy, residual nematicity, evaluation FLOPs,
and execution time over N_repeats = 5 per configuration.
Exports to: data/processed/benchmarks/m1_7/bloch_supercell_scaling_dense.dat
"""

import time
from pathlib import Path

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_supercell
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.observables import compute_sublattice_magnetizations
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def main() -> None:
    out_dir = Path("data/processed/benchmarks/m1_7")
    out_dir.mkdir(parents=True, exist_ok=True)
    backend = TorchBackend(device="cpu")
    h_val = 0.20
    h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
    params = ModelParameters(K=-1.0, h=h_vec)

    # Supercell geometries: (n1, n2) such that n_c = 2 * n1 * n2
    # {2: (1,1), 4: (2,1), 6: (3,1), 8: (2,2), 12: (3,2), 16: (4,2), 24: (4,3), 32: (4,4)}
    geometries = [
        (1, 1, 2),
        (2, 1, 4),
        (3, 1, 6),
        (2, 2, 8),
        (3, 2, 12),
        (4, 2, 16),
        (4, 3, 24),
        (4, 4, 32),
    ]

    total_k_target = 1152
    config_base = SolverConfig(
        max_iterations=40,
        tol_correlator=1e-9,
        tol_gauge=1e-9,
        method="diis",
        mixing_alpha=0.6,
        skew_mode="retract_step",
    )

    rows = []
    print("=== Generating Dense Bloch Supercell Scaling Benchmark ===")

    e_primitive = None

    for n1, n2, n_c in geometries:
        cell = make_honeycomb_supercell(n1, n2)
        n_k1 = max(4, int(np.round(np.sqrt(total_k_target / n_c))))
        n_k2 = max(4, int(np.round((total_k_target / n_c) / n_k1)))
        n_k = n_k1 * n_k2
        cfg = SolverConfig(
            k_grid_shape=(n_k1, n_k2),
            max_iterations=config_base.max_iterations,
            tol_correlator=config_base.tol_correlator,
            tol_gauge=config_base.tol_gauge,
            method=config_base.method,
            mixing_alpha=config_base.mixing_alpha,
            skew_mode=config_base.skew_mode,
        )

        M = 4 * cell.num_sublattices
        dim_x = 16 * len(cell.bonds) + 16 * cell.num_sublattices + 3 * cell.num_sublattices

        times = []
        n_repeats = 5
        final_st = None
        final_met = None

        for _rep in range(n_repeats):
            seed = make_qsl_seed(cell)
            t0 = time.perf_counter()
            st, met = solve_self_consistent_mean_field(cell, params, seed, cfg, backend)
            t_ms = (time.perf_counter() - t0) * 1000.0
            times.append(t_ms)
            final_st = st
            final_met = met

        assert final_st is not None and final_met is not None
        mean_time = float(np.mean(times))
        std_time = float(np.std(times))
        e_breakdown = compute_gaussian_energy(final_st, cell, params)
        e_site = e_breakdown.E_total / cell.num_sublattices
        if e_primitive is None:
            e_primitive = e_site
        e_diff = abs(e_site - e_primitive)

        # Nematicity
        _, spin_mag = compute_sublattice_magnetizations(final_st, cell)
        nematicities = [
            abs(spin_mag[s, 0] - spin_mag[s, 1]) + abs(spin_mag[s, 1] - spin_mag[s, 2])
            for s in range(cell.num_sublattices)
        ]
        max_nem = float(np.max(nematicities))

        # Flops per evaluation: 16/3 * M^3 * n_k
        flops_eval = (16.0 / 3.0) * (M**3) * n_k
        mflops_eval = flops_eval / 1e6

        time_per_eval = mean_time / final_met.map_evaluations
        std_eval = std_time / final_met.map_evaluations
        throughput_gflops = (flops_eval / (time_per_eval * 1e-3)) / 1e9

        rows.append(
            (
                n_c,
                M,
                dim_x,
                n_k,
                final_met.map_evaluations,
                mean_time,
                std_time,
                time_per_eval,
                std_eval,
                mflops_eval,
                throughput_gflops,
                e_site,
                e_diff,
                max_nem,
            )
        )
        print(
            f"n_c={n_c:2d} (M={M:3d}, Nk={n_k:3d}) | Time={mean_time:7.2f} ± {std_time:5.2f} ms | Eval={time_per_eval:6.2f} ms | GFLOP/s={throughput_gflops:5.2f} | dE={e_diff:.2e} | Nem={max_nem:.2e}"
        )

    out_file = out_dir / "bloch_supercell_scaling_dense.dat"
    with open(out_file, "w") as f:
        f.write("# Dense Bloch Supercell Scaling Benchmark (h = 0.20, N_repeats = 5)\n")
        f.write(
            "# Columns: N_C  M_DIM  X_DIM  N_K  ITERS  TIME_MS  TIME_STD_MS  MS_PER_EVAL  MS_EVAL_STD  MFLOPS_EVAL  GFLOPS  E_SITE  E_DIFF  NEMATICITY\n"
        )
        for r in rows:
            f.write(
                f"{r[0]:<3d} {r[1]:<5d} {r[2]:<6d} {r[3]:<4d} {r[4]:<5d} {r[5]:<9.2f} {r[6]:<9.2f} {r[7]:<11.3f} {r[8]:<11.3f} {r[9]:<12.3f} {r[10]:<7.2f} {r[11]:<14.10f} {r[12]:<11.3e} {r[13]:<11.3e}\n"
            )
    print(f"Exported dense supercell scaling to {out_file}")


if __name__ == "__main__":
    main()
