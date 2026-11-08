#!/usr/bin/env python3
"""Run T=0 real-space observable atlas on honeycomb torus clusters.

Phase 2, Milestone M2.1, Task T2.1.1.
Uses Track B Newton-Schulz purification with damped L-Broyden accelerator (m=10, alpha=0.5).
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


def run_atlas(
    L: int = 12,
    h_min: float = 0.0,
    h_max: float = 0.50,
    num_h: int = 26,
    max_iter: int = 65,
    tol: float = 1e-6,
    output_path: Path | None = None,
) -> None:
    print(
        f"=== Running Real-Space Observable Atlas: L={L} (N={2 * L * L} sites, M={8 * L * L} Majoranas) ==="
    )
    cluster = RealSpaceTorusCluster(L=L)
    cell2 = make_honeycomb_2site()
    seed_state = make_qsl_seed(cell2)

    h_values = np.linspace(h_min, h_max, num_h)
    records = []

    curr_state = seed_state
    for i, h_mag in enumerate(h_values):
        h_vec = (h_mag / np.sqrt(3.0), h_mag / np.sqrt(3.0), h_mag / np.sqrt(3.0))
        params = ModelParameters(K=-1.0, h=h_vec)

        t0 = time.perf_counter()
        omega, final_state, met = cluster.solve_real_space_purification(
            params=params,
            initial_state=curr_state,
            max_iter=max_iter,
            tol=tol,
            accelerator="lbroyden",
            mixing_alpha=0.5,
            lbroyden_history_size=10,
        )
        elapsed = time.perf_counter() - t0

        conv = met["converged"]
        iters = met["iterations"]
        # Exact variational Hartree-Fock ground state energy per site (including Wick subtraction)
        breakdown = compute_gaussian_energy(final_state, cell2, params)
        e_site = breakdown.per_site

        # Observables from converged state
        # V has shape (2, 4, 4) for A and B sublattices
        # spin m_gamma = V_{0, gamma}
        m_A = np.array([final_state.V[0, 0, 1], final_state.V[0, 0, 2], final_state.V[0, 0, 3]])
        m_B = np.array([final_state.V[1, 0, 1], final_state.V[1, 0, 2], final_state.V[1, 0, 3]])
        m_avg = 0.5 * (m_A + m_B)
        m_111 = float(np.sum(m_avg) / np.sqrt(3.0))

        # Matter and gauge bond correlators
        # U has shape (3, 4, 4)
        # In make_honeycomb_2site():
        #   bond 0: z-bond (bond_type = 3) -> gauge component is gamma = 3
        #   bond 1: x-bond (bond_type = 1) -> gauge component is gamma = 1
        #   bond 2: y-bond (bond_type = 2) -> gauge component is gamma = 2
        u_00 = float(np.mean([final_state.U[b, 0, 0] for b in range(3)]))
        u_z = final_state.U[0, 3, 3]
        u_x = final_state.U[1, 1, 1]
        u_y = final_state.U[2, 2, 2]
        u_gg = float(np.mean([u_x, u_y, u_z]))

        # Plaquette flux W_p: product of 6 bond correlators u_b = U_b^{\gamma\gamma}
        # For symmetric 2-site unit cell, W_p = (u_x * u_y * u_z)^2
        w_p = float((u_x * u_y * u_z) ** 2)

        # Lagrange multiplier norm
        lambda_norm = float(np.linalg.norm(final_state.lambdas))

        print(
            f"[{i + 1:02d}/{num_h}] h={h_mag:.4f} | E/site={e_site:.6f} | "
            f"m_111={m_111:.5f} | W_p={w_p:.5f} | ||lambda||={lambda_norm:.4f} | "
            f"conv={conv} (iter={iters}, {elapsed:.2f}s)"
        )

        records.append(
            {
                "h": h_mag,
                "e_site": e_site,
                "m_111": m_111,
                "m_x": m_avg[0],
                "m_y": m_avg[1],
                "m_z": m_avg[2],
                "u_00": u_00,
                "u_gg": u_gg,
                "w_p": w_p,
                "lambda_norm": lambda_norm,
                "conv": int(conv),
                "iters": iters,
            }
        )

        # Continuous continuation along the branch
        curr_state = final_state

    # Compute numerical susceptibility chi_111 = d m_111 / d h
    for j in range(len(records)):
        if j == 0:
            dh = records[1]["h"] - records[0]["h"]
            chi = (records[1]["m_111"] - records[0]["m_111"]) / dh
        elif j == len(records) - 1:
            dh = records[-1]["h"] - records[-2]["h"]
            chi = (records[-1]["m_111"] - records[-2]["m_111"]) / dh
        else:
            dh = records[j + 1]["h"] - records[j - 1]["h"]
            chi = (records[j + 1]["m_111"] - records[j - 1]["m_111"]) / dh
        records[j]["chi_111"] = float(chi)

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("# Real-space observable atlas for pure FM Kitaev (K=-1.0) under h || [111]\n")
            f.write(f"# Lattice: Honeycomb Torus L={L} (N={2 * L * L} sites)\n")
            f.write(
                "# Solver: Track B Newton-Schulz Purification + Damped L-Broyden (m=10, alpha=0.5)\n"
            )
            f.write("# Generated on: 2026-11-07\n")
            f.write(
                "# h\te_site\tm_111\tchi_111\tm_x\tm_y\tm_z\tu_00\tu_gg\tw_p\tlambda_norm\tconv\titers\n"
            )
            for r in records:
                f.write(
                    f"{r['h']:.6e}\t{r['e_site']:.8e}\t{r['m_111']:.8e}\t{r['chi_111']:.8e}\t"
                    f"{r['m_x']:.8e}\t{r['m_y']:.8e}\t{r['m_z']:.8e}\t"
                    f"{r['u_00']:.8e}\t{r['u_gg']:.8e}\t{r['w_p']:.8e}\t{r['lambda_norm']:.8e}\t"
                    f"{r['conv']}\t{r['iters']}\n"
                )
        print(f"\n[Atlas] Saved {len(records)} points to {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Real-space observable atlas for FM Kitaev [111]")
    parser.add_argument("--L", type=int, default=12, help="Linear cluster dimension")
    parser.add_argument("--h-min", type=float, default=0.0, help="Minimum field")
    parser.add_argument("--h-max", type=float, default=0.50, help="Maximum field")
    parser.add_argument("--num-h", type=int, default=26, help="Number of field points")
    parser.add_argument(
        "--output",
        type=str,
        default="data/processed/benchmarks/m2_1/real_space_observable_atlas.dat",
    )
    args = parser.parse_args()

    out_p = Path(args.output)
    run_atlas(
        L=args.L,
        h_min=args.h_min,
        h_max=args.h_max,
        num_h=args.num_h,
        output_path=out_p,
    )


if __name__ == "__main__":
    main()
