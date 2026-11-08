#!/usr/bin/env python3
# ruff: noqa: E402
# Repository-path bootstrap must precede local imports in this standalone entrypoint.
"""Benchmark and multi-theory comparison script for Phase 4:
Compares 4 theories for the Kitaev FM [111] model on the 3x3 torus (N=18):
1. Exact Diagonalization (ED) in the physical spin space
2. Perturbative Kitaev (2006) 3rd order asymptotic predictions
3. Semi-Analytical Mean-Field (Schrieffer-Wolff algebraic prefactor C_MF)
4. Numerical Self-Consistent Mean-Field (fully relaxed D_sym=5 DIIS solver)

Outputs ACID benchmark table to data/processed/benchmarks/phase_4/four_theory_observables_comparison.dat
"""

import sys
import time
from pathlib import Path

# Ensure workspace root is in python path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from kitaev_mf.services.phase_4.comprehensive_observables import (
    compute_comprehensive_ed,
    compute_comprehensive_mf,
)
from kitaev_mf.services.phase_4.exact_diagonalization import HoneycombCluster


def main() -> None:
    print("=" * 80)
    print("PHASE 4: FOUR-THEORY COMPREHENSIVE OBSERVABLES BENCHMARK (N=18 TORUS)")
    print("=" * 80)

    # 3x3 Torus (N=18)
    cluster = HoneycombCluster.make_torus(3, 3)
    print(f"Cluster: {cluster.n1}x{cluster.n2} torus ({cluster.num_spins} spins)")
    print(f"Total Hilbert space dimension: 2^{cluster.num_spins} = {1 << cluster.num_spins}")

    h_values = [0.000, 0.010, 0.020, 0.030, 0.040, 0.050]
    out_dir = REPO_ROOT / "data" / "processed" / "benchmarks" / "phase_4"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / "four_theory_observables_comparison.dat"

    # Perturbative Kitaev constants (canonical units: K=1.0)
    # E_0/N = -0.78729862
    # Delta_Kitaev = 28.6806 * h^3
    # chi_m = 0.1884
    c_kitaev = 28.6806
    c_mf_semi = 4.0 * 2.6634  # 10.6536

    rows = []

    print(f"\nEvaluating {len(h_values)} magnetic field points across 4 frameworks...")
    t0_all = time.time()

    prev_ed_state = None

    for h in h_values:
        t0 = time.time()
        print(f"\n--- Field h = {h:.4f} ---", flush=True)

        # 1. Exact Diagonalization
        ed, prev_ed_state = compute_comprehensive_ed(cluster, h=h, K=1.0, prev_state=prev_ed_state)
        print(
            f"  [ED]  E/N = {ed.energy_per_spin:.6f}, Gap = {ed.gap_e1_e0:.6e}, "
            f"m111 = {ed.magnetization_111:.6f}, <W_p> = {ed.plaquette_flux:.6f}, "
            f"S_vN = {ed.entanglement_entropy_bipartite:.4f}, chi_F = {ed.fidelity_susceptibility:.4e}",
            flush=True,
        )

        # 2. Perturbative Kitaev
        e_per_spin_kit = -0.78729862 - 0.5 * 0.1884 * (h**2)
        gap_kit = c_kitaev * (h**3)
        m111_kit = 0.1884 * h

        # 3. Semi-Analytical MF
        gap_semi = c_mf_semi * (h**3)

        # 4. Numerical Self-Consistent MF
        mf = compute_comprehensive_mf(h=h, K=1.0, nk=24)
        print(
            f"  [MF]  E/N = {mf.energy_per_spin:.6f}, Gap = {mf.dirac_gap_canonical:.6e}, "
            f"m111 = {mf.magnetization_111:.6f}, mu = {mf.contraction_factor_mu:.4f}, "
            f"Gamma = {mf.relaxation_rate_gamma:.4f}, Iters = {mf.iterations_to_converge}",
            flush=True,
        )

        dt = time.time() - t0
        print(f"  Point completed in {dt:.2f} s", flush=True)

        rows.append(
            (
                h,
                ed.energy_per_spin,
                mf.energy_per_spin,
                e_per_spin_kit,
                ed.gap_e1_e0,
                gap_kit,
                gap_semi,
                mf.dirac_gap_canonical,
                ed.magnetization_111,
                mf.magnetization_111,
                m111_kit,
                ed.susceptibility_chi_m,
                mf.susceptibility_chi_m,
                ed.corr_nn_parallel,
                ed.corr_nn_transverse,
                ed.three_spin_chiral,
                ed.plaquette_flux,
                ed.entanglement_entropy_bipartite,
                ed.fidelity_susceptibility,
                mf.contraction_factor_mu,
                mf.relaxation_rate_gamma,
            )
        )

    total_time = time.time() - t0_all
    print(f"\nAll calculations completed in {total_time:.2f} s.")

    # Write ACID dat file
    with open(out_file, "w") as f:
        f.write("# ACID BENCHMARK: PHASE 4 FOUR-THEORY COMPREHENSIVE COMPARISON\n")
        f.write("# Model: Kitaev Ferromagnetic [111] (K=1.0, Pauli sigma in {-1, 1})\n")
        f.write("# Geometry: 3x3 Torus (N=18 spins, PBC, contains Dirac K/K' points)\n")
        f.write("# Theories compared:\n")
        f.write("#   1. Exact Diagonalization (Spin Space)\n")
        f.write("#   2. Perturbative Kitaev (2006) 3rd order asymptotic (Delta = 28.68 h^3)\n")
        f.write("#   3. Semi-Analytical MF (Schrieffer-Wolff resolvent, Delta = 10.6536 h^3)\n")
        f.write("#   4. Numerical Self-Consistent MF (symmetrized D_sym=5 DIIS solver)\n")
        f.write("# Columns:\n")
        f.write(
            "# 1:h 2:E_ED 3:E_MF 4:E_Kit 5:Gap_ED 6:Gap_Kit 7:Gap_Semi 8:Gap_MF "
            "9:m111_ED 10:m111_MF 11:m111_Kit 12:chi_ED 13:chi_MF "
            "14:C_nn_par 15:C_nn_perp 16:C_chi 17:W_p 18:S_vN 19:chi_F 20:mu_MF 21:Gamma_MF\n"
        )
        for r in rows:
            f.write(" ".join(f"{x:16.8e}" for x in r) + "\n")

    print(f"Benchmark table written to: {out_file}")


if __name__ == "__main__":
    main()
