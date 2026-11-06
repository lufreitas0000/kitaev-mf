"""Generate dense real-space cluster benchmarks across L in {4, 5, 6, 7, 8, 9, 10, 11, 12, 14}.

Measures:
1. Matrix condition number kappa(A) highlighting the L = 0 mod 3 Dirac resonance.
2. Complex eigensolver (zheevd) wall time and standard deviation over N_repeats = 5.
3. Real skew eigensolver (dsyevd on -A^2) wall time, condition number kappa(S) = kappa(A)^2, and mantissa loss.
4. Newton-Schulz matrix sign purification wall time, throughput (GFLOP/s), and machine-precision accuracy.
Exports to: data/processed/benchmarks/m1_8/real_space_dense_scaling_with_errors.dat
"""

import time
from pathlib import Path

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.purification import purify_newton_schulz
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


def main() -> None:
    out_dir = Path("data/processed/benchmarks/m1_8")
    out_dir.mkdir(parents=True, exist_ok=True)

    h_val = 0.20
    h_vec = (h_val / np.sqrt(3.0), h_val / np.sqrt(3.0), h_val / np.sqrt(3.0))
    params = ModelParameters(K=-1.0, h=h_vec)
    cell2 = make_honeycomb_2site()
    state = make_qsl_seed(cell2)

    sizes = [4, 5, 6, 7, 8, 9, 10, 11, 12, 14]
    n_repeats = 5
    rows = []

    print("=== Generating Dense Real-Space Scaling with Statistical Error Bars ===")

    for L in sizes:
        cluster = RealSpaceTorusCluster(L=L)
        M = cluster.num_majoranas
        n_sites = cluster.num_sites
        A = cluster.assemble_skew_hamiltonian(params, state)

        # Spectral analysis: singular values of A (or eigenvalues of i*A)
        # s_max / s_min gives condition number
        H_c = 1j * A.astype(np.complex128)
        w, _ = np.linalg.eigh(H_c)
        abs_w = np.sort(np.abs(w))
        s_min = float(abs_w[0])
        s_max = float(abs_w[-1])
        kappa_A = s_max / max(s_min, 1e-15)
        is_mod3 = L % 3 == 0

        # 1. Complex eigensolve (zheevd) repeats
        t_c_list = []
        for _ in range(n_repeats):
            t0 = time.perf_counter()
            omega_c, e_c = cluster.solve_covariance_complex_eigh(A)
            t_c_list.append((time.perf_counter() - t0) * 1000.0)
        mean_t_c = float(np.mean(t_c_list))
        std_t_c = float(np.std(t_c_list))

        # 2. Real skew eigensolve (dsyevd on -A^2) repeats
        t_r_list = []
        for _ in range(n_repeats):
            t0 = time.perf_counter()
            omega_r, e_r = cluster.solve_covariance_real_eigh(A)
            t_r_list.append((time.perf_counter() - t0) * 1000.0)
        mean_t_r = float(np.mean(t_r_list))
        std_t_r = float(np.std(t_r_list))

        diff_frob_tr_a = float(np.linalg.norm(omega_c - omega_r)) / np.sqrt(M)

        # 3. Track B: Newton-Schulz Purification repeats
        t_p_list = []
        purif_res = None
        for _ in range(n_repeats):
            t0 = time.perf_counter()
            omega_p, met_p = purify_newton_schulz(A, tol=1e-11, max_iter=30)
            t_p_list.append((time.perf_counter() - t0) * 1000.0)
            purif_res = (omega_p, met_p)
        assert purif_res is not None
        mean_t_p = float(np.mean(t_p_list))
        std_t_p = float(np.std(t_p_list))

        diff_frob_tr_b = float(np.linalg.norm(purif_res[0] - omega_c)) / np.sqrt(M)
        speedup = mean_t_c / max(mean_t_p, 1e-4)

        # Purification throughput: 4 * iters * M^3 FLOPs
        flops_purif = 4.0 * purif_res[1].iterations * (M**3)
        throughput_purif = (flops_purif / (mean_t_p * 1e-3)) / 1e9

        rows.append(
            (
                L,
                n_sites,
                M,
                int(is_mod3),
                s_min,
                kappa_A,
                mean_t_c,
                std_t_c,
                mean_t_r,
                std_t_r,
                diff_frob_tr_a,
                mean_t_p,
                std_t_p,
                speedup,
                purif_res[1].iterations,
                throughput_purif,
                diff_frob_tr_b,
            )
        )
        print(
            f"L={L:2d} (M={M:4d}, mod3={is_mod3}) | s_min={s_min:.2e}, kappa={kappa_A:.2e} | "
            f"zheevd: {mean_t_c:7.2f}±{std_t_c:5.2f} ms | dsyevd: {mean_t_r:7.2f}±{std_t_r:5.2f} ms (frob={diff_frob_tr_a:.2e}) | "
            f"Purif: {mean_t_p:6.2f}±{std_t_p:4.2f} ms ({speedup:5.1f}x, frob={diff_frob_tr_b:.2e})"
        )

    out_file = out_dir / "real_space_dense_scaling_with_errors.dat"
    with open(out_file, "w") as f:
        f.write("# Dense Real-Space Scaling with Statistical Errors (h = 0.20, N_repeats = 5)\n")
        f.write(
            "# Columns: L  N_SITES  M  IS_MOD3  S_MIN  KAPPA_A  TIME_C_MS  STD_C_MS  TIME_R_MS  STD_R_MS  DIFF_TR_A  TIME_P_MS  STD_P_MS  SPEEDUP  ITERS  GFLOPS_P  DIFF_TR_B\n"
        )
        for r in rows:
            f.write(
                f"{r[0]:<3d} {r[1]:<7d} {r[2]:<6d} {r[3]:<8d} {r[4]:<11.4e} {r[5]:<11.4e} {r[6]:<10.2f} {r[7]:<9.2f} {r[8]:<10.2f} {r[9]:<9.2f} {r[10]:<12.3e} {r[11]:<10.2f} {r[12]:<9.2f} {r[13]:<8.2f} {r[14]:<6d} {r[15]:<9.2f} {r[16]:<12.3e}\n"
            )
    print(f"Exported dense real space scaling to {out_file}")


if __name__ == "__main__":
    main()
