"""ACID benchmark exporter for Wave-2B (UOW-03 / HC-02) Gutzwiller SU(2)_R & flux stability."""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import subprocess
from math import sqrt
from pathlib import Path

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.services.gutzwiller_projection import (
    _exact_spin_observables,
    assemble_torus_majorana_matrix,
    build_thouless_data,
    exact_wavefunction_from_thouless,
)
from kitaev_mf.services.gutzwiller_stability import (
    FluxSectorName,
    evaluate_flux_sector_gutzwiller,
    evaluate_su2_gauge_stability,
    make_staggered_su2_patterns,
)
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def export_wave2b_gutzwiller_stability_benchmark(
    repo_root: Path,
    *,
    n_chains: int = 50,
    burn_in_sweeps: int = 15,
    samples_per_chain: int = 50,
) -> tuple[Path, Path]:
    """Generate ACID-tracked Wave-2B HC-02 Gutzwiller SU(2)_R & flux-sector stability benchmark."""
    cell = make_honeycomb_2site()
    cfg = SolverConfig(k_grid_shape=(24, 24), tol_correlator=1e-10, tol_gauge=1e-10)
    st0, met0 = solve_self_consistent_mean_field(
        cell, ModelParameters(K=-1.0), make_qsl_seed(cell), cfg
    )
    if not met0.converged:
        raise RuntimeError("h=0 scMFT solve failed to converge.")
    h25 = 0.25
    c25 = h25 / sqrt(3.0)
    st25, met25 = solve_self_consistent_mean_field(
        cell, ModelParameters(K=-1.0, h=(c25, c25, c25)), st0, cfg
    )
    if not met25.converged:
        raise RuntimeError("h=0.25 scMFT solve failed to converge.")

    sha = subprocess.check_output(
        ["git", "-C", str(repo_root), "rev-parse", "HEAD"], text=True
    ).strip()
    out_dir = repo_root / "data" / "processed" / "benchmarks" / "wave2"
    dat_path = out_dir / "hc02_gutzwiller_su2_flux_stability.dat"
    meta_path = dat_path.with_suffix(".meta.json")

    lines = [
        f"# Benchmark: {dat_path.name}",
        f"# Generated: {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
        f"# Git Commit: {sha[:7]} ({sha})",
        "# Generator: kitaev_mf.services.gutzwiller_benchmark.export_wave2b_gutzwiller_stability_benchmark",
        "# Description: Wave-2B (UOW-03 / HC-02) Gutzwiller D_j=+1 VMC stability across Z_2 flux sectors and staggered SU(2)_R gauges",
        "# Columns: kind  N_sites  nx  ny  h  label  bc_x  bc_y  qp_exc  cond  rot_norm  exact_fid  probe_eloc_diff  e_exact  e_vmc  ci_low  ci_high  e_loc_std  wp_exact  g_sokal  tau_wolff  g_ips  n_eff_ips  split_rhat  acc1  acc2",
        "#" + "-" * 130,
    ]

    flux_cases: list[tuple[int, int, FluxSectorName]] = [
        (2, 2, "zero_flux"),
        (2, 2, "two_vison"),
        (2, 2, "pi_flux"),
        (2, 3, "zero_flux"),
        (2, 3, "two_vison"),
        (2, 3, "pi_flux"),
        (3, 3, "zero_flux"),
        (3, 3, "two_vison"),
    ]
    for idx_f, (nx, ny, sector) in enumerate(flux_cases):
        res_f = evaluate_flux_sector_gutzwiller(
            nx,
            ny,
            0.0,
            st0,
            sector,
            n_chains=n_chains,
            burn_in_sweeps=burn_in_sweeps,
            samples_per_chain=samples_per_chain,
            seed_base=4200 + 100 * idx_f,
        )
        e_ex_s = (
            f"{res_f.e_proj_exact:15.10f}" if res_f.e_proj_exact is not None else f"{'nan':>15s}"
        )
        wp_s = (
            f"{res_f.plaquette_flux_exact:+9.5f}"
            if res_f.plaquette_flux_exact is not None
            else f"{'nan':>9s}"
        )
        d = res_f.diagnostics
        lines.append(
            f"{'flux':6s} {res_f.n_sites:3d} {nx:2d} {ny:2d} {0.0:6.3f} {sector:20s} "
            f"{res_f.bc_twist[0]:3d} {res_f.bc_twist[1]:3d} {res_f.n_excited_qp:2d} "
            f"{res_f.ref_condition_number:8.3f} {0.0:9.4f} {1.0:10.7f} {0.0:10.3e} {e_ex_s} "
            f"{res_f.vmc_e_mean:15.10f} {res_f.vmc_e_ci_low:15.10f} {res_f.vmc_e_ci_high:15.10f} "
            f"{res_f.vmc_e_loc_std:10.3e} {wp_s} {d.g_sokal:7.3f} {d.tau_int_wolff:7.3f} "
            f"{d.g_ips_geyer:7.3f} {d.n_eff_ips:8.1f} {d.split_rhat:7.4f} "
            f"{d.acc_rate_1flip:6.4f} {d.acc_rate_2flip:6.4f}"
        )

    su2_cases: list[tuple[int, int, tuple[str, ...]]] = [
        (
            2,
            2,
            (
                "identity",
                "sublattice_staggered",
                "pi_flips_xyz",
                "z2_2pi_staggered",
                "random_inhomogeneous",
            ),
        ),
        (2, 3, ("identity", "sublattice_staggered", "pi_flips_xyz", "random_inhomogeneous")),
        (3, 3, ("identity", "sublattice_staggered", "pi_flips_xyz")),
    ]
    for idx_g, (nx, ny, pat_names) in enumerate(su2_cases):
        n_sites = 2 * nx * ny
        all_pats = make_staggered_su2_patterns(n_sites, seed=2026 + idx_g)
        sub_pats = {k: all_pats[k] for k in pat_names}
        rows_g = evaluate_su2_gauge_stability(
            nx,
            ny,
            h25,
            st25,
            patterns=sub_pats,
            n_chains=n_chains,
            burn_in_sweeps=burn_in_sweeps,
            samples_per_chain=samples_per_chain,
            seed_base=5100 + 200 * idx_g,
            independent_vmc_seeds=True,
        )
        base_obs_e: float | None = None
        base_obs_wp: float | None = None
        if n_sites <= 12:
            a_b = assemble_torus_majorana_matrix(
                nx,
                ny,
                ModelParameters(K=-1.0, h=(c25, c25, c25)),
                st25,
                bc_twist=rows_g[0].bc_twist,
            )
            td_b = build_thouless_data(a_b, n_sites)
            assert td_b is not None
            psi_b = exact_wavefunction_from_thouless(td_b[0], td_b[1])
            assert psi_b is not None
            base_obs_e, _, base_obs_wp, _ = _exact_spin_observables(nx, ny, h25, psi_b)

        for r in rows_g:
            fid_s = f"{r.exact_fidelity:10.7f}" if r.exact_fidelity is not None else f"{'nan':>10s}"
            e_ex_s = f"{base_obs_e:15.10f}" if base_obs_e is not None else f"{'nan':>15s}"
            wp_s = f"{base_obs_wp:+9.5f}" if base_obs_wp is not None else f"{'nan':>9s}"
            d = r.diagnostics
            lines.append(
                f"{'su2_r':6s} {r.n_sites:3d} {nx:2d} {ny:2d} {h25:6.3f} {r.pattern_name:20s} "
                f"{r.bc_twist[0]:3d} {r.bc_twist[1]:3d} {0:2d} "
                f"{r.cond_rot:8.3f} {r.matrix_rotation_norm:9.4f} {fid_s} "
                f"{r.max_probe_eloc_diff:10.3e} {e_ex_s} "
                f"{r.vmc_e_mean:15.10f} {r.vmc_e_ci_low:15.10f} {r.vmc_e_ci_high:15.10f} "
                f"{r.vmc_e_loc_std:10.3e} {wp_s} {d.g_sokal:7.3f} {d.tau_int_wolff:7.3f} "
                f"{d.g_ips_geyer:7.3f} {d.n_eff_ips:8.1f} {d.split_rhat:7.4f} "
                f"{d.acc_rate_1flip:6.4f} {d.acc_rate_2flip:6.4f}"
            )

    _atomic_write(dat_path, "\n".join(lines) + "\n")
    meta = {
        "schema_version": 1,
        "dataset": "hc02_gutzwiller_su2_flux_stability",
        "source_repo": "kitaev-mf",
        "git_sha": sha,
        "sha256": hashlib.sha256(dat_path.read_bytes()).hexdigest(),
        "unit_convention": "pauli_k_minus_one_field_111",
        "geometry": "Honeycomb tori 2x2 (N=8), 2x3 (N=12), 3x3 (N=18)",
        "parameters": {
            "uow_id": "UOW-03",
            "vmc_chains": n_chains,
            "burn_in_sweeps": burn_in_sweeps,
            "samples_per_chain": samples_per_chain,
            "bootstrap_method": "BCa 95% CI (B=2000)",
            "autocorrelation_estimators": [
                "Sokal_g",
                "Wolff_tau_int=g/2",
                "Geyer_IPS_grand_mean",
                "Split_Rhat",
            ],
            "flux_sectors_h0": ["zero_flux", "two_vison", "pi_flux"],
            "su2_r_patterns_h025": list(su2_cases[0][2]),
        },
    }
    _atomic_write(meta_path, json.dumps(meta, indent=2, sort_keys=True) + "\n")
    return dat_path, meta_path


__all__ = ["export_wave2b_gutzwiller_stability_benchmark"]
