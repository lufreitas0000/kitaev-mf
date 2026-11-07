#!/usr/bin/env python3
"""Benchmark Script for Task T3.W1.3: Thermodynamic Spectral Sum Rules up to L=37.

Exports ACID-compliant benchmark table with metadata provenance:
- Table: data/processed/benchmarks/phase_3/mzm_thermodynamic_spectral_benchmarks.dat
"""

import datetime
import os
import subprocess
from pathlib import Path

from kitaev_mf.services.phase_3.thermodynamic_spectral_benchmarks import (
    benchmark_thermodynamic_spectral_response,
)


def get_git_revision() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii").strip()
    except Exception:
        return "unknown"


def main() -> None:
    git_rev = get_git_revision()
    now_iso = datetime.datetime.now(datetime.timezone.utc).isoformat()
    out_dir = Path(__file__).resolve().parents[2] / "data/processed/benchmarks/phase_3"
    os.makedirs(out_dir, exist_ok=True)

    print("Benchmarking thermodynamic spectral response up to L=37...")
    L_values = [6, 12, 18, 24, 30, 36, 37]
    data = benchmark_thermodynamic_spectral_response(L_values=L_values)

    path = os.path.join(out_dir, "mzm_thermodynamic_spectral_benchmarks.dat")
    with open(path, "w") as f:
        f.write("# ACID DATA PROVENANCE HEADER\n")
        f.write("# Task: T3.W1.3 Thermodynamic Spectral Sum Rules and Demarcation Gap\n")
        f.write(f"# Date: {now_iso}\n")
        f.write(f"# Git-Commit: {git_rev}\n")
        f.write("# Parameters: K=1.0, kappa=0.15, boundary=obc\n")
        f.write(f"# Evaluated L values: {L_values}\n")
        f.write(
            "# Schema: L n_sites is_sing_core is_sing_bulk core_sum_rule bulk_sum_rule omega_MZM bulk_min_energy demarcation_gap quad_err\n"
        )
        for d in data:
            f.write(
                f"{d.L:4d} {d.n_sites:7d} {str(d.is_singular_core):>12s} {str(d.is_singular_bulk):>12s} "
                f"{d.core_sum_rule:14.6f} {d.bulk_sum_rule:14.6f} {d.omega_MZM:10.6f} "
                f"{d.bulk_min_energy:15.6f} {d.demarcation_gap:15.6f} {d.quadrature_error:12.4e}\n"
            )
    print(f"Exported {path}")


if __name__ == "__main__":
    main()
