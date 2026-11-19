#!/usr/bin/env python3
"""Wave 1 HC-01 runner: Feshbach-Schur cubic Dirac gap decomposition and continuum limit."""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import subprocess
from pathlib import Path

from kitaev_mf.services.cubic_gap_perturbation import extrapolate_cubic_prefactor_continuum


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def _git_sha(root: Path) -> str:
    return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()


def main() -> int:
    root = Path(__file__).resolve().parents[2]
    out_dir = root / "data" / "processed" / "benchmarks" / "wave1"
    dat_path = out_dir / "hc01_cubic_gap_reconciliation.dat"
    sha = _git_sha(root)

    rows, summary = extrapolate_cubic_prefactor_continuum()
    lines = [
        f"# Benchmark: {dat_path.name}",
        f"# Generated: {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
        f"# Git Commit: {sha[:7]} ({sha})",
        "# Description: HC-01 Feshbach-Schur cubic gap decomposition (canonical 4*H_B units, K=-1, [111])",
        f"# Summary: c_MF_continuum={summary.c_mf_continuum:.6f}, ||lambda*||/h_cont={summary.lambda_over_h_continuum:.6f}, slope={summary.exponent_fit:.6f}",
        "# Columns: h  nk  conv  norm_lam_over_h  gauge_dimer_gap_can  delta_can_direct  delta_schur_full  delta_schur_bare_zeeman  delta_3rd_dressed  delta_3rd_qq_lambda_only  delta_3rd_zeeman_uv_no_lambda  delta_3rd_lambda_only  delta_direct_over_h3  norm_dv_pp",
        "#" + "-" * 78,
    ]
    for r in rows:
        lines.append(
            f"{r.h:8.5f} {r.nk:4d} {int(r.converged):1d} {r.norm_lambda_over_h:12.6f} "
            f"{r.gauge_dimer_gap_can:12.6f} {r.delta_can_direct:16.8e} {r.delta_schur_full:16.8e} "
            f"{r.delta_schur_bare_zeeman:14.6e} {r.delta_3rd_dressed:16.8e} "
            f"{r.delta_3rd_qq_lambda_only:16.8e} {r.delta_3rd_zeeman_uv_no_lambda:14.6e} "
            f"{r.delta_3rd_lambda_only:16.8e} {r.delta_direct_over_h3:12.6f} {r.norm_dv_pp:12.4e}"
        )

    _atomic_write(dat_path, "\n".join(lines) + "\n")
    meta = {
        "schema_version": 1,
        "dataset": "hc01_cubic_gap_reconciliation",
        "source_repo": "kitaev-mf",
        "git_sha": sha,
        "sha256": hashlib.sha256(dat_path.read_bytes()).hexdigest(),
        "unit_convention": "pauli_k_minus_one_canonical_4HB",
        "geometry": "Honeycomb 2-site BZ grid Nk x Nk",
        "parameters": {
            "nk_grid": list(summary.nk_grid),
            "c_mf_by_nk": {str(k): v for k, v in summary.c_mf_by_nk.items()},
            "c_mf_continuum": summary.c_mf_continuum,
            "lambda_over_h_continuum": summary.lambda_over_h_continuum,
            "exponent_fit": summary.exponent_fit,
        },
    }
    _atomic_write(
        dat_path.with_suffix(".meta.json"), json.dumps(meta, indent=2, sort_keys=True) + "\n"
    )
    print(
        f"Exported HC-01 benchmark -> {dat_path} (c_cont={summary.c_mf_continuum:.6f}, p={summary.exponent_fit:.6f})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
