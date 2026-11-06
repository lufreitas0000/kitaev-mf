"""Shared utilities, data schemas, and file I/O for M1.6 exploratory benchmarks."""

from __future__ import annotations

import datetime
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import torch

BENCHMARK_DIR = Path("data/processed/benchmarks/m1_6")
BUILD_DIR = Path("build/rnd/m1_6")
FIGURES_DIR = Path("build/rnd/m1_6/figures")


def setup_cpu_environment(num_threads: int = 1) -> None:
    """Enforce single-threaded CPU execution for process-level parallel sweeps."""
    torch.set_num_threads(num_threads)
    torch.set_num_interop_threads(num_threads)


def get_git_commit() -> str:
    """Retrieve short git commit SHA."""
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def write_dat_benchmark(
    filename: str,
    header_lines: list[str],
    data_lines: list[str],
) -> Path:
    """Write structured tabular benchmark data file with provenance header."""
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


def append_jsonl_record(log_name: str, record: dict[str, Any]) -> None:
    """Append execution telemetry record to an ephemeral JSONL file in build/."""
    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    out_path = BUILD_DIR / log_name

    def default_serializer(obj: Any) -> Any:
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.floating, np.integer)):
            return obj.item()
        return str(obj)

    line = json.dumps(record, default=default_serializer)
    with out_path.open("a", encoding="utf-8") as f:
        f.write(line + "\n")
