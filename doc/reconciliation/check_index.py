"""Read-only verification of pinned catalog coverage and absorbed reconciliation inputs."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

ABSORBED_ARCHIVE_SHA256 = {
    "data/processed/benchmarks/provisional_phase3/mzm_thermodynamic_spectral_benchmarks.dat": "2820b8d191257674e544abb1b3bb86448f07d727372d678df144173c3f3c00c7",
}

REQUIRED_ABSORBED_FILES = [
    "doc/audits/w001_acceptance.yaml",
    "doc/audits/w001_foundation_review.md",
    "doc/audits/w001_perturbation_recovery.md",
    "doc/rnd/phase_3_dynamical_spin_correlations/theory/t3_w1_3_thermodynamic_sum_rules_and_spectral_functions.md",
    "doc/reconciliation/CODE_SALVAGE_QUEUE.md",
]


def _ref_exists(root: Path, ref: str) -> bool:
    return (
        subprocess.run(
            ["git", "-C", str(root), "cat-file", "-e", f"{ref}^{{commit}}"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode
        == 0
    )


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    inventory = json.loads((root / "doc/reconciliation/INVENTORY.json").read_text(encoding="utf-8"))

    for label in ("main_input", "w001_input"):
        ref = inventory[label]
        if not _ref_exists(root, ref):
            continue
        output = subprocess.check_output(["git", "-C", str(root), "ls-tree", "-rl", ref], text=True)
        actual = {}
        for line in output.splitlines():
            metadata, path = line.split("\t", 1)
            _, kind, blob, size = metadata.split()
            if kind == "blob":
                actual[path] = (blob, int(size))
        indexed = {
            entry["path"]: (entry["git_blob"], entry["bytes"])
            for entry in inventory["catalog"]
            if entry["input"] == label
        }
        if indexed != actual:
            raise SystemExit(f"FAIL: catalog coverage or blob mismatch for {label}")

    for rel in REQUIRED_ABSORBED_FILES:
        if not (root / rel).is_file():
            raise SystemExit(f"FAIL: missing absorbed reconciliation artifact: {rel}")

    for rel, expected_sha in ABSORBED_ARCHIVE_SHA256.items():
        path = root / rel
        if not path.is_file():
            raise SystemExit(f"FAIL: missing absorbed benchmark artifact: {rel}")
        actual_sha = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual_sha != expected_sha:
            raise SystemExit(f"FAIL: SHA-256 mismatch for {rel}: {actual_sha} != {expected_sha}")

    print("PASS: absorbed reconciliation artifacts and pinned hashes verified")


if __name__ == "__main__":
    main()
