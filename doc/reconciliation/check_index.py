"""Read-only verification of pinned catalog coverage and original untracked inputs."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    inventory = json.loads((root / "doc/reconciliation/INVENTORY.json").read_text())
    for label in ("main_input", "w001_input"):
        output = subprocess.check_output(
            ["git", "-C", str(root), "ls-tree", "-rl", inventory[label]], text=True
        )
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
    for worktree in inventory["worktrees"]:
        # This task's own worktree is a changing candidate, not a preservation input.
        if worktree["branch"] == "refs/heads/doc/reconciliation-foundation-index-20261116":
            continue
        for entry in worktree["untracked_files"]:
            path = Path(worktree["worktree"]) / entry["path"]
            if not path.is_file():
                raise SystemExit(f"UNKNOWN: preservation input moved or unavailable: {path}")
            if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
                raise SystemExit(f"FAIL: preservation input differs from dated observation: {path}")
    print("PASS: exact pinned catalog coverage and original untracked input hashes")


if __name__ == "__main__":
    main()
