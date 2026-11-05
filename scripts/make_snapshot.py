#!/usr/bin/env python3
"""Generate a Markdown context snapshot of the repository."""

from __future__ import annotations

import datetime
import subprocess
from pathlib import Path

SNAPSHOT_FILENAME = "REPO_SNAPSHOT.md"

ROOT_FILES = [
    "AGENTS.md",
    "PROJECT_STATE.yaml",
    "README.md",
    "CHANGELOG.md",
    "pyproject.toml",
    "Makefile",
]

DOC_PATTERNS = [
    "doc/*.md",
    "doc/governance/*.md",
    ".agents/rules/*.md",
    "doc/adr/*.md",
    "doc/rnd/*.md",
]

CODE_PATTERNS = [
    "scripts/*.py",
    "src/kitaev_mf/**/*.py",
    "test/**/*.py",
]

EXCLUDED_PARTS = {
    ".git",
    ".venv",
    "venv",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    "build",
    "dist",
    "src/kitaev_mf.egg-info",
    "doc/reference",
}


def get_git_info() -> tuple[str, str]:
    try:
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"], text=True
        ).strip()
        commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        branch, commit = "unknown", "unknown"
    return branch, commit


def get_directory_tree(root: Path) -> list[str]:
    files: list[str] = []
    for p in sorted(root.rglob("*")):
        if p.name == SNAPSHOT_FILENAME or "repo_snapshot" in p.name:
            continue
        rel_parts = p.relative_to(root).parts
        if any(ex in rel_parts for ex in EXCLUDED_PARTS):
            continue
        if p.is_file():
            files.append(p.relative_to(root).as_posix())
    return files


def collect_files(root: Path, patterns: list[str]) -> list[Path]:
    collected: set[Path] = set()
    for pattern in patterns:
        for p in sorted(root.glob(pattern)):
            if not p.is_file():
                continue
            rel_parts = p.relative_to(root).parts
            if any(ex in rel_parts for ex in EXCLUDED_PARTS):
                continue
            collected.add(p)
    return sorted(collected)


def main() -> None:
    root = Path.cwd()
    branch, commit = get_git_info()
    now_utc = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines: list[str] = [
        f"# Repository Context Snapshot: {root.name}",
        "",
        f"- **Generated:** {now_utc}",
        f"- **Git Branch:** {branch}",
        f"- **Git Commit:** {commit}",
        "",
        "## 1. Directory Tree",
        "```text",
    ]

    for f_path in get_directory_tree(root):
        lines.append(f_path)
    lines.extend(["```", "", "## 2. Governance, State & Documentation", ""])

    # Root files
    for r_name in ROOT_FILES:
        p = root / r_name
        if p.is_file():
            lines.extend(
                [
                    f"### File: `{r_name}`",
                    "```markdown"
                    if p.suffix == ".md"
                    else "```yaml"
                    if p.suffix in (".yaml", ".yml")
                    else "```toml"
                    if p.suffix == ".toml"
                    else "```makefile"
                    if p.name == "Makefile"
                    else "```text",
                    p.read_text(encoding="utf-8").strip(),
                    "```",
                    "",
                ]
            )

    # Doc files
    for p in collect_files(root, DOC_PATTERNS):
        rel = p.relative_to(root).as_posix()
        if rel in ROOT_FILES:
            continue
        lines.extend(
            [f"### File: `{rel}`", "```markdown", p.read_text(encoding="utf-8").strip(), "```", ""]
        )

    lines.extend(["## 3. Source Code & Tests", ""])
    for p in collect_files(root, CODE_PATTERNS):
        rel = p.relative_to(root).as_posix()
        lines.extend(
            [f"### File: `{rel}`", "```python", p.read_text(encoding="utf-8").strip(), "```", ""]
        )

    output = "\n".join(lines) + "\n"
    out_path = root / SNAPSHOT_FILENAME
    out_path.write_text(output, encoding="utf-8")
    print(f"Snapshot written to {out_path} ({len(output.splitlines())} lines)")


if __name__ == "__main__":
    main()
