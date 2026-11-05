#!/usr/bin/env python3
"""Generate a lean, context-optimized Markdown snapshot of the repository.

By default, includes active milestone specifications, core governance,
production code, and tests. Use --full to embed all historical notes.
"""

from __future__ import annotations

import argparse
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

# Core governance and active milestone R&D notes
ACTIVE_DOC_PATTERNS = [
    "doc/*.md",
    "doc/governance/*.md",
    ".agents/rules/*.md",
    "doc/adr/*.md",
    "doc/rnd/note_04f_advanced_invariants_and_sum_rules_backlog.md",
    "doc/rnd/note_04g_variational_kkt_and_convergence_asymptotics.md",
    "doc/rnd/note_05_quasifree_states_and_variational_stability.md",
]

FULL_DOC_PATTERNS = [
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

EXCLUDED_DIR_NAMES = {
    ".git",
    ".venv",
    "venv",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    "__pycache__",
    "build",
    "dist",
    "src/kitaev_mf.egg-info",
    "doc/reference",
}

EXCLUDED_EXTENSIONS = {
    ".pyc",
    ".pdf",
    ".png",
    ".jpg",
    ".jpeg",
    ".svg",
}


def get_git_info() -> tuple[str, str]:
    try:
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"], text=True
        ).strip()
        commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True
        ).strip()
    except Exception:
        branch, commit = "unknown", "unknown"
    return branch, commit


def get_directory_tree(root: Path) -> list[str]:
    files: list[str] = []
    for p in sorted(root.rglob("*")):
        if p.name == SNAPSHOT_FILENAME or "repo_snapshot" in p.name:
            continue
        if p.name.startswith(".coverage"):
            continue
        if p.suffix in EXCLUDED_EXTENSIONS:
            continue
        rel_parts = p.relative_to(root).parts
        if any(part in EXCLUDED_DIR_NAMES for part in rel_parts):
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
            if p.suffix in EXCLUDED_EXTENSIONS:
                continue
            rel_parts = p.relative_to(root).parts
            if any(part in EXCLUDED_DIR_NAMES for part in rel_parts):
                continue
            collected.add(p)
    return sorted(collected)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate repository context snapshot")
    parser.add_argument(
        "--full",
        action="store_true",
        help="Embed all historical and scratchpad R&D notes instead of active set",
    )
    args = parser.parse_args()

    root = Path.cwd()
    branch, commit = get_git_info()
    now_utc = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    doc_patterns = FULL_DOC_PATTERNS if args.full else ACTIVE_DOC_PATTERNS

    lines: list[str] = [
        f"# Repository Context Snapshot: {root.name}",
        "",
        f"- **Generated:** {now_utc}",
        f"- **Git Branch:** {branch}",
        f"- **Git Commit:** {commit}",
        f"- **Scope:** {'Full Archive' if args.full else 'Active Milestone (M1.4.5) Context'}",
        "",
        "## 1. Directory Tree (Filtered)",
        "```text",
    ]

    for f_path in get_directory_tree(root):
        lines.append(f_path)
    lines.extend(["```", "", "## 2. Governance, State & Active Documentation", ""])

    # Root files
    for r_name in ROOT_FILES:
        p = root / r_name
        if p.is_file():
            fence = (
                "```markdown"
                if p.suffix == ".md"
                else "```yaml"
                if p.suffix in (".yaml", ".yml")
                else "```toml"
                if p.suffix == ".toml"
                else "```makefile"
                if p.name == "Makefile"
                else "```text"
            )
            lines.extend([f"### File: `{r_name}`", fence, p.read_text(encoding="utf-8").strip(), "```", ""])

    # Selected Doc files
    for p in collect_files(root, doc_patterns):
        rel = p.relative_to(root).as_posix()
        if rel in ROOT_FILES:
            continue
        lines.extend([f"### File: `{rel}`", "```markdown", p.read_text(encoding="utf-8").strip(), "```", ""])

    lines.extend(["## 3. Source Code & Tests", ""])
    for p in collect_files(root, CODE_PATTERNS):
        rel = p.relative_to(root).as_posix()
        lines.extend([f"### File: `{rel}`", "```python", p.read_text(encoding="utf-8").strip(), "```", ""])

    output = "\n".join(lines) + "\n"
    out_path = root / SNAPSHOT_FILENAME
    out_path.write_text(output, encoding="utf-8")
    print(f"Snapshot written to {out_path} ({len(output.splitlines())} lines)")


if __name__ == "__main__":
    main()
