"""Generate a single uncommitted Markdown snapshot of the repository for AI context ingestion."""

from __future__ import annotations

import datetime
import subprocess
from pathlib import Path

OUTPUT_FILE = Path("REPO_SNAPSHOT.md")

EXCLUDE_DIRS = {
    ".git",
    ".venv",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    "__pycache__",
    "build",
    "dist",
    "doc/reference/thesis",
}

DOC_PATTERNS = [
    "AGENTS.md",
    "PROJECT_STATE.yaml",
    "pyproject.toml",
    "Makefile",
    ".agents/rules/*.md",
    "doc/adr/*.md",
    "doc/rnd/*.md",
]


def _get_git_info() -> tuple[str, str]:
    try:
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
        commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
        return branch, commit
    except Exception:
        return "unknown", "unknown"


def _should_exclude(path: Path) -> bool:
    path_str = str(path).replace("\\", "/")
    return any(ex in path_str for ex in EXCLUDE_DIRS) or path.suffix == ".egg-info"


def _get_syntax_fence(path: Path) -> str:
    ext_map = {
        ".py": "python",
        ".md": "markdown",
        ".yaml": "yaml",
        ".yml": "yaml",
        ".toml": "toml",
        ".txt": "text",
        ".in": "text",
    }
    if path.name == "Makefile":
        return "makefile"
    return ext_map.get(path.suffix, "text")


def build_snapshot() -> None:
    root = Path(".")
    branch, commit = _get_git_info()
    now_utc = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines: list[str] = [
        "# Repository Context Snapshot: kitaev-mf",
        "",
        f"- **Generated:** {now_utc}",
        f"- **Git Branch:** {branch}",
        f"- **Git Commit:** {commit}",
        "",
        "## 1. Directory Tree",
        "```text",
    ]

    # Build filtered directory tree
    tree_paths: list[str] = []
    for p in sorted(root.rglob("*")):
        if (
            p.is_file()
            and not _should_exclude(p)
            and p.name not in (OUTPUT_FILE.name, "repo_snapshot.txt")
        ):
            tree_paths.append(str(p.as_posix()))

    lines.extend(tree_paths)
    lines.extend(["```", "", "## 2. Governance, State & Documentation", ""])

    included_paths: set[Path] = set()

    for pattern in DOC_PATTERNS:
        for p in sorted(root.glob(pattern)):
            if p.is_file() and not _should_exclude(p):
                included_paths.add(p)
                fence = _get_syntax_fence(p)
                lines.append(f"### File: `{p.as_posix()}`")
                lines.append(f"```{fence}")
                lines.append(p.read_text(encoding="utf-8").rstrip())
                lines.append("```")
                lines.append("")

    lines.extend(["## 3. Source Code & Tests", ""])

    for sub_dir in ("src", "test"):
        for p in sorted((root / sub_dir).rglob("*.py")):
            if p.is_file() and not _should_exclude(p):
                included_paths.add(p)
                fence = _get_syntax_fence(p)
                lines.append(f"### File: `{p.as_posix()}`")
                lines.append(f"```{fence}")
                lines.append(p.read_text(encoding="utf-8").rstrip())
                lines.append("```")
                lines.append("")

    OUTPUT_FILE.write_text("\n".join(lines) + "\n", encoding="utf-8")
    line_count = len(OUTPUT_FILE.read_text(encoding="utf-8").splitlines())
    print(f"Snapshot written to {OUTPUT_FILE} ({line_count} lines, {len(included_paths)} files).")


if __name__ == "__main__":
    build_snapshot()
