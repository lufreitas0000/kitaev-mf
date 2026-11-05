#!/usr/bin/env python3
"""Generate a context-optimized Markdown snapshot of a repository for Web/Local LLMs.

Supports three context modes:
- compact: For local LLMs (Qwen/Ollama) or tight token budgets (governance + state + active task targets)
- active:  Default for Web LLMs (ChatGPT, Claude, Gemini, Qwen Web, DeepSeek)
- full:    Complete repository dump for long-context audits

Optionally prepends a role-specific task prompt header (--role rnd|labbook|manuscript|referee).
"""

from __future__ import annotations

import argparse
import datetime
import subprocess
from pathlib import Path
from typing import Any

import yaml

SNAPSHOT_FILENAME = "REPO_SNAPSHOT.md"

ROLE_MAP = {
    "rnd": "role_rnd_coder.md",
    "labbook": "role_labbook_synthesizer.md",
    "manuscript": "role_manuscript_editor.md",
    "referee": "role_adversarial_referee.md",
}

ROOT_FILES = [
    "AGENTS.md",
    "PROJECT_STATE.yaml",
    "LOG_DISCOVERY_AND_PIVOTS.md",
    "README.md",
    "CHANGELOG.md",
    "pyproject.toml",
    "lakefile.toml",
    "Project.toml",
    "Makefile",
]

ACTIVE_DOC_PATTERNS = [
    "doc/*.md",
    "doc/governance/*.md",
    ".agents/rules/*.md",
    "doc/adr/*.md",
]

FULL_DOC_PATTERNS = [
    "doc/*.md",
    "doc/governance/*.md",
    ".agents/rules/*.md",
    "doc/adr/*.md",
    "doc/rnd/**/*.md",
    "doc/shared/*.tex",
]

CODE_PATTERNS = [
    "scripts/*.py",
    "src/**/*.py",
    "src/**/*.jl",
    "test/**/*.py",
    "tests/**/*.py",
    "**/*.lean",
]

EXCLUDED_DIR_NAMES = {
    ".git",
    ".venv",
    "venv",
    ".lake",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
    "__pycache__",
    "build",
    "cache_tikz",
    "dist",
    "doc/reference",
}

EXCLUDED_EXTENSIONS = {
    ".pyc",
    ".pdf",
    ".png",
    ".jpg",
    ".jpeg",
    ".svg",
    ".h5",
    ".npz",
    ".pt",
    ".olean",
}


def get_git_info(root: Path) -> tuple[str, str]:
    try:
        branch = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "--abbrev-ref", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        commit = subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "--short", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        branch, commit = "unknown", "unknown"
    return branch, commit


def is_excluded(path: Path, root: Path) -> bool:
    if path.name == SNAPSHOT_FILENAME or "repo_snapshot" in path.name.lower():
        return True
    if path.name.startswith(".coverage"):
        return True
    if path.suffix in EXCLUDED_EXTENSIONS:
        return True
    rel_parts = path.relative_to(root).parts
    return any(part in EXCLUDED_DIR_NAMES or part.endswith(".egg-info") for part in rel_parts)


def get_directory_tree(root: Path) -> list[str]:
    files: list[str] = []
    for p in sorted(root.rglob("*")):
        if p.is_file() and not is_excluded(p, root):
            files.append(p.relative_to(root).as_posix())
    return files


def collect_files(root: Path, patterns: list[str]) -> list[Path]:
    collected: set[Path] = set()
    for pattern in patterns:
        for p in sorted(root.glob(pattern)):
            if p.is_file() and not is_excluded(p, root):
                collected.add(p)
    return sorted(collected, key=lambda x: x.relative_to(root).as_posix())


def get_active_milestone_targets(root: Path) -> list[Path]:
    """Read PROJECT_STATE.yaml to discover active milestone R&D notes and target files."""
    state_file = root / "PROJECT_STATE.yaml"
    if not state_file.is_file():
        return []
    try:
        data: dict[str, Any] = yaml.safe_load(state_file.read_text(encoding="utf-8")) or {}
    except Exception:
        return []

    paths: set[Path] = set()
    active_ms = data.get("active_milestone") or {}
    ms_id = str(active_ms.get("id", "")).lower().replace(".", "_")
    if ms_id:
        for p in root.glob(f"doc/rnd/**/{ms_id}_*.md"):
            if p.is_file() and not is_excluded(p, root):
                paths.add(p)

    for task in active_ms.get("active_tasks") or []:
        for target_rel in task.get("target_files") or []:
            candidate = root / target_rel
            if candidate.is_file() and not is_excluded(candidate, root):
                paths.add(candidate)

    return sorted(paths, key=lambda x: x.relative_to(root).as_posix())


def load_role_prompt(root: Path, role: str | None) -> str:
    if not role:
        return ""
    filename = ROLE_MAP.get(role)
    if not filename:
        return ""

    candidates = [
        root / ".agents" / "prompts" / "roles" / filename,
        root / "prompts" / "roles" / filename,
        Path(__file__).resolve().parent.parent / "prompts" / "roles" / filename,
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.read_text(encoding="utf-8").strip()
    return ""


def lang_for_suffix(suffix: str) -> str:
    return {
        ".py": "python",
        ".md": "markdown",
        ".yaml": "yaml",
        ".yml": "yaml",
        ".toml": "toml",
        ".tex": "latex",
        ".lean": "lean",
        ".jl": "julia",
        ".dat": "text",
    }.get(suffix, "text")


def build_snapshot(root: Path, mode: str, role: str | None, output_path: Path) -> None:
    branch, commit = get_git_info(root)
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    root_files = [root / name for name in ROOT_FILES if (root / name).is_file()]
    active_targets = get_active_milestone_targets(root)

    if mode == "compact":
        selected_files = sorted(
            set(root_files) | set(active_targets),
            key=lambda x: x.relative_to(root).as_posix(),
        )
    elif mode == "active":
        docs = collect_files(root, ACTIVE_DOC_PATTERNS)
        code = collect_files(root, CODE_PATTERNS)
        selected_files = sorted(
            set(root_files) | set(docs) | set(active_targets) | set(code),
            key=lambda x: x.relative_to(root).as_posix(),
        )
    else:
        docs = collect_files(root, FULL_DOC_PATTERNS)
        code = collect_files(root, CODE_PATTERNS)
        selected_files = sorted(
            set(root_files) | set(docs) | set(code),
            key=lambda x: x.relative_to(root).as_posix(),
        )

    lines: list[str] = []
    role_header = load_role_prompt(root, role)
    if role_header:
        lines.append(role_header)
        lines.append("\n---\n")

    lines.append(f"# Repository Snapshot: `{root.name}`")
    lines.append(
        f"- **Generated:** `{now}`\n"
        f"- **Branch / Commit:** `{branch}` (`{commit}`)\n"
        f"- **Snapshot Mode:** `{mode}`" + (f" (Role: `{role}`)" if role else "")
    )
    lines.append("\n## 1. Repository File Tree\n```text")
    lines.extend(get_directory_tree(root))
    lines.append("```\n")

    lines.append("## 2. Bundled File Contents\n")
    for file_path in selected_files:
        rel = file_path.relative_to(root).as_posix()
        lang = lang_for_suffix(file_path.suffix)
        try:
            content = file_path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        lines.append(f"### File: `{rel}`\n````{lang}\n{content.rstrip()}\n````\n")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines), encoding="utf-8")
    print(
        f"✓ Wrote {output_path} ({len(selected_files)} files, mode={mode}"
        + (f", role={role}" if role else "")
        + ")"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate LLM context snapshot.")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Repository root directory (default: current directory)",
    )
    parser.add_argument(
        "--mode",
        choices=["compact", "active", "full"],
        default="active",
        help="Context bundle budget: compact (local LLM), active (web LLM), full (audit)",
    )
    parser.add_argument(
        "--role",
        choices=list(ROLE_MAP.keys()),
        default=None,
        help="Optional role prompt header to prepend (rnd, labbook, manuscript, referee)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output path (default: <root>/REPO_SNAPSHOT.md)",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    out = args.output.resolve() if args.output else (root / SNAPSHOT_FILENAME)
    build_snapshot(root, args.mode, args.role, out)


if __name__ == "__main__":
    main()
