#!/usr/bin/env python3
"""Sanitize file trailing whitespace and EOF blank lines."""

from __future__ import annotations

import sys
from pathlib import Path

EXCLUDED_DIRS = {
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

TEXT_EXTENSIONS = {
    ".py",
    ".md",
    ".yaml",
    ".yml",
    ".toml",
    ".txt",
    ".in",
    ".sh",
    "Makefile",
}


def is_text_file(path: Path) -> bool:
    if path.name == "Makefile" or path.suffix in TEXT_EXTENSIONS:
        return True
    return False


def is_excluded(path: Path, root: Path) -> bool:
    rel_parts = path.relative_to(root).parts
    for excluded in EXCLUDED_DIRS:
        ex_parts = Path(excluded).parts
        if len(rel_parts) >= len(ex_parts) and rel_parts[: len(ex_parts)] == ex_parts:
            return True
    return False


def fix_whitespace_in_file(path: Path) -> bool:
    try:
        content = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return False

    lines = [line.rstrip(" \t\r") for line in content.splitlines()]

    while lines and lines[-1] == "":
        lines.pop()

    new_content = "\n".join(lines) + "\n" if lines else ""

    if new_content != content:
        path.write_text(new_content, encoding="utf-8")
        return True
    return False


def main() -> int:
    root = Path.cwd()
    changed_count = 0

    if len(sys.argv) > 1:
        targets = [Path(arg).resolve() for arg in sys.argv[1:]]
    else:
        targets = [
            p
            for p in root.rglob("*")
            if p.is_file() and is_text_file(p) and not is_excluded(p, root)
        ]

    for path in targets:
        if path.is_file() and not is_excluded(path, root):
            if fix_whitespace_in_file(path):
                print(f"Fixed whitespace: {path.relative_to(root)}")
                changed_count += 1

    print(f"Total files updated: {changed_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
