#!/usr/bin/env python3
"""Ingest and verify external cross-repository benchmark datasets under data/external/.

Enforces ADR-003 cross-repository provenance:
- Quarantine datasets under data/external/historical_unverified/ are excluded from active
  producer checks provided the quarantine README exists.
- Every active .dat file under data/external/<producer>/ MUST have a sibling .meta.json
  sidecar with valid source_repo, 40-hex git_sha, matching sha256, unit_convention,
  geometry, and parameters, plus a leading '#' provenance comment header in the .dat file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path

REQUIRED_META_KEYS = (
    "source_repo",
    "git_sha",
    "sha256",
    "unit_convention",
    "geometry",
    "parameters",
)
HEX40_RE = re.compile(r"^[0-9a-f]{40}$")
HEX64_RE = re.compile(r"^[0-9a-f]{64}$")
QUARANTINE_DIR = "historical_unverified"


def _compute_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_unique_json(path: Path) -> object:
    def _no_dups(pairs: list[tuple[str, object]]) -> dict[str, object]:
        obj: dict[str, object] = {}
        for key, value in pairs:
            if key in obj:
                raise ValueError(f"Duplicate JSON key '{key}' in {path}")
            obj[key] = value
        return obj

    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_no_dups)


def _has_comment_header(dat_path: Path) -> bool:
    for line in dat_path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        return stripped.startswith("#")
    return False


def _normalize_repo_dir(source_repo: str) -> str:
    normalized = source_repo.strip().lower().replace("-", "_")
    if not normalized or "/" in normalized or "\\" in normalized or ".." in normalized:
        raise ValueError(f"Invalid source_repo identifier: {source_repo!r}")
    return normalized


def verify_external_data(root: Path) -> tuple[bool, list[str], int, int]:
    """Verify all external datasets under root/data/external/.

    Returns:
        Tuple of (ok, errors, verified_active_count, quarantined_count).
    """
    ext_root = root / "data" / "external"
    errors: list[str] = []

    if not ext_root.is_dir():
        return False, [f"Missing external data directory: {ext_root}"], 0, 0

    for required_readme in (ext_root / "README.md", ext_root / QUARANTINE_DIR / "README.md"):
        if not required_readme.is_file():
            errors.append(
                f"Missing required documentation file: {required_readme.relative_to(root)}"
            )

    quarantine_root = ext_root / QUARANTINE_DIR
    quarantined_count = len(list(quarantine_root.rglob("*.dat"))) if quarantine_root.is_dir() else 0

    dat_files = sorted(ext_root.rglob("*.dat"))
    meta_files = sorted(ext_root.rglob("*.meta.json"))

    active_dats = [p for p in dat_files if QUARANTINE_DIR not in p.relative_to(ext_root).parts]
    active_metas = [p for p in meta_files if QUARANTINE_DIR not in p.relative_to(ext_root).parts]

    expected_metas = {p.with_suffix(".meta.json") for p in active_dats}
    for meta_path in active_metas:
        if meta_path not in expected_metas:
            errors.append(
                f"Orphan sidecar without matching .dat file: {meta_path.relative_to(root)}"
            )

    verified_active = 0
    for dat_path in active_dats:
        rel_dat = dat_path.relative_to(ext_root)
        errors_before = len(errors)
        if len(rel_dat.parts) < 2:
            errors.append(
                f"External .dat must reside in data/external/<source_repo>/: data/external/{rel_dat}"
            )
            continue

        producer_dir = rel_dat.parts[0]
        if not _has_comment_header(dat_path):
            errors.append(
                f"Active external .dat lacks '#' provenance header: data/external/{rel_dat}"
            )

        meta_path = dat_path.with_suffix(".meta.json")
        if not meta_path.is_file():
            errors.append(f"Missing sibling metadata sidecar: {meta_path.relative_to(root)}")
            continue

        try:
            raw_meta = _parse_unique_json(meta_path)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            errors.append(f"Invalid JSON in {meta_path.relative_to(root)}: {exc}")
            continue

        if not isinstance(raw_meta, dict):
            errors.append(f"Metadata sidecar must be a JSON object: {meta_path.relative_to(root)}")
            continue

        missing = [k for k in REQUIRED_META_KEYS if k not in raw_meta]
        if missing:
            errors.append(
                f"Missing required metadata keys {missing} in {meta_path.relative_to(root)}"
            )
            continue

        source_repo = raw_meta["source_repo"]
        if not isinstance(source_repo, str) or not source_repo.strip():
            errors.append(
                f"source_repo must be a non-empty string in {meta_path.relative_to(root)}"
            )
        else:
            try:
                expected_dir = _normalize_repo_dir(source_repo)
                if expected_dir != producer_dir:
                    errors.append(
                        f"Directory {producer_dir!r} does not match normalized source_repo "
                        f"{source_repo!r} in {meta_path.relative_to(root)}"
                    )
            except ValueError as exc:
                errors.append(str(exc))

        git_sha = raw_meta["git_sha"]
        if not isinstance(git_sha, str) or not HEX40_RE.match(git_sha):
            errors.append(f"git_sha must be 40-char lowercase hex in {meta_path.relative_to(root)}")

        recorded_sha = raw_meta["sha256"]
        actual_sha = _compute_sha256(dat_path)
        if not isinstance(recorded_sha, str) or not HEX64_RE.match(recorded_sha):
            errors.append(f"sha256 must be 64-char lowercase hex in {meta_path.relative_to(root)}")
        elif recorded_sha != actual_sha:
            errors.append(
                f"SHA-256 mismatch for data/external/{rel_dat}: "
                f"recorded={recorded_sha} actual={actual_sha}"
            )

        for str_key in ("unit_convention", "geometry"):
            val = raw_meta[str_key]
            if not isinstance(val, str) or not val.strip():
                errors.append(
                    f"{str_key} must be a non-empty string in {meta_path.relative_to(root)}"
                )

        if not isinstance(raw_meta["parameters"], dict):
            errors.append(f"parameters must be a JSON object in {meta_path.relative_to(root)}")

        if len(errors) == errors_before:
            verified_active += 1

    return len(errors) == 0, errors, verified_active, quarantined_count


def _parse_params(raw_params: list[str]) -> dict[str, object]:
    parsed: dict[str, object] = {}
    for item in raw_params:
        if "=" not in item:
            raise ValueError(f"Expected key=value for --param, got {item!r}")
        key, val = item.split("=", 1)
        key = key.strip()
        if not key:
            raise ValueError(f"Empty parameter key in {item!r}")
        try:
            parsed[key] = json.loads(val)
        except json.JSONDecodeError:
            parsed[key] = val
    return parsed


def ingest_external_data(
    root: Path,
    dat_source: Path,
    source_repo: str,
    git_sha: str,
    unit_convention: str,
    geometry: str,
    parameters: dict[str, object],
) -> tuple[Path, Path]:
    """Copy a validated external .dat file into data/external/<source_repo>/ with .meta.json."""
    if not dat_source.is_file():
        raise FileNotFoundError(f"Input .dat file not found: {dat_source}")
    if dat_source.suffix != ".dat":
        raise ValueError(f"Input file must have .dat extension: {dat_source}")
    if not _has_comment_header(dat_source):
        raise ValueError(f"Input .dat must begin with '#' provenance comment lines: {dat_source}")
    if not HEX40_RE.match(git_sha):
        raise ValueError(f"git_sha must be 40-char lowercase hex, got {git_sha!r}")
    if not unit_convention.strip() or not geometry.strip():
        raise ValueError("unit_convention and geometry must be non-empty strings")

    producer_dir = _normalize_repo_dir(source_repo)
    if producer_dir == QUARANTINE_DIR:
        raise ValueError(f"Cannot ingest active data into quarantine directory {QUARANTINE_DIR}")

    dest_dir = root / "data" / "external" / producer_dir
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest_dat = dest_dir / dat_source.name
    shutil.copy2(dat_source, dest_dat)

    meta: dict[str, object] = {
        "source_repo": source_repo.strip(),
        "git_sha": git_sha,
        "sha256": _compute_sha256(dest_dat),
        "unit_convention": unit_convention.strip(),
        "geometry": geometry.strip(),
        "parameters": parameters,
    }
    dest_meta = dest_dat.with_suffix(".meta.json")
    dest_meta.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return dest_dat, dest_meta


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify or ingest cross-repository external benchmark datasets."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Repository root directory (defaults to kitaev-mf root).",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Verify data/external/ READMEs and all active .dat + .meta.json provenance sidecars.",
    )
    parser.add_argument(
        "--ingest",
        type=Path,
        default=None,
        help="Path to an upstream .dat dataset to ingest into data/external/<source_repo>/.",
    )
    parser.add_argument("--source-repo", default="", help="Producer repository name.")
    parser.add_argument("--git-sha", default="", help="Producer 40-char lowercase Git commit SHA.")
    parser.add_argument("--unit-convention", default="", help="Energy/coupling unit convention.")
    parser.add_argument("--geometry", default="", help="Lattice/cluster geometry specification.")
    parser.add_argument(
        "--param",
        action="append",
        default=[],
        help="Optional metadata parameter key=value (repeatable).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    root: Path = args.root.resolve()

    if args.ingest is not None:
        if (
            not args.source_repo
            or not args.git_sha
            or not args.unit_convention
            or not args.geometry
        ):
            parser.error(
                "--ingest requires --source-repo, --git-sha, --unit-convention, and --geometry"
            )
        try:
            params = _parse_params(list(args.param))
            dest_dat, dest_meta = ingest_external_data(
                root=root,
                dat_source=args.ingest,
                source_repo=str(args.source_repo),
                git_sha=str(args.git_sha),
                unit_convention=str(args.unit_convention),
                geometry=str(args.geometry),
                parameters=params,
            )
        except (OSError, ValueError) as exc:
            print(f"FAIL: ingest error: {exc}", file=sys.stderr)
            return 1
        print(f"INGESTED: {dest_dat.relative_to(root)} and {dest_meta.relative_to(root)}")
        if not args.verify:
            return 0

    if args.verify or args.ingest is None:
        ok, errors, active_count, quarantined_count = verify_external_data(root)
        if not ok:
            for err in errors:
                print(f"FAIL: {err}", file=sys.stderr)
            return 1
        print(
            f"PASS: verified {active_count} active external dataset(s); "
            f"{quarantined_count} legacy dataset(s) quarantined in {QUARANTINE_DIR}/"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
