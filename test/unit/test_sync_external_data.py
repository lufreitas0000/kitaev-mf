"""Unit tests for scripts/sync_external_data.py external dataset provenance contract."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_sync_module() -> ModuleType:
    module_path = REPO_ROOT / "scripts" / "sync_external_data.py"
    spec = importlib.util.spec_from_file_location("sync_external_data", module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["sync_external_data"] = module
    spec.loader.exec_module(module)
    return module


_SYNC_MOD = _load_sync_module()
ingest_external_data = _SYNC_MOD.ingest_external_data
main = _SYNC_MOD.main
verify_external_data = _SYNC_MOD.verify_external_data


def _setup_external_tree(tmp_path: Path) -> Path:
    ext = tmp_path / "data" / "external"
    (ext / "historical_unverified").mkdir(parents=True)
    (ext / "kitaev_numerics").mkdir(parents=True)
    (ext / "README.md").write_text("# External Data\n", encoding="utf-8")
    (ext / "historical_unverified" / "README.md").write_text("# Quarantine\n", encoding="utf-8")
    return ext


def test_live_repo_external_data_passes() -> None:
    ok, errors, active_count, quarantined_count = verify_external_data(REPO_ROOT)
    assert ok, errors
    assert active_count >= 0
    assert quarantined_count == 4


def test_ingest_and_verify_roundtrip(tmp_path: Path) -> None:
    _setup_external_tree(tmp_path)
    src_dat = tmp_path / "ed_24_site.dat"
    src_dat.write_text("# source: kitaev-numerics\nh e0\n0.0 -0.5\n", encoding="utf-8")

    sha = "0123456789abcdef0123456789abcdef01234567"
    dest_dat, dest_meta = ingest_external_data(
        root=tmp_path,
        dat_source=src_dat,
        source_repo="kitaev-numerics",
        git_sha=sha,
        unit_convention="K=-1",
        geometry="honeycomb_24_site_torus",
        parameters={"N_sites": 24, "h_direction": "[111]"},
    )
    assert dest_dat.is_file()
    assert dest_meta.is_file()

    ok, errors, active_count, _ = verify_external_data(tmp_path)
    assert ok, errors
    assert active_count == 1
    assert main(["--root", str(tmp_path), "--verify"]) == 0


def test_rejects_sha256_mismatch(tmp_path: Path) -> None:
    ext = _setup_external_tree(tmp_path)
    dat = ext / "kitaev_numerics" / "gap_curve.dat"
    dat.write_text("# header\n0.1 0.002\n", encoding="utf-8")
    meta = {
        "source_repo": "kitaev-numerics",
        "git_sha": "a" * 40,
        "sha256": "0" * 64,
        "unit_convention": "K=-1",
        "geometry": "torus_24",
        "parameters": {},
    }
    dat.with_suffix(".meta.json").write_text(json.dumps(meta), encoding="utf-8")

    ok, errors, _, _ = verify_external_data(tmp_path)
    assert not ok
    assert any("SHA-256 mismatch" in msg for msg in errors)


def test_rejects_missing_sidecar_and_invalid_git_sha(tmp_path: Path) -> None:
    ext = _setup_external_tree(tmp_path)
    dat = ext / "kitaev_numerics" / "missing_sidecar.dat"
    dat.write_text("# header\n1 2\n", encoding="utf-8")

    ok, errors, _, _ = verify_external_data(tmp_path)
    assert not ok
    assert any("Missing sibling metadata sidecar" in msg for msg in errors)

    bad_sha = hashlib.sha256(dat.read_bytes()).hexdigest()
    dat.with_suffix(".meta.json").write_text(
        json.dumps(
            {
                "source_repo": "kitaev-numerics",
                "git_sha": "not-a-sha",
                "sha256": bad_sha,
                "unit_convention": "K=-1",
                "geometry": "torus_16",
                "parameters": {},
            }
        ),
        encoding="utf-8",
    )
    ok2, errors2, _, _ = verify_external_data(tmp_path)
    assert not ok2
    assert any("git_sha" in msg for msg in errors2)


def test_quarantine_dat_is_skipped_when_documented(tmp_path: Path) -> None:
    ext = _setup_external_tree(tmp_path)
    legacy = ext / "historical_unverified" / "old_table.dat"
    legacy.write_text("0.0 1.0\n", encoding="utf-8")

    ok, errors, active_count, quarantined_count = verify_external_data(tmp_path)
    assert ok, errors
    assert active_count == 0
    assert quarantined_count == 1
