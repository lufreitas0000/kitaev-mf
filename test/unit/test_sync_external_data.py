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


def _valid_meta(dat_path: Path, source_repo: str = "kitaev-numerics") -> dict[str, object]:
    return {
        "source_repo": source_repo,
        "git_sha": "a" * 40,
        "sha256": hashlib.sha256(dat_path.read_bytes()).hexdigest(),
        "unit_convention": "K=-1",
        "geometry": "torus_24",
        "parameters": {},
    }


def test_live_repo_external_data_passes() -> None:
    ok, errors, active_count, quarantined_count = verify_external_data(REPO_ROOT)
    assert ok, errors
    assert active_count >= 3
    assert quarantined_count >= 1


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


def test_rejects_sha256_mismatch_without_counting_verified(tmp_path: Path) -> None:
    ext = _setup_external_tree(tmp_path)
    dat = ext / "kitaev_numerics" / "gap_curve.dat"
    dat.write_text("# header\n0.1 0.002\n", encoding="utf-8")
    meta = _valid_meta(dat)
    meta["sha256"] = "0" * 64
    dat.with_suffix(".meta.json").write_text(json.dumps(meta), encoding="utf-8")

    ok, errors, verified_count, _ = verify_external_data(tmp_path)
    assert not ok
    assert verified_count == 0
    assert any("SHA-256 mismatch" in msg for msg in errors)


def test_rejects_missing_sidecar_and_invalid_git_sha(tmp_path: Path) -> None:
    ext = _setup_external_tree(tmp_path)
    dat = ext / "kitaev_numerics" / "missing_sidecar.dat"
    dat.write_text("# header\n1 2\n", encoding="utf-8")

    ok, errors, verified_count, _ = verify_external_data(tmp_path)
    assert not ok
    assert verified_count == 0
    assert any("Missing sibling metadata sidecar" in msg for msg in errors)

    meta = _valid_meta(dat)
    meta["git_sha"] = "not-a-sha"
    dat.with_suffix(".meta.json").write_text(json.dumps(meta), encoding="utf-8")
    ok2, errors2, verified_count2, _ = verify_external_data(tmp_path)
    assert not ok2
    assert verified_count2 == 0
    assert any("git_sha" in msg for msg in errors2)


def test_rejects_headerless_dat_dir_mismatch_duplicate_json_and_orphan_sidecar(
    tmp_path: Path,
) -> None:
    ext = _setup_external_tree(tmp_path)

    headerless = ext / "kitaev_numerics" / "no_header.dat"
    headerless.write_text("0.0 -0.5\n", encoding="utf-8")
    headerless.with_suffix(".meta.json").write_text(
        json.dumps(_valid_meta(headerless)), encoding="utf-8"
    )
    ok_hdr, err_hdr, ver_hdr, _ = verify_external_data(tmp_path)
    assert not ok_hdr
    assert ver_hdr == 0
    assert any("provenance header" in msg for msg in err_hdr)

    headerless.unlink()
    headerless.with_suffix(".meta.json").unlink()

    mismatch = ext / "kitaev_numerics" / "wrong_repo.dat"
    mismatch.write_text("# header\n0.0 -0.5\n", encoding="utf-8")
    mismatch.with_suffix(".meta.json").write_text(
        json.dumps(_valid_meta(mismatch, source_repo="kitaev-walk-geometry")), encoding="utf-8"
    )
    ok_dir, err_dir, ver_dir, _ = verify_external_data(tmp_path)
    assert not ok_dir
    assert ver_dir == 0
    assert any("does not match normalized source_repo" in msg for msg in err_dir)

    mismatch.unlink()
    mismatch.with_suffix(".meta.json").unlink()

    dup_dat = ext / "kitaev_numerics" / "dup_keys.dat"
    dup_dat.write_text("# header\n0.0 -0.5\n", encoding="utf-8")
    dup_dat.with_suffix(".meta.json").write_text(
        '{"source_repo": "kitaev-numerics", "source_repo": "kitaev-numerics"}',
        encoding="utf-8",
    )
    ok_dup, err_dup, ver_dup, _ = verify_external_data(tmp_path)
    assert not ok_dup
    assert ver_dup == 0
    assert any("Duplicate JSON key" in msg for msg in err_dup)

    dup_dat.unlink()
    dup_dat.with_suffix(".meta.json").unlink()

    orphan = ext / "kitaev_numerics" / "orphan.meta.json"
    orphan.write_text("{}", encoding="utf-8")
    ok_orp, err_orp, ver_orp, _ = verify_external_data(tmp_path)
    assert not ok_orp
    assert ver_orp == 0
    assert any("Orphan sidecar" in msg for msg in err_orp)


def test_quarantine_dat_is_skipped_when_documented(tmp_path: Path) -> None:
    ext = _setup_external_tree(tmp_path)
    legacy = ext / "historical_unverified" / "old_table.dat"
    legacy.write_text("0.0 1.0\n", encoding="utf-8")

    ok, errors, active_count, quarantined_count = verify_external_data(tmp_path)
    assert ok, errors
    assert active_count == 0
    assert quarantined_count == 1
