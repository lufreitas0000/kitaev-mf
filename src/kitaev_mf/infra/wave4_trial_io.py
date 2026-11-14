"""Create-only, atomic JSONL storage for Wave 4 raw trial records."""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import asdict, fields
from pathlib import Path
from typing import Any, Sequence

from kitaev_mf.domain.wave4_trials import TrialRecord

SCHEMA_VERSION = 1
_RECORD_FIELDS = {field.name for field in fields(TrialRecord)}


def _reject_duplicate_keys(items: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in items:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _encode_record(record: TrialRecord) -> str:
    return json.dumps(
        {"schema_version": SCHEMA_VERSION, "record": asdict(record)},
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def write_trial_records(path: Path, records: Sequence[TrialRecord]) -> None:
    """Atomically create a new JSONL file; never replace prior evidence.

    Does not assert statistical approval or scientific validity of the records.
    """
    path = Path(path)
    if path.suffix != ".jsonl":
        raise ValueError("raw trial output must have .jsonl extension")
    if not records:
        raise ValueError("refusing to write an empty trial file")
    identities = {
        (r.protocol_id, r.stratum_id, r.cohort, r.method, r.case_id) for r in records
    }
    if len(identities) != len(records):
        raise ValueError("duplicate method/case trial identity")
    path.parent.mkdir(parents=True, exist_ok=True)
    temp_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="\n", dir=path.parent, delete=False
        ) as temporary:
            temp_path = Path(temporary.name)
            for record in records:
                temporary.write(_encode_record(record) + "\n")
            temporary.flush()
            os.fsync(temporary.fileno())
        # Atomic creation; unlike os.replace, hard-linking cannot overwrite a file.
        os.link(temp_path, path)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)


def _decode_record(line: str) -> TrialRecord:
    try:
        document = json.loads(line, object_pairs_hook=_reject_duplicate_keys)
    except (json.JSONDecodeError, TypeError) as exc:
        raise ValueError("malformed JSONL record") from exc
    if not isinstance(document, dict) or set(document) != {"schema_version", "record"}:
        raise ValueError("unexpected JSONL envelope")
    if type(document["schema_version"]) is not int or document["schema_version"] != SCHEMA_VERSION:
        raise ValueError("unsupported Wave 4 JSONL schema")
    payload = document["record"]
    if not isinstance(payload, dict) or set(payload) != _RECORD_FIELDS:
        raise ValueError("missing or unexpected trial record fields")
    if type(payload["converged"]) is not bool:
        raise ValueError("converged must be a boolean")
    for name in ("iterations", "map_evaluations"):
        if type(payload[name]) is not int:
            raise ValueError(f"{name} must be an integer")
    if type(payload["wall_seconds"]) not in (int, float):
        raise ValueError("wall_seconds must be numeric")
    if not isinstance(payload["residual_history"], list) or any(
        type(x) not in (int, float) for x in payload["residual_history"]
    ):
        raise ValueError("residual_history must be a numeric list")
    payload["residual_history"] = tuple(payload["residual_history"])
    try:
        return TrialRecord(**payload)
    except (TypeError, ValueError) as exc:
        raise ValueError("invalid trial record") from exc


def read_trial_records(path: Path) -> list[TrialRecord]:
    """Read and validate complete records; reject truncation or duplicate trials."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines or any(not line.strip() for line in lines):
        raise ValueError("empty line or empty trial file")
    records = [_decode_record(line) for line in lines]
    identities = {
        (r.protocol_id, r.stratum_id, r.cohort, r.method, r.case_id) for r in records
    }
    if len(identities) != len(records):
        raise ValueError("duplicate method/case trial identity")
    return records
