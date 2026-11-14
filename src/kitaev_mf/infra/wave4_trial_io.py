"""Create-only, checksum-verified case archives and raw Wave 4 trial records."""

from __future__ import annotations

import base64
import binascii
import json
import os
import tempfile
from collections.abc import Sequence
from dataclasses import asdict, fields
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.domain.wave4_contract import ExperimentContract
from kitaev_mf.domain.wave4_trials import PairedCase, TrialRecord

SCHEMA_VERSION = 2
_RECORD_FIELDS = {field.name for field in fields(TrialRecord)}
_CASE_FIELDS = {
    "protocol_id", "case_id", "stratum_id", "cohort", "measure_ref", "population_id",
    "experiment_json", "experiment_digest", "sample_seed", "initial_state_hash", "arrays",
}
_ARRAY_NAMES = ("U", "V", "lambdas")


def _reject_duplicate_keys(items: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in items:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _encode(document: dict[str, Any]) -> str:
    return json.dumps(document, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _read_jsonl(path: Path) -> list[dict[str, Any]]:
    raw = Path(path).read_bytes()
    if not raw or not raw.endswith(b"\n"):
        raise ValueError("empty or incomplete JSONL file: final newline required")
    try:
        lines = raw.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise ValueError("malformed UTF-8 archive") from exc
    if not lines or any(not line.strip() for line in lines):
        raise ValueError("empty line or empty archive")
    decoded: list[dict[str, Any]] = []
    for line in lines:
        try:
            document = json.loads(line, object_pairs_hook=_reject_duplicate_keys)
        except (json.JSONDecodeError, TypeError) as exc:
            raise ValueError("malformed JSONL record") from exc
        if not isinstance(document, dict):
            raise ValueError("JSONL entries must be objects")
        decoded.append(document)
    return decoded


def _create_jsonl(path: Path, documents: Sequence[dict[str, Any]]) -> None:
    path = Path(path)
    if path.suffix != ".jsonl":
        raise ValueError("Wave 4 archive must use the .jsonl extension")
    if not documents:
        raise ValueError("refusing to write an empty archive")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="\n", dir=path.parent, delete=False
        ) as temporary:
            temporary_path = Path(temporary.name)
            for document in documents:
                temporary.write(_encode(document) + "\n")
            temporary.flush()
            os.fsync(temporary.fileno())
        # Create-only publication: a preexisting file is never overwritten.
        os.link(temporary_path, path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _case_key(case: PairedCase) -> tuple[str, str, str, str]:
    return (case.protocol_id, case.stratum_id, case.cohort, case.case_id)


def _case_document(case: PairedCase) -> dict[str, Any]:
    state = case.materialize_state()
    reconstructed = PairedCase.from_state(
        protocol_id=case.protocol_id,
        case_id=case.case_id,
        stratum_id=case.stratum_id,
        cohort=case.cohort,
        measure_ref=case.measure_ref,
        experiment=case.experiment,
        sample_seed=case.sample_seed,
        state=state,
    )
    if reconstructed.initial_state_hash != case.initial_state_hash:
        raise ValueError("case snapshot checksum mismatch")
    arrays = {}
    for name, array in (("U", state.U), ("V", state.V), ("lambdas", state.lambdas)):
        fixed = np.asarray(array, dtype=np.dtype("<f8"), order="C")
        arrays[name] = {
            "dtype": "<f8",
            "shape": list(fixed.shape),
            "data_base64": base64.b64encode(fixed.tobytes(order="C")).decode("ascii"),
        }
    return {
        "protocol_id": case.protocol_id,
        "case_id": case.case_id,
        "stratum_id": case.stratum_id,
        "cohort": case.cohort,
        "measure_ref": case.measure_ref,
        "population_id": case.population_id,
        "experiment_json": case.experiment.document_json,
        "experiment_digest": case.experiment.digest,
        "sample_seed": case.sample_seed,
        "initial_state_hash": case.initial_state_hash,
        "arrays": arrays,
    }


def write_paired_cases(path: Path, cases: Sequence[PairedCase]) -> None:
    """Persist exact float64 initial arrays, with full contract and verified hashes."""
    if len({_case_key(case) for case in cases}) != len(cases):
        raise ValueError("duplicate case identity")
    _create_jsonl(
        path, [{"schema_version": SCHEMA_VERSION, "case": _case_document(case)} for case in cases]
    )


def _decode_case(document: dict[str, Any]) -> PairedCase:
    if set(document) != {"schema_version", "case"} or (
        type(document["schema_version"]) is not int or document["schema_version"] != SCHEMA_VERSION
    ):
        raise ValueError("invalid case archive version/envelope")
    payload = document["case"]
    if not isinstance(payload, dict) or set(payload) != _CASE_FIELDS:
        raise ValueError("missing or unexpected case fields")
    experiment = ExperimentContract(payload["experiment_json"])
    if payload["experiment_digest"] != experiment.digest:
        raise ValueError("case experiment checksum mismatch")
    encoded_arrays = payload["arrays"]
    if not isinstance(encoded_arrays, dict) or set(encoded_arrays) != set(_ARRAY_NAMES):
        raise ValueError("missing or unexpected initial arrays")
    arrays: dict[str, NDArray[np.float64]] = {}
    for name in _ARRAY_NAMES:
        item = encoded_arrays[name]
        if not isinstance(item, dict) or set(item) != {"dtype", "shape", "data_base64"}:
            raise ValueError("invalid array metadata")
        if item["dtype"] != "<f8":
            raise ValueError("only little-endian float64 snapshot arrays are supported")
        shape = item["shape"]
        if not isinstance(shape, list) or any(type(n) is not int or n < 0 for n in shape):
            raise ValueError("invalid snapshot shape")
        try:
            decoded = base64.b64decode(item["data_base64"], validate=True)
        except (ValueError, binascii.Error, TypeError) as exc:
            raise ValueError("invalid snapshot encoding") from exc
        if len(decoded) != 8 * int(np.prod(shape, dtype=np.int64)):
            raise ValueError("truncated or oversized snapshot array")
        arrays[name] = np.frombuffer(decoded, dtype=np.dtype("<f8")).reshape(shape).copy()
    state = MeanFieldState(U=arrays["U"], V=arrays["V"], lambdas=arrays["lambdas"])
    case = PairedCase.from_state(
        protocol_id=payload["protocol_id"],
        case_id=payload["case_id"],
        stratum_id=payload["stratum_id"],
        cohort=payload["cohort"],
        measure_ref=payload["measure_ref"],
        experiment=experiment,
        sample_seed=payload["sample_seed"],
        state=state,
    )
    if case.population_id != payload["population_id"] or (
        case.initial_state_hash != payload["initial_state_hash"]
    ):
        raise ValueError("initial-state snapshot checksum mismatch")
    return case


def read_paired_cases(path: Path) -> list[PairedCase]:
    """Reconstruct exact initial states; reject mismatched checksums or partial archives."""
    cases = [_decode_case(document) for document in _read_jsonl(path)]
    if len({_case_key(case) for case in cases}) != len(cases):
        raise ValueError("duplicate case identity")
    return cases


def _record_key(record: TrialRecord) -> tuple[str, str, str, str, str]:
    return (
        record.protocol_id, record.stratum_id, record.cohort, record.method, record.case_id
    )


def _verify_case_links(records: Sequence[TrialRecord], cases: Sequence[PairedCase]) -> None:
    case_map = {_case_key(case): case for case in cases}
    if not case_map or len(case_map) != len(cases):
        raise ValueError("invalid or duplicate case archive")
    for record in records:
        key = (record.protocol_id, record.stratum_id, record.cohort, record.case_id)
        if key not in case_map:
            raise ValueError("trial lacks a persisted initial-state snapshot")
        case = case_map[key]
        if (
            record.initial_state_hash != case.initial_state_hash
            or record.experiment_json != case.experiment.document_json
            or record.experiment_digest != case.experiment.digest
            or record.measure_ref != case.measure_ref
            or record.population_id != case.population_id
            or record.sample_seed != case.sample_seed
        ):
            raise ValueError("trial and persisted initial-state contract mismatch")


def write_trial_records(
    path: Path, records: Sequence[TrialRecord], *, cases_path: Path
) -> None:
    """Atomically create trial telemetry only after verifying saved initial states."""
    if len({_record_key(r) for r in records}) != len(records):
        raise ValueError("duplicate method/case trial identity")
    _verify_case_links(records, read_paired_cases(cases_path))
    _create_jsonl(
        path, [{"schema_version": SCHEMA_VERSION, "record": asdict(r)} for r in records]
    )


def _decode_record(document: dict[str, Any]) -> TrialRecord:
    if set(document) != {"schema_version", "record"} or (
        type(document["schema_version"]) is not int or document["schema_version"] != SCHEMA_VERSION
    ):
        raise ValueError("unsupported Wave 4 trial schema")
    payload = document["record"]
    if not isinstance(payload, dict) or set(payload) != _RECORD_FIELDS:
        raise ValueError("missing or unexpected trial record fields")
    if type(payload["converged"]) is not bool:
        raise ValueError("converged must be a boolean")
    for name in ("iterations", "map_evaluations"):
        if type(payload[name]) is not int:
            raise ValueError(f"{name} must be an integer")
    if payload["sample_seed"] is not None and type(payload["sample_seed"]) is not int:
        raise ValueError("sample_seed must be an integer or null")
    for name in (
        "wall_seconds", "final_correlator_residual", "final_onsite_residual",
        "final_gauge_residual",
    ):
        if payload[name] is not None and type(payload[name]) not in (int, float):
            raise ValueError(f"{name} must be numeric or null")
    for name in ("residual_history", "warnings"):
        if not isinstance(payload[name], list):
            raise ValueError(f"{name} must be a list")
    if any(type(x) not in (int, float) for x in payload["residual_history"]):
        raise ValueError("residual_history must contain numeric values")
    if any(not isinstance(x, str) for x in payload["warnings"]):
        raise ValueError("warnings must contain strings")
    payload["residual_history"] = tuple(payload["residual_history"])
    payload["warnings"] = tuple(payload["warnings"])
    try:
        return TrialRecord(**payload)
    except (TypeError, ValueError) as exc:
        raise ValueError("invalid trial record") from exc


def read_trial_records(path: Path, *, cases_path: Path) -> list[TrialRecord]:
    """Read complete records only with corresponding checksum-verified case archive."""
    records = [_decode_record(document) for document in _read_jsonl(path)]
    if len({_record_key(r) for r in records}) != len(records):
        raise ValueError("duplicate method/case trial identity")
    _verify_case_links(records, read_paired_cases(cases_path))
    return records
