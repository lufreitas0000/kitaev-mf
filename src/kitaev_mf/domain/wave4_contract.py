"""Canonical immutable metadata contracts for Wave 4 numerical attempts.

These records capture declared settings; they do not certify their physical validity.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any


def _object(value: Any, name: str) -> dict[str, Any]:
    if not isinstance(value, dict) or not all(isinstance(k, str) for k in value):
        raise ValueError(f"{name} must be a JSON object with string keys")
    return value


def _keys(value: dict[str, Any], keys: set[str], name: str) -> None:
    if set(value) != keys:
        raise ValueError(f"{name} must contain precisely {sorted(keys)}")


def _text(value: Any, name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")


def _positive_int(value: Any, name: str) -> None:
    if type(value) is not int or value <= 0:
        raise ValueError(f"{name} must be a positive integer")


def _number(value: Any, name: str, *, positive: bool = False) -> None:
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f"{name} must be finite and numeric")
    if positive and value <= 0:
        raise ValueError(f"{name} must be positive")


def _mapping_of_numbers(value: Any, name: str, *, positive: bool = False) -> None:
    obj = _object(value, name)
    if not obj:
        raise ValueError(f"{name} cannot be empty")
    for key, item in obj.items():
        _text(key, name)
        _number(item, f"{name}.{key}", positive=positive)


def _validate_experiment(document: dict[str, Any]) -> None:
    _keys(document, {"model", "geometry", "numerics", "sampling", "provenance"}, "experiment")
    model = _object(document["model"], "model")
    _keys(model, {"name", "couplings", "field_xyz"}, "model")
    _text(model["name"], "model.name")
    _mapping_of_numbers(model["couplings"], "model.couplings")
    field = model["field_xyz"]
    if not isinstance(field, list) or len(field) != 3:
        raise ValueError("model.field_xyz must have three components")
    for item in field:
        _number(item, "model.field_xyz")

    geometry = _object(document["geometry"], "geometry")
    _keys(geometry, {"cell", "boundary_conditions", "k_grid_shape"}, "geometry")
    if not _object(geometry["cell"], "geometry.cell"):
        raise ValueError("geometry.cell must describe the cell rather than name a file")
    _text(geometry["boundary_conditions"], "geometry.boundary_conditions")
    shape = geometry["k_grid_shape"]
    if not isinstance(shape, list) or len(shape) != 2:
        raise ValueError("geometry.k_grid_shape must contain two dimensions")
    for item in shape:
        _positive_int(item, "geometry.k_grid_shape")

    numerics = _object(document["numerics"], "numerics")
    _keys(
        numerics,
        {"precision", "stopping_rule", "residual_tolerances", "max_map_evaluations"},
        "numerics",
    )
    _text(numerics["precision"], "numerics.precision")
    _text(numerics["stopping_rule"], "numerics.stopping_rule")
    _mapping_of_numbers(
        numerics["residual_tolerances"], "numerics.residual_tolerances", positive=True
    )
    _positive_int(numerics["max_map_evaluations"], "numerics.max_map_evaluations")

    sampling = _object(document["sampling"], "sampling")
    _keys(
        sampling,
        {"population_id", "generator", "version", "parameters", "measure_definition"},
        "sampling",
    )
    if sampling["population_id"] != "algebraic_unverified":
        raise ValueError("physical Gaussian or symmetry population has not been certified")
    for key in ("generator", "version", "measure_definition"):
        _text(sampling[key], f"sampling.{key}")
    _object(sampling["parameters"], "sampling.parameters")

    provenance = _object(document["provenance"], "provenance")
    _keys(
        provenance,
        {
            "git_sha",
            "git_dirty",
            "python_version",
            "numpy_version",
            "torch_version",
            "platform",
            "hardware",
            "blas_backend",
            "threads",
        },
        "provenance",
    )
    git_sha = provenance["git_sha"]
    if (
        not isinstance(git_sha, str)
        or len(git_sha) != 40
        or any(c not in "0123456789abcdef" for c in git_sha)
    ):
        raise ValueError("provenance.git_sha must be a full lower-case Git SHA")
    if type(provenance["git_dirty"]) is not bool:
        raise ValueError("provenance.git_dirty must be a boolean")
    for key in (
        "python_version",
        "numpy_version",
        "torch_version",
        "platform",
        "hardware",
        "blas_backend",
    ):
        _text(provenance[key], f"provenance.{key}")
    _positive_int(provenance["threads"], "provenance.threads")


def _validate_solver(document: dict[str, Any]) -> None:
    _keys(document, {"method", "implementation", "parameters"}, "solver")
    _text(document["method"], "solver.method")
    _text(document["implementation"], "solver.implementation")
    _object(document["parameters"], "solver.parameters")


def _canonical(document: Mapping[str, Any]) -> str:
    try:
        return json.dumps(document, sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError("configuration must be finite, JSON-serializable data") from exc


@dataclass(frozen=True)
class ExperimentContract:
    """Canonical self-contained model, sampling, numerics and runtime provenance."""

    document_json: str

    def __post_init__(self) -> None:
        try:
            document = _object(json.loads(self.document_json), "experiment")
        except (json.JSONDecodeError, TypeError) as exc:
            raise ValueError("malformed experiment contract") from exc
        _validate_experiment(document)
        if _canonical(document) != self.document_json:
            raise ValueError("experiment contract must use canonical JSON")

    @classmethod
    def from_mapping(cls, document: Mapping[str, Any]) -> ExperimentContract:
        return cls(_canonical(document))

    @property
    def digest(self) -> str:
        return hashlib.sha256(self.document_json.encode("utf-8")).hexdigest()

    @property
    def document(self) -> dict[str, Any]:
        return _object(json.loads(self.document_json), "experiment")

    @property
    def source_commit(self) -> str:
        result: str = self.document["provenance"]["git_sha"]
        return result

    @property
    def population_id(self) -> str:
        result: str = self.document["sampling"]["population_id"]
        return result

    @property
    def measure_definition(self) -> str:
        result: str = self.document["sampling"]["measure_definition"]
        return result


@dataclass(frozen=True)
class SolverConfiguration:
    """Concrete per-method settings; distinct algorithms may have distinct values."""

    document_json: str

    def __post_init__(self) -> None:
        try:
            document = _object(json.loads(self.document_json), "solver")
        except (json.JSONDecodeError, TypeError) as exc:
            raise ValueError("malformed solver configuration") from exc
        _validate_solver(document)
        if _canonical(document) != self.document_json:
            raise ValueError("solver configuration must use canonical JSON")

    @classmethod
    def from_mapping(cls, document: Mapping[str, Any]) -> SolverConfiguration:
        return cls(_canonical(document))

    @property
    def method(self) -> str:
        result: str = _object(json.loads(self.document_json), "solver")["method"]
        return result

    @property
    def digest(self) -> str:
        return hashlib.sha256(self.document_json.encode("utf-8")).hexdigest()
