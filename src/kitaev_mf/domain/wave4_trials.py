"""Wave 4 paired-case identities and descriptive trial estimands.

No physical sampler, numerical threshold, or branch oracle is implemented here.
"""

from __future__ import annotations

import hashlib
import math
from collections.abc import Sequence
from dataclasses import dataclass
from statistics import mean
from typing import Literal, Protocol

import numpy as np
from numpy.typing import NDArray

from kitaev_mf.domain.state import MeanFieldState

BranchLabel = Literal["target", "other", "unknown"]
Cohort = Literal["synthetic", "pilot", "confirmatory"]


def _require_nonempty(value: str, field_name: str) -> None:
    if not value or not value.strip():
        raise ValueError(f"{field_name} must be nonempty")


@dataclass(frozen=True)
class PairedCase:
    """Immutable snapshot of an algebraically checked initial mean-field state.

    This constructor does not certify Gaussian representability or C3 invariance.
    Every method obtains a fresh copy from materialize_state().
    """

    protocol_id: str
    case_id: str
    stratum_id: str
    cohort: Cohort
    measure_ref: str
    initial_state_hash: str
    population_id: str
    _shapes: tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]
    _payload: tuple[bytes, bytes, bytes]

    @classmethod
    def from_state(
        cls,
        *,
        protocol_id: str,
        case_id: str,
        stratum_id: str,
        cohort: Cohort,
        measure_ref: str,
        state: MeanFieldState,
    ) -> PairedCase:
        for field_name, field_value in (
            ("protocol_id", protocol_id),
            ("case_id", case_id),
            ("stratum_id", stratum_id),
            ("measure_ref", measure_ref),
        ):
            _require_nonempty(field_value, field_name)
        if cohort not in ("synthetic", "pilot", "confirmatory"):
            raise ValueError("invalid cohort")

        arrays: list[NDArray[np.float64]] = []
        digest = hashlib.sha256(b"wave4-algebraic-initial-state-v1")
        shapes: list[tuple[int, ...]] = []
        payload: list[bytes] = []
        for name, value in (("U", state.U), ("V", state.V), ("lambdas", state.lambdas)):
            array = np.asarray(value)
            if array.dtype != np.dtype("float64"):
                raise ValueError(f"{name} must use float64 (no implicit conversion)")
            if not np.isfinite(array).all():
                raise ValueError(f"{name} contains nonfinite values")
            frozen = np.array(array, dtype=np.dtype("<f8"), order="C", copy=True)
            shape = tuple(int(dim) for dim in frozen.shape)
            raw = frozen.tobytes(order="C")
            digest.update(name.encode("ascii"))
            digest.update(len(shape).to_bytes(1, "little"))
            for dim in shape:
                digest.update(dim.to_bytes(8, "little"))
            digest.update(raw)
            arrays.append(frozen)
            shapes.append(shape)
            payload.append(raw)

        # Validate a snapshot, not just the possibly mutated caller-owned arrays.
        MeanFieldState(U=arrays[0], V=arrays[1], lambdas=arrays[2])
        return cls(
            protocol_id=protocol_id,
            case_id=case_id,
            stratum_id=stratum_id,
            cohort=cohort,
            measure_ref=measure_ref,
            initial_state_hash=digest.hexdigest(),
            population_id="algebraic_unverified",
            _shapes=(shapes[0], shapes[1], shapes[2]),
            _payload=(payload[0], payload[1], payload[2]),
        )

    def materialize_state(self) -> MeanFieldState:
        arrays = [
            np.frombuffer(raw, dtype=np.dtype("<f8")).reshape(shape).copy()
            for raw, shape in zip(self._payload, self._shapes, strict=True)
        ]
        return MeanFieldState(U=arrays[0], V=arrays[1], lambdas=arrays[2])


@dataclass(frozen=True)
class BranchDecision:
    """A decision supplied by an independently specified branch oracle."""

    label: BranchLabel
    oracle_ref: str | None = None

    def __post_init__(self) -> None:
        if self.label not in ("target", "other", "unknown"):
            raise ValueError("invalid branch label")
        if self.label != "unknown":
            if self.oracle_ref is None:
                raise ValueError("known branch requires independent oracle provenance")
            _require_nonempty(self.oracle_ref, "oracle_ref")


class BranchOracle(Protocol):
    """Interface only; neither a reference branch nor its equivalence is assumed."""

    def classify(self, state: MeanFieldState, case: PairedCase) -> BranchDecision: ...


@dataclass(frozen=True)
class TrialObservation:
    """Telemetry supplied by an external runner, including work on failures."""

    converged: bool
    iterations: int
    map_evaluations: int
    wall_seconds: float
    stop_reason: str
    residual_history: tuple[float, ...]
    convergence_check_ref: str | None
    final_state: MeanFieldState | None = None


@dataclass(frozen=True)
class TrialRecord:
    """A single method/case outcome; unknown is never silently classified as other."""

    protocol_id: str
    case_id: str
    stratum_id: str
    cohort: Cohort
    measure_ref: str
    population_id: str
    initial_state_hash: str
    method: str
    config_ref: str
    source_commit: str
    converged: bool
    branch: BranchLabel
    oracle_ref: str | None
    iterations: int
    map_evaluations: int
    wall_seconds: float
    stop_reason: str
    residual_history: tuple[float, ...]
    convergence_check_ref: str | None

    def __post_init__(self) -> None:
        for name in (
            "protocol_id",
            "case_id",
            "stratum_id",
            "measure_ref",
            "population_id",
            "initial_state_hash",
            "method",
            "config_ref",
            "source_commit",
            "stop_reason",
        ):
            _require_nonempty(getattr(self, name), name)
        if self.cohort not in ("synthetic", "pilot", "confirmatory"):
            raise ValueError("invalid cohort")
        if self.population_id != "algebraic_unverified":
            raise ValueError("no symmetry/Gaussian certification is implemented")
        if self.branch not in ("target", "other", "unknown"):
            raise ValueError("invalid branch label")
        if len(self.initial_state_hash) != 64 or any(
            c not in "0123456789abcdef" for c in self.initial_state_hash
        ):
            raise ValueError("initial_state_hash must be a SHA256 hex digest")
        if self.iterations < 0 or self.map_evaluations < 0:
            raise ValueError("actual iteration and map counts cannot be negative")
        if not math.isfinite(self.wall_seconds) or self.wall_seconds < 0:
            raise ValueError("wall_seconds must be finite and nonnegative")
        if any(not math.isfinite(r) or r < 0 for r in self.residual_history):
            raise ValueError("residual_history must contain finite nonnegative values")
        if not self.converged and self.branch != "unknown":
            raise ValueError("nonconverged run cannot have an identified branch")
        if self.converged and (
            self.convergence_check_ref is None or not self.convergence_check_ref.strip()
        ):
            raise ValueError("convergence requires a validation reference")
        if self.branch != "unknown" and (self.oracle_ref is None or not self.oracle_ref.strip()):
            raise ValueError("classified branch requires an oracle reference")


def make_trial_record(
    *,
    case: PairedCase,
    method: str,
    config_ref: str,
    source_commit: str,
    observation: TrialObservation,
    oracle: BranchOracle | None = None,
) -> TrialRecord:
    """Build a record without executing a solver or choosing physical tolerances."""
    decision = BranchDecision("unknown")
    if observation.converged and oracle is not None and observation.final_state is not None:
        decision = oracle.classify(observation.final_state, case)
    return TrialRecord(
        protocol_id=case.protocol_id,
        case_id=case.case_id,
        stratum_id=case.stratum_id,
        cohort=case.cohort,
        measure_ref=case.measure_ref,
        population_id=case.population_id,
        initial_state_hash=case.initial_state_hash,
        method=method,
        config_ref=config_ref,
        source_commit=source_commit,
        converged=observation.converged,
        branch=decision.label if observation.converged else "unknown",
        oracle_ref=decision.oracle_ref if observation.converged else None,
        iterations=observation.iterations,
        map_evaluations=observation.map_evaluations,
        wall_seconds=observation.wall_seconds,
        stop_reason=observation.stop_reason,
        residual_history=observation.residual_history,
        convergence_check_ref=observation.convergence_check_ref,
    )


@dataclass(frozen=True)
class DescriptiveSummary:
    n: int
    n_converged: int
    n_target: int
    n_other: int
    n_unknown_converged: int
    p_converged: float
    p_target_lower: float
    p_target_upper: float
    p_same_given_converged_lower: float | None
    p_same_given_converged_upper: float | None
    unknown_given_converged: float | None
    mean_actual_map_evaluations: float
    mean_actual_wall_seconds: float


def summarize_trials(records: Sequence[TrialRecord]) -> DescriptiveSummary:
    """Descriptive finite-sample fractions and identification bounds, not confidence."""
    if not records:
        raise ValueError("cannot summarize an empty cohort")
    keys = {
        (r.protocol_id, r.stratum_id, r.cohort, r.measure_ref, r.population_id, r.method)
        for r in records
    }
    if len(keys) != 1:
        raise ValueError("summary cannot mix strata, cohorts, populations or methods")
    if len({r.case_id for r in records}) != len(records):
        raise ValueError("duplicate case_id in a method arm")
    n = len(records)
    n_conv = sum(r.converged for r in records)
    n_target = sum(r.converged and r.branch == "target" for r in records)
    n_other = sum(r.converged and r.branch == "other" for r in records)
    n_unknown = n_conv - n_target - n_other
    return DescriptiveSummary(
        n=n,
        n_converged=n_conv,
        n_target=n_target,
        n_other=n_other,
        n_unknown_converged=n_unknown,
        p_converged=n_conv / n,
        p_target_lower=n_target / n,
        p_target_upper=(n_target + n_unknown) / n,
        p_same_given_converged_lower=n_target / n_conv if n_conv else None,
        p_same_given_converged_upper=(n_target + n_unknown) / n_conv if n_conv else None,
        unknown_given_converged=n_unknown / n_conv if n_conv else None,
        mean_actual_map_evaluations=mean(r.map_evaluations for r in records),
        mean_actual_wall_seconds=mean(r.wall_seconds for r in records),
    )


@dataclass(frozen=True)
class PairedComparison:
    n_pairs: int
    convergence_difference_a_minus_b: float
    mean_map_cost_difference_a_minus_b: float
    a_only_converged: int
    b_only_converged: int


def compare_paired_trials(
    arm_a: Sequence[TrialRecord], arm_b: Sequence[TrialRecord]
) -> PairedComparison:
    """Match exact initial snapshots, without pretending paired arms are iid."""
    if not arm_a or not arm_b:
        raise ValueError("paired comparison requires two nonempty arms")
    if len({r.case_id for r in arm_a}) != len(arm_a):
        raise ValueError("duplicate case_id in first arm")
    if len({r.case_id for r in arm_b}) != len(arm_b):
        raise ValueError("duplicate case_id in second arm")
    a_by_id = {r.case_id: r for r in arm_a}
    b_by_id = {r.case_id: r for r in arm_b}
    if a_by_id.keys() != b_by_id.keys():
        raise ValueError("arms have different case sets")
    if len({r.method for r in arm_a}) != 1 or len({r.method for r in arm_b}) != 1:
        raise ValueError("each arm must contain one method")
    if arm_a[0].method == arm_b[0].method:
        raise ValueError("paired comparison requires distinct methods")
    pairs = [(a_by_id[k], b_by_id[k]) for k in sorted(a_by_id)]
    for a, b in pairs:
        if (
            a.initial_state_hash != b.initial_state_hash
            or a.protocol_id != b.protocol_id
            or a.stratum_id != b.stratum_id
            or a.cohort != b.cohort
            or a.measure_ref != b.measure_ref
            or a.population_id != b.population_id
            or a.method != arm_a[0].method
            or b.method != arm_b[0].method
        ):
            raise ValueError("pairing/protocol mismatch")
    return PairedComparison(
        n_pairs=len(pairs),
        convergence_difference_a_minus_b=mean(
            int(a.converged) - int(b.converged) for a, b in pairs
        ),
        mean_map_cost_difference_a_minus_b=mean(
            a.map_evaluations - b.map_evaluations for a, b in pairs
        ),
        a_only_converged=sum(a.converged and not b.converged for a, b in pairs),
        b_only_converged=sum(b.converged and not a.converged for a, b in pairs),
    )
