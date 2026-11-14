"""Synthetic-only regression tests for Wave 4 paired benchmark infrastructure."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.domain.wave4_trials import (
    BranchDecision,
    PairedCase,
    TrialObservation,
    TrialRecord,
    compare_paired_trials,
    make_trial_record,
    summarize_trials,
)
from kitaev_mf.infra.wave4_trial_io import read_trial_records, write_trial_records

SOURCE_SHA = "a" * 40


def _state(value: float = 0.0) -> MeanFieldState:
    u = np.zeros((3, 4, 4), dtype=np.float64)
    u[0, 0, 0] = value
    return MeanFieldState(
        U=u,
        V=np.zeros((2, 4, 4), dtype=np.float64),
        lambdas=np.zeros((2, 3), dtype=np.float64),
    )


def _case(case_id: str = "c0", value: float = 0.0) -> PairedCase:
    return PairedCase.from_state(
        protocol_id="synthetic-protocol",
        case_id=case_id,
        stratum_id="synthetic-stratum",
        cohort="synthetic",
        measure_ref="synthetic-measure-only",
        state=_state(value),
    )


class _Oracle:
    def __init__(self, label: str) -> None:
        self.label = label

    def classify(self, state: MeanFieldState, case: PairedCase) -> BranchDecision:
        assert isinstance(state, MeanFieldState)
        assert case.population_id == "algebraic_unverified"
        if self.label == "target":
            return BranchDecision("target", "synthetic-oracle-v1")
        if self.label == "other":
            return BranchDecision("other", "synthetic-oracle-v1")
        return BranchDecision("unknown")


def _trial(
    case: PairedCase,
    method: str = "diis",
    *,
    converged: bool = True,
    label: str | None = None,
    cost: int = 5,
) -> TrialRecord:
    result: TrialRecord = make_trial_record(
        case=case,
        method=method,
        config_ref="synthetic-config",
        source_commit=SOURCE_SHA,
        observation=TrialObservation(
            converged=converged,
            iterations=cost,
            map_evaluations=cost,
            wall_seconds=0.1 * cost,
            stop_reason="synthetic-terminal",
            residual_history=(0.2, 0.1),
            convergence_check_ref="synthetic-convergence-gate" if converged else None,
            final_state=case.materialize_state() if converged else None,
        ),
        oracle=_Oracle(label) if label is not None else None,
    )
    return result


def test_identical_snapshot_for_each_algorithm_and_caller_mutation() -> None:
    source = _state(0.125)
    case = PairedCase.from_state(
        protocol_id="synthetic-protocol",
        case_id="seed-1",
        stratum_id="test",
        cohort="synthetic",
        measure_ref="test-measure",
        state=source,
    )
    first = case.materialize_state()
    second = case.materialize_state()
    np.testing.assert_array_equal(first.U, second.U)
    np.testing.assert_array_equal(first.V, second.V)
    np.testing.assert_array_equal(first.lambdas, second.lambdas)
    source.U[0, 0, 0] = 999.0
    first.U[0, 0, 0] = -999.0
    assert second.U[0, 0, 0] == 0.125
    assert case.materialize_state().U[0, 0, 0] == 0.125
    assert _case("x", 0.1).initial_state_hash != _case("x", 0.2).initial_state_hash
    assert _case("x", 0.1).initial_state_hash == _case("y", 0.1).initial_state_hash


def test_reject_nonfinite_dtype_and_mutated_skew() -> None:
    bad = _state()
    bad.U[0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="nonfinite"):
        PairedCase.from_state(
            protocol_id="p", case_id="a", stratum_id="s",
            cohort="synthetic", measure_ref="m", state=bad
        )
    source = _state()
    bad = MeanFieldState(U=source.U.astype(np.float32), V=source.V, lambdas=source.lambdas)
    with pytest.raises(ValueError, match="float64"):
        PairedCase.from_state(
            protocol_id="p", case_id="a", stratum_id="s",
            cohort="synthetic", measure_ref="m", state=bad
        )
    bad = _state()
    bad.V[0, 0, 1] = 1.0
    with pytest.raises(ValueError, match="skew-symmetric"):
        PairedCase.from_state(
            protocol_id="p", case_id="a", stratum_id="s",
            cohort="synthetic", measure_ref="m", state=bad
        )


def test_three_way_classification_and_identification_bounds() -> None:
    records = [
        _trial(_case("target"), label="target", cost=2),
        _trial(_case("other"), label="other", cost=6),
        _trial(_case("unresolved"), label=None, cost=11),
        _trial(_case("failed"), converged=False, label="target", cost=17),
    ]
    summary = summarize_trials(records)
    assert (summary.n, summary.n_converged, summary.n_target) == (4, 3, 1)
    assert (summary.n_other, summary.n_unknown_converged) == (1, 1)
    assert summary.p_converged == 0.75
    assert summary.p_target_lower == 0.25
    assert summary.p_target_upper == 0.50
    assert summary.p_same_given_converged_lower == pytest.approx(1 / 3)
    assert summary.p_same_given_converged_upper == pytest.approx(2 / 3)
    assert summary.unknown_given_converged == pytest.approx(1 / 3)
    assert summary.mean_actual_map_evaluations == 9.0
    assert summary.mean_actual_wall_seconds == pytest.approx(0.9)
    assert records[3].branch == "unknown"


def test_no_convergence_keeps_conditional_quantities_undefined() -> None:
    records = [_trial(_case("a"), converged=False, cost=12)]
    summary = summarize_trials(records)
    assert summary.p_converged == 0.0
    assert summary.p_target_lower == summary.p_target_upper == 0.0
    assert summary.p_same_given_converged_lower is None
    assert summary.p_same_given_converged_upper is None
    assert summary.unknown_given_converged is None
    assert summary.mean_actual_map_evaluations == 12.0


def test_nonconverged_and_unproven_classification_rejected() -> None:
    record = _trial(_case("a"), label="target")
    with pytest.raises(ValueError, match="nonconverged"):
        replace(record, converged=False)
    with pytest.raises(ValueError, match="oracle reference"):
        replace(record, oracle_ref=None)
    with pytest.raises(ValueError, match="convergence"):
        replace(record, convergence_check_ref=None)
    with pytest.raises(ValueError, match="population"):
        replace(record, population_id="gaussian_certified")
    with pytest.raises(ValueError, match="negative"):
        replace(record, map_evaluations=-1)
    with pytest.raises(ValueError, match="nonnegative"):
        replace(record, wall_seconds=float("nan"))


def test_paired_convergence_and_actual_failure_cost() -> None:
    cases = [_case(f"c{i}", value=float(i)) for i in range(3)]
    a = [
        _trial(cases[0], "fp", converged=True, cost=12),
        _trial(cases[1], "fp", converged=True, cost=14),
        _trial(cases[2], "fp", converged=False, cost=30),
    ]
    b = [
        _trial(cases[2], "jfnk", converged=True, cost=35),
        _trial(cases[0], "jfnk", converged=False, cost=50),
        _trial(cases[1], "jfnk", converged=True, cost=16),
    ]
    comparison = compare_paired_trials(a, b)
    assert comparison.n_pairs == 3
    assert comparison.convergence_difference_a_minus_b == 0.0
    assert comparison.a_only_converged == 1
    assert comparison.b_only_converged == 1
    assert comparison.mean_map_cost_difference_a_minus_b == pytest.approx(-15.0)


def test_unmatched_initial_state_and_duplicate_cases_are_rejected() -> None:
    x = _trial(_case("same", 0.0), "fp")
    y = _trial(_case("same", 1.0), "jfnk")
    with pytest.raises(ValueError, match="pairing"):
        compare_paired_trials([x], [y])
    with pytest.raises(ValueError, match="duplicate"):
        compare_paired_trials([x, x], [y, y])
    with pytest.raises(ValueError, match="different case sets"):
        compare_paired_trials([x], [_trial(_case("other"), "jfnk")])
    with pytest.raises(ValueError, match="duplicate"):
        summarize_trials([x, x])


def test_jsonl_roundtrip_and_create_only(tmp_path: Path) -> None:
    records = [
        _trial(_case("case-a"), method="fp", label="target", cost=4),
        _trial(_case("case-a"), method="jfnk", converged=False, cost=21),
    ]
    output = tmp_path / "wave4-synthetic.jsonl"
    write_trial_records(output, records)
    assert read_trial_records(output) == records
    original = output.read_bytes()
    with pytest.raises(FileExistsError):
        write_trial_records(output, records)
    assert output.read_bytes() == original
    with pytest.raises(ValueError, match="extension"):
        write_trial_records(tmp_path / "legacy-scorecard.dat", records)


def test_jsonl_rejects_partial_bad_schema_duplicates_and_incorrect_types(tmp_path: Path) -> None:
    record = _trial(_case("one"), label="target")
    output = tmp_path / "scratch.jsonl"
    write_trial_records(output, [record])
    original_line = output.read_text(encoding="utf-8").strip()
    document = json.loads(original_line)

    document["record"].pop("oracle_ref")
    output.write_text(json.dumps(document) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing"):
        read_trial_records(output)

    document = json.loads(original_line)
    document["schema_version"] = 999
    output.write_text(json.dumps(document) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="schema"):
        read_trial_records(output)

    document = json.loads(original_line)
    document["record"]["converged"] = "false"
    output.write_text(json.dumps(document) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="boolean"):
        read_trial_records(output)

    output.write_text(original_line + "\n" + original_line + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate"):
        read_trial_records(output)

    output.write_text(original_line[:50] + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="malformed"):
        read_trial_records(output)

    output.write_text(original_line + "\n\n", encoding="utf-8")
    with pytest.raises(ValueError, match="empty line"):
        read_trial_records(output)


def test_undefined_branch_oracle_does_not_invent_target() -> None:
    case = _case("a")
    unknown = _trial(case)
    assert unknown.branch == "unknown"
    assert unknown.oracle_ref is None
    resolved = _trial(case, label="other")
    assert resolved.branch == "other"
    assert resolved.oracle_ref == "synthetic-oracle-v1"
