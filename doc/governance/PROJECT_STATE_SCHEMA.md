# PROJECT_STATE.yaml Schema and Field Semantics

## Purpose

PROJECT_STATE.yaml is a compact, machine-readable snapshot of durable project state and the currently active milestone execution. It is not a Git/CI dashboard and is not a replacement for the milestone roadmap, R&D notes, or test evidence.

The canonical planning policy is .agents/rules/05_project_management.md.

## Top-level fields

| Field | Required | Meaning |
|---|---:|---|
| schema_version | yes | Version of this YAML structure. Start at integer 1. |
| project_name | yes | Stable repository/project identifier. |
| active_phase | yes | Current macro phase, using the roadmap's phase label. |
| current_milestone | yes | ID of the active milestone, e.g. M1.4.5. |
| milestone_history | yes | Completed milestone summaries with evidence references. |
| active_milestone | yes | Objective, status, sprint goal, and current execution tasks for the active milestone. |
| hardware_targets | optional | Relevant known development/compute targets; not a machine-status claim. |
| verified_physics_invariants | optional | Invariants with evidence provenance and scope. |
| open_blockers | yes | Current unresolved blockers; use an empty list when none are known. |

Package version remains authoritative in pyproject.toml; do not duplicate it here.

## Milestone history

Each milestone_history entry should contain:
- id, name, and status: completed;
- a concise acceptance summary;
- evidence paths to tests, source, analytical notes, or audit records;
- validation provenance where relevant (analytical, automated, numerical, or user-reported).

Do not retain task-level lists in closed milestone records. Preserve evidence references, not obsolete execution bookkeeping.

## Active milestone

The active_milestone object must have:
- id, matching current_milestone;
- description and objective;
- status (normally active, blocked, or awaiting_local_validation);
- sprint_goal;
- active_tasks, a list of task objects.

Each active task has:
- id: T1, T2, ...; unique within the active milestone only;
- title;
- status: planned, in_progress, blocked, awaiting_local_validation, or completed;
- target_files: exact paths, with planned-to-create paths identified as such;
- contracts: stable contract IDs and concise descriptions;
- validation: exact commands or a clearly marked command to be added with the target implementation.

A task marked completed must have evidence. A task marked awaiting_local_validation must not be reported as fully verified.

## Prohibited transient fields

Do not store the current Git branch, working-tree cleanliness, current HEAD, current CI status, or unqualified "last test passed" claims. These change independently of durable project state. Record them only in a dated, commit-pinned audit report when historically relevant.

## Invariants and scientific provenance

Each invariant should identify its mathematical statement, status, assumptions, and provenance. Distinguish:
- verified_analytical: a derivation/proof is recorded and reviewed;
- verified_pytest: automated tests cover the stated finite computational scope;
- analytical_formula_with_finite_mesh_regression: an analytical expression is paired with finite-grid numerical evidence;
- user_reported_local_validation: a local run was reported by the user, with commands/results when available.

Do not infer analytical proof from a passing test. Do not infer thermodynamic-limit accuracy from a finite mesh without a controlled convergence argument.

## Consistency checks

Before committing a state update, verify:
1. current_milestone == active_milestone.id.
2. The active milestone appears in the roadmap with matching phase and status.
3. Active task IDs are unique and sequential from T1 unless a documented reason requires otherwise.
4. Closed milestones have no active task list.
5. Every completed claim has a resolvable evidence path.
6. No transient Git/CI dashboard fields remain.
7. YAML parses successfully and git diff --check is clean.
