# Rule 05: Project Management — Phases, Milestones, Tasks, and Sprints

## 1. Purpose and planning hierarchy

Use a lightweight, Scrum-inspired process to coordinate implementation and research. The repository is not required to reproduce the full Scrum framework; ceremonies, roles, and artifacts should remain proportional to a small research-software project.

The planning hierarchy is:

Phase → Milestone → Task

- **Phase** (Phase N): a macro-horizon describing a physical or engineering regime. A phase groups related strategic outcomes and may span multiple releases.
- **Milestone** (M<phase>.<index>): a persistent strategic backlog item and acceptance gate. A milestone states its objective, scope, assumptions, dependencies, target deliverables, and acceptance criteria. A nested milestone may use an additional suffix, e.g. M1.4.5, when explicitly subordinate to a parent milestone; do not use extra numeric components to encode tasks or batches.
- **Task** (T<index>): a short-lived execution unit belonging only to the active milestone's current sprint. Task numbering restarts at T1 when a new milestone becomes active. A task must have a concrete outcome, exact target files or interfaces, analytical/engineering contracts, and verification commands.

Milestones persist in the roadmap. Tasks are not a second strategic backlog: do not pre-allocate task IDs or granular implementation plans to future milestones.

## 2. Roadmap and project-state responsibilities

- doc/MILESTONE_ROADMAP.md is the persistent product/research backlog. It records phases, milestone objectives, dependencies, and acceptance gates.
- PROJECT_STATE.yaml is the concise current-state record. It identifies the active phase and milestone, the active milestone's current tasks, completed milestone evidence, verified invariants, and blockers.
- R&D notes contain derivations, exploratory calculations, failed approaches, and scientific qualifications. A roadmap entry is not evidence that a scientific claim has been proved.
- Tests and reproducibility artifacts provide evidence only for their actual scope. Keep analytical proof, numerical tests, finite-size evidence, and approximations distinct.

Do not duplicate full task descriptions in both files. The roadmap should identify the active milestone and link to the state record for its current task breakdown. The state record should not reproduce the entire future roadmap.

## 3. Milestone contract

Every roadmap milestone must define:

1. Stable ID, phase, title, and status.
2. Objective and physical/software scope.
3. Assumptions, conventions, dependencies, and non-goals where relevant.
4. Target deliverables or repository areas. Mark paths that do not yet exist as planned targets.
5. Acceptance criteria that are observable and independently reviewable.
6. Verification evidence when closed, including exact commands, relevant tests, analytical records, and validation provenance.

Use statuses consistently: backlog, ready, active, blocked, and completed. A milestone is completed only after its acceptance gate is reviewed. Do not promote a planned numerical benchmark or conjectural identity into a verified invariant merely by listing it as a criterion.

Future milestones must remain at milestone granularity. Do not add numbered tasks, sub-batches, or implementation schedules until that milestone becomes active.

## 4. Active sprint and task contract

The active milestone is the current sprint-level execution focus. It may require more than one time-boxed work cycle, but only the currently selected execution tasks are itemized in PROJECT_STATE.yaml.

Each active task must include:
- ID (T1, T2, ...), concise title, and status;
- exact target files or interfaces, distinguishing existing files from files to be created;
- contracts or questions to resolve, with stable contract IDs where useful;
- validation commands and expected evidence;
- dependencies and material assumptions.

Keep tasks independently reviewable and small enough to produce a verifiable outcome. Split a task when it mixes unrelated scientific contracts or has no coherent acceptance test. Do not encode performance budgets as hard pass/fail criteria unless the hardware, precision, workload, and measurement method are specified.

When a milestone closes:
1. Move its acceptance evidence into the milestone's closure record in the roadmap and/or a linked audit/R&D note.
2. Remove its active task list from PROJECT_STATE.yaml.
3. Reset task numbering for the next active milestone.
4. Preserve completed milestone IDs and evidence; do not preserve obsolete task IDs as a second backlog.

## 5. Lightweight Scrum-inspired operating model

- The **product backlog** is the ordered milestone roadmap. The project owner determines priorities and accepts scope changes.
- **Sprint planning** selects the active milestone goal and a small set of tasks with explicit acceptance criteria.
- The **sprint goal** is the milestone-level outcome, not a list of unrelated activities.
- **Execution review** checks task evidence, test results, blockers, and scope changes before tasks are marked complete.
- **Retrospective** records process changes only when they are likely to improve future work; do not create ceremony documents by default.
- The **Scrum Master function** is facilitation: keep the backlog coherent, expose dependencies and blockers, protect scope, and ensure evidence is recorded. It does not override the project owner's scientific or prioritization decisions.

For a solo workflow, these are responsibilities rather than mandatory separate human roles. Do not claim Scrum compliance merely because the terminology is used.

## 6. Change discipline and branch ownership

- Update the active task status and project-state record in the same commit as a material implementation, test, scope, or verification change on the milestone branch.
- Keep PROJECT_STATE.yaml and doc/MILESTONE_ROADMAP.md mutually consistent when a milestone opens, changes scope, or closes.
- Preserve unrelated user changes and do not report unrun local checks as passed.
- Do not copy an in-progress milestone state from one development branch to another. State and roadmap updates belong to the branch that owns that milestone until reviewed integration.
- Governance-policy changes belong on governance/project-governance; ordinary implementation and research milestones use their dedicated development branches and integrate through reviewed pull requests.

For multi-file state/roadmap updates, create one Git commit containing all related files. Inspect the branch head immediately before advancing its ref and use an expected-head check; do not split a logically atomic update into independent file commits.
