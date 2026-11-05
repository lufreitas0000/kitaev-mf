# Rule 05: Project Management — Phases, Milestones, Tasks, and Sprints

## 1. Purpose and planning hierarchy

Use a lightweight, Scrum-inspired process to coordinate implementation and research. The repository is not required to reproduce the full Scrum framework; ceremonies, roles, and artifacts should remain proportional to a small research-software project.

The planning hierarchy is:

Phase → Milestone → Task

- **Phase** (Phase N): A macro-horizon describing a physical or engineering regime (e.g. Phase 1: Analytical Rigor & Domain Foundations).
- **Milestone** (M<phase>.<index>): A persistent strategic deliverable gate. A milestone defines its objective, physical/software scope, assumptions, dependencies, target deliverables, and acceptance criteria.
- **Task** (T<index>): A short-lived operational execution unit active exclusively within the current milestone sprint. Task numbering restarts at T1 when a new milestone opens.

## 2. Roadmap and project-state responsibilities

- `doc/MILESTONE_ROADMAP.md` is the persistent backlog. It records phases, milestones, dependencies, and acceptance criteria.
- `PROJECT_STATE.yaml` is the concise current-state record. It identifies the active phase, current milestone, active sprint tasks, completed milestone evidence, verified invariants, and blockers.
- `doc/rnd/phase_<N>/` contains the Tier-1 machine-first discovery records.
- `doc/labbook/` contains the Tier-2 human-readable living lab book.

## 3. Milestone completion contract

A milestone cannot be transitioned to `completed` in `PROJECT_STATE.yaml` and `doc/MILESTONE_ROADMAP.md` until:
1. All planned task deliverables exist and relevant automated tests pass 100% under `make verify`.
2. All corresponding R&D notes in `doc/rnd/phase_<N>/` are filed under deterministic filenames (`m<phase>_<milestone>_t<task>_<slug>.md`) and carry valid YAML metadata headers.
3. The corresponding chapter in the living labbook (`doc/labbook/chapters/ch_m<phase>_<milestone>_*.tex`) is updated with derivations, failure modes, and embedded `.dat` benchmark tables.
4. `make labbook` compiles `doc/labbook/main.pdf` cleanly with zero TeX errors.
5. Numerical benchmarks are exported to `data/processed/benchmarks/*.dat`.
