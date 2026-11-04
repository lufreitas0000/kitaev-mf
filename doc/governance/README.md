# Repository Governance Architecture

**Status:** Initial policy baseline  
**Scope:** Repository organization, document authority, and governance workflow  
**Applies to:** Human contributors and automated agents

This document defines how repository guidance, executable configuration, scientific records, and project planning relate to one another. It establishes document responsibilities; it does not certify the correctness of any scientific result or imply that planned directories and workflows already exist.

## 1. Governance principles

1. **One authoritative home per kind of information.** Avoid maintaining the same policy, configuration value, or project status independently in several files.
2. **Executable configuration takes precedence for tool behavior.** The actual behavior of CI, linters, type checkers, tests, packaging, and build commands is determined by their checked-in configuration and workflows, not by prose summaries.
3. **Rules describe required practice.** A rule must be actionable, testable where practical, and scoped to the repository's current needs.
4. **Decisions preserve rationale.** Record durable architectural choices and their trade-offs separately from operational instructions.
5. **State records evidence, not aspiration.** A completed milestone must be supported by identifiable artifacts and verification; plans and hypotheses must not be presented as completed facts.
6. **Exploration remains distinguishable from accepted results.** Research notes may contain incomplete derivations, failed approaches, numerical observations, and conjectures.
7. **Governance should remain proportional.** Add a document, workflow, or abstraction only when it resolves a real coordination or maintenance problem.

## 2. Document authority and precedence

Different repository artifacts answer different questions. They should not be treated as interchangeable.

| Artifact | Authoritative for | Not authoritative for |
|---|---|---|
| GitHub repository settings and branch protections | Server-enforced permissions and merge restrictions | Scientific correctness or local workflow preferences |
| CI workflows and tool configuration | Checks and automated behavior actually executed | Requirements not encoded or checked |
| `pyproject.toml`, requirements inputs and lockfiles | Package metadata, dependency constraints, and configured Python tooling | Scientific derivations or project priorities |
| `Makefile` and maintained scripts | Named local development commands and their concrete behavior | Requirements that the commands do not implement |
| `AGENTS.md` and `.agents/rules/` | Contributor and agent operating policy | Whether a scientific claim has been proved |
| Accepted ADRs in `doc/adr/` | The rationale and consequences of architectural decisions | Current implementation status when code has since changed |
| `PROJECT_STATE.yaml` | A concise index of current durable project state and verified milestones | Raw evidence, transient Git status, or a substitute for tests |
| Roadmaps and backlogs, when established | Intended future work and prioritization | Completed or verified work |
| R&D notes in `doc/rnd/` | Exploratory derivations, calculations, experiments, and research history | Production guarantees or publication-ready claims by default |
| Source code and tests | Implemented behavior and the computational checks currently encoded | A general mathematical proof beyond the tested scope |
| `README.md` | Project orientation and entry points | A second, independently maintained copy of detailed policy |
| `CHANGELOG.md` | User-visible changes associated with package development and releases | A complete commit or research diary |

### 2.1 Resolving inconsistencies

When two artifacts disagree:

1. Identify the type of question: server policy, executable behavior, contributor policy, architectural rationale, scientific result, or project status.
2. Consult the authoritative artifact for that question using the table above.
3. Treat an inconsistency as a maintenance issue; do not silently choose whichever statement is more convenient.
4. Update the stale or incorrectly scoped artifact in a reviewed change.
5. If the inconsistency exposes an unresolved design decision, record that decision before changing implementation.

A prose rule must not claim that a check is enforced when no corresponding CI job, tool configuration, or documented manual review process exists. Conversely, executable configuration must not be changed merely to make it conform to stale prose.

Scientific claims require their own evidence standard. An analytical proof, an exact finite-size calculation, a numerical test, and an extrapolation are different forms of support. A passing test demonstrates only the property and parameter domain that the test actually covers.

## 3. Repository areas and responsibilities

The following map describes the repository at the time this baseline was written. Existing paths are distinguished from proposed locations.

### 3.1 Existing canonical areas

| Path | Responsibility |
|---|---|
| `AGENTS.md` | Entry point for contributor and agent instructions; links to the detailed rules |
| `.agents/rules/` | Modular normative policies |
| `src/kitaev_mf/domain/` | Mathematical domain objects and interfaces, subject to the implemented package architecture |
| `src/kitaev_mf/services/` | Physics and computational services |
| `src/kitaev_mf/infra/` | Concrete numerical, hardware, storage, and I/O adapters |
| `test/` | Automated verification; tests should be organized by scope and responsibility |
| `doc/adr/` | Architecture Decision Records |
| `doc/rnd/` | Exploratory theoretical and numerical research notes |
| `doc/reference/` | Reference material, including the thesis source snapshot |
| `PROJECT_STATE.yaml` | Compact durable project-state index |
| `pyproject.toml` | Package metadata and Python tool configuration |
| `requirements.in`, `requirements-dev.in` | Direct runtime and development dependency specifications |
| `requirements.txt`, `requirements-dev.txt` | Resolved dependency lockfiles |
| `.github/workflows/` | GitHub Actions automation |
| `Makefile` | Local development commands |
| `README.md` | Project overview and developer entry point |
| `CHANGELOG.md` | User-visible package change history |

The table is a responsibility map, not a guarantee that every listed subdirectory is populated or that every planned architectural layer is complete. In particular, do not create empty directories merely to make the tree match this document.

### 3.2 Proposed governance documentation area

The following locations are reserved as possible homes for governance planning and audit records. They are **not established requirements to create all of these files immediately**:

- `doc/governance/README.md`: governance architecture and navigation.
- `doc/governance/ROADMAP.md`: ordered governance and repository-improvement milestones.
- `doc/governance/BACKLOG.md`: actionable, prioritized governance tasks.
- `doc/governance/DEVELOPMENT_HISTORY.md`: concise record of material development and governance transitions.
- `doc/governance/audits/`: dated or commit-pinned audit reports.

Create these artifacts incrementally when the corresponding workflow is adopted. Avoid copying transient branch status, test output, or entire commit histories into durable documents.

### 3.3 Boundaries between documentation areas

- Put reusable contributor requirements in `.agents/rules/`, not in an individual research note.
- Put the reasoning behind a durable architectural choice in `doc/adr/`; link to it from the relevant rule when useful.
- Put derivations and research experiments in `doc/rnd/`, not in an ADR merely because they motivate a design.
- Keep external or upstream reference material in `doc/reference/`; record provenance and avoid silently modifying imported sources.
- Put planned work in a roadmap or backlog, not in `PROJECT_STATE.yaml` as though it were implemented.
- Keep `AGENTS.md` as an entry point and concise operating summary. Detailed rules belong in their dedicated files.
- Keep tool settings in their native configuration files. Documentation may explain those settings but should not become a second configuration source.

## 4. Architecture Decision Records

ADRs document decisions with consequences beyond a single implementation detail.

An ADR should normally include:

- A stable identifier and descriptive title.
- A status such as **Proposed**, **Accepted**, **Superseded**, or **Rejected**.
- The context and concrete problem.
- The decision and its scope.
- Important alternatives considered and why they were not selected.
- Consequences, including costs, limitations, and follow-up work.
- References to relevant code, tests, issues, or research notes.

An accepted ADR records the rationale that was accepted at the time. Do not silently rewrite that rationale when the architecture evolves. Preserve the historical record, update its status when appropriate, and create a new ADR for a materially different decision. Small editorial corrections that do not change meaning need not create a new decision record.

An ADR is not a proof. Mathematical statements used to justify an architectural decision must be supported by an appropriate derivation or primary reference.

## 5. Project state, planning, and evidence

These records have deliberately separate roles:

- **Project state:** a concise snapshot of durable facts, such as the active milestone, completed milestones, verified invariants, and known blockers.
- **Roadmap:** the intended sequence of larger objectives.
- **Backlog:** actionable work items with scope, dependencies, and acceptance criteria.
- **R&D notes:** the technical record of exploration, including unsuccessful approaches and unresolved questions.
- **Tests and reproducibility artifacts:** executable evidence for the properties they actually check.

Do not use `PROJECT_STATE.yaml` to assert that a working tree is clean, that a branch is current, or that CI is passing as a timeless fact. Such information is transient. If a dated audit needs to record it, identify the commit and observation time, and treat it as a historical snapshot.

A milestone should be marked complete only when its acceptance criteria are explicit and the corresponding checks or evidence have been reviewed. Record the scope of verification: a unit test, a finite-volume benchmark, a convergence study, and an analytical proof are not interchangeable.

## 6. Governance change workflow

Changes to governance are repository changes and follow the same review discipline as other changes.

1. Identify the policy or organizational problem and its affected files.
2. Check for existing rules, ADRs, configuration, and workflows before adding a new document or mechanism.
3. Make the smallest coherent change. Avoid mixing governance changes with unrelated physics implementation.
4. Update cross-references and examples affected by the change.
5. Validate Markdown links and paths, configuration syntax, and any automated checks relevant to the modified files.
6. Open a pull request against `main` describing the motivation, scope, compatibility implications, and verification.
7. After integration, remove obsolete or superseded guidance only when its replacement is clear and traceable.

The governance branch is a long-lived working branch for developing and reviewing repository policy. It is not the parent branch for ordinary feature or research work. Governance changes reach `main` through reviewed pull requests; development branches should use `main` as their integration baseline unless a task explicitly requires otherwise.

## 7. Maintenance of this document

This document defines the current governance architecture, not an immutable directory specification. Update it when repository structure or document responsibilities materially change. Keep the existing-versus-proposed distinction accurate, and do not expand the governance system without a demonstrated need.

