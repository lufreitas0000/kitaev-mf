# Operating Protocol for AI Agents and Collaborators

This repository develops Majorana parton mean-field theory for extended Kitaev models and maintains a separate research track on exact perturbation theory of the pure Kitaev spin liquid.

All contributors must follow the modular rules in `.agents/rules/`:

- `.agents/rules/00_workspace_hygiene.md`: file granularity, directory taxonomy, and generated-artifact hygiene.
- `.agents/rules/01_repository_safety.md`: branch safety and simulation-data integrity.
- `.agents/rules/02_tdd_protocol.md`: test-first implementation and numerical verification.
- `.agents/rules/03_editorial_rigor.md`: mathematical exposition, claim provenance, notation, and bibliography.
- `.agents/rules/04_versioning_and_releases.md`: SemVer, changelog maintenance, and release procedure.

If a task appears to conflict with a rule, identify the conflict and resolve it explicitly before implementation. Do not silently bypass a rule.

---

## 1. Engineering and architecture

1. **No zero-shot implementation.** Before non-trivial code, state the mathematical contract, dependencies, interfaces, and verification strategy. Write focused failing tests before implementation.
2. **Package layout.** This research branch currently retains the legacy `src/domain/` and `src/infra/` layout. The production architecture on `main` is under `src/kitaev_mf/`, with `domain/`, `services/`, `infra/`, and `app/` subpackages. Do not refer to the `main` package paths as if they already exist on this branch. Architecture synchronization is a separate integration task; until then, follow the imports and interfaces actually present in the checked-out branch.
3. **Dependency direction.** Domain code must not depend on services, infrastructure, or application orchestration. Physics services should depend on domain interfaces, not concrete backend implementations.
4. **Pure logic and side effects.** Keep mathematical transformations deterministic and side-effect free where practical. Isolate I/O, device allocation, persistence, and workflow orchestration behind explicit interfaces.
5. **Backend discipline.** Route backend-specific linear algebra through the established backend protocol. Do not leak Torch/device-specific allocation into domain objects or pure physics algorithms.
6. **Test-driven development.** Follow the Red-Green-Refactor process in Rule 02. Add tests for analytical identities, limiting cases, invariances, and failure modes. Do not claim a milestone is complete without running its acceptance checks.
7. **Repository safety.** Never use destructive Git operations or force-push without explicit authorization. Work on the requested branch and inspect the working tree before changes.

## 2. Session and project-state protocol

At the beginning of a substantial task:

1. Read `PROJECT_STATE.yaml` and identify the current production milestone, verified benchmarks, and recorded research tracks.
2. Inspect relevant architecture decision records under `doc/adr/`.
3. Read the applicable R&D roadmap, backlog, and notes before extending an existing research task.
4. Confirm the current branch and compare against `main` when shared interfaces or conventions may have changed.

At the end of a substantial task:

1. Report changed files, commits, checks run, and unresolved issues.
2. Update project state or R&D status only when the recorded milestone, benchmark, or research status actually changes; do not manufacture status updates for routine edits.
3. Record failed approaches and material convention decisions in the appropriate R&D note or ADR.

## 3. Scientific standards

1. Distinguish exact algebraic results, finite-volume results, numerical evidence, approximations, and conjectures.
2. State Hilbert space, boundary conditions, graph orientation, Fourier convention, gauge choice, and normalization whenever they affect a result.
3. Do not infer thermodynamic statements from finite-size calculations without an explicit limiting argument or controlled extrapolation.
4. Verify bibliographic metadata and numerical benchmarks against primary sources. Never invent identifiers, coefficients, or citations.
5. Keep exploratory derivations in `doc/rnd/`. Do not present exploratory notes as a finished publication or production specification.
6. Follow Rule 03 for claim provenance and references. Existing notes that predate the rule require an explicit migration task; do not retroactively claim compliance.

## 4. Separate research track: exact pure-Kitaev perturbation theory

The research track lives in `doc/rnd/exact-perturbation-theory/` and is intentionally distinct from the production mean-field solver.

1. **Analytical-first workflow.** Derive and agree on conventions, operator identities, physical constraints, and numerical contracts before implementing solver modules.
2. **Physical projection is mandatory.** Distinguish the enlarged Majorana Fock space from the physical spin Hilbert space. Derive the gauge constraint and matter-parity selection rule before interpreting Gaussian flux-sector energies as physical thresholds.
3. **Normalization is not assumed.** The production mean-field Bloch convention and the research notes' Pauli-spin/Majorana convention have not yet been reconciled. Do not identify their coupling parameters or dispersion prefactors by inspection. Complete backlog items B-001 and B-002 before using production spectra as quantitative validation.
4. **Notation is scope-dependent.** Production mean-field code uses its established four-flavor convention. The exact Kitaev derivation may use the standard `b_i^\alpha,c_i` representation. Do not conflate symbols or transfer prefactors between these representations without an explicit operator-level map.
5. **Independent validation.** Use transparent finite-volume calculations and, where feasible, spin-space exact diagonalization as independent checks. Numerical agreement is evidence for the tested finite systems, not by itself a proof of a thermodynamic claim.
6. **No premature production integration.** Keep exploratory code isolated until mathematical contracts, tests, and compatibility with the production domain/backend interfaces have been reviewed.
7. **Generated artifacts.** Place untracked arrays, logs, plots, and build products under `build/rnd/exact-perturbation-theory/`. Commit only compact, reviewed results with reproducibility metadata.

## 5. Versioning and release governance

1. The package version has one source of truth: `[project].version` in `pyproject.toml`.
2. Follow Rule 04 for every version bump, changelog entry, pre-release, and release tag.
3. Only publish releases from `main`. Research milestones and R&D note identifiers do not trigger package version changes.
4. Do not claim a release exists until its matching tag and GitHub release have been published.

## 6. Conventions and implementation changes

- The production convention documented in the current architecture and Rule 02 is authoritative for existing production mean-field code.
- The exact perturbation track retains its explicitly stated Hamiltonian and Majorana conventions until B-001/B-002 derive the mapping.
- Any change to a shared convention must include: operator-level derivation, affected files/interfaces, regression tests, and an update to the relevant documentation.
- Do not change a physical prefactor merely to make two existing formulas agree. Identify the Hamiltonian normalization, bond counting, Majorana anticommutation convention, and quadratic-form convention first.
