# Operating Protocol for AI Agents and Collaborators

Welcome to `kitaev-mf`, a research codebase for Majorana parton mean-field theory on extended Kitaev models.

All contributors must follow the modular rules in `.agents/rules/`:

- `.agents/rules/00_workspace_hygiene.md`: file granularity, directory taxonomy, and generated-artifact hygiene.
- `.agents/rules/01_repository_safety.md`: branch discipline, reviewed pull requests, and simulation-data integrity.
- `.agents/rules/02_tdd_protocol.md`: test-first implementation and numerical verification.
- `.agents/rules/03_editorial_rigor.md`: mathematical exposition, claim provenance, notation, and bibliography.
- `.agents/rules/04_versioning_and_releases.md`: SemVer, changelog maintenance, and release procedure.

If a task conflicts with a rule, identify the conflict and resolve it explicitly before implementation. Do not silently bypass a rule.

## 1. Engineering and architecture

1. **No zero-shot implementation.** Before non-trivial code, state the mathematical contract, dependencies, interfaces, and verification strategy. Write focused failing tests before implementation.
2. **Hexagonal architecture.**
   - `src/kitaev_mf/domain/`: pure mathematical entities and port protocols; no dependencies on services or infrastructure.
   - `src/kitaev_mf/services/`: physics algorithms depending on domain entities and interfaces.
   - `src/kitaev_mf/infra/`: concrete hardware, storage, and I/O adapters.
   - `src/kitaev_mf/app/`: CLI and workflow orchestration.
3. **Dependency direction.** Domain code must not depend on services, infrastructure, or application orchestration. Physics services should depend on domain interfaces, not concrete backend implementations.
4. **Pure logic and side effects.** Keep mathematical transformations deterministic and side-effect free where practical. Isolate I/O, device allocation, persistence, and workflow orchestration behind explicit interfaces.
5. **Backend discipline.** Route backend-specific linear algebra through the established backend protocol. Do not leak Torch/device-specific allocation into domain objects or pure physics algorithms.
6. **Test-driven development.** Follow Rule 02. Add tests for analytical identities, limiting cases, invariances, and failure modes. Do not claim a milestone is complete without running its acceptance checks.

## 2. Session and project-state protocol

At the beginning of a substantial task:

1. Read `PROJECT_STATE.yaml` and identify the current production milestone, verified benchmarks, and recorded research tracks.
2. Inspect relevant architecture decision records under `doc/adr/`.
3. Read the applicable R&D roadmap, backlog, and notes before extending an existing research task.
4. Confirm the current branch and compare against `main` when shared interfaces or conventions may have changed.

At the end of a substantial task:

1. Report changed files, commits, checks run, and unresolved issues.
2. Update project state or R&D status only when the recorded milestone, benchmark, or research status actually changes.
3. Record failed approaches and material convention decisions in the appropriate R&D note or ADR.

## 3. Scientific standards

1. Distinguish exact algebraic results, finite-volume results, numerical evidence, approximations, and conjectures.
2. State Hilbert space, boundary conditions, graph orientation, Fourier convention, gauge choice, and normalization whenever they affect a result.
3. Do not infer thermodynamic statements from finite-size calculations without an explicit limiting argument or controlled extrapolation.
4. Verify bibliographic metadata and numerical benchmarks against primary sources. Never invent identifiers, coefficients, or citations.
5. Keep exploratory derivations in `doc/rnd/`. Do not present exploratory notes as a finished publication or production specification.
6. Follow Rule 03 for claim provenance and references. Existing notes that predate the rule require an explicit migration task; do not retroactively claim compliance.

## 4. Versioning and release governance

1. The package version has one source of truth: `[project].version` in `pyproject.toml`.
2. Follow Rule 04 for every version bump, changelog entry, pre-release, and release tag.
3. Only publish releases from `main`. Research milestones and R&D note identifiers do not trigger package version changes.
4. Do not claim a release exists until its matching tag and GitHub release have been published.

## 5. Separate research track: exact pure-Kitaev perturbation theory

The exact perturbation research track is maintained under `doc/rnd/exact-perturbation-theory/` and is distinct from the production mean-field solver.

1. **Analytical-first workflow.** Derive and agree on conventions, operator identities, physical constraints, and numerical contracts before implementing solver modules.
2. **Physical projection is mandatory.** Distinguish the enlarged Majorana Fock space from the physical spin Hilbert space. Derive the gauge constraint and matter-parity selection rule before interpreting Gaussian flux-sector energies as physical thresholds.
3. **Normalization is not assumed.** The production mean-field Bloch convention and the research notes' Pauli-spin/Majorana convention have not yet been reconciled. Do not identify their coupling parameters or dispersion prefactors by inspection. Complete backlog items B-001 and B-002 before using production spectra as quantitative validation.
4. **Notation is scope-dependent.** Production mean-field code uses its established four-flavor convention. The exact Kitaev derivation may use the standard `b_i^\\alpha,c_i` representation. Do not conflate symbols or transfer prefactors between these representations without an explicit operator-level map.
5. **Independent validation.** Use transparent finite-volume calculations and, where feasible, spin-space exact diagonalization as independent checks. Numerical agreement is evidence for the tested finite systems, not by itself a proof of a thermodynamic claim.
6. **No premature production integration.** Keep exploratory code isolated until mathematical contracts, tests, and compatibility with production domain/backend interfaces have been reviewed.
7. **Generated artifacts.** Place untracked arrays, logs, plots, and build products under `build/rnd/exact-perturbation-theory/`. Commit only compact, reviewed results with reproducibility metadata.
