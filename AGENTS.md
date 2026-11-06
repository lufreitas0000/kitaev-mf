# Agent and Contributor Protocol

This file is the repository entry point for human contributors and automated agents. It defines the operating sequence and points to the detailed policies; it is not a duplicate of those policies.

---

## Strategic Project Vision: Multi-Phase Architecture

The long-term objective of `kitaev-mf` is structured into dedicated, model-family phases:

1. **Phase 1: Pure Ferromagnetic Kitaev Honeycomb ($K = -1.0$) under $\mathbf{h} \parallel [111]$ (Numerical Quantification & Algorithmic Foundations)**
   * *Exclusive Scope:* Exhaustive numerical characterization of all solver algorithms (Picard, DIIS, Full Broyden, Limited-Memory Broyden, JFNK, Matrix Sign Purification) on the real orthogonal Grassmannian $\Gr(4, 8)$ and real-space torus clusters.
   * *Deliverable:* Complete operational mastery of every numerical nuance: convergence basins, critical slowing down at spinodal transitions ($h_{c2}^+ = 0.400$), step and memory complexity, precision loss on Dirac resonances ($3 \mid L$), and explicit hardware resource/FLOP calibration on CPU/GPU.
   * *Status:* **Completed & Certified** (11 chapters in Living Lab Book, 20 verified invariants `K-001`--`K-020`, 150/150 passing tests).

2. **Phase 2: Physics of the Ferromagnetic Kitaev $[111]$ Model (Ground State, Vortices, Thermodynamics & Dynamical Response)**
   * *Exclusive Scope:* Direct physical continuation of Phase 1 in the *identical* model family ($K = -1.0$, $\mathbf{h} \parallel [111]$), selecting the optimal solver discovered in Phase 1 (Joint DIIS and warm-started matrix sign purification) to map the complete physical phase diagram and compute all physical observables.
   * *Core Physical Targets:*
     - Ground-state phase diagram and critical exponents ($h \in [0, 0.6]$).
     - Localized topological flux defects: $T=0$ states with isolated pairs of $\mathbb{Z}_2$ visons ($W_p = -1$), interaction potential $V(R)$, bound Majorana zero modes, and emergent electric polarization textures.
     - Finite-temperature thermodynamic states: thermal flux proliferation, vison deconfinement crossover, and specific heat peaks.
     - Static and dynamical response: two-spin dynamical correlators $S^{\alpha\beta}(\mathbf{q}, \omega)$, optical conductivity, and inelastic neutron scattering cross-sections.

3. **Subsequent Phases (Dedicated Model Families):**
   * Subsequent phases will tackle each distinct physical model family sequentially (e.g. Antiferromagnetic Kitaev $K > 0$, Extended $J$--$K$--$\Gamma$--$\Gamma'$--$D$ models on honeycomb, and 3D Kitaev lattices).

---

## 1. Required guidance

Read the rules relevant to the task before changing files:

- `.agents/rules/00_workspace_hygiene.md` — workspace hygiene, file granularity, and repository layout.
- `.agents/rules/01_repository_safety.md` — Git discipline, review, and simulation-data integrity.
- `.agents/rules/02_tdd_protocol.md` — testing and numerical verification.
- `.agents/rules/03_editorial_rigor.md` — scientific exposition, claim provenance, notation, and references.
- `.agents/rules/04_versioning_and_releases.md` — package versioning, changelog, and releases.
- `.agents/rules/05_project_management.md` — phase, milestone, task, and sprint planning.
- `.agents/rules/06_milestone_progress.md` — development-branch state updates, verification evidence, and task handover.

Also consult the relevant accepted ADRs in `doc/adr/`, research records in `doc/rnd/`, and project-state records when they materially affect the task. The governance architecture in `doc/governance/README.md` defines the responsibilities and precedence of repository documents.

If guidance conflicts with executable configuration or the current repository structure, identify the discrepancy. Do not silently claim that a planned directory, check, or workflow exists. Resolve policy conflicts explicitly before proceeding.

## 2. Session start

For every substantial task:

1. Inspect the current branch, working-tree status, and existing diff before editing. Preserve unrelated user changes.
2. Read this file and the applicable rules; do not assume that prior-session branch or CI status is still current.
3. Inspect the relevant source, tests, configuration, ADRs, and R&D notes. Read `PROJECT_STATE.yaml` when milestone status, conventions, or verified invariants matter.
4. Confirm the task's scope, dependencies, acceptance criteria, and verification plan. For work on a long-lived branch, compare with `main` when intervening changes may affect the task.

For a small, self-contained task, keep this inspection proportional; do not perform unrelated repository-wide refactoring.

## 3. Branch and change discipline

- Use a dedicated branch for routine development and integrate through a reviewed pull request to `main`. Follow Rule 01 for branch naming and history protection.
- Ordinary feature and research branches should use `main` as their integration baseline unless a task explicitly requires otherwise.
- `governance/project-governance` is reserved for developing repository governance. It is not a parent branch for ordinary feature or research work.
- Do not push routine development changes directly to `main`.
- Do not reset, force-push, clean, delete branches, or discard user changes without explicit authorization.
- Keep each change set coherent and commits focused. Avoid mixing physics implementation, governance, dependency changes, and unrelated cleanup.
- Before editing an existing file, inspect its current contents and preserve established conventions unless the task explicitly changes them.

## 4. Implementation protocol

### Mathematical and scientific contracts

Before non-trivial physics implementation, state the mathematical contract: objects and domains, assumptions, conventions, boundary conditions, expected invariances or limiting cases, and relevant failure modes. Identify the intended interfaces and dependencies before selecting an implementation.

Distinguish:
- exact algebraic statements and analytical proofs;
- finite-volume or finite-dimensional results;
- numerical tests and benchmarks;
- approximations, extrapolations, and conjectures.

Do not promote finite-size numerical agreement into a thermodynamic claim without a controlled limiting argument. State Fourier, gauge, graph-orientation, normalization, and boundary conventions whenever they affect the result. Follow Rule 03 for scientific documentation and provenance.

### Software architecture

The intended dependency direction is:

- `domain/`: mathematical entities and interfaces; independent of services and infrastructure.
- `services/`: physics and computational algorithms depending on domain abstractions.
- `infra/`: concrete numerical backends, hardware, storage, and I/O adapters.
- `app/`: application and workflow orchestration, if and when this layer is established.

Do not create missing architectural layers merely to satisfy this description. Keep pure mathematical logic separate from side effects. Route backend-specific allocation and operations through established interfaces; do not leak device-specific behavior into domain objects.

Prefer the smallest implementation that satisfies the contract. Reuse existing abstractions, avoid speculative generality, and keep interfaces explicit and testable.

### Tests and numerical verification

Follow Rule 02 and the repository's actual test/tool configuration. Add focused tests for analytical identities, limiting cases, invariances, and failure modes. Run the relevant tests and verify that a new test fails for the intended reason before relying on it as a regression check when feasible.

Do not use a universal numerical tolerance by habit. Justify tolerances in relation to scale, conditioning, precision, algorithm, and the expected numerical error. Separate exact symbolic identities from floating-point acceptance criteria.

## 5. Verification and reporting

Before declaring a task complete:

1. Run checks relevant to the changed files and affected interfaces.
2. Inspect the final diff for scope, accidental generated files, whitespace errors, and stale references.
3. Report the exact commands run and their outcomes. Distinguish passed, failed, skipped, and not-run checks.
4. Do not claim CI passed unless the corresponding remote workflow actually ran and its result was inspected.
5. Identify unresolved issues and limitations rather than silently broadening the acceptance claim.

A test passing establishes only the tested property and parameter domain. It is not, by itself, a proof of a general mathematical statement.

## 6. Session handover and project records

At the end of a substantial task, provide a concise handover containing:

- the branch and commit/PR state;
- the files changed and the purpose of each change;
- verification commands and results;
- remaining work, blockers, and decisions requiring user input.

Update `PROJECT_STATE.yaml` only when durable project facts, milestone status, verified invariants, or blockers change. Do not use it as a live Git/CI dashboard: branch names, clean-tree claims, and CI status become stale immediately. Put dated audit observations in an explicitly time- and commit-scoped audit record when needed.

Update R&D notes when a material derivation, failed approach, numerical result, or research decision should be preserved. Use an ADR for a durable architectural decision, not for routine implementation detail. Keep roadmap/backlog planning separate from completed project state.

## 7. Project-specific research track: exact pure-Kitaev perturbation theory

The exact pure-Kitaev perturbation work is scientifically distinct from the production mean-field solver. Its intended research-note location is `doc/rnd/exact-perturbation-theory/`; create or use that area only when the corresponding work is active. Keep exploratory code and generated outputs isolated from production modules.

For this track:

1. Derive operator identities, conventions, physical constraints, and numerical contracts before implementing solver modules.
2. Distinguish the enlarged Majorana Fock space from the physical spin Hilbert space. Derive the gauge constraint and matter-parity selection rule before interpreting Gaussian flux-sector energies as physical thresholds.
3. Do not identify production mean-field coupling parameters or dispersion prefactors with those in the Pauli-spin/Majorana derivation by inspection. Reconcile conventions through an explicit operator-level map before quantitative comparison; consult the relevant R&D backlog and notes when present.
4. The standard `b_i^\alpha,c_i` representation and the production four-flavor convention are scope-dependent. Do not transfer symbols or prefactors between them without an explicit map.
5. Use transparent finite-volume calculations and, where feasible, spin-space exact diagonalization as independent checks. Numerical agreement is evidence for the tested finite systems, not a thermodynamic proof.
6. Do not integrate exploratory code into production until its mathematical contract, tests, and compatibility with domain/backend interfaces have been reviewed.
7. Keep generated arrays, logs, plots, and build products under the designated ignored build area (currently `build/rnd/exact-perturbation-theory/`); commit only compact, reviewed results with reproducibility metadata.
