# Project Roadmap and Milestone Backlog: kitaev-mf

## Planning hierarchy

- **Phase**: macro physical and engineering horizon.
- **Milestone**: persistent strategic backlog item and acceptance gate. IDs use M<phase>.<index>; a nested milestone may add a subordinate suffix, as in M1.4.5.
- **Task**: short-lived execution unit, numbered T1, T2, ... only within the active milestone. Task IDs restart when the active milestone changes.

The roadmap contains persistent milestone objectives and acceptance gates. Current execution tasks are maintained in PROJECT_STATE.yaml. Future milestones are not decomposed into task lists before activation.

## Phase 1: Analytical Rigor and Domain Foundations

### Completed milestones

#### M1.1 — Local Majorana Algebra and Exact Parton Representation

- **Status:** completed.
- **Scope:** local Clifford CAR algebra, spin/gauge representation, physical projection, Pauli algebra, and exact quartic exchange Hamiltonian.
- **Acceptance:** contracts A01–A20 verified across 25 unit tests.
- **Evidence:** R&D notes 00a–00e; unit tests for local algebra, representations, physical projection, and exchange Hamiltonian. See PROJECT_STATE.yaml for exact paths and provenance.

#### M1.2 — Generalized Decoupling and Momentum-Space Bloch Matrices

- **Status:** completed.
- **Scope:** Bloch-Majorana Hamiltonian construction for two-site and four-site unit cells.
- **Acceptance:** pure-Kitaev analytical dispersion and supercell Brillouin-zone folding parity verified, including general-coupling folding tests.
- **Evidence:** src/kitaev_mf/services/hamiltonian.py, test/unit/test_hamiltonian.py, and doc/rnd/note_02_momentum_bloch_formalism.md. K-004 and K-005 are recorded in PROJECT_STATE.yaml.

#### M1.3 — Reciprocal Spectral Projection and Self-Consistent Solver

- **Status:** completed.
- **Scope:** zero-temperature spectral projection, Brillouin-zone integration, correlator reconstruction, and joint Anderson acceleration (DIIS).
- **Acceptance:** correlator and solver regression suites pass; the pure-Kitaev matter correlator has the analytical Brillouin-zone expression and finite-mesh regression evidence.
- **Evidence:** src/kitaev_mf/services/correlators.py, src/kitaev_mf/services/solver.py, test/unit/test_correlators.py, test/unit/test_solver.py, and doc/rnd/note_03_correlators_saddle_point.md. The matter correlator is negative, approximately -0.5248657458506, in the convention recorded as K-006; the gauge-link correlator is +1.

#### M1.4 — Production Observables, Plaquette Wilson Loops, and Variational Energetics

- **Status:** completed.
- **Scope:** sublattice magnetizations, bond spin correlators, plaquette Wilson loops, Gaussian energy decomposition, spectral ground-state energy, and energy-consistency checks.
- **Acceptance:** production service tests and the M1 end-to-end system-convergence test are present; local WSL validation was reported by the project owner.
- **Evidence:** src/kitaev_mf/services/observables.py, src/kitaev_mf/services/energetics.py, test/unit/test_production_observables.py, test/unit/test_production_energetics.py, test/integration/test_m1_system_convergence.py, and doc/rnd/note_04e_m1_4_integration_and_audit.md.

### Active milestone

#### M1.4.5 — Consolidated Verification and Invariant Hardening

- **Status:** active.
- **Objective:** establish independent, layered verification of kinematic sum rules, Gaussian-state representability, solvable limits, symmetry selection rules, variational response, dualities, and backend parity.
- **Scope:** production services and test suites under src/, test/, scripts/, and Makefile.
- **Research specification:** doc/rnd/note_04f_advanced_invariants_and_sum_rules_backlog.md. This note is a contract proposal; each identity's assumptions, conventions, and numerical realization must be audited before it is treated as a regression gate.
- **Acceptance gate:** all four active tasks T1–T4 completed; tests and architecture checks pass; CPU verification is recorded; CUDA parity is recorded when CUDA is available. Any skipped hardware check must be reported explicitly.
- **Execution detail:** see active_tasks in PROJECT_STATE.yaml. Tasks are not copied into this roadmap as a second source of truth.

The previous package labels M1.4.5a–M1.4.5d are retired as task identifiers. Their scopes are represented by T1–T4. These tasks are planned; no implementation or verification success is implied.

### Strategic backlog

#### M1.5 — Standalone Research Monograph

- **Status:** backlog.
- **Objective:** deliver a publication-grade LaTeX monograph synthesizing the theoretical and computational results from Parts I–V.
- **Acceptance:** coherent manuscript and reproducible PDF build; equations, claims, figures, and numerical results trace to source derivations or tests using stable CLM-M1 identifiers.
- **Planned deliverable:** doc/report/main.pdf. The doc/report/ directory is a proposed target and is not asserted to exist.

## Phase 2: Homogeneous Phase Diagrams and Magnetic Ordering

#### M2.1 — Multi-Sublattice Order-Parameter Framework

- **Status:** backlog.
- **Objective:** extract and classify Néel, ferromagnetic, zigzag, and stripy order parameters across supported unit cells.
- **Acceptance:** symmetry-aware order parameters are defined, tested on limiting states, and reported with cell/normalization conventions.

#### M2.2 — Automated Phase-Boundary Tracking

- **Status:** backlog.
- **Objective:** implement pseudo-arclength continuation and gap-closing detection for phase-boundary tracking.
- **Acceptance:** continuation reproduces benchmark boundaries, detects documented gap closures, and reports convergence/failure diagnostics.

#### M2.3 — Anisotropic Critical-Field Mapping

- **Status:** backlog.
- **Objective:** map critical fields as a function of field orientation.
- **Acceptance:** angular scans are reproducible and compared against explicitly identified DMRG or other independent benchmarks.

## Phase 3: Real-Space Inhomogeneous Vortices and Bound States

#### M3.1 — Large-Scale Real-Space Lattice Engine

- **Status:** backlog.
- **Objective:** support real-space open and periodic clusters, with N up to 3200 sites as a target to be benchmarked.
- **Acceptance:** boundary conditions, graph conventions, scaling, and numerical accuracy are validated against smaller exactly tractable clusters.

#### M3.2 — Metastable Vortex Configurations

- **Status:** backlog.
- **Objective:** stabilize and characterize localized vison-pair configurations and associated string branch cuts.
- **Acceptance:** configurations are reproducible, satisfy gauge/flux constraints, and have controlled energy and stability diagnostics.

#### M3.3 — Induced Polarization and Quadrupole Interactions

- **Status:** backlog.
- **Objective:** compute local charge-imbalance profiles and effective vortex interactions.
- **Acceptance:** polarization profiles and interaction energies are validated against symmetry constraints and finite-size/boundary-condition checks.

## Phase 4: Spectroscopic and Dynamical Observables

#### M4.1 — Scanning Tunneling Spectroscopy

- **Status:** backlog.
- **Objective:** compute local Majorana Green functions and tunneling/cotunneling observables relevant to zero-bias features.
- **Acceptance:** spectral functions obey normalization and causality checks and reproduce appropriate limiting cases.

#### M4.2 — Dynamical Spin Structure Factor

- **Status:** backlog.
- **Objective:** compute S(q, omega) and characterize the inelastic response continuum.
- **Acceptance:** frequency/momentum resolution, sum rules, and comparison benchmarks are documented.

#### M4.3 — Thermal Hall Conductance and Chiral Edge States

- **Status:** backlog.
- **Objective:** compute Berry-curvature/Chern diagnostics, ribbon spectra, chiral edge modes, and thermal Hall response.
- **Acceptance:** bulk topological invariants and edge spectra are mutually consistent in validated regimes.

## Phase 5: Cloud Production Pipelines and Open-Science Archiving

#### M5.1 — GPU Cluster Orchestration

- **Status:** backlog.
- **Objective:** automate reproducible parameter sweeps on cloud GPU resources, with L4/A100-class hardware as candidate targets rather than guaranteed infrastructure.
- **Acceptance:** jobs are restartable, resource-aware, and emit structured logs and provenance metadata.

#### M5.2 — Open-Science Provenance and Archiving

- **Status:** backlog.
- **Objective:** establish versioned data pipelines, atomic HDF5 outputs, and archival release workflows.
- **Acceptance:** datasets are reproducible and traceable to code/configuration; DVC and Zenodo DOI workflows are evaluated and adopted where appropriate.
