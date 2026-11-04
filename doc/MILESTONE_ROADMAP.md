# Milestone Roadmap: kitaev-mf

## Milestone M1: Analytical Foundations, Bloch Formalism, and Self-Consistent Solvers

### Milestone M1.1 — Local Algebra and Exact Parton Representation (Closed)
- **Status:** `completed`
- **Contracts:** A01–A20 verified across 25 unit tests.
- **Evidence:** `test/unit/test_local_majorana_algebra.py`, `test_local_majorana_representations.py`, `test_physical_projection.py`, `test_exchange_hamiltonian.py`.

### Milestone M1.2 — Bloch Formalism and Folding Parity (Closed)
- **Status:** `completed`
- **Contracts:** Exact Kitaev dispersion (K-004) and 2-site to 4-site supercell folding parity (K-005) verified to `< 1e-14`.
- **Evidence:** `src/kitaev_mf/services/hamiltonian.py`, `test/unit/test_hamiltonian.py`.

### Milestone M1.3 — Reciprocal Projection and Self-Consistent Solver (Closed)
- **Status:** `completed`
- **Contracts:** Zero-temperature spectral projection, real-space BZ integration, and Joint Anderson Acceleration (DIIS).
- **Evidence:** `src/kitaev_mf/services/correlators.py`, `src/kitaev_mf/services/solver.py`, `test/unit/test_correlators.py`, `test/unit/test_solver.py`.

### Milestone M1.4 — Physical Observables, Plaquette Fluxes, and Exact Energetics

#### M1.4.1: Observable and Energetics Specification (Closed)
- **Status:** `completed`
- **Evidence:** `doc/rnd/note_04_observables_energetics.md`.

#### M1.4.2: Independent Finite-Fock Observable Tests (Closed)
- **Status:** `completed`
- **Evidence:** `test/unit/test_observables.py`, `doc/rnd/note_04a_observable_test_specification.md`.

#### M1.4.3: Gaussian Energetics and Quadrature Provenance (Closed)
- **Status:** `completed`
- **Validation:** User-reported WSL verification passed (`make verify` green across 70 tests).
- **Reconciliation:** Reconciled `T04b-03` reference value to certified benchmark `-0.5248657458506`.
- **Evidence:** `test/unit/test_energetics.py`, `src/kitaev_mf/services/quadrature.py`, `scripts/k006_quadrature.py`, `doc/rnd/note_04c_energy_normalization_audit.md`, `doc/rnd/note_04d_k006_quadrature_provenance.md`.

#### M1.4.4: Production Observables and Energetics Services (Active)
- **Objective:** Deploy production services wrapping verified Fock oracles, implement plaquette Wilson loops, and build system integration tests.
- **Batches:**
  - `M1.4.4a`: Define contracts and test specifications for production observables and energetics.
  - `M1.4.4b`: Implement `src/kitaev_mf/services/observables.py` (magnetization and bond spin-correlators).
  - `M1.4.4c`: Define `Plaquette` domain value object and implement `compute_plaquette_fluxes`.
  - `M1.4.4d`: Define `EnergyBreakdown` domain value object and write normalization tests.
  - `M1.4.4e`: Implement `src/kitaev_mf/services/energetics.py` (variational Gaussian functional).
  - `M1.4.4f`: Integrate quadratic spectral ground-state energy ($E_{\text{quad}, 0} - C$) and assert consistency.
  - `M1.4.4g`: System integration test suite across solver, observables, and energetics.

---

### Subsequent Milestones
- **M1.4.5:** Consolidated verification (K-001–K-006, CPU/CUDA parity, import-cycle checks).
- **M1.5:** Standalone LaTeX research monograph (`doc/report/` / `main.pdf`).
