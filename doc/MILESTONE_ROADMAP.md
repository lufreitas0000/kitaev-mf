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
- **Reconciliation:** Reconciled `T04b-03` reference value to certified benchmark `-0.5248657458506`.
- **Evidence:** `test/unit/test_energetics.py`, `src/kitaev_mf/services/quadrature.py`, `scripts/k006_quadrature.py`, `doc/rnd/note_04c_energy_normalization_audit.md`, `doc/rnd/note_04d_k006_quadrature_provenance.md`.

#### M1.4.4: Production Observables and Energetics Services (Closed)
- **Status:** `completed`
- **Deliverables:**
  - `src/kitaev_mf/domain/lattice.py`: Plaquette topology.
  - `src/kitaev_mf/domain/state.py`: EnergyBreakdown normalization.
  - `src/kitaev_mf/services/observables.py`: Sublattice magnetization, bond spin-correlators, Wilson loop fluxes.
  - `src/kitaev_mf/services/energetics.py`: Variational Gaussian functional, spectral energy, and consistency check.
  - `test/unit/test_production_observables.py`: Unit contracts.
  - `test/unit/test_production_energetics.py`: Unit contracts.
  - `test/integration/test_m1_system_convergence.py`: End-to-end integration suite.

#### M1.4.5: Consolidated Verification and Quality Gate (Active)
- **Objective:** Run full test suite across unit and integration targets, verify clean imports, and test CPU/CUDA execution parity.

---

### Subsequent Milestones
- **M1.5:** Standalone LaTeX research monograph (`doc/report/` / `main.pdf`).
