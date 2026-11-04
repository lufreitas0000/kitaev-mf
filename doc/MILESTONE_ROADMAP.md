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

### Milestone M1.4 — Physical Observables, Plaquette Fluxes, and Exact Energetics (Closed)
- **Status:** `completed`
- **Evidence:**
  - `src/kitaev_mf/services/observables.py`: Sublattice magnetizations, bond spin-correlators, plaquette fluxes.
  - `src/kitaev_mf/services/energetics.py`: Variational energy, spectral ground-state energy, consistency check.
  - `test/unit/test_production_observables.py`, `test_production_energetics.py`, `test/integration/test_m1_system_convergence.py`.
  - Audit record: `doc/rnd/note_04e_m1_4_integration_and_audit.md`.

---

### Milestone M1.4.5 — Consolidated Verification and Invariant Hardening (Active)

- **Objective:** Deploy the full verification suite testing kinematic sum rules, semidefinite representability constraints, solvable reductions, discrete selection rules, and canonical dualities across the codebase.
- **Parent R&D Specification:** `doc/rnd/note_04f_advanced_invariants_and_sum_rules_backlog.md`.

#### Batches & Execution Schedule:

1. **M1.4.5a: Package 1 — Kinematic and Semidefinite Gates (Tier 1)**
   - Target Module: `test/unit/test_kinematic_and_semidefinite.py`.
   - Contracts:
     - `TST-KIN-01`: Parseval–Frobenius second spectral moment identity (`CLM-M1-INV-001`).
     - `TST-KIN-02`: Pointwise projector-Hamiltonian commutation $[\mathcal{H}(\mathbf{k}), P(\mathbf{k})] \equiv 0$ (`CLM-M1-INV-002`).
     - `TST-KIN-03`: Reciprocal projector idempotency $P^2 = 2P$ and reflection sum rule $P(-\mathbf{k}) + P^T(\mathbf{k}) = 2\mathbb{1}$ (`CLM-M1-INV-002`).
     - `TST-SDP-01`: Local two-site density matrix positivity LMI $\lambda_{\min}(\Gamma_{ij}) \ge 0$ (`CLM-M1-INV-003`).
     - `TST-SDP-02`: Covariance operator norm bounds $\sigma_{\max}(\mathrm{U}_b) \le 1.0$ and $\max \vert{}\nu_k(\mathrm{V}_s)\vert{} \le 1.0$ (`CLM-M1-INV-004`).
     - `TST-SDP-03`: Physical dipole envelope $\Vert{}\langle \boldsymbol{\sigma}_s \rangle\Vert{}_2 \le 1.0$ and parity bounds $\vert{}\operatorname{Pf}(\mathrm{V}_s)\vert{} \le 1.0$ (`CLM-M1-INV-004`).
   - Performance: $< 15\,\mathrm{ms}$.

2. **M1.4.5b: Package 2 — Solvable Limits, Symmetries & Asymptotics (Tier 2)**
   - Target Modules: `test/unit/test_solvable_limits.py`, `test/unit/test_symmetry_selection_rules.py`.
   - Contracts:
     - `TST-RED-01`: Decoupled 1D zigzag chain limit ($K_z = 0 \implies \bar{U} = 4/(3\pi)$) (`CLM-M1-INV-005`).
     - `TST-RED-02`: Isolated dimer flat-band limit ($\vert{}K_z\vert{} \gg \vert{}K_x\vert{}, \vert{}K_y\vert{}$) (`CLM-M1-INV-006`).
     - `TST-SYM-01`: Time-reversal odd correlator vanishing ($\mathrm{U}_b^{0\gamma} \equiv 0$ at $\mathbf{h}=\mathbf{0}, D=0$) (`CLM-M1-INV-007`).
     - `TST-SYM-02`: Sublattice inversion parity in uniform fields ($\langle \boldsymbol{\sigma}_A \rangle = \langle \boldsymbol{\sigma}_B \rangle$).
     - `TST-SYM-03`: Permutation symmetry on the isotropic $C_3$ line ($E_{b,x} = E_{b,y} = E_{b,z}$).
     - `TST-ASY-01`: High-field polarized limit ($\vert{}\mathbf{h}\vert{}/\vert{}K\vert{} \sim 25 \implies \Vert{}\mathrm{U}_b\Vert{} \to 0$).
   - Performance: $\sim 60\,\mathrm{ms}$.

3. **M1.4.5c: Package 3 — Variational Response & Canonical Dualities (Tier 3)**
   - Target Modules: `test/integration/test_variational_response.py`, `test/integration/test_dualities_and_asymptotics.py`.
   - Contracts:
     - `TST-VAR-01`: Numerical Hellmann–Feynman derivative $\frac{\mathrm{d}E}{\mathrm{d}h^\gamma} = -\sum \mathrm{V}_s^{0\gamma}$ (`CLM-M1-INV-008`).
     - `TST-VAR-02`: Thermodynamic energy concavity $\chi^{\gamma\gamma} \ge 0$ (`CLM-M1-INV-008`).
     - `TST-RED-03`: Bipartite sign inversion duality $E(K) \equiv E(-K)$ at $\mathbf{h} = \mathbf{0}$.
     - `TST-DUA-01`: Four-sublattice Klein duality in the 4-site cell ($E(J, K, \Gamma) \equiv E(-J, K+2J, -\Gamma)$) (`CLM-M1-INV-009`).
   - Performance: $\sim 800\,\mathrm{ms}$.

4. **M1.4.5d: Hardware Parity & Final M1 Quality Gate**
   - Consolidated execution of all unit and integration test suites on CPU and CUDA.
   - Zero-import cycle verification.
   - 100% green compliance under `make verify`.

---

### Subsequent Milestone

- **M1.5: Standalone Research Monograph Delivery** (`doc/report/` / `main.pdf`)
  - Integration of Part V (Observables, Variational Energetics, and Plaquette Dynamics).
  - Traceability mapping between LaTeX theorems and pytest assertion targets via `CLM-M1-...` tags.
