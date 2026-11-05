# Project Roadmap & Milestone Backlog: kitaev-mf

## Taxonomy
- **Phase**: Macro engineering and physical regimes.
- **Milestone (`M<Phase>.<Index>`)**: Strategic deliverable gates and persistent backlog units.
- **Task (`T<Index>`)**: Immediate operational execution units active exclusively within the current milestone sprint.

---

## Phase 1: Analytical Rigor & Domain Foundations

### Completed Milestones
- **M1.1: Local Algebra & Exact Parton Representation**
  - *Contracts*: A01–A20 verified across 25 unit tests (`K-001`, `K-002`, `K-003`).
- **M1.2: Bloch-Majorana Formalism & Folding Parity**
  - *Contracts*: Exact Kitaev dispersion (`K-004`) and 2-site to 4-site supercell folding parity (`K-005`).
- **M1.3: Reciprocal Spectral Projection & DIIS Solver**
  - *Contracts*: Zero-temperature spectral projector, BZ integration, and Joint Anderson Acceleration. High-resolution reference quadrature $\bar{U} \approx 0.5248657458506$ (`K-006`).
- **M1.4: Production Observables & Gaussian Energetics**
  - *Contracts*: Physical spin observables, plaquette Wilson loops, extensive and density energetics, scalar double-counting subtraction $C = E_{\text{exchange}}$, and end-to-end system convergence.
  - *Audit Record*: `doc/rnd/note_04e_m1_4_integration_and_audit.md`.

---

### Active Milestone: M1.4.5 — Consolidated Verification & Invariant Hardening
- **Objective**: Execute multi-tiered verification covering kinematic sum rules, semidefinite representability constraints, solvable reductions, discrete selection rules, variational response, and canonical dualities.
- **Parent R&D Specifications**:
  - `doc/rnd/note_04f_advanced_invariants_and_sum_rules_backlog.md`
  - `doc/rnd/note_04g_variational_kkt_and_convergence_asymptotics.md`
  - `doc/rnd/note_05_quasifree_states_and_variational_stability.md`

#### Active Sprint Tasks
- **T1: Kinematic & Semidefinite Gates (Tier 1)** — **COMPLETED**
  - *Target*: `test/unit/test_kinematic_and_semidefinite.py`
  - *Contracts Verified*: Parseval–Frobenius spectral moment identity (`TST-KIN-01`), pointwise $[\mathcal{H}(\mathbf{k}), P(\mathbf{k})] \equiv 0$ (`TST-KIN-02`), projector idempotency (`TST-KIN-03`), local two-site LMI positivity $\lambda_{\min}(\mathcal{C}_{ij}) \ge 0$ (`TST-SDP-01`), operator norm bounds (`TST-SDP-02`), and dipole/parity envelopes (`TST-SDP-03`).

- **T2: Solvable Limits & Symmetry Selection Rules (Tier 2)** — **COMPLETED**
  - *Targets*: `test/unit/test_solvable_limits.py`, `test/unit/test_symmetry_selection_rules.py`
  - *Contracts Verified*: 1D zigzag chain limit $U^{00} = -2/\pi$ (`TST-RED-01`), isolated dimer flat-band saturation with gauge zero modes (`TST-RED-02`), high-field classical polarization and energy convergence (`TST-ASY-01`), time-reversal-odd correlator vanishing (`TST-SYM-01`), sublattice inversion parity (`TST-SYM-02`), and $C_3$ ray permutation invariance (`TST-SYM-03`).

- **T3: Variational Response & Canonical Dualities (Tier 3)** — **IN PROGRESS**
  - *Targets*: `src/kitaev_mf/services/response.py`, `test/integration/test_variational_response.py`, `test/integration/test_dualities_and_asymptotics.py`
  - *Contracts*: Variational Hellmann–Feynman derivative $\frac{\mathrm{d}E}{\mathrm{d}h^\gamma} = -\langle \sigma_{\mathrm{tot}}^\gamma \rangle$ with $\mathcal{O}(\delta^2)$ scaling (`TST-VAR-01`), thermodynamic concavity $\chi^{\gamma\gamma} > 0$ (`TST-VAR-02`), bipartite sign inversion duality (`TST-RED-03`), and symmetric Klein duality with cubic boundary suppression (`TST-DUA-01`).

- **T4: Hardware Parity & Verification Gate** — **PLANNED**
  - *Scope*: Full test suite execution on CPU and CUDA through `TorchBackend`, zero-import-cycle validation, and 100% green pass under `make verify`.

---

### Strategic Backlog (Future Milestones)

#### Phase 1: Analytical Rigor & Domain Foundations (Closing Milestone)
- **M1.5: Standalone Research Monograph**
  - Deliver publication-grade LaTeX deliverable (`doc/report/main.pdf`) synthesizing Parts I–V with end-to-end provenance tags (`CLM-M1-...`).
