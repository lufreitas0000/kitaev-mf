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
  - *Contracts*: Physical spin observables, plaquette Wilson loops, extensive and normalized energy densities, double-counting subtraction $C = E_{\text{exchange}}$, and end-to-end system convergence.
  - *Audit Record*: `doc/rnd/note_04e_m1_4_integration_and_audit.md`.

---

### Active Milestone: M1.4.5 — Consolidated Verification & Invariant Hardening
- **Objective**: Execute multi-tiered verification covering kinematic sum rules, semidefinite representability constraints, solvable reductions, discrete selection rules, and canonical dualities.
- **Parent R&D Specifications**:
  - `doc/rnd/note_04f_advanced_invariants_and_sum_rules_backlog.md`
  - `doc/rnd/note_05_quasifree_states_and_variational_stability.md`

#### Active Sprint Tasks
- **T1: Kinematic & Semidefinite Gates (Tier 1)**
  - *Target*: `test/unit/test_kinematic_and_semidefinite.py`
  - *Contracts*: Parseval–Frobenius spectral moment identity (`TST-KIN-01`), pointwise $[\mathcal{H}(\mathbf{k}), P(\mathbf{k})] \equiv 0$ (`TST-KIN-02`), projector idempotency (`TST-KIN-03`), local two-site LMI positivity $\lambda_{\min}(\mathcal{C}_{ij}) \ge 0$ (`TST-SDP-01`), operator norm bounds (`TST-SDP-02`), and dipole/parity envelopes (`TST-SDP-03`).
  - *Budget*: Zero-iteration linear algebra ($< 15\,\mathrm{ms}$).

- **T2: Solvable Limits & Symmetry Selection Rules (Tier 2)**
  - *Targets*: `test/unit/test_solvable_limits.py`, `test/unit/test_symmetry_selection_rules.py`
  - *Contracts*: 1D zigzag chain limit $\bar{U}_{\mathrm{1D}} = 4/(3\pi)$ (`TST-RED-01`), isolated dimer flat-band saturation (`TST-RED-02`), time-reversal-odd correlator vanishing (`TST-SYM-01`), sublattice inversion parity (`TST-SYM-02`), $C_3$ ray permutation invariance (`TST-SYM-03`), and high-field classical energy convergence (`TST-ASY-01`).
  - *Budget*: Single-shot BZ integrations on coarse grids ($\sim 60\,\mathrm{ms}$).

- **T3: Variational Response & Canonical Dualities (Tier 3)**
  - *Targets*: `test/integration/test_variational_response.py`, `test/integration/test_dualities_and_asymptotics.py`
  - *Contracts*: Numerical Hellmann–Feynman derivative $\frac{\mathrm{d}E}{\mathrm{d}h^\gamma} = -\sum \mathrm{V}_s^{0\gamma}$ (`TST-VAR-01`), thermodynamic concavity $\chi^{\gamma\gamma} \ge 0$ (`TST-VAR-02`), bipartite sign inversion duality (`TST-RED-03`), and four-sublattice Klein duality in the 4-site cell (`TST-DUA-01`).
  - *Budget*: Iterative solver convergence loops ($\sim 800\,\mathrm{ms}$).

- **T4: Hardware Parity & Verification Gate**
  - *Scope*: Full test suite execution on CPU and CUDA through `TorchBackend`, zero-import-cycle validation, and 100% green pass under `make verify`.

---

### Strategic Backlog (Future Milestones)

#### Phase 1: Analytical Rigor & Domain Foundations (Closing Milestone)
- **M1.5: Standalone Research Monograph**
  - Deliver publication-grade LaTeX deliverable (`doc/report/main.pdf`) synthesizing Parts I–V with end-to-end provenance tags (`CLM-M1-...`).

#### Phase 2: Homogeneous Phase Diagrams & Magnetic Ordering
- **M2.1: Multi-Sublattice Order Parameter Framework** (Néel, Ferromagnetic, Zigzag, Stripy order parameter extraction).
- **M2.2: Automated Phase Boundary Tracking** (Pseudo-arclength continuation and gap-closing detection).
- **M2.3: Anisotropic Critical Field Mapping** (Field angle dependence and comparison with DMRG benchmarks).

#### Phase 3: Real-Space Inhomogeneous Vortices & Bound States
- **M3.1: Large-Scale Real-Space Lattice Engine** (Open/periodic clusters up to $N = 3200$ sites).
- **M3.2: Metastable Vortex Configurations** (Localized vison pair stabilization and string branch cuts).
- **M3.3: Induced Polarization & Quadrupole Interactions** (Charge imbalance profiles $\delta n_l$ and vortex interactions).

#### Phase 4: Spectroscopic & Dynamical Observables
- **M4.1: Scanning Tunneling Spectroscopy (STS)** (Local Majorana Green's functions and cotunneling zero-bias peaks).
- **M4.2: Dynamical Spin Structure Factor $S(\mathbf{q}, \omega)$** (Inelastic neutron scattering continuum).
- **M4.3: Thermal Hall Conductance & Chiral Edge States** (Berry curvature Chern numbers and ribbon edge modes).

#### Phase 5: Cloud Production Pipelines & Open Science Archiving
- **M5.1: GPU Cluster Orchestration** (Automated sweeps on NVIDIA L4/A100 instances).
- **M5.2: Open Science Provenance** (DVC pipelines, atomic HDF5 datasets, Zenodo release DOIs).
