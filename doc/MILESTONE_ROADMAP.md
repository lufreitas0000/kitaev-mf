# Project Roadmap & Milestone Backlog: kitaev-mf

## Taxonomy
- **Phase**: Macro physical and engineering regimes.
- **Milestone (`M<Phase>.<Index>`)**: Strategic deliverable gates and persistent backlog units.
- **Task (`T<Index>`)**: Operational execution units active exclusively within the current milestone sprint.

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
- **M1.4.5: Consolidated Verification & Invariant Hardening**
  - *Contracts*: Verified TST-KIN-01..03, TST-SDP-01..03, TST-RED-01..03, TST-SYM-01..03, TST-ASY-01, TST-VAR-01..02, TST-DUA-01 across 103 unit and integration tests.
  - *Artifacts*: `data/processed/benchmarks/k006_quadrature.dat`, `hellmann_feynman_scaling.dat`, `susceptibility_curve.dat`, `klein_duality_splitting.dat`.

---

### Active Milestone: M1.5 — Dual-Track Documentation & Living Lab Book
- **Objective**: Standardize machine-readable R&D notes (`doc/rnd/phase_1/`) with mandatory YAML frontmatter, scaffold the human-first living LaTeX labbook (`doc/labbook/`) using `memoir` with soft-beige/charcoal aesthetics, and establish automated benchmark-to-table pipelines.
- **Active Sprint Tasks**:
  - **T1: R&D Note Reorganization & YAML Metadata Injection**
    - Target: `doc/rnd/phase_1/m1_*_t*_<slug>.md`
    - Scope: Rename unstructured notes to deterministic phase/milestone/task filenames; inject standardized YAML headers.
  - **T2: Living Lab Book Scaffolding (Memoir, Beige/Charcoal Theme)**
    - Target: `doc/labbook/main.tex`, `preamble.tex`, `Makefile`
    - Scope: Configure `memoir` layout, `#FAF6EE` background, `#2B2B2B` charcoal font, and `pgfplotstable` integration.
  - **T3: Cumulative Synthesis of Completed Milestones (Chapters M1.1 to M1.4.5)**
    - Target: `doc/labbook/chapters/ch_m1_1_*.tex` through `ch_m1_4_5_*.tex`
    - Scope: Synthesize derivations, failure modes, and benchmark plots directly from `data/processed/benchmarks/*.dat`.

---

### Concrete Sprint Backlog (Next Milestones)

- **M1.6: Fixed-Point Solver Algorithms & Basin Multistability (2-Site Cell, FM Kitaev + [111] Field)**
  - *Physical Model*: Ferromagnetic Kitaev model $K = -1.0, J = \Gamma = \Gamma' = D = 0$, external Zeeman field $\mathbf{h} = \frac{h}{\sqrt{3}}(1, 1, 1)$.
  - *Scope*:
    1. Implement and benchmark four convergence algorithms on the Grassmannian manifold: Damped Picard, Joint DIIS, Broyden's Quasi-Newton, and Jacobian-Free Newton-Krylov (JFNK).
    2. Map basin of attraction across distinct seeds: QSL zero-flux seed, field-polarized seed, random Grassmannian states $\Omega_{\text{rand}} = Q\Omega_0 Q^T$, and adiabatic continuation sweeps ($h: 0 \to 1$ vs. $h: 1 \to 0$).
    3. Document fixed-point contraction rates, spectral gap closing $\Delta_g(h)$, and algorithmic efficiency in `doc/labbook/chapters/ch_m1_6_solvers_and_basins.tex`.

- **M1.7: Phase Diagram of the [111]-Field FM Kitaev Model (2-Site $\to$ 4-Site $\to$ 8-Site Cells)**
  - *Physical Model*: $K = -1.0, \mathbf{h} \parallel [111]$, scanning $h/|K| \in [0, 0.8]$ via adaptive bisection based on susceptibility curvature $\chi^{[111]}(h)$.
  - *Scope*:
    1. Map energy density $E_{\text{site}}(h)$, magnetization $\langle \sigma^{[111]} \rangle(h)$, susceptibility $\chi^{[111]}(h)$, Wilson loop $W_p(h)$, and hybrid hoppings $U_b^{0\gamma}(h)$ on the 2-site cell.
    2. Scale to 4-sublattice rectangular and 8-sublattice cells to test for spontaneous translational symmetry breaking or multi-sublattice magnetic ordering near the critical field $h_c$.
    3. Record verified boundaries, phase transition order, and figures in `doc/labbook/chapters/ch_m1_7_fm_kitaev_field_111.tex`.

---

### Adaptive Research Backlog (Post-M1.7 Horizons)

The following trajectories probe the parameter space $\mathbf{p} = (J, K, \Gamma, \Gamma', D, \mathbf{h}) \in \mathbb{R}^8$ along physically motivated, literature-benchmarked manifolds:

- **M1.8.1: Antiferromagnetic Kitaev Model ($K = +1.0$) under [111] Field**
  - Investigate the intermediate gapless/chiral phase prior to polarization.
- **M1.8.2: Kitaev-Heisenberg $J\text{-}K$ Trajectory ($\theta \in [0, 2\pi)$)**
  - Benchmark Néel, Stripy, Zigzag, and QSL boundaries against literature benchmarks (Chaloupka, Jackeli, Khaliullin 2010; Rau et al. 2014).
- **M1.8.3: Off-Diagonal Exchange $\Gamma$ and Dzyaloshinskii-Moriya $D$ Interactions**
  - Canted zigzag domains, magnetic torque response, and non-coplanar structures.
- **M1.8.4: Plaquette Flux-Sector Inversion ($\mathbb{Z}_2$ Vortex Solver)**
  - Direct simulation of vison-pair and localized vortex states ($W_p = -1$) on supercells.
- **M1.8.5: Arbitrary Field Directions & Multi-$\mathbf{Q}$ Star Orders**
  - Trajectories along [001] and [110]; triple-$\mathbf{Q}$ states in 8-sublattice geometries.
- **M1.9: Final Phase 1 Research Monograph (`doc/report/main.pdf`)**
  - Publication-grade synthesis distilling Phase 1 into formal definitions, propositions, and certified phase diagrams with persistent claim tags (`CLM-M1-*`).
