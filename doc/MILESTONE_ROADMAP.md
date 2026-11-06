# Project Roadmap & Milestone Backlog: kitaev-mf

## Taxonomy
- **Phase**: Macro physical and engineering regimes.
- **Milestone (`M<Phase>.<Index>`)**: Strategic deliverable gates and persistent backlog units.
- **Task (`T<Index>`)**: Operational execution units active exclusively within the current milestone sprint.

---

## Phase 1: Analytical Foundations, Algorithmic Scaling & Pure Kitaev Benchmarking

Phase 1 establishes the mathematical foundations, verified algebraic invariants, rigorous solver benchmarking, and numerical scaling laws for pure ferromagnetic Kitaev honeycomb magnets ($K = -1.0$) under magnetic field $\mathbf{h} \parallel [111]$.

### Completed Foundations (M1.1 – M1.6 Pilot)
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
- **M1.5: Dual-Track Documentation & Living Lab Book**
  - *Contracts*: Standardized machine-readable R&D notes (`doc/rnd/phase_1/`) with YAML frontmatter, compiled Chapters 1–5 in `memoir` parchment format, and established automated benchmark-to-table pipelines.
- **M1.6 Pilot: Fixed-Point Solvers & Multi-Basin Landscape (Chapter 6)**
  - *Contracts*: Evaluated Picard, DIIS, Broyden, JFNK; mapped first-order crossing $h_c^* = 0.3325$ and forward spinodal $h_{c2}^+ = 0.400$ (`K-011`); documented power-law Majorana gap fit $\Delta_K(h) \approx 2.592 h^{2.994}$ (`K-010`); authored Living Lab Book Chapter 6.

---

### Active Sprints & Phase 1 Execution Gates (Overnight Sequence)

- **M1.6: Formal Statistical Scorecards, Contraction Theory & Concrete FLOP Profiler (Chapter 7)**
  - *Focus*: Rigorous Computer Science and statistical characterization of all solver paradigms on FM Kitaev $[111]$.
  - *Deliverables*:
    1. Pre-experiment Front-1 notes: `m1_6_t6` (statistical metrics & Wilson CIs), `m1_6_t7` (solver mechanics & contraction theory), `m1_6_t8` (concrete FLOP & memory methodology), `m1_6_t9` (saddle trapping vs. critical slowing down).
    2. Infrastructure: `src/kitaev_mf/infra/profiler.py` (BLAS/LAPACK FLOP counter) and `LimitedMemoryBroydenSolver` (`L-Broyden`, $\mathcal{O}(mD)$ space/FLOPs).
    3. Hardware micro-calibration: `scripts/m1_6/calibrate_hardware_flops.py` on Intel i5-13500 CPU.
    4. Statistical benchmark: `scripts/m1_6/statistical_solver_benchmark.py` ($N_{\mathrm{trials}}=40$ $C_3$/sublattice-symmetric seeds $\times$ 8 perturbation amplitudes $\delta \in [0.005, 0.50]$ $\times$ 6 solvers at $h=0.40$ and $h=0.20$).
    5. Post-experiment evaluation: `m1_6_t10_statistical_benchmark_evaluation.md`.
    6. Living Lab Book Chapter 7 (`doc/labbook/chapters/ch_m1_7_solver_primer_and_scorecard.tex`).

- **M1.7: Bloch Solvers on $N_c$-Site Supercells (Chapter 8)**
  - *Focus*: Scaling and symmetry robustness in momentum-space supercells.
  - *Deliverables*:
    1. Front-1 pre-experiment note: `m1_7_t1_supercell_bloch_scaling_contracts.md`.
    2. General supercell generator: `make_honeycomb_supercell(n1, n2)` in `src/kitaev_mf/domain/lattice.py` supporting $n_c = 2 n_1 n_2 \in \{2, 4, 8, 16, 32\}$.
    3. Benchmark script: `scripts/m1_7/benchmark_bloch_supercells.py` verifying $\mathcal{O}(N_k n_c^3)$ scaling, intensive energy parity, and unconstrained generic random seed symmetry breaking in 4-site and 8-site cells.
    4. Evaluation note: `m1_7_t2_supercell_benchmark_evaluation.md`.
    5. Living Lab Book Chapter 8 (`doc/labbook/chapters/ch_m1_8_supercell_bloch_scaling.tex`).

- **M1.8: Real-Space Torus Clusters — Track A (Dense Real/Complex) & Track B (Eigenvector-Free Purification) (Chapter 9)**
  - *Focus*: Large real-space clusters ($L \in \{6, 12, 18, 24\}$, $N_{\mathrm{sites}} = 2L^2 \in \{72, 288, 648, 1152\}$, $M = 8L^2 \in \{288, 1152, 2592, 4608\}$) comparing dense eigensolvers vs. matrix sign purification.
  - *Track A (Original Baseline)*: Real skew-symmetric algebra (`float64`) vs. complex Hermitian (`complex128` `zheevd`), testing the $4\times$ FLOP and $2\times$ memory reduction.
  - *Track B (Parallel Purification Track)*: Evaluating $\Phi$ without eigenvectors via matrix sign purification $\Omega = \mathrm{i}\operatorname{sgn}(\mathrm{i} A)$ using cold- and warm-started Newton–Schulz / McWeeny ($Y_{k+1} = \frac{1}{2}Y_k(3I + Y_k^2)$) and SP2 iterations.
  - *Deliverables*:
    1. Front-1 theory notes: `m1_8_t1` (Track A real arithmetic) and `m1_8_t2` (Track B matrix sign purification).
    2. Implementations: `src/kitaev_mf/services/real_space.py` and `src/kitaev_mf/services/purification.py`.
    3. Benchmark script: `scripts/m1_8/benchmark_real_space_tracks.py`.
    4. Evaluation note: `m1_8_t3_real_space_tracks_evaluation.md`.
    5. Living Lab Book Chapter 9 (`doc/labbook/chapters/ch_m1_9_real_space_and_purification.tex`).

- **M1.9: Front 3 "Shop Window" Executive Manuscript & Phase 1 Certification**
  - *Focus*: High-impact, persuasive publication-grade paper (`doc/manuscript/main.pdf`) communicating the core breakthroughs of Phase 1 to external researchers.
  - *Deliverables*:
    1. Front-3 modular sections in `doc/manuscript/sections/`: `sec01_introduction.tex`, `sec02_algebraic_foundations_and_invariants.tex`, `sec03_field_response_and_multistability.tex`, `sec04_grassmannian_solvers_and_real_space_scaling.tex`, `sec05_conclusion_and_outlook.tex`.
    2. Compiled publication artifact: `doc/manuscript/main.pdf`.
    3. Full verification and governance sign-off: `make verify && make labbook && make manuscript`.

---

## Future Strategic Phases & Milestone Backlog

Following the completion of Phase 1, each subsequent phase focuses systematically on specific physical models and phenomena.

### Phase 2: AFM Kitaev & Intermediate Quantum Spin Liquids (M2)
- **M2.1: Antiferromagnetic Kitaev Model ($K = +1.0$) under [111] Field**
  - *Scope*: Detailed phase diagram scan for intermediate gapless / nematic spin liquids, vison condensation, and emergent Fermi surfaces.
- **M2.2: Gutzwiller Projection Gap & Finite-Cluster Pfaffians**
  - *Scope*: Pfaffian Monte Carlo on small clusters ($N = 6, 8, 12, 16$) comparing exact diagonalization (ED), Gaussian mean-field $\langle \boldsymbol{\sigma} \rangle_{\mathrm{G}}$, and projected expectations $\langle \boldsymbol{\sigma} \rangle_{\mathrm{phys}}$. Determine the quantitative projection gap $\Delta_{\mathrm{proj}}(\boldsymbol{\lambda})$.
- **M2.3: Variational 2-RDM Semidefinite Programming Lower Bounds**
  - *Scope*: Implement 2-body reduced density matrix (2-RDM) semidefinite programming for rigorous energy lower bounds $E_{\mathrm{SDP}} \le E_{\mathrm{exact}}$.

### Phase 3: Topological Defects, Visons & Real-Space Textures (M3)
- **M3.1: Isolated Vison Pairs & Triangular Lattice Green's Functions**
  - *Scope*: Real-space insertion of localized vison pairs ($W_p = -1$) on large torus clusters, computing the vison bound-state gap $\Delta_{\mathrm{vison}}$ and static gauge susceptibility.
- **M3.2: Field-Induced Skyrmionic and Dimer Textures in Extended Cells**
  - *Scope*: Systematic study of broken-translation and broken-rotational ground states under tilted fields $\mathbf{h} \parallel [001], [110]$.

### Phase 4: Non-Integrable Materials & Multi-Parametric Horizons (M4)
- **M4.1: Kitaev-Heisenberg $J\text{-}K\text{-}\Gamma\text{-}\Gamma'$ Global Phase Diagram**
  - *Scope*: Benchmark the full parameter space $\mathbf{p} = (J, K, \Gamma, \Gamma', D, \mathbf{h}) \in \mathbb{R}^8$, mapping Néel, Stripy, Zigzag, and canted magnetic domains against literature standards.
- **M4.2: Quantized Thermal Hall Conductivity & Finite-Temperature Majorana Transport**
  - *Scope*: Finite-temperature Fermi-Dirac occupation and calculation of the half-integer quantized thermal Hall effect $\kappa_{xy}(T) / T \to \frac{\pi}{12} \frac{k_B^2}{\hbar} \mathcal{C}$.
- **M4.3: Dynamical Spin Structure Factor $S(\mathbf{q}, \omega)$**
  - *Scope*: Compute dynamical response functions via random phase approximation (RPA) and two-particle bubble convolutions to simulate inelastic neutron scattering (INS).
