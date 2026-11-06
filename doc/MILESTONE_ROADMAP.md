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
- **M1.5: Dual-Track Documentation & Living Lab Book**
  - *Contracts*: Standardized machine-readable R&D notes (`doc/rnd/phase_1/`) with YAML frontmatter, compiled Chapters 1–5 in `memoir` parchment format, and established automated benchmark-to-table pipelines.
- **M1.6: Fixed-Point Solver Algorithms & Multi-Basin Landscape (2-Site & 4-Site Cells)**
  - *Contracts*: Evaluated Picard, DIIS, Broyden, JFNK; verified skew-symmetric retraction A/B stability; mapped first-order crossing $h_c^* = 0.3325$ and forward spinodal $h_{c2}^+ = 0.400$ (`K-011`); proved cubic Majorana gap scaling $\Delta_K(h) \approx 2.592 h^3$ (`K-010`); verified 2-site vs 4-site scale parity ($< 2.2 \times 10^{-4}$); authored Living Lab Book Chapter 6 (42 pages).

---

### Active Milestone Backlog (Immediate Horizon)

- **M1.7: Thermodynamic Response & Phase Diagram of the [111]-Field FM Kitaev Model**
  - *Physical Model*: $K = -1.0, \mathbf{h} \parallel [111]$, scanning $h/|K| \in [0, 0.8]$.
  - *Scope*:
    1. Magnetic susceptibility curve $\chi^{[111]}(h) \equiv \partial \langle \sigma^{[111]} \rangle / \partial h$ via central finite differences across the $h_c^* = 0.3325$ first-order transition.
    2. Topological flux susceptibility $\chi_W(h) \equiv -\partial W_p / \partial h$ tracking vison creation and flux collapse.
    3. Direct Berry curvature integration and Chern number evaluation $\mathcal{C}(h) = \pm 1 \to 0$ of occupied Bogoliubov bands across the transition.
    4. Author Chapter 7 of Living Lab Book (`doc/labbook/chapters/ch_m1_7_fm_kitaev_field_111.tex`).

- **M1.8: Antiferromagnetic Kitaev Model ($K = +1.0$) under [111] Field**
  - *Physical Context*: Unlike the FM model which jumps directly from QSL to polarized state, the AFM Kitaev model under $[111]$ field is widely conjectured to host an exotic intermediate gapless / nematic $U(1)$ spin liquid phase.
  - *Scope*:
    1. Map the phase diagram for $K = +1.0$ under $\mathbf{h} \parallel [111]$.
    2. Test for intermediate gapless regimes, vison condensation, and emergent Fermi surfaces.
    3. Document in Living Lab Book Chapter 8 (`ch_m1_8_afm_kitaev_field_111.tex`).

---

## Phase 2: Correlated Partons, Projection Corrections & Gauge Fluctuations (M2)

- **M2.1: Gutzwiller Projection Gap & Finite-Cluster Pfaffians (BL-01)**
  - *Scope*: Implement Pfaffian Monte Carlo on small clusters ($N = 6, 8, 12, 16$) comparing exact diagonalization (ED), Gaussian mean-field $\langle \boldsymbol{\sigma} \rangle_{\mathrm{G}}$, and projected expectations $\langle \boldsymbol{\sigma} \rangle_{\mathrm{phys}}$. Determine the quantitative projection gap $\Delta_{\mathrm{proj}}(\boldsymbol{\lambda})$.
- **M2.2: Localized $\mathbb{Z}_2$ Vortex / Vison Dynamics & Green's Functions (BL-02, BL-05)**
  - *Scope*: Direct simulation of localized vison pairs ($W_p = -1$) on supercells, calculation of the vison excitation gap $\Delta_{\mathrm{vison}}$ via triangular lattice Green's functions, and static gauge susceptibility.
- **M2.3: Variational 2-RDM Semidefinite Programming Lower Bounds (BL-06)**
  - *Scope*: Implement 2-body reduced density matrix (2-RDM) semidefinite programming to compute rigorous, variational-independent energy lower bounds $E_{\mathrm{SDP}} \le E_{\mathrm{exact}}$.

---

## Phase 3: Thermal Transport, Dynamics & Multi-Parametric Horizons (M3)

- **M3.1: Kitaev-Heisenberg $J\text{-}K\text{-}\Gamma\text{-}\Gamma'$ Global Phase Diagram**
  - *Scope*: Benchmark the full parameter space $\mathbf{p} = (J, K, \Gamma, \Gamma', D, \mathbf{h}) \in \mathbb{R}^8$, mapping Néel, Stripy, Zigzag, and canted magnetic domains against literature standards.
- **M3.2: Quantized Thermal Hall Conductivity & Finite-Temperature Majorana Transport**
  - *Scope*: Finite-temperature Fermi-Dirac occupation and calculation of the half-integer quantized thermal Hall effect $\kappa_{xy}(T) / T \to \frac{\pi}{12} \frac{k_B^2}{\hbar} \mathcal{C}$.
- **M3.3: Dynamical Spin Structure Factor $S(\mathbf{q}, \omega)$**
  - *Scope*: Compute dynamical response functions via random phase approximation (RPA) and two-particle bubble convolutions to simulate inelastic neutron scattering (INS).

---

*(Note: Milestone M1.9 [Publication Manuscript] has been deferred indefinitely per research steering decision, prioritizing core physical breakthroughs across M1.7, M1.8, and Phase 2).*
