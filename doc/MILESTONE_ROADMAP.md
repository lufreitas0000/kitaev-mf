# Project Roadmap & Milestone Backlog: kitaev-mf

## Taxonomy
- **Phase**: Macro model-family and methodological regimes.
- **Milestone (`M<Phase>.<Index>`)**: Strategic deliverable gates and persistent backlog units.
- **Task (`T<Index>`)**: Operational execution units active exclusively within the current milestone sprint.

---

## Phase 1: Pure Ferromagnetic Kitaev ($K = -1.0$) in $[111]$ Field — Numerical Quantification & Algorithmic Foundations (COMPLETED & CERTIFIED)

Phase 1 established the mathematical foundations, verified algebraic invariants (`K-001`–`K-023`), rigorous solver benchmarking, failure-mode taxonomy, and hardware resource scaling laws for the pure ferromagnetic Kitaev honeycomb model ($K = -1.0$) under magnetic field $\mathbf{h} \parallel [111]$. All 11 milestones below are completed and certified across 153 automated tests.

### Completed Milestones (M1.1 – M1.11)
- **M1.1: Local Algebra & Exact Parton Representation (Chapters 1–2)**
  - *Contracts*: A01–A20 verified across 25 unit tests (`K-001`, `K-002`, `K-003`).
- **M1.2: Bloch-Majorana Formalism & Folding Parity (Chapter 3)**
  - *Contracts*: Exact Kitaev dispersion (`K-004`) and 2-site to 4-site supercell folding parity (`K-005`).
- **M1.3: Reciprocal Spectral Projection & DIIS Solver (Chapter 4)**
  - *Contracts*: Zero-temperature spectral projector, BZ integration, and Joint Anderson Acceleration. High-resolution reference quadrature $\bar{U} \approx 0.5248657458506$ (`K-006`).
- **M1.4: Production Observables & Gaussian Energetics (Chapter 5)**
  - *Contracts*: Physical spin observables, plaquette Wilson loops, extensive and density energetics, scalar double-counting subtraction $C = E_{\text{exchange}}$ (`K-007`), and end-to-end system convergence.
- **M1.4.5: Consolidated Verification & Invariant Hardening (Chapter 6)**
  - *Contracts*: Verified `TST-KIN-01..03`, `TST-SDP-01..03`, `TST-RED-01..03`, `TST-SYM-01..03`, `TST-ASY-01`, `TST-VAR-01..02`, `TST-DUA-01` (`K-008`, `K-009`).
- **M1.5: Dual-Track Documentation & Living Lab Book**
  - *Contracts*: Standardized machine-readable R&D notes (`doc/rnd/phase_1/`) with YAML frontmatter, `memoir` parchment Living Lab Book (`doc/labbook/labbook_kitaev_mf.pdf`), and automated `.dat` $\to$ PGFPlots figure pipelines.
- **M1.6: Fixed-Point Solvers, Statistical Scorecards & Multi-Basin Landscape (Chapter 7)**
  - *Contracts*: Evaluated Picard, Joint DIIS, Full Broyden, L-Broyden, and JFNK across $N_{\mathrm{trials}}=30$ symmetric perturbations with 95% Wilson CIs (`K-012`, `K-013`); calibrated host hardware FLOP roofline (`K-014`); mapped first-order crossing $h_c^* = 0.3325$ and forward spinodal $h_{c2}^+ = 0.400$ (`K-011`); proved cubic Majorana gap law $\Delta_K(h) = 2.5921 h^3$ (`K-010`) and reconciled it with Kitaev's 3rd-order perturbation theory (`K-019`); audited post-spinodal KKT saddle excursions and coexistence domain-wall tension $\sigma_{\mathrm{DW}} = +0.042|K| > 0$ (`K-023`).
- **M1.7: Bloch Solvers on $N_c$-Site Supercells (Chapter 8)**
  - *Contracts*: General supercell generator `make_honeycomb_supercell(n1, n2)` for $n_c \in \{2, 4, 8, 16, 32\}$; verified $\mathcal{O}(N_k n_c^3)$ scaling, exact intensive energy parity $\le 5.55 \times 10^{-16}$ (`K-015`), and universal isotropic ground-state recovery $P_{\mathrm{iso}} = 1.000$ (`K-016`).
- **M1.8: Real-Space Torus Clusters — Track A (Real Skew) & Track B (Eigenvector-Free Purification) (Chapter 9)**
  - *Contracts*: Track A real skew-symmetric algebra ($2\times$ memory reduction, `K-017`) and diagnosis of $3 \mid L$ Dirac resonance squaring ($\kappa(\mathbf{S}) \approx 3.6 \times 10^{31}$); Track B Newton–Schulz matrix sign purification $\Omega = \mathrm{i}\operatorname{sgn}(\mathrm{i}\mathbf{A})$ achieving $63.53\times$ speedup over eigensolvers and $104.5\text{--}114.6\text{ GFLOP/s}$ sustained `dgemm` throughput (`K-018`).
- **M1.9: Publication Overhaul, 3-Agent Adversarial Peer Review & Hypothesis Contracts (Chapter 11)**
  - *Contracts*: Authored SciPost executive monograph (`doc/manuscript/scipost_kitaev_mf.pdf`), passed 3-round adversarial peer review, and established the centralized Hypothesis Contracts Registry (`doc/HYPOTHESIS_CONTRACTS.md` and Lab Book Chapter 11).
- **M1.10: Real-Space Vortices, Vison Pairs & In-Gap Majorana Bound States (Chapter 10)**
  - *Contracts*: Dual triangular lattice minimal branch-cut generator for isolated $\mathbb{Z}_2$ vortex pairs ($W_{p_1} = W_{p_2} = -1$), monotonic attractive two-vison potential $V(R) \in [0.114, 0.168]|K|$, in-gap Majorana bound modes ($\epsilon_1 / \Delta_{\mathrm{bulk}} \approx 13\%$), and self-consistent defect trapping (`K-020`).
- **M1.11: Wave 3 Advanced Solvers, Disordered Defect Clusters & Master Solver Resource Table**
  - *Contracts*: Self-consistent purification mechanics with outer DIIS and damped L-Broyden ($\alpha = 0.5$, `K-021`), quenched bond disorder and site vacancy zero-mode conditioning ($\kappa(\mathbf{A}) \approx 3.82 \times 10^4$, `K-022`), inhomogeneous gauge screening ($\xi_{\mathrm{gauge}} \approx 1.1 a$), and the **Phase 1 Master Solver Selection & Hardware Resource Table** (`HC-P1-01..03`).

---

## Phase 2: Pure Ferromagnetic Kitaev ($K = -1.0$) in $[111]$ Field — Low-Energy Physics, Phase Diagram, Vortices, Thermodynamics & Dynamical Correlations (ACTIVE PREPARATION)

Phase 2 is the direct physical continuation of Phase 1 in the **identical model regime** ($K = -1.0$, $\mathbf{h} \parallel [111]$). Deploying the certified optimal solvers from Phase 1 (**Joint DIIS** for reciprocal-space Bloch cells and **Newton–Schulz purification + damped L-Broyden $\alpha=0.5$** for large real-space clusters), Phase 2 focuses entirely on extracting the complete low-energy physics and observables across four physical settings (`HC-P2-01` through `HC-P2-05`):

- **M2.1: $T = 0$ Bulk Ground-State Phase Diagram, Complete Observable Atlas & Variational Bounds (`HC-P2-01`, `HC-P2-02`)**
  - *Scope*: High-resolution mapping of the $T=0$ phase diagram across $h \in [0, 0.60]$ using Joint DIIS ($N_k = 36^2\text{--}64^2$):
    1. Complete observable atlas: sublattice magnetization $\mathbf{m}_s(h)$, magnetic susceptibility tensor $\chi^{\alpha\beta}(h) = \partial m^\alpha / \partial h^\beta$, matter and gauge bond expectations $\mathrm{U}_b^{\mu\nu}(h)$, spectral gap $\Delta(\mathbf{k}, h)$, and bulk Chern number $\nu \in \{-1, 0\}$ via Fukui–Hatsugai–Suzuki lattice Berry curvature.
    2. Comparison of factored 6-link Wilson loop $W_p(h)$ versus exact 12-Majorana Pfaffian $\mathcal{W}_p(h) = \operatorname{Pf}(\Omega|_{\partial p})$ (`H2` / `HC-P2-01`).
    3. Gutzwiller projection norm overlap $\mathcal{N}_G(h) = \operatorname{Tr}(\rho_{\mathrm{MF}} P_G)$ on finite ED clusters ($N \le 24$) and global 2-RDM semidefinite programming ($D, Q, G$) lower bounds $E_{\mathrm{SDP}}^{(2)} \le E_0 \le E_{\mathrm{MF}}$ (`H1`, `H3` / `HC-P2-02`).

- **M2.2: $T = 0$ States with Isolated Pairs of $\mathbb{Z}_2$ Vortices (`HC-P2-03`)**
  - *Scope*: Large real-space torus simulations ($L \in \{12, 16, 24\}$, $M \in \{1152, 2048, 4608\}$) using Track B purification + damped L-Broyden ($\alpha=0.5$):
    1. Full self-consistent two-vison interaction potential $V(R, h)$ across separations $R \in [1, L/2]$ and fields $h \in [0, h_c^*)$.
    2. Non-Abelian Majorana zero-mode hybridization splitting $\delta\epsilon_0(R) \propto \exp(-R/\xi(h)) \cos(k_F R + \phi)$ and comparison against the bulk coherence length $\xi_{\mathrm{bulk}}(h) = v_F / \Delta_K(h)$.
    3. Spatial profile of local spin magnetization $\delta\mathbf{m}_i$, bond energy density $\delta E_b$, and inhomogeneous gauge multiplier screening cloud $\delta\boldsymbol{\lambda}_i$ ($\xi_{\mathrm{gauge}} \approx 1.1 a$) around vortex cores (`H11` / `HC-P2-03`).

- **M2.3: Finite-Temperature ($T > 0$) States & Thermal Crossover Thermodynamics (`HC-P2-04`)**
  - *Scope*: Self-consistent finite-temperature Fermi–Dirac covariance solver $\Omega_T = \mathrm{i}\tanh(\mathrm{i}\beta\mathbf{A}/2)$ minimizing the variational Helmholtz free energy $\mathcal{F}(\Omega_T, \boldsymbol{\lambda}; T) = \mathcal{E}(\Omega_T, \boldsymbol{\lambda}) - T S_{\mathrm{FD}}(\Omega_T)$:
    1. Two-stage specific heat $C_v(T)$ and entropy release $S(T)$ from $S(\infty) = \ln 2$ per spin across the high-temperature itinerant Majorana crossover ($T_H \approx 0.375|K|$) and low-temperature gauge/vison scale ($T_L$).
    2. Thermal melting of the first-order spinodal hysteresis loop $[h_{c1}^-(T), h_{c2}^+(T)]$ terminating at a finite-temperature critical end point $(h_{\mathrm{CEP}}, T_{\mathrm{CEP}})$.

- **M2.4: Static and Dynamical Spin Correlations $S^{\alpha\beta}(\mathbf{q})$ and $S^{\alpha\beta}(\mathbf{q}, \omega)$ (`HC-P2-05`)**
  - *Scope*: Connected four-Majorana Wick contraction and Lehmann spectral representation across the $(h, T, \mathbf{q}, \omega)$ manifold:
    1. Equal-time static structure factor $S^{\alpha\beta}(\mathbf{q})$ across the Brillouin zone in the $T=0$ ground state, two-vortex state, and $T>0$ thermal ensemble.
    2. Dynamical spin structure factor $S^{\alpha\beta}(\mathbf{q}, \omega)$, verifying the exact sum rule $\int \frac{\mathrm{d}\omega}{2\pi} S^{\alpha\alpha}(\mathbf{q}, \omega) = S^{\alpha\alpha}(\mathbf{q})$, the sharp single-vison + matter-fermion excitation gap threshold at $h=0$, and field-induced magnon-like collective modes in the polarized phase ($h > h_c^*$).

---

## Subsequent Strategic Phases (Organized by Model Family)

Once the pure ferromagnetic Kitaev $[111]$ physics program of Phase 2 is complete, subsequent phases will be developed sequentially—one dedicated phase per model family:

### Phase 3: Antiferromagnetic Kitaev Honeycomb Model ($K = +1.0$) under Magnetic Field
- **M3.1:** Intermediate field-induced gapless quantum spin liquid phase ($h_{c1} < h < h_{c2}$), emergent neutral spinon Fermi surfaces, and U(1) gauge fluctuations.
- **M3.2:** Vison energetics, chiral gap scaling, and dynamical structure factors in the AFM Kitaev regime.

### Phase 4: Extended Non-Integrable Honeycomb Models ($J\text{--}K\text{--}\Gamma\text{--}\Gamma'\text{--}D$)
- **M4.1:** Competing zigzag, stripy, Néel, and 120° magnetic orders on $n_c \in \{4, 8, 16, 32\}$ supercells (`H5`).
- **M4.2:** Field-induced quantum spin liquid stability windows and inelastic neutron scattering modeling for candidate materials ($\alpha\text{-RuCl}_3$, $\text{Na}_2\text{IrO}_3$).

### Phase 5: Three-Dimensional Tricoordinated Kitaev Lattices
- **M5.1:** Hyperhoneycomb and stripyhoneycomb 3D Kitaev models with Majorana nodal loops and Fermi surfaces using Track B matrix sign purification.
