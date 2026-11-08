# Changelog

All notable user-visible changes to `kitaev-mf` are documented here.

The package version is maintained in `pyproject.toml`. See
[Rule 04: Versioning, Changelog, and Releases](.agents/rules/04_versioning_and_releases.md)
for the SemVer policy and release procedure.

## [Unreleased]

### Added
- Separate finite chiral-matter transport with gauge, spectral, overlap, and
  closure diagnostics; retain the legacy mean-field braiding API.
- Isolated, explicitly unaccepted Gaussian quadrature prototype and byte-faithful
  source manifest for the PI-supplied perturbation bundle.

### Fixed
- Include the already-declared matplotlib runtime dependency in CI lockfiles so
  plotting tests can be collected in a clean environment.

### Changed
- Retire the unsupported mean-field enclosing-loop guarantee from required CI by
  PI decision, preserving its unchanged failing assertions as an opt-in research
  check. This does not accept the old physical claim.

### Changed
- Audited perturbative foundation claims against the original contracts and pinned
  primary source; retained conditional finite algebra and recorded bounded
  diagnostic mismatches without repairing or recertifying production routines.

- Collected historical research on the provisional exploration line with explicit
  acceptance status; isolated document publication from ordinary verification.

- Added a concrete Git/PR/data workflow and preservation checkpoints; verified
  original perturbation notes remain unchanged at compared refs and corrected the
  earlier concern from possible alteration to neglect of the original contracts.

- Established model-agnostic main/exploration PR routing, wave/task planning,
  mandatory task intake, multiple classification tags, and claim-level acceptance.
- Recorded PI disapproval of historical Phase 4 calculations and scheduled
  independent recovery of earlier exact-perturbation work.
- Specified one canonical labbook and optional incremental review notes;
  neutral PRL manuscript migration remains a separate task.

## [1.0.0] - 2026-11-06

### Phase 1: Analytical Rigor & Algorithmic Foundations (Certified Closure)

#### Added
- **Algebraic Foundations & CAR Algebra (M1.1):** Verified local 4-Majorana CAR Clifford algebra, Casimir invariants, Dirac fermion representation, and single-generator parity $\mathcal{D} = -c^0 c^1 c^2 c^3$ across 25 unit tests (`K-001`--`K-003`).
- **Bloch Momentum & Supercell Folding (M1.2):** Explicit 2-site and 4-site reciprocal Bloch Hamiltonian matrices $\mathcal{H}(\mathbf{k})$, exact pure-Kitaev matter and flat gauge dispersions (`K-004`), and machine-precision unitary folding equivalence matrix $W(\mathbf{q})$ (`K-005`).
- **Spectral Projector & DIIS Fixed-Point Solver (M1.3):** Zero-temperature negative-energy spectral projector $P(\mathbf{k}) = 2\Pi^-(\mathbf{k})$, real-space correlator quadrature, and Joint Anderson Acceleration (DIIS) with relative scale-invariant Tikhonov regularization $\epsilon_{\mathrm{reg}} = 10^{-12}$.
- **Reference Quadrature Invariant (K-006):** Analytical proof of matter bond correlator $3\bar{U} = W_3(1) \approx 1.5745972375506$ via 3-step Pearson planar random walk and Richardson-extrapolated BZ quadrature.
- **Production Observables & Energetics (M1.4):** Local spin magnetization $m_s^a = \mathrm{V}_s^{0a}$, two-site spin correlators, factored plaquette Wilson loops $W_p = \prod u_b$, and scalar double-counting subtraction theorem $E_{\mathrm{total}} = \langle H_{\mathrm{quad}} \rangle_{\mathrm{G}} - E_{\mathrm{exchange}}$.
- **Kinematic Invariants & Semidefinite Gates (M1.4.5):** Parseval-Frobenius spectral sum rule, two-site semidefinite representability LMI $\mathcal{C}_{ij} \succeq 0 \iff \|\mathrm{U}_{ij}\|_2 \le 1.0$, and cubic Klein duality boundary splitting law $\delta E \propto \epsilon^3$.
- **Statistical Solver Scorecards & Geometry (M1.6):** Nonlinear fixed-point formulation on Grassmannian fiber bundle $\Gr_{\mathbb{C}}(2N_{\mathrm{sub}}, 4N_{\mathrm{sub}})$, first-order crossing $h_c^* = 0.3325$, forward spinodal $h_{c2}^+ = 0.400$ with critical slowing down $\mu_k \to 1^-$, and analytical proof of topological Majorana Dirac gap $\Delta_K \propto h^3$.
- **Limited-Memory Broyden (L-Broyden):** Implemented $\mathcal{O}(mD)$ limited-memory secant solver with vector accumulation and adaptive damping $\alpha = 0.5$.
- **Generalized Bloch Supercells (M1.7):** Scalable supercell generator `make_honeycomb_supercell(n1, n2)` across $n_c \in \{2, 4, 8, 16, 32\}$, verifying intensive energy scale parity to $5.55 \times 10^{-16}$ and absence of spontaneous nematicity ($N_{\mathrm{nem}} < 7.3 \times 10^{-8}$).
- **Real-Space Torus Clusters & Pure Real Arithmetic (M1.8 Track A):** Exact real skew-symmetric algebra $\Omega = -\mathbf{A}(-\mathbf{A}^2)^{-1/2}$, cutting memory by $2.0\times$ ($8 M^2$ bytes vs $16 M^2$ bytes) with $\frac{26}{3}M^3$ total real FLOPs (`K-017`).
- **Eigenvector-Free Matrix Sign Purification (M1.8 Track B):** Cubic Newton-Schulz matrix sign purification $\Omega = \ii \operatorname{sgn}(\ii\mathbf{A})$ in pure BLAS-3 `dgemm`, delivering a **$2.31\times$ mean speedup** ($N_{\mathrm{repeats}} = 5$) over multi-threaded divide-and-conquer eigensolvers at $M=1152$ ($613.0\text{ ms}$ vs $1418.6\text{ ms}$, and up to $63.53\times$ against worst-case deflation stalls) with $109.73\text{--}114.64\text{ GFLOP/s}$ sustained throughput and machine-precision parity ($4.01 \times 10^{-15}$, `K-018`).
- **Real-Space Vortex Geometry & Bound States (M1.10):** Dual-triangular lattice branch-cut generator for isolated $\mathbb{Z}_2$ vortex pairs, two-vison attractive potential $V(R)$, sub-gap Majorana bound modes ($\epsilon_0 / \Delta_{\mathrm{bulk}} \approx 13\%$), and self-consistent flux defect trapping (`K-020`).
- **Self-Consistent Purification, Disorder Conditioning & Post-Spinodal Audit (M1.11):** Outer DIIS and damped L-Broyden ($\alpha=0.5$) real-space purification solver (`K-021`), bond-disorder and vacancy zero-mode conditioning ($\kappa(\mathbf{A}) \approx 3.82 \times 10^4$, `K-022`), post-spinodal KKT saddle excursion audit and coexistence domain-wall tension $\sigma_{\mathrm{DW}} = +0.042|K|$ (`K-023`), and the **Phase 1 Master Solver Selection & Hardware Resource Table** (`HC-P1-01`--`03`).
- **Centralized Hypothesis Contracts Registry (`doc/HYPOTHESIS_CONTRACTS.md` & Chapter 11):** Consolidated all 22 foundational (`PH-01`--`06`, `AH-01`--`05`), closure (`HC-P1-01`--`03`, `C-M1.9-T1`--`T6`), and Phase 2 (`H1`--`H11` / `HC-P2-01`--`05`) hypothesis contracts with pass/refuted verdicts and mathematical post-mortems.
- **Living Lab Book & SciPost Executive Monograph:** Authored and compiled the 74-page R&D Living Lab Book (`doc/labbook/labbook_kitaev_mf.pdf`, 11 chapters, beige parchment theme) and the 9-page SciPost Executive Monograph (`doc/manuscript/scipost_kitaev_mf.pdf`) with automated GitHub Actions release publishing.

### Changed
- Refactored project governance to strict multi-phase architecture: Phase 1 (Numerical Quantification & Algorithmic Foundations of FM Kitaev $[111]$) $\to$ Phase 2 (Ground State, Isolated Vortex Pairs, Finite Temperature & Dynamical Correlations of FM Kitaev $[111]$) $\to$ Phase 3+ (Subsequent Model Families).
