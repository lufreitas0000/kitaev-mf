# Changelog

All notable user-visible changes to `kitaev-mf` are documented here.

The package version is maintained in `pyproject.toml`. See
[Rule 04: Versioning, Changelog, and Releases](.agents/rules/04_versioning_and_releases.md)
for the SemVer policy and release procedure.

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
- **Real-Space Torus Clusters & Pure Real Arithmetic (M1.8 Track A):** Exact real skew-symmetric algebra $\Omega = -\mathbf{A}(-\mathbf{A}^2)^{-1/2}$, cutting memory by $2.0\times$ (8 bytes vs 16 bytes) and FLOPs by $4.0\times$.
- **Eigenvector-Free Matrix Sign Purification (M1.8 Track B):** Cubic Newton-Schulz matrix sign purification $\Omega = \ii \operatorname{sgn}(\ii\mathbf{A})$ in pure BLAS-3 `dgemm`, delivering a **$63.53\times$ speedup** over divide-and-conquer eigensolvers at $M=1152$ ($702\text{ ms}$ vs $44.6\text{ s}$) with machine-precision parity ($4.22 \times 10^{-15}$).
- **Real-Space Vortex Geometry & Bound States (M1.9):** Dual-triangular lattice branch-cut generator for isolated $\mathbb{Z}_2$ vortex pairs, two-vison attractive potential $V(R)$, sub-gap Majorana bound modes ($\epsilon_0 / \Delta_{\mathrm{bulk}} \approx 13\%$), and self-consistent flux defect trapping.
- **Scientific Hypothesis Contracts (Chapter 11):** Formulated 11 falsifiable scientific and algorithmic hypothesis contracts (`H1`--`H11`) governing Gutzwiller projection gaps, 12-Majorana Pfaffians, 2-RDM SDP bounds, and real-space vortex polarization clouds.
- **Living Lab Book Monograph:** Authored and compiled 72-page comprehensive R&D monograph `doc/labbook/labbook_kitaev_mf.pdf` (11 chapters, beige parchment theme, 0 warnings/errors).

### Changed
- Refactored project governance to strict multi-phase architecture: Phase 1 (Numerical Quantification & Algorithmic Foundations) $\to$ Phase 2 (Ground State, Vortices, Thermodynamics & Dynamics of FM Kitaev $[111]$).
