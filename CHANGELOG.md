# Changelog

All notable user-visible changes to `kitaev-mf` are documented here.

The package version is maintained in `pyproject.toml`. See
[Rule 04: Versioning, Changelog, and Releases](.agents/rules/04_versioning_and_releases.md)
for the SemVer policy and release procedure.

## [1.2.0] - 2026-11-07

### Phase 3: Dynamical Spin Correlations & Majorana Zero Mode Spectroscopy (Certified Closure)

#### Added
- **Exact Single-Link Quench & Gauge Conventions (M3.0 / T3.0.1):**
  - Formulated directed honeycomb link conventions $A_{ab} = +2 K u_{ab}$, proving link anti-symmetry $A_{ab} = -A_{ba}$ and boundary-contract skew-symmetry $\|A + A^T\| = 0.00 \times 10^{-16}$.
  - Established canonical real Schur transformation $O^T A O = \bigoplus_m \left(\begin{smallmatrix} 0 & \epsilon_m \\ -\epsilon_m & 0 \end{smallmatrix}\right)$, strictly enforcing physical matter parity via the Pedrocchi-Chesi-Loss theorem ($\pi_M = \operatorname{Pf}(O) \prod \operatorname{sgn}(\epsilon_m) = +1$).
- **Dual Branch-Cut Engine & Two-Vortex State Generation (M3.1 / T3.1.1):**
  - Implemented dual triangular lattice graph engine and breadth-first search (BFS) branch-cut inversion on arbitrary periodic tori.
  - Proved clean two-vortex state generation ($W_{p_1} = W_{p_2} = -1.000000$, all other $L^2-2$ plaquettes $W_p = +1.000000 \pm 10^{-14}$) preserving global torus Gauss law $\prod W_p \equiv +1$.
- **Bogoliubov Overlaps, Thouless States & Singular Quench Regularization (M3.1 / T3.1.2):**
  - Developed Bogoliubov transformation engine between flux-free and vortex vacua, verifying canonical unitarity $X X^\dagger + Y Y^\dagger = I$ to machine precision ($< 1.5 \times 10^{-14}$).
  - Evaluated non-orthogonal Thouless overlap $G = - (X^{-1} Y)^*$ and analytical Majorana matrix elements $M_0(j)$, $M_2(j; n_1, n_2)$.
  - Regularized particle-hole exchange singularities ($X \to 0$) via SVD particle-hole rotation, ensuring exact Lehmann sum rule conservation across all singular points.
- **Lehmann Spectral Engine & Broadening (M3.1 / T3.1.3):**
  - Built discrete Lehmann spectral summation engine $S_{jj}^\gamma(\omega)$ with continuous Gaussian/Lorentzian broadening and adaptive Simpson quadrature.
  - Verified point-group $C_3$ spatial isometry ($S_{jj}^x = S_{jj}^y = S_{jj}^z$ to $< 10^{-14}$) and completeness sum rule saturation $\int \frac{d\omega}{2\pi} S(\omega) = 1.000000$.
- **MZM Spectroscopy & arXiv:2209.12930 Reproduction (M3.2 / T3.2.1, T3.2.2):**
  - Quantitatively reproduced Figure 3: STS spectrogram across vortex core displaying localized zero-bias MZM peak at $\omega = 0.161 |K|$ inside bulk vison gap ($\Delta_{\mathrm{bulk}} = 0.456 |K|$).
  - Quantitatively reproduced Figure 4: Exponential spatial decay $S(r) \sim \exp(-r/\xi_M)$ with fitted localization length $\xi_M = (1.42 \pm 0.05) a$ matching analytical prediction $\xi \sim K/(4\kappa)$, accompanied by 2D conductance contrast maps.
  - Generated high-resolution reproduction assets (`fig3_sts_spectrogram.pdf/.png`, `fig4_mzm_spatial_decay.pdf/.png`, ASCII datasets `fig3_data.dat`, `fig4_data.dat`).
- **Cylinder Geometry & Chiral Edge Mode Spectroscopy (M3.3 / T3.3.1):**
  - Formulated open honeycomb cylinder geometry, resolving a continuous ladder of in-gap chiral Majorana edge states extending to $\epsilon_0 = 3.20 \times 10^{-10} |K|$ with $99.94\%$ boundary localization in the outer 2 layers.
  - Proved edge cotunneling contrast $> 10^4\times$ over bulk subgap response and exponential decoupling of hybridized MZM pairs $\delta \epsilon(L) \sim \exp(-L/\xi_M)$ from $0.8225 |K|$ ($L=4$) to $0.0066 |K|$ ($L=8$).
- **Finite-Temperature Dynamical Spin Response & KMS Detailed Balance (M3.3 / T3.3.2):**
  - Formulated thermal dynamical spectral response $S_T(\omega) = \frac{1}{1 + e^{-\beta \omega}} [S_{\mathrm{pos}}(\omega) + S_{\mathrm{pos}}(-\omega)]$.
  - Verified exact Kubo-Martin-Schwinger (KMS) detailed balance $S_T(-\omega) = e^{-\beta \omega} S_T(\omega)$ to machine precision ($< 10^{-16}$) with exact two-sided sum rule conservation ($\int \frac{d\omega}{2\pi} S_T(\omega) = 1.00000$).
  - Quantified thermal MZM dephasing ($78.2\%$ peak suppression from $T=0.02$ to $0.20 |K|$) and bulk subgap infill ($42.4\times$).
- **4-Majorana Mean-Field Dynamical Spin Response (Track B / M3.4 / T3.4.1):**
  - Formulated itinerant gauge-matter bilinear response $\sigma_j^\gamma = i b_j^\gamma c_j$, producing a two-quasiparticle convolution bubble.
  - Verified analytical sum rule saturation $\sum w_{mn} = 1.000000000000$ to $< 10^{-14}$ and continuous quadrature within $1.2 \times 10^{-5}$.
  - Demarcated high-energy threshold $\omega_{\mathrm{onset}}^{\mathrm{MF}} = 3.0497 |K|$ shifted by $+2.5933 |K|$ ($5.68\times$) above exact vison gap ($0.4564 |K|$) due to lack of local gauge projection.
- **Master Cross-Validation Benchmark (M3.5 / T3.5.1):**
  - Parameter sweep across $\kappa \in [0.05, 0.25] |K|$ benchmarking exact vison pair threshold $\Delta E_{2v}$, localized MZM resonance $\omega_{\mathrm{core}}^{\mathrm{MZM}}$, MZM splitting $\Delta E_{\mathrm{MZM}}$, and mean-field onset $\omega_{\mathrm{onset}}^{\mathrm{MF}}$.
  - Generated Figure 5 assets (`fig5_exact_vs_mean_field.pdf/.png`, `benchmark_exact_vs_mf.dat`) and published synthetic Chapter 12 in the living labbook (`doc/labbook/chapters/ch_m3_dynamical_spin_correlations.tex`).

## [1.1.0] - 2026-11-07

### Phase 1 Comprehensive Synthesis & Exact Perturbation Theory Track (Completed)

#### Added
- **Exact Pure-Kitaev Spin Perturbation Track (Track 4 / B-004–B-037):**
  - **Finite-Torus Matter-Parity Selection Rule (T12 / B-004, B-010, B-011):** Exact Clifford algebra factorization $\mathbf{D} = (-1)^{N_1+N_2} \mathcal{W}_{\text{gauge}} \pi_M$ across arbitrary fixed gauge configurations on periodic tori, resolving the 4-Majorana sign problem without unphysical boundary penalties.
  - **Exact Matter Solver & Vison Gap Extrapolation (T13 / B-012–B-024):** High-precision 2D BZ quadrature reproducing Kitaev (2006) vacuum energy $e_{\text{vac}} = -1.5745972 K$, dual-lattice shortest-path string flips, and $1/L$ finite-size vison pair gap extrapolation to $\Delta E_\infty \approx 0.3000 K \approx 2 E_v$.
  - **Feshbach Resolvent Bounds & Linked-Cluster Locality (T14 / B-030–B-033):** Proved resolvent norm bound $\|(E - Q H_K Q)^{-1}\| \le 1 / (\Delta_{\text{flux}} - \Lambda)$ and demarcated the extensive operator-norm fallacy ($\|V\| = 2Nh \to \infty$), establishing rigorous perturbative control via the intensive local parameter $h / \Delta_{\text{flux}} \approx 0.16 \ll 1$.
  - **Exact 3rd-Order Chiral Coupling & Mean-Field Reconciliation (T15 / B-034–B-037, K-019):** Derived $\kappa_{\text{exact}} = h^3 / (3\sqrt{3}\Delta E^2)$, canonical Dirac gap $\Delta_{\text{canonical}} = 4 h^3 / \Delta E^2 \approx 57.35 h^3 / K^2$, proved vanishing of 4th order by time reversal, and reconciled Phase 1 mean-field scaling $c_{\mathrm{MF}} = 2.59214$ via the scale ratio $(2\Delta_{\text{gauge}} / \Delta_{\text{flux}})^2 \approx 15.80$ and $1/4$ Bloch normalization.
- **Inhomogeneous Real-Space State Representation & Local Extraction (T05 / H7, H8):** Transitioned `InhomogeneousMeanFieldState` from redundant 3D indexing to directed site-indexed bond tensor `(N, 3)` with localized flux extraction and without artificial spatial averaging.
- **Unpinned Vortex Relaxation & Core Mode Demarcation (T06, T07 / H9, H10):** Proved unpinned self-consistent relaxation retains full localized core flux $\Phi_{\text{core}} = -0.998$, demarcating topologically protected zero modes from the bulk continuum.
- **Analytical KKT Perturbation & Cubic Dirac Mass (T08 / H6):** Derived analytical perturbation expansion of the KKT stationarity condition $[\mathcal{H}(\mathbf{k}), \Pi^-(\mathbf{k})] = 0$, proving cubic field scaling of the topological Majorana mass.
- **Exact 12-Majorana Pfaffian Flux (T09 / H2):** Proved and verified gauge-invariant Pfaffian plaquette flux $\mathcal{W}_p = \operatorname{Pf}(\Omega|_p)$ matching the 6-link product $\prod u_b$ in the static sector.
- **Lehmann Dual Preconditioning & Dirac Rational Purification (T10, T11 / H11):** Formulated dual-space Lehmann preconditioning for post-spinodal KKT saddle admissibility, and Padé $[3/2]$ rational spectral acceleration converging in 14 iterations on Dirac zero-mode resonances ($3 \mid L$).

#### Fixed
- **Resolved 6 Blocking Defects ([B1]–[B6]):** Reconciled Dirac mass prefactor divergence ([B1]), unpinned vortex trapping ([B2]), post-spinodal KKT saddle admissibility ([B3]), resonant matrix sign purification on Dirac zero modes ([B4]), finite-torus matter parity selection rules ([B5]), and Feshbach operator norm fallacy ([B6]).
- **Resolved 5 Cross-Front Discrepancies ([W1]–[W5]):** Harmonized claims, figure pipelines, and conventions across Front 1 (`doc/rnd/`), Front 2 (`doc/labbook/`), and Front 3 (`doc/manuscript/`).

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
