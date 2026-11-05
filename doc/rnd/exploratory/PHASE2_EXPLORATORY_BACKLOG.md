---
id: RND-EXP-BACKLOG-001
title: "Phase 2 Exploratory Backlog: Unresolved Foundational Questions, Projection Gaps, and Advanced Numerical Manifolds"
phase: 2
status: planned
backlog_items:
  - id: BL-01
    title: "Gutzwiller Projection Gap and Exponential Orthogonality Catastrophe"
    first_referenced: "ch_m1_1b, ch_m1_4"
    status: exploratory
  - id: BL-02
    title: "Factored vs. Full Multi-Link Plaquette Wilson Loop Operator"
    first_referenced: "ch_m1_1b, ch_m1_4"
    status: exploratory
  - id: BL-03
    title: "Incommensurate Spiral Magnetic Ordering and Artificial Supercell Pitch Locking"
    first_referenced: "ch_m1_2, ch_m1_4_5"
    status: exploratory
  - id: BL-04
    title: "Coordinate System Shear and Passive Geometry Transformation"
    first_referenced: "ch_m1_2"
    status: exploratory
  - id: BL-05
    title: "Infrared Divergence of Static Gauge Susceptibility at Gapless Dirac Nodes"
    first_referenced: "ch_m1_3"
    status: exploratory
  - id: BL-06
    title: "Semidefinite Programming Lower Bounds via Variational 2-RDM Hierarchy"
    first_referenced: "ch_m1_1b, ch_m1_4_5"
    status: exploratory
  - id: BL-07
    title: "Non-Gaussianity Quantum Relative Entropy Metric as an A Posteriori Error Bound"
    first_referenced: "ch_m1_1b"
    status: exploratory
---

# Phase 2 Exploratory Backlog: Foundational Questions & Future Roadmaps

This document serves as the formal epistemological repository for deep physical and mathematical questions identified during the Phase 1 multi-agent adversarial audit of `kitaev-mf`. In accordance with workspace scientific governance principles, complex phenomena without immediate, verified, closed-form solutions are isolated here rather than obscured by ad-hoc approximations.

---

## BL-01: Gutzwiller Projection Gap and Exponential Orthogonality Catastrophe

- **Location & First Appearance:** Lab Book Chapters `ch_m1_1b` (§7) and `ch_m1_4` (§2); Note `m1_1_t3`.
- **Physical Context & Mathematical Statement:**
  Mean-field states are defined as pure quasi-free states $|\Psi_0\rangle$ on the enlarged $4N$-Majorana Fock space $\mathcal{F}$. Physical spin observables are evaluated as unprojected Gaussian expectations:
  $$\langle \mathcal{O} \rangle_{\mathrm{G}} \equiv \langle \Psi_0 | \mathcal{O} | \Psi_0 \rangle.$$
  However, the physical spin state is the Gutzwiller-projected state:
  $$|\Psi_{\mathrm{phys}}\rangle = \frac{1}{\sqrt{\mathcal{Z}}} P_G |\Psi_0\rangle, \quad P_G = \prod_{i=1}^N \frac{\Id + \mathcal{D}_i}{2}, \quad \mathcal{D}_i = -c_i^0 c_i^1 c_i^2 c_i^3.$$
  Because each local projector has eigenvalue $1$ on a 2D subspace of the local 4D Fock space, the unprojected state norm satisfies the exponential decay bound:
  $$\mathcal{Z} = \langle \Psi_0 | P_G | \Psi_0 \rangle \sim 2^{-N}.$$
- **Why It Cannot Be Solved Simultaneously in Phase 1:**
  Exact Gutzwiller projection destroys the quasi-free property of the state. Evaluating $\langle \Psi_{\mathrm{phys}} | \mathcal{O} | \Psi_{\mathrm{phys}} \rangle$ analytically in the thermodynamic limit requires Pfaffian determinant Monte Carlo (VMC) or non-Abelian tensor network contractions (PEPS), which lie outside quadratic mean-field theory.
- **Minimal Experiment / Attack Vector:**
  Implement a finite-cluster Pfaffian evaluation tool on small honeycomb clusters ($N = 6, 8, 12, 16$) comparing exact diagonalization (ED), unprojected mean-field expectations $\langle \sigma_i^\alpha \sigma_j^\beta \rangle_{\mathrm{G}}$, and Gutzwiller-projected expectations $\langle \sigma_i^\alpha \sigma_j^\beta \rangle_{\mathrm{phys}}$.
- **Acceptance Signal:**
  Quantitative calibration curve of the projection gap $\Delta_{\mathrm{proj}} \equiv \|\langle \boldsymbol{\sigma} \rangle_{\mathrm{phys}} - \langle \boldsymbol{\sigma} \rangle_{\mathrm{G}}\|$ as a function of the Lagrange multipliers $\boldsymbol{\lambda}$.

---

## BL-02: Factored vs. Full Multi-Link Plaquette Wilson Loop Operator

- **Location & First Appearance:** Lab Book Chapters `ch_m1_1b` (§6) and `ch_m1_4` (§3); Note `m1_4_t1`.
- **Physical Context & Mathematical Statement:**
  The hexagonal plaquette Wilson loop operator on plaquette $p$ is defined as the product of directed bond gauge operators around boundary cycle $\partial p$:
  $$\hat{W}_p = \prod_{b \in \partial p} \hat{u}_b = \prod_{b=(i,j) \in \partial p} (-i c_i^{\gamma_b} c_j^{\gamma_b}).$$
  In Phase 1 production services (`observables.py`), the observable is evaluated under the **factored approximation**:
  $$W_p^{\mathrm{factored}} \equiv \prod_{b \in \partial p} \langle \hat{u}_b \rangle_{\mathrm{G}} = \prod_{b \in \partial p} \mathrm{U}_b^{\gamma_b \gamma_b}.$$
  By Wick's theorem, the true Gaussian expectation $\langle \hat{W}_p \rangle_{\mathrm{G}}$ is the Pfaffian of a $12 \times 12$ skew-symmetric covariance sub-block, which contains connected 4-link and 6-link contractions:
  $$\langle \hat{W}_p \rangle_{\mathrm{G}} = \operatorname{Pf}\left(\Omega_{\partial p}\right) = W_p^{\mathrm{factored}} + \sum_{\text{cross-contractions}} \mathrm{U}_{b_1} \mathrm{U}_{b_2} \cdots.$$
- **Why It Cannot Be Solved Simultaneously in Phase 1:**
  The pure Kitaev zero-flux sector has exact $\langle \hat{u}_b \rangle = \pm 1$, causing cross-contractions to vanish identically in that solvable limit. However, in the presence of $\Gamma, \Gamma', h$, off-diagonal cross-contractions become active, and their quantitative magnitude across phase boundaries has not been mapped.
- **Minimal Experiment / Attack Vector:**
  Add a service function `compute_exact_plaquette_pfaffian(state, p)` in `services/observables.py` and benchmark the discrepancy $|W_p^{\mathrm{full}} - W_p^{\mathrm{factored}}|$ across a sweep of $h/|K|$ in $[0, 1.0]$.
- **Acceptance Signal:**
  Publication figure and `.dat` benchmark establishing the parameter domain where the factored approximation remains within $10^{-4}$ of the true operator expectation.

---

## BL-03: Incommensurate Spiral Magnetic Ordering and Supercell Pitch Locking

- **Location & First Appearance:** Lab Book Chapters `ch_m1_2` (§4) and `ch_m1_4_5` (§3); Note `multisublattice_saddle_survey.md`.
- **Physical Context & Mathematical Statement:**
  Under extended Heisenberg $J$ and off-diagonal exchange $\Gamma, \Gamma'$, classical Monte Carlo and infinite DMRG find incommensurate spiral and non-coplanar phases whose magnetic pitch vector $\mathbf{Q}_{\mathrm{spiral}}$ varies continuously with exchange parameters. The commensurate 2-sublattice primitive cell and 4-sublattice orthorhombic supercell artificially constrain allowed momentum vectors to discrete reciprocal fractions ($0, \pi, \pi/2$).
- **Why It Cannot Be Solved Simultaneously in Phase 1:**
  Resolving an incommensurate pitch $\mathbf{Q} \notin \frac{1}{L}\mathbb{Z}^2$ requires either generalized spiral boundary conditions, continuous Fourier variational ansätze, or large-scale supercells ($L \times L \times 2$ with $L \ge 24$).
- **Minimal Experiment / Attack Vector:**
  Develop a continuous pitch vector optimization module $\mathbf{Q} \in \mathrm{BZ}$ parameterizing twisted spatial boundary conditions:
  $$c_{i + \mathbf{R}}^\mu = e^{i \mathbf{Q} \cdot \mathbf{R}} c_i^\mu.$$
- **Acceptance Signal:**
  Continuous energy minimization over $\mathbf{Q}$ reproducing the known classical incommensurate spiral wavevector for the $J\text{-}K\text{-}\Gamma$ model.

---

## BL-04: Coordinate System Shear and Passive Geometry Transformation

- **Location & First Appearance:** Lab Book Chapter `ch_m1_2` (§2); Note `multisublattice_saddle_survey.md`.
- **Physical Context & Mathematical Statement:**
  The exploratory note `multisublattice_saddle_survey.md` adopts a coordinate system where zigzag chains are oriented along the Cartesian $x$-axis with primitive translation vectors:
  $$\mathbf{a}_1 = (a, 0), \quad \mathbf{a}_2 = (a/2, \sqrt{3}a/2).$$
  Conversely, `src/kitaev_mf/domain/lattice.py` and standard Kitaev literature use coordinates rotated by $30^\circ$ or oriented along the armchair direction.
- **Why It Cannot Be Solved Simultaneously in Phase 1:**
  Changing coordinate bases in production code would invalidate existing unit test assertions (`test_lattice.py`, `test_hamiltonian.py`).
- **Minimal Experiment / Attack Vector:**
  Formulate and document the explicit $\mathrm{SO}(2)$ coordinate rotation matrix and corresponding bond permutation map $\Pi \in S_3$ mapping the saddle survey frame to `lattice.py`.
- **Acceptance Signal:**
  Automated translation test verifying that dispersion relations and high-symmetry paths match between frames within machine precision.

---

## BL-05: Infrared Divergence of Static Gauge Susceptibility at Gapless Dirac Nodes

- **Location & First Appearance:** Lab Book Chapter `ch_m1_3` (§4); Note `car_gauge_formalism_survey.md`.
- **Physical Context & Mathematical Statement:**
  At the gapless Kitaev point, conical Dirac crossings occur at $K$ and $K'$. The static gauge susceptibility tensor:
  $$\chi^{\mathrm{gauge}}_{b b'}(\omega=0) \sim \int_{\mathrm{BZ}} \frac{\dd^2\mathbf{k}}{(2\pi)^2} \sum_{m,n} \frac{|\langle \psi_n(\mathbf{k}) | \hat{u}_b | \psi_m(\mathbf{k}) \rangle|^2}{\mathcal{E}_n(\mathbf{k}) + \mathcal{E}_m(\mathbf{k})},$$
  exhibits infrared sensitivity due to the vanishing excitation gap $\Delta \to 0$. In self-consistent DIIS and Picard loops, this induces slow algebraic convergence $O(1/\sqrt{\ell})$ or oscillations along gauge directions.
- **Why It Cannot Be Solved Simultaneously in Phase 1:**
  Addressing this requires formulating an infrared-regularized / dual-preconditioned Grassmannian Newton solver, which is the primary subject of Milestone M1.6.
- **Minimal Experiment / Attack Vector:**
  Implement a spectral preconditioning operator $P_{\mathrm{reg}}(\mathbf{k}) = (\mathcal{H}(\mathbf{k})^2 + \epsilon_{\mathrm{IR}}^2 \Id)^{-1/2}$ inside the Jacobian-Free Newton-Krylov (JFNK) inner loop in M1.6.
- **Acceptance Signal:**
  Convergence of the mean-field residual $\|R\| < 10^{-10}$ in fewer than 20 iterations at the gapless point.

---

## BL-06: Semidefinite Programming Lower Bounds via Variational 2-RDM Hierarchy

- **Location & First Appearance:** Lab Book Chapter `ch_m1_1b` (§8) and `ch_m1_4_5` (§2).
- **Physical Context & Mathematical Statement:**
  Mean-field theory on the Grassmannian provides a strict variational **upper bound** to the true ground-state energy:
  $$E_0 \le \min_{M \in \Gr} E(M).$$
  Conversely, the fermionic $N$-representability problem yields a rigorous sequence of convex **lower bounds** via the semidefinite programming (SDP) moment-matrix hierarchy on 2-body reduced density matrices (2-RDMs):
  $$E_{\mathrm{SDP}}^{(2)} = \min_{\rho^{(2)} \ge 0} \Tr[H \rho^{(2)}] \quad \text{s.t.} \quad D \ge 0, \quad Q \ge 0, \quad G \ge 0.$$
  Sandwiching $E_0$ between the SDP lower bound and the quasi-free upper bound yields a mathematically unassailable error certificate:
  $$E_{\mathrm{SDP}}^{(2)} \le E_0 \le E_{\mathrm{MF}}.$$
- **Why It Cannot Be Solved Simultaneously in Phase 1:**
  Requires integrating a high-performance interior-point SDP solver (such as MOSEK, SCS, or Clarabel via `cvxpy`) for fermionic 2-RDM constraints on honeycomb clusters.
- **Minimal Experiment / Attack Vector:**
  Formulate the 2-site and 4-site reduced Hamiltonian SDP relaxation in Python and compute the bounds for the pure Kitaev and extended Kitaev models.
- **Acceptance Signal:**
  Automated certificate script confirming $E_{\mathrm{SDP}}^{(2)} \le E_{\mathrm{exact}} \le E_{\mathrm{MF}}$ across the phase diagram.

---

## BL-07: Non-Gaussianity Quantum Relative Entropy Metric as an A Posteriori Error Bound

- **Location & First Appearance:** Lab Book Chapter `ch_m1_1b` (§9).
- **Physical Context & Mathematical Statement:**
  For any state $\rho$, its non-Gaussianity can be rigorously quantified by the quantum relative entropy to its Gaussian surrogate $\omega_{\mathrm{G}}[\rho]$ having the identical 1-RDM:
  $$D_{\mathrm{NG}}(\rho) \equiv S(\rho \| \omega_{\mathrm{G}}[\rho]) = S(\omega_{\mathrm{G}}[\rho]) - S(\rho) \ge 0.$$
  By Klein's inequality and Pinsker's theorem, this relative entropy upper-bounds the trace distance to the quasi-free manifold:
  $$\|\rho - \omega_{\mathrm{G}}\|_1 \le \sqrt{2 \ln 2 \cdot D_{\mathrm{NG}}(\rho)}.$$
- **Why It Cannot Be Solved Simultaneously in Phase 1:**
  Evaluating $D_{\mathrm{NG}}(\rho)$ requires access to the exact many-body ground-state density matrix $\rho$ from DMRG or exact diagonalization.
- **Minimal Experiment / Attack Vector:**
  Extract ground states of small clusters from QuSpin or TenPy and compute $D_{\mathrm{NG}}(\rho)$ across the transition from the QSL to the polarized phase.
- **Acceptance Signal:**
  Correlate $D_{\mathrm{NG}}(\rho)$ with mean-field observable errors, proving that mean-field accuracy is governed by non-Gaussianity.
