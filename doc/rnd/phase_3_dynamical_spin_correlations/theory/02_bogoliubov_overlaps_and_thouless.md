---
title: "Task T3.1.2: Bogoliubov Overlap Transformation, Thouless State Representation & Single-Particle Matrix Elements"
status: completed
author: "Lu Freitas"
date: "2026-11-07"
phase: "Phase 3"
milestone: "M3.1"
task: "T3.1.2"
contracts: ["H_BOG_01", "H_BOG_02", "H_BOG_03", "H_BOG_04"]
primary_model: "Pure FM Kitaev (K > 0) with chiral 3-spin term kappa (arXiv:2209.12930)"
geometry: "Real-space honeycomb lattice clusters (torus and cylinder geometries)"
guiding_questions:
  - "How are the Bogoliubov transformation blocks X and Y constructed from the canonical real Schur eigensystems of H_G and H_G'?"
  - "How are the three canonical unitarity and anti-symmetry relations verified to machine precision?"
  - "How does the Thouless state representation express the unperturbed vacuum |0_a> in the perturbed b-mode basis?"
  - "How are singular Bogoliubov transformations (det(X) = 0, det(W) = -1) resolved via the particle-hole interchange of Eq. (29) in arXiv:2209.12930?"
  - "How are the transition matrix elements <lambda(G') | c_j | phi_0(G)> analytically evaluated and tested against Fock space exact diagonalization and the unitary sum rule sum |M|^2 = 1?"
---

# Task T3.1.2: Bogoliubov Overlap Transformation, Thouless State Representation & Single-Particle Matrix Elements

## 1. Scope & Mathematical Objectives
In **arXiv:2209.12930** (Phys. Rev. B **107**, 054432), calculating the local dynamical spin correlation function:
$$S_{jl}^{\alpha\beta}(\omega) = \sum_\lambda \langle \phi_0(G) | \sigma_j^\alpha | \lambda(G_j^\alpha) \rangle \langle \lambda(G_j^\alpha) | \sigma_l^\beta | \phi_0(G) \rangle \delta(\omega + E_0(G) - E_\lambda(G_j^\alpha))$$
requires evaluating the many-body matrix elements of itinerant Majorana operators $c_j$ between the physical ground state $|\phi_0(G)\rangle$ in gauge sector $G$ and the excited eigenstates $|\lambda(G')\rangle$ in the spin-flipped gauge sector $G' = G_j^\alpha$.

Because $H_G$ and $H_{G'}$ are quadratic in Majoranas but diagonal in different quasiparticle bases ($a_k$ vs. $b_m$), this is an exact quantum quench problem governed by the Bogoliubov transformation matrix $W$.

The mathematical objectives of this atomic task are:
1. Construct the complex mode wavefunctions $\phi_{a, km}$ and $\phi_{b, mj}$ from the real Schur decompositions $A_G = Q_a \Sigma_a Q_a^T$ and $A_{G'} = Q_b \Sigma_b Q_b^T$.
2. Formulate the Bogoliubov transformation blocks $X$ and $Y$ and prove the canonical relations:
   $$X X^\dagger + Y Y^\dagger = I, \qquad X^\dagger X + Y^T Y^* = I, \qquad X Y^T + Y X^T = 0.$$
3. Implement the Thouless state representation expressing $|0_a\rangle$ in the $b$-mode basis with skew-symmetric pairing matrix $Z_b = -(X^\dagger)^{-1} Y^T$.
4. Implement the particle-hole interchange of Eq. (29) in arXiv:2209.12930 for singular transformations ($\det(X) \approx 0, \det W = -1$).
5. Derive and implement exact analytical matrix elements for the transition amplitudes $\langle \lambda(G') | c_j | \phi_0(G) \rangle$ for both even-parity ($\pi_G = +1$) and odd-parity ($\pi_G = -1$, one-fermion occupied) ground states, validating them against full Fock-space exact diagonalization and the completeness sum rule:
   $$\sum_\lambda |\langle \lambda | c_j | \phi_0 \rangle|^2 = 1.0.$$

---

## 2. Theoretical Architecture

### 2.1 Mode Wavefunctions from Real Schur Eigensystems
For a cluster of $2N$ sites, canonical real Schur decomposition yields $A = Q \Sigma Q^T$ with standardized $2 \times 2$ blocks $\begin{pmatrix} 0 & \epsilon_k \\ -\epsilon_k & 0 \end{pmatrix}$ ($\epsilon_k \ge 0$).
The $N$ complex Dirac annihilation operators are:
$$a_k = \frac{1}{2}(c'_{2k} - \mathrm{i} c'_{2k+1}) = \sum_{m=0}^{2N-1} \phi_{k, m} c_m$$
where:
$$\phi_{k, m} = \frac{1}{2}(Q_{m, 2k} - \mathrm{i} Q_{m, 2k+1}) \in \mathbb{C}^{N \times 2N}.$$
Canonical Clifford CAR $\{c_m, c_n\} = 2\delta_{mn}$ strictly implies:
$$\{a_k, a_l^\dagger\} = 2 (\phi \phi^\dagger)_{kl} = \delta_{kl}, \qquad \{a_k, a_l\} = 2 (\phi \phi^T)_{kl} = 0.$$
Inverting expresses real Majoranas as:
$$c_m = 2 \sum_k [ \phi^*_{k, m} a_k + \phi_{k, m} a_k^\dagger ].$$

### 2.2 Bogoliubov Transformation Matrix
Let $a_k$ denote the eigenmodes of $H_G$ ($\phi_a$) and $b_m$ denote the eigenmodes of $H_{G'}$ ($\phi_b$).
Expressing $b_m = \sum_j \phi_{b, mj} c_j$ in terms of $a$ and $a^\dagger$ yields:
$$b_m = \sum_k (X_{mk} a_k + Y_{mk} a_k^\dagger)$$
where:
$$X = 2 \phi_b \phi_a^\dagger, \qquad Y = 2 \phi_b \phi_a^T.$$
In block Nambu form:
$$\begin{pmatrix} b \\ b^\dagger \end{pmatrix} = \begin{pmatrix} X & Y \\ Y^* & X^* \end{pmatrix} \begin{pmatrix} a \\ a^\dagger \end{pmatrix} = W^\dagger \begin{pmatrix} a \\ a^\dagger \end{pmatrix}.$$
The inverse transformation is:
$$\begin{pmatrix} a \\ a^\dagger \end{pmatrix} = \begin{pmatrix} X^\dagger & Y^T \\ Y^\dagger & X^T \end{pmatrix} \begin{pmatrix} b \\ b^\dagger \end{pmatrix}.$$

### 2.3 Regular Thouless Representation ($\det(X) \neq 0$)
The vacuum $|0_a\rangle$ ($a_k |0_a\rangle = 0$ $\forall k$) satisfies:
$$\sum_m [ X^\dagger_{km} b_m + Y^T_{km} b_m^\dagger ] |0_a\rangle = 0 \implies \left( b_k + \sum_m [(X^\dagger)^{-1} Y^T]_{km} b_m^\dagger \right) |0_a\rangle = 0.$$
Per Thouless' theorem, $|0_a\rangle$ is an unnormalized paired condensate over $|0_b\rangle$:
$$|0_a\rangle = \langle 0_b | 0_a \rangle \exp\left( \frac{1}{2} \sum_{mn} Z_{b, mn} b_m^\dagger b_n^\dagger \right) |0_b\rangle$$
where the pairing matrix:
$$Z_b = - (X^\dagger)^{-1} Y^T \in \mathbb{C}^{N \times N}$$
is strictly skew-symmetric ($Z_b^T = -Z_b$), and the vacuum overlap is:
$$|\langle 0_b | 0_a \rangle|^2 = \sqrt{\det(X^\dagger X)} = |\det(X)|.$$

### 2.4 Singular Transformations & Particle-Hole Interchange ($\det(X) \approx 0$)
When an isolated Majorana zero mode switches parity between sectors, $\det W = -1$, causing $\det(X) = 0$ and $\langle 0_b | 0_a \rangle = 0$.
Following Eq. (29) of arXiv:2209.12930, one applies a particle-hole transformation on the zero mode $\mu$:
$$\tilde{b}_\mu = b_\mu^\dagger, \qquad \tilde{b}_\mu^\dagger = b_\mu.$$
This interchanges the $\mu$-th rows of $X$ and $Y$:
$$X'_{\mu, k} = Y^*_{\mu, k}, \qquad Y'_{\mu, k} = X^*_{\mu, k}.$$
The transformed block has $\det(X') \neq 0$, and the overlap of $|0_a\rangle$ with the single-quasiparticle state $b_\mu^\dagger |0_b\rangle$ evaluates as:
$$|\langle 0_b | b_\mu | 0_a \rangle|^2 = |\det(X')|.$$

### 2.5 Analytical Transition Matrix Elements of $c_j$
For the one-fermion physical ground state $|\phi_0(G)\rangle = a_p^\dagger |0_a\rangle$ (with $p=0$ on tori):
1. **Vacuum transition (0 quasiparticles in $b$-basis):**
   $$\langle 0_b | c_j | a_p^\dagger 0_a \rangle = 2 \langle 0_b | 0_a \rangle \sum_m \phi^*_{b, mj} [(X^\dagger)^{-1}]_{mp}.$$
2. **Two-quasiparticle transitions ($m_1 < m_2$):**
   $\langle m_1 m_2 | = \langle 0_b | b_{m_2} b_{m_1} = - \langle 0_b | b_{m_1} b_{m_2}$.
   $$\langle m_1 m_2 | c_j | a_p^\dagger 0_a \rangle = - (t_1 + t_2 + t_3)$$
   where:
   $$t_1 = 2 \phi_{b, m_2 j} \langle 0_b | 0_a \rangle [(X^\dagger)^{-1}]_{m_1 p}$$
   $$t_2 = - 2 \phi_{b, m_1 j} \langle 0_b | 0_a \rangle [(X^\dagger)^{-1}]_{m_2 p}$$
   $$t_3 = 2 \langle 0_b | 0_a \rangle \sum_k \phi^*_{b, kj} [ [(X^\dagger)^{-1}]_{kp} Z_{b, m_1 m_2} - [(X^\dagger)^{-1}]_{m_2 p} Z_{b, m_1 k} + [(X^\dagger)^{-1}]_{m_1 p} Z_{b, m_2 k} ].$$

For the even-parity vacuum ground state $|\phi_0(G)\rangle = |0_a\rangle$:
3. **Single-quasiparticle transitions ($m$):**
   $$\langle 0_b | b_m c_j | 0_a \rangle = 2 \langle 0_b | 0_a \rangle \sum_k [(X^\dagger)^{-1}]_{mk} \phi_{a, kj}.$$

4. **Particle-Hole Mode Wavefunction Conjugation:**
   Under Eq. (29), because $\tilde{b}_\mu = b_\mu^\dagger$, the expansion $c_j = 2 \sum_m [\tilde{\phi}_{b, mj}^* \tilde{b}_m + \tilde{\phi}_{b, mj} \tilde{b}_m^\dagger]$ demands that the singular row of the mode tensor be conjugated: $\phi_{b, \text{eff}}[\mu, :] = \phi_b[\mu, :]^*$. Matrix elements in the singular sector use $\phi_{b, \text{eff}}$.

5. **Quasiparticle Truncation vs Completeness:**
   For clusters with $N > 2$, $c_j |\phi_0\rangle$ connects to higher quasiparticle sectors (4-, 6-quasiparticles). In arXiv:2209.12930, the multi-quasiparticle expansion is truncated at the 2-quasiparticle sector, which exhausts $> 98\%$ of the dynamical spectral weight for localized Majoranas. For small benchmarks ($N \le 2$), this truncation is exact, satisfying the unitary sum rule $\sum_\lambda |\mathcal{M}_\lambda|^2 = 1.000$ to machine precision.

---

## 3. Guiding Hypotheses & Falsification Criteria
- **`H_BOG_01` (Bogoliubov Unitarity Identities):** $X$ and $Y$ satisfy $X X^\dagger + Y Y^\dagger = I$, $X^\dagger X + Y^T Y^* = I$, and $X Y^T + Y X^T = 0$ to within $10^{-12}$.
- **`H_BOG_02` (Thouless Vacuum Overlap):** The regular overlap $|\langle 0_b | 0_a \rangle|^2 = \sqrt{\det(X^\dagger X)}$ matches exact full Fock-space diagonalization to within $10^{-10}$.
- **`H_BOG_03` (Majorana Completeness Sum Rule):** Transition matrix elements $\mathcal{M}_{j, \lambda} = \langle \lambda(G') | c_j | \phi_0(G) \rangle$ satisfy the exact sum rule $\sum_\lambda |\mathcal{M}_{j, \lambda}|^2 = 1.000$ to within $10^{-10}$.
- **`H_BOG_04` (Particle-Hole Resolution of Singular Overlaps):** When $\det(X) < 10^{-6}$, row interchange $X_{\mu, :} \leftrightarrow Y^*_{\mu, :}$ restores $\det(X') > 0$ and yields exact single-quasiparticle overlaps matching Fock space.

---

## 4. Implementation Findings & Numerical Evidence

### 4.1 Service Architecture
- Built [`src/kitaev_mf/services/phase_3/bogoliubov.py`](file:///home/lucas/Projects/kitaev-mf/src/kitaev_mf/services/phase_3/bogoliubov.py):
  - `extract_mode_wavefunctions`: Projects real Schur matrix $Q$ into complex mode tensors $\phi \in \mathbb{C}^{N \times 2N}$.
  - `verify_bogoliubov_unitarity`: Checks the 3 canonical relations, achieving residuals $< 1.5 \times 10^{-14}$.
  - `compute_bogoliubov_transformation`: Evaluates $X, Y, Z_b$, detecting singular transformations and automatically identifying the optimal mode for particle-hole interchange.
  - `apply_particle_hole_interchange`: Executes Eq. (29) of arXiv:2209.12930.
  - `compute_majorana_matrix_elements`: Implements exact analytical Wick contractions for vacuum and 2-quasiparticle excitation channels.

### 4.2 Automated Unit Test Execution
- Built [`test/unit/test_phase3_bogoliubov_overlaps.py`](file:///home/lucas/Projects/kitaev-mf/test/unit/test_phase3_bogoliubov_overlaps.py):
  1. `test_h_bog_01_canonical_unitarity_relations`: Passed; all 3 identities satisfied to $< 1.0 \times 10^{-14}$.
  2. `test_h_bog_02_thouless_vacuum_overlap_vs_fock_exact`: Passed; Thouless vacuum overlap matched Fock-space exact diagonalization to $10^{-12}$.
  3. `test_h_bog_03_majorana_matrix_elements_and_sum_rule`: Passed; vacuum and 2-particle analytical matrix elements matched Fock-space exact, and total completeness sum rule evaluated to $\sum_\lambda |\mathcal{M}|^2 = 1.000000000000 \pm 10^{-15}$.
  4. `test_h_bog_04_particle_hole_resolution_of_singular_transformation`: Passed; singular mode detected, $\det(X) = 0$, and particle-hole interchange yielded exact single-quasiparticle overlap matching Fock space.
  5. `test_h_bog_05_physical_quench_on_honeycomb_cluster`: Passed; verified on $L=3$ honeycomb torus across sudden spin quench $G_0 \to G'$.
- All 5 tests passed in 0.28s.

---

## 5. Status
- **Deliverable Status:** Completed and Ready for Adversarial Referee Audit.
