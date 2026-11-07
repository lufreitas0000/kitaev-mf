---
title: "Task T3.0.1: Theoretical Foundations, Hamiltonian Conventions & Physical Parity Formulation"
status: completed
author: "Lu Freitas"
date: "2026-11-07"
phase: "Phase 3"
milestone: "M3.0"
task: "T3.0.1"
contracts: ["H_DYN_01", "H_DYN_02", "H_DYN_03"]
primary_model: "Pure FM Kitaev (K > 0) with chiral 3-spin term kappa (arXiv:2209.12930)"
geometry: "Real-space honeycomb lattice clusters (torus and cylinder geometries)"
guiding_questions:
  - "What are the exact sign conventions and prefactors connecting spin operators, Majorana CAR, and the skew-symmetric matrix A?"
  - "How is the next-nearest-neighbor chiral interaction traversed, and what is its analytical Dirac mass gap?"
  - "How is the physical matter-fermion parity formulated on periodic tori per Pedrocchi et al. (2011) and Zschocke & Vojta (2015)?"
  - "How does a single-spin gauge quench perturb the Z2 flux background?"
---

# Task T3.0.1: Theoretical Foundations, Hamiltonian Conventions & Physical Parity Formulation

## 1. Scope & Mathematical Objectives
This atomic task establishes the rigorous theoretical foundation, exact sign conventions, operator representations, and physical matter-fermion parity constraints required to reproduce the dynamical spin correlation functions and scanning tunneling spectroscopy of Majorana zero modes from **arXiv:2209.12930** (*Scanning tunneling spectroscopy of Majorana zero modes in a Kitaev spin liquid*, Phys. Rev. B **107**, 054432).

### 1.1 Model Hamiltonian
The physical spin Hamiltonian studied in the target paper is:
$$H_K = -K \sum_{\langle j l \rangle_\alpha} \sigma_j^\alpha \sigma_l^\alpha - \kappa \sum_{\langle j k \rangle_\alpha, \langle k l \rangle_\beta} \sigma_j^\alpha \sigma_k^\gamma \sigma_l^\beta$$
where $\alpha, \beta, \gamma \in \{x, y, z\}$ are mutually distinct bond types sharing site $k$, $K > 0$ (ferromagnetic Kitaev exchange, e.g. $K = 1.0$), and $\kappa$ is the effective time-reversal-breaking 3-spin chiral coupling induced by magnetic field $\mathbf{h} \parallel [111]$:
$$\kappa \approx \frac{h_x h_y h_z}{\Delta_{\mathrm{flux}}^2}.$$

### 1.2 Exact Majorana & Gauge Representation
On each site $j$, the spin operator is represented in terms of four Majorana operators $(c_j, b_j^x, b_j^y, b_j^z)$ satisfying the canonical Clifford CAR algebra:
$$\{c_j, c_l\} = 2\delta_{jl}, \quad \{b_j^\alpha, b_l^\beta\} = 2\delta_{jl}\delta_{\alpha\beta}, \quad \{c_j, b_l^\alpha\} = 0.$$
The spin operators are:
$$\sigma_j^\alpha = \mathrm{i} b_j^\alpha c_j.$$
The bond operators define static $\mathbb{Z}_2$ link variables on nearest-neighbor $\alpha$-bonds connecting site $j \in \mathcal{A}$ to site $l \in \mathcal{B}$:
$$\hat{u}_{jl}^\alpha = \mathrm{i} b_j^\alpha b_l^\alpha = \pm 1.$$
Under directed bond reversal, $\hat{u}_{lj}^\alpha = \mathrm{i} b_l^\alpha b_j^\alpha = -\mathrm{i} b_j^\alpha b_l^\alpha = -\hat{u}_{jl}^\alpha$.

### 1.3 Exact Hamiltonian Sign Derivation (`[B1]`)
Let $j \in \mathcal{A}$ and $l \in \mathcal{B}$. The two-spin exchange operator evaluates as:
$$\sigma_j^\alpha \sigma_l^\alpha = (\mathrm{i} b_j^\alpha c_j)(\mathrm{i} b_l^\alpha c_l) = b_j^\alpha b_l^\alpha c_j c_l = (-\mathrm{i} \hat{u}_{jl}^\alpha) c_j c_l = -\mathrm{i} \hat{u}_{jl}^\alpha c_j c_l.$$
Consequently, the ferromagnetic exchange interaction term with $K > 0$ yields:
$$-K \sigma_j^\alpha \sigma_l^\alpha = -K (-\mathrm{i} \hat{u}_{jl}^\alpha c_j c_l) = +\mathrm{i} K \hat{u}_{jl}^\alpha c_j c_l.$$
Writing the matter Hamiltonian as a quadratic form in itinerant Majoranas:
$$H_G = \frac{\mathrm{i}}{4} \sum_{m, n} A_{mn}(G) c_m c_n$$
the contribution of the pair $(j, l)$ is:
$$\frac{\mathrm{i}}{4} [ A_{jl} c_j c_l + A_{lj} c_l c_j ] = \frac{\mathrm{i}}{2} A_{jl} c_j c_l.$$
Matching this to $+\mathrm{i} K u_{jl}^\alpha c_j c_l$ requires:
$$A_{jl} = +2 K u_{jl}^\alpha, \qquad A_{lj} = -2 K u_{jl}^\alpha \quad (j \in \mathcal{A}, l \in \mathcal{B}).$$
This matches Eq. (265) of arXiv:2209.12930.

### 1.4 Directed Traversal of the Chiral 3-Spin Term (`[B2]`)
For the chiral term $-\kappa \sigma_j^\alpha \sigma_k^\gamma \sigma_l^\beta$ connecting $j \in \mathcal{A} \to k \in \mathcal{B} \to l \in \mathcal{A}$ with right-handed bond orientation $\epsilon^{\alpha \gamma \beta} = +1$:
$$\sigma_j^\alpha \sigma_k^\gamma \sigma_l^\beta = (\mathrm{i} b_j^\alpha c_j)(\mathrm{i} b_k^\gamma c_k)(\mathrm{i} b_l^\beta c_l) = -\mathrm{i} u_{jk}^\alpha u_{kl}^\beta c_j c_l.$$
Because the bond between $k \in \mathcal{B}$ and $l \in \mathcal{A}$ is canonically defined with orientation $l \in \mathcal{A} \to k \in \mathcal{B}$, traversing backwards yields:
$$u_{kl}^\beta = -u_{lk}^\beta.$$
Thus, in terms of standard links $u_{jk}^\alpha$ and $u_{lk}^\beta$:
$$u_{jk}^\alpha u_{kl}^\beta = - u_{jk}^\alpha u_{lk}^\beta.$$
In the uniform vortex-free gauge where $u_{AB} = +1$, this product is $-1$.
Matching $-\mathrm{i}\kappa u_{jk}^\alpha u_{kl}^\beta c_j c_l$ to $\frac{\mathrm{i}}{2} A_{jl} c_j c_l$ produces:
$$A_{jl} += -2\kappa (u_{jk}^\alpha u_{kl}^\beta) \delta_{\langle\langle jl \rangle\rangle}.$$
At the Dirac points $\mathbf{K}$ and $\mathbf{K}'$, this opens the exact topological bulk Dirac mass gap:
$$\epsilon_{\mathrm{Dirac}} = 3\sqrt{3}\kappa \implies \Delta_{\mathrm{bulk}} = 6\sqrt{3}\kappa.$$

### 1.5 Physical Matter-Fermion Parity (`[B3]`)
The physical states $|\Psi\rangle$ are projectors under $D_j = b_j^x b_j^y b_j^z c_j$:
$$\mathcal{P}_0 = \prod_{j=1}^{2N} \frac{1 + D_j}{2} \implies \hat{D} = \prod_j D_j \overset{!}{=} +1.$$
Per Pedrocchi, Chesi, and Loss (Phys. Rev. B **84**, 165414, 2011) and Zschocke & Vojta (Phys. Rev. B **92**, 014403, 2015):
$$\hat{D} = (-1)^\theta (-1)^{N_\chi} \det(Q^u) \hat{\pi}$$
where:
- $\hat{\pi} = (-1)^{N_a}$ is the matter fermion parity of the BdG quasiparticle eigenmodes.
- $(-1)^{N_\chi} = \prod_b u_b$ tracks flipped gauge links relative to the uniform gauge.
- $\det(Q^u) \in \{-1, +1\}$ is the determinant of the real Schur canonical transformation $A = Q \Sigma Q^T$.
- $(-1)^\theta$ is the topological phase from reordering Majorana operators ($(-1)^\theta = +1$ for $L_1 = L_2 = L$).

**Zschocke & Vojta Universal Theorem:** For any clean translation-invariant honeycomb torus in the gapless or chiral phase:
$$(-1)^\theta \det(Q^u) = -1 \quad \implies \quad \pi_{G_0} = -1.$$
All physical states in the vortex-free sector have **odd** matter fermion parity ($\hat{\pi} = -1$). The naive unprojected vacuum $|0\rangle$ is unphysical; the true physical ground state contains one matter fermion excitation in the lowest single-particle mode $\epsilon_1$:
$$E_{G_0, 0}^{\mathrm{phys}} = E_{G_0, 0} + \epsilon_1.$$

---

## 2. Guiding Hypotheses & Falsification Tests
- **`H_DYN_01` (Gauge Quench Locality):** Toggling link $u_{jl}^\alpha \to -u_{jl}^\alpha$ flips the gauge flux $W_p \to -W_p$ on exactly two adjacent plaquettes and leaves all other plaquettes undisturbed.
- **`H_DYN_02` (Chiral Dirac Bulk Gap):** For $K=1.0$, the chiral 3-spin term opens an exact topological Dirac mass gap $\epsilon_{\mathrm{min}} = 3\sqrt{3}\kappa$, giving bulk band gap $\Delta_{\mathrm{bulk}} = 6\sqrt{3}\kappa$.
- **`H_DYN_03` (Pedrocchi–Zschocke Parity Theorem):** On a periodic torus in the vortex-free sector, $\det(Q^u) = -1$, yielding odd physical matter parity $\pi_G = -1$ and shifting the physical ground state energy to $E_0 + \epsilon_1$.

---

## 3. Implementation Summary
- **Module:** [`src/kitaev_mf/services/phase_3/gauge_theory.py`](file:///home/lucas/Projects/kitaev-mf/src/kitaev_mf/services/phase_3/gauge_theory.py)
  - `DirectedBond`: Strict bipartite $A \to B$ mapping avoiding undirected key ambiguities (`[B4]`).
  - `GaugeConfiguration`: Directed link storage with exact anti-symmetry $u_{ba} = -u_{ab}$.
  - `assemble_matter_skew_matrix`: Implements $+2 K u_{ab}$ and directed chiral traversal (`[B1]`, `[B2]`).
  - `diagonalize_matter_hamiltonian`: Implements canonical real Schur decomposition, computes $\det(Q^u)$, and extracts the exact physical matter parity $\pi_G$ (`[B3]`).
  - `build_honeycomb_torus_cluster`: Generates $L \times L$ periodic honeycomb clusters with verified triplet chiralities.
- **Unit Test Suite:** [`test/unit/test_phase3_theory_conventions.py`](file:///home/lucas/Projects/kitaev-mf/test/unit/test_phase3_theory_conventions.py)
  - 4 automated unit tests verifying locality, exact chiral gap $3\sqrt{3}\kappa$, Zschocke-Vojta parity theorem, and directed link anti-symmetry (`[B5]`).
  - All tests passed in 0.26s.

---

## 4. Status & Audit Outcome
- **Deliverable Status:** Completed, Reconciled, and Verified.
