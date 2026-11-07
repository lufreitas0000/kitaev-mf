---
title: "Task T3.1.1: Static Flux Configurations, Dual Branch Cuts & Real-Space Two-Vortex State Generation"
status: completed
author: "Lu Freitas"
date: "2026-11-07"
phase: "Phase 3"
milestone: "M3.1"
task: "T3.1.1"
contracts: ["H_VOR_01", "H_VOR_02", "H_VOR_03", "H_VOR_04"]
primary_model: "Pure FM Kitaev (K > 0) with chiral 3-spin term kappa (arXiv:2209.12930)"
geometry: "Real-space honeycomb lattice torus of L x L unit cells"
guiding_questions:
  - "How are two isolated Z2 vortices constructed on a periodic honeycomb torus via minimal-length dual branch cuts?"
  - "How do we verify that the Wilson loop W_p = -1 on exactly the two target plaquettes and +1 everywhere else?"
  - "How does the single-particle BdG spectrum respond to two well-separated vortices in the chiral phase (kappa > 0), and does an isolated Majorana zero mode (MZM) emerge?"
  - "How does a local spin operator sigma_j^alpha act on the two-vortex reference state: vortex translation vs. additional vortex-pair creation?"
---

# Task T3.1.1: Static Flux Configurations, Dual Branch Cuts & Real-Space Two-Vortex State Generation

## 1. Scope & Mathematical Objectives
In **arXiv:2209.12930** (*Scanning tunneling spectroscopy of Majorana zero modes in a Kitaev spin liquid*, Phys. Rev. B **107**, 054432), the central calculation of Majorana zero mode (MZM) scanning tunneling spectroscopy is performed in a reference gauge configuration $G$ containing **two well-separated $\mathbb{Z}_2$ vortices (visons)** on a periodic $L \times L$ honeycomb lattice.

The objective of this atomic task is to develop and verify the exact geometric and gauge-theoretic engine for:
1. Generating minimal-length branch cut strings on the dual triangular lattice connecting any two target plaquettes $p_1 = (m_1, m_2)$ and $p_2 = (n_1, n_2)$.
2. Constructing the static two-vortex gauge background $G(p_1, p_2)$ with $W_{p_1} = -1$, $W_{p_2} = -1$, and $W_p = +1$ for all other $L^2 - 2$ plaquettes.
3. Quantifying the isolated in-gap Majorana zero mode (MZM) hybridization $\epsilon_0(d_{12}) \sim \mathrm{e}^{-d_{12}/\xi_M}$ in the chiral phase ($\kappa > 0$).
4. Categorizing the action of local spin operators $\sigma_j^\alpha$:
   - **Vortex translation:** If bond $\langle jl \rangle_\alpha$ is adjacent to $p_1$, toggling $u_{jl}^\alpha$ shifts the vortex from $p_1$ to the neighboring plaquette sharing that bond, maintaining a 2-vortex sector ($W_{p_1} \to +1, W_{p_1'} \to -1$).
   - **Vortex-pair creation:** If bond $\langle jl \rangle_\alpha$ is in the bulk away from $p_1$ and $p_2$, toggling $u_{jl}^\alpha$ creates an additional pair of vortices, transitioning to a 4-vortex sector ($W_{p_1} = -1, W_{p_2} = -1, W_{p_3} = -1, W_{p_4} = -1$).

---

## 2. Theoretical Architecture & Plaquette Geometry

### 2.1 Hexagonal Plaquette Boundary
On an $L \times L$ unit-cell torus with periodic boundary conditions, each unit cell $(m_1, m_2)$ contains one hexagonal plaquette $p = (m_1, m_2)$.
The boundary $\partial p$ consists of 6 directed bonds:
1. $z$-bond at cell $(m_1, m_2)$: connecting $s_A(m_1, m_2) \to s_B(m_1, m_2)$.
2. $x$-bond at cell $(m_1, m_2 + 1)$: connecting $s_A(m_1, m_2 + 1) \to s_B(m_1, m_2)$.
3. $y$-bond at cell $(m_1, m_2 + 1)$: connecting $s_A(m_1, m_2 + 1) \to s_B(m_1 - 1, m_2 + 1)$.
4. $z$-bond at cell $(m_1 - 1, m_2 + 1)$: connecting $s_A(m_1 - 1, m_2 + 1) \to s_B(m_1 - 1, m_2 + 1)$.
5. $x$-bond at cell $(m_1 - 1, m_2 + 1)$: connecting $s_A(m_1 - 1, m_2 + 1) \to s_B(m_1 - 1, m_2)$.
6. $y$-bond at cell $(m_1, m_2)$: connecting $s_A(m_1, m_2) \to s_B(m_1 - 1, m_2)$.

The gauge-invariant Wilson loop operator is:
$$W_p = \prod_{b \in \partial p} u_b \in \{-1, +1\}$$
where $u_b = u_{AB} \in \{-1, +1\}$ is the link variable on bond $b$.

### 2.2 Dual Graph & Minimal Branch Cut
Every bond $b$ is shared by exactly two adjacent plaquettes $p_A, p_B$. The dual graph whose vertices are plaquettes and whose edges are shared bonds forms a 2D triangular lattice with coordination number $z = 6$.
Given two target plaquettes $p_1$ and $p_2$, the minimal branch cut is the shortest path on the dual triangular graph, computed via Breadth-First Search (BFS) under torus periodicity.
Inverting the link variables along this path:
$$u_b \longrightarrow -u_b \quad \forall b \in \mathrm{Cut}(p_1, p_2)$$
guarantees by telescoping cancellation:
$$W_{p_1} \longrightarrow -W_{p_1}, \qquad W_{p_2} \longrightarrow -W_{p_2}, \qquad W_p \longrightarrow W_p \quad (\forall p \neq p_1, p_2).$$

### 2.3 Majorana Zero Modes & In-Gap Spectrum
In the presence of broken time-reversal symmetry ($\kappa > 0$), each isolated $\mathbb{Z}_2$ vortex binds a single localized Majorana bound state:
$$\gamma_v = \sum_j \psi_v(j) c_j.$$
For two well-separated vortices at distance $d_{12} \gg \xi_M$, the two Majoranas $\gamma_1, \gamma_2$ hybridize to form one Dirac mode:
$$a_0 = \frac{1}{2}(\gamma_1 - \mathrm{i}\gamma_2)$$
with energy splitting decaying exponentially:
$$\epsilon_0(d_{12}) \propto \mathrm{e}^{-d_{12}/\xi_M}, \qquad \xi_M \approx \frac{K}{\Delta_{\mathrm{bulk}}} = \frac{1}{6\sqrt{3}\kappa}$$
where $\Delta_{\mathrm{bulk}} = 6\sqrt{3}\kappa$ is the full bulk band gap, and the single-particle BdG Dirac threshold is $\epsilon_{\mathrm{Dirac}} = 3\sqrt{3}\kappa = \frac{1}{2}\Delta_{\mathrm{bulk}}$.

---

## 3. Hypothesis Contracts & Falsification Criteria
- **`H_VOR_01` (Dual Branch Cut Exactness):** Toggling links along $\mathrm{Cut}(p_1, p_2)$ inverts flux on $p_1$ and $p_2$ exclusively: $W_{p_1} = -1, W_{p_2} = -1$, and $W_p = +1$ on all remaining $L^2 - 2$ plaquettes.
- **`H_VOR_02` (MZM In-Gap Isolation):** In the chiral phase ($\kappa > 0$), inserting two well-separated vortices creates a sub-gap eigenvalue $\epsilon_0 \ll \epsilon_{\mathrm{bulk}}$, with sub-gap ratio $\epsilon_0 / \epsilon_{\mathrm{clean}} < 0.25$.
- **`H_VOR_03` (Vortex Translation Selection Rule):** Toggling a bond link $u_b$ adjacent to a vortex plaquette $p_1$ shifts the vortex to the adjacent plaquette $p_1'$, preserving the total number of vortices $N_v = 2$.
- **`H_VOR_04` (Pair Creation Selection Rule):** Toggling a bulk bond link unshared with $p_1$ or $p_2$ creates two additional vortices, transitioning the system to $N_v = 4$.

---

## 4. Implementation Findings & Numerical Evidence

### 4.1 Implementation Architecture
- Built [`src/kitaev_mf/services/phase_3/vortex_engine.py`](file:///home/lucas/Projects/kitaev-mf/src/kitaev_mf/services/phase_3/vortex_engine.py):
  - `HoneycombGeometry`: Implements exact dual triangular graph adjacency, shortest-path BFS branch cuts, and torus Euclidean distances.
  - `create_two_vortex_gauge_configuration`: Generates two-vortex configurations via minimal branch-cut inversion.
  - `evaluate_all_plaquette_fluxes` & `identify_vortices`: Evaluates all $L^2$ Wilson loops $W_p = \prod_{b \in \partial p} u_b$.
  - `apply_local_spin_quench` & `categorize_spin_quench`: Formulates sudden quench $\sigma_j^\alpha$ into vortex translation ($N_v = 2 \to 2$) or pair creation ($N_v = 2 \to 4$).

### 4.2 Automated Unit Test Execution
- Built [`test/unit/test_phase3_vortex_configurations.py`](file:///home/lucas/Projects/kitaev-mf/test/unit/test_phase3_vortex_configurations.py):
  1. `test_h_vor_01_dual_branch_cut_exactness`: Verified on $L=4$ (16 plaquettes) that flipping links along the BFS cut produces $W_{(0, 0)} = -1, W_{(2, 2)} = -1$, while the remaining 14 plaquettes satisfy $W_p = +1$.
  2. `test_h_vor_02_mzm_in_gap_isolation`: On $L=4$ with $K=1.0, \kappa=0.15$, clean vacuum mode is $\epsilon_{\mathrm{clean}} \approx 1.400$. Two vortices bind an in-gap mode $\epsilon_0 \approx 0.2946$ (sub-gap ratio $0.210 < 0.25$).
  3. `test_h_vor_03_vortex_translation_vs_pair_creation`: Quenching an adjacent bond shifted vortex $(0, 0) \to (0, 1)$ ($N_v = 2 \to 2$), while quenching a distant bulk bond created an additional pair ($N_v = 2 \to 4$).
  4. `test_h_vor_04_torus_distance_and_geometry`: Verified distance symmetry, positive-definiteness, and periodic image equivalence.
- All 4 tests passed in 0.80s.

---

## 5. Status
- **Deliverable Status:** Completed and Ready for Adversarial Referee Audit.
