---
title: "Task T2.1.3: Variational Projection Bounds — Gutzwiller Overlap N_G(h) on ED Clusters & 2-RDM Semidefinite Lower Bounds"
task_id: "T2.1.3"
milestone: "M2.1"
author: "Lu Freitas"
date: "2026-11-07"
status: "in_progress"
model_selection:
  candidate_model: "M1: Gaussian-projected Gutzwiller overlap N_G(h) = N_0 * exp(-kappa * N * h^2)"
  null_model: "M0: Linear Gutzwiller relaxation N_G(h) = N_0 * (1 - c * h)"
  target_delta_bic: 10.0
  target_bayes_factor: 148.0
energy_bracket_contracts:
  strict_inequality: "E_SDP_2(h) <= E_exact(h) <= E_proj(h)"
  max_duality_gap: 0.15
  ed_cluster_geometries:
    - name: "hexagon_6"
      num_spins: 6
      pbc: false
    - name: "honeycomb_torus_8"
      num_spins: 8
      pbc: true
popperian_falsification:
  support_if:
    - "E_SDP_2(h) <= E_exact(h) + 1e-12 for all h in [0.0, 0.40]"
    - "E_exact(h) <= E_proj(h) + 1e-12 for all h in [0.0, 0.40]"
    - "Delta_gap(h) = (E_proj - E_SDP_2) / |E_exact| <= 0.15 across the intermediate field regime h <= 0.35"
    - "Delta_BIC_01 > 10.0 decisively favoring quadratic exponential decay M1 over linear M0"
  refute_if:
    - "E_proj(h) < E_exact(h) - 1e-10 (Rayleigh-Ritz variational violation)"
    - "E_SDP_2(h) > E_exact(h) + 1e-10 (Semidefinite relaxation bound violation)"
    - "Delta_BIC_01 <= 0.0 or Bayes factor B_10 < 1.0"
---

# Task T2.1.3: Variational Projection Bounds — Gutzwiller Overlap $\mathcal{N}_G(h)$ on ED Clusters & 2-RDM Semidefinite Lower Bounds

## 1. Executive Summary & Context

Variational mean-field theories for quantum spin liquids express physical spin-1/2 operators in an enlarged Majorana Fock space:
$$\sigma_j^\alpha = -\frac{i}{2} \epsilon_{\alpha\beta\gamma} c_j^\beta c_j^\gamma = i c_j^0 c_j^\alpha$$
subject to the local gauge constraint:
$$\mathcal{D}_j = -c_j^0 c_j^x c_j^y c_j^z = +1 \quad \forall j \in \{1, \dots, N\}.$$
The bare 4-Majorana saddle point represents the ground state by a Gaussian fermionic state $|\Psi_{\mathrm{MF}}\rangle$. Because the local gauge constraint is enforced only on average via Lagrange multipliers $\langle \mathcal{D}_j \rangle = 1$, the unprojected state $|\Psi_{\mathrm{MF}}\rangle$ contains unphysical states belonging to the $\mathcal{D}_j = -1$ sector.

The physical wave function is recovered by the Gutzwiller projector:
$$\mathcal{P}_G = \prod_{j=1}^N \left( \frac{1 + \mathcal{D}_j}{2} \right), \quad |\Psi_{\mathrm{proj}}\rangle = \frac{\mathcal{P}_G |\Psi_{\mathrm{MF}}\rangle}{\sqrt{\mathcal{N}_G}},$$
where $\mathcal{N}_G = \langle \Psi_{\mathrm{MF}} | \mathcal{P}_G | \Psi_{\mathrm{MF}} \rangle$ is the Gutzwiller norm overlap.

In this task, we establish the rigorous mathematical foundation of the mean-field variational quality by computing:
1. **Exact Spin-Space Diagonalization (ED):** Ground state energy $E_{\mathrm{exact}}(h)$ and wave function $|\Psi_{\mathrm{exact}}\rangle$ on $N=6$ (isolated hexagon) and $N=8$ ($2 \times 2$ periodic honeycomb torus) clusters under field $\mathbf{h} \parallel [111]$.
2. **Gutzwiller Projection:** Exact Fock-space embedding of $|\Psi_{\mathrm{MF}}\rangle$, evaluation of norm overlap $\mathcal{N}_G(h)$, and projected physical energy $E_{\mathrm{proj}}(h) = \langle \Psi_{\mathrm{proj}} | H_{\mathrm{spin}} | \Psi_{\mathrm{proj}} \rangle$.
3. **2-Body Reduced Density Matrix (2-RDM) Semidefinite Lower Bounds ($E_{\mathrm{SDP}}^{(2)}$):** A convex semidefinite relaxation bounding the exact ground state energy strictly from below.
4. **Variational Duality Bracket:**
   $$E_{\mathrm{SDP}}^{(2)}(h) \le E_{\mathrm{exact}}(h) \le E_{\mathrm{proj}}(h).$$
5. **Bayesian Model Selection on Gutzwiller Overlap:** Testing whether field-induced projection decay is quadratic-exponential $\mathcal{N}_G(h) = \mathcal{N}_0 \mathrm{e}^{-\kappa N h^2}$ ($\mathcal{M}_1$) or linear relaxation $\mathcal{N}_G(h) = \mathcal{N}_0 (1 - c h)$ ($\mathcal{M}_0$).

---

## 2. Mathematical Formalism

### 2.1 Fock-Space Construction and Local Constraint
Each site $j$ hosts 4 Majorana fermions $\{c_j^0, c_j^1, c_j^2, c_j^3\}$, which can be paired into 2 complex Dirac fermions:
$$f_{j,1} = \frac{1}{2}(c_j^0 + i c_j^1), \quad f_{j,2} = \frac{1}{2}(c_j^2 + i c_j^3).$$
The local fermionic Fock space $\mathcal{H}_j$ has dimension $2^2 = 4$, with basis:
$$\{|00\rangle_j, |10\rangle_j, |01\rangle_j, |11\rangle_j\}.$$
The local constraint operator $\mathcal{D}_j = -c_j^0 c_j^1 c_j^2 c_j^3$ acts on this basis as:
$$\mathcal{D}_j = (-1)^{n_{j,1} + n_{j,2}} = \begin{pmatrix} +1 & 0 & 0 & 0 \\ 0 & -1 & 0 & 0 \\ 0 & 0 & -1 & 0 \\ 0 & 0 & 0 & +1 \end{pmatrix}.$$
Therefore, the physical subspace $\mathcal{D}_j = +1$ is spanned by the even-parity states $\{|00\rangle_j, |11\rangle_j\}$, corresponding to spin-1/2 $|\!\uparrow\rangle_j$ and $|\!\downarrow\rangle_j$:
$$|\!\uparrow\rangle_j \equiv |00\rangle_j, \quad |\!\downarrow\rangle_j \equiv |11\rangle_j.$$
The odd-parity states $\{|10\rangle_j, |01\rangle_j\}$ have $\mathcal{D}_j = -1$ and are unphysical gauge artifacts.

### 2.2 Semidefinite Programming (SDP) 2-RDM Lower Bound
The exact spin Hamiltonian is:
$$H_{\mathrm{spin}} = -K \sum_{\langle ij \rangle_\gamma} \sigma_i^\gamma \sigma_j^\gamma - \sum_{i} \mathbf{h} \cdot \boldsymbol{\sigma}_i.$$
Let $\chi$ be the two-body spin correlation matrix:
$$\chi_{(i\alpha), (j\beta)} = \langle \sigma_i^\alpha \sigma_j^\beta \rangle, \quad m_i^\alpha = \langle \sigma_i^\alpha \rangle.$$
For any physical state, the block matrix:
$$\mathcal{M}_1 = \begin{pmatrix} 1 & \mathbf{m}^T \\ \mathbf{m} & \chi \end{pmatrix} \succeq 0.$$
Furthermore, for every pair of sites $(i, j)$, the two-site density matrix $\rho_{ij}$ in $\mathbb{C}^{4 \times 4}$:
$$\rho_{ij} = \frac{1}{4} \left[ I_4 + \sum_\alpha m_i^\alpha (\sigma^\alpha \otimes I_2) + \sum_\beta m_j^\beta (I_2 \otimes \sigma^\beta) + \sum_{\alpha\beta} \chi_{(i\alpha), (j\beta)} (\sigma^\alpha \otimes \sigma^\beta) \right]$$
must satisfy the partial trace and positivity conditions:
$$\rho_{ij} \succeq 0, \quad \operatorname{Tr}(\rho_{ij}) = 1.$$
Minimizing $\langle H_{\mathrm{spin}} \rangle = \sum_{\langle ij \rangle_\gamma} -K \chi_{(i\gamma),(j\gamma)} - \sum_{i\alpha} h_\alpha m_i^\alpha$ subject to these linear matrix inequalities (LMIs) yields the 2-RDM lower bound:
$$E_{\mathrm{SDP}}^{(2)} = \min_{\mathbf{m}, \chi} \operatorname{Tr}(H \rho) \quad \text{s.t.} \quad \mathcal{M}_1 \succeq 0, \; \rho_{ij} \succeq 0 \; \forall \langle ij \rangle, \; \chi_{(i\alpha),(i\alpha)} = 1.$$
Since every physical density matrix $\rho$ satisfies these conditions, the SDP minimum is a rigorous lower bound:
$$E_{\mathrm{SDP}}^{(2)} \le E_{\mathrm{exact}}.$$

---

## 3. Pre-Registered Hypotheses & Falsification Criteria

1. **Rayleigh-Ritz Variational Safety:**
   $$E_{\mathrm{proj}}(h) \ge E_{\mathrm{exact}}(h) - 10^{-12}.$$
2. **2-RDM Semidefinite Lower Bound:**
   $$E_{\mathrm{SDP}}^{(2)}(h) \le E_{\mathrm{exact}}(h) + 10^{-12}.$$
3. **Duality Gap Bound:**
   $$\Delta_{\mathrm{gap}}(h) = \frac{E_{\mathrm{proj}}(h) - E_{\mathrm{SDP}}^{(2)}(h)}{|E_{\mathrm{exact}}(h)|} \le 0.15 \quad \forall h \in [0.0, 0.35].$$
4. **Bayesian Model Selection on Gutzwiller Overlap:**
   Quadratic model $\mathcal{M}_1: \ln \mathcal{N}_G(h) = \ln \mathcal{N}_0 - \kappa N h^2$ versus Linear null model $\mathcal{M}_0: \mathcal{N}_G(h) = \mathcal{N}_0(1 - c h)$.
   Decisive selection requires:
   $$\Delta\mathrm{BIC}_{01} = \mathrm{BIC}_0 - \mathrm{BIC}_1 > 10.0, \quad B_{10} = \exp(\Delta\mathrm{BIC}_{01} / 2) > 148.0.$$
