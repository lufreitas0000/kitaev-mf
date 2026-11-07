---
id: t05_inhomogeneous_state_and_local_extraction
title: Inhomogeneous Real-Space Mean-Field State Representation and Non-Averaging Correlator Extraction
status: completed
date: 2026-11-07
questions:
  - "How to represent broken translation invariance in real space without sacrificing cache locality and BLAS alignment?"
  - "How to index directed bonds by their source site A on the honeycomb torus so that bond tensors map bijectively to physical coordinates?"
  - "How to extract local Majorana correlators U_{(m1, m2), gamma} and V_{(m1, m2), s} from covariance Omega without spatial averaging?"
  - "What is the exact variational dimension D of the unconstrained inhomogeneous manifold on an L x L torus?"
depends_on:
  - doc/rnd/phase_1/m1_9_t1_matrix_sign_purification_real_space_torus.md
  - doc/rnd/phase_1/t04_purification_and_lbroyden_safeguards.md
contracts:
  - K-020: Exact local Majorana covariance extraction without spatial averaging
  - K-021: Inhomogeneous Hamiltonian assembly preserving local bond and on-site variations
artifacts:
  - src/kitaev_mf/domain/state.py
  - src/kitaev_mf/services/real_space.py
  - test/unit/test_real_space_and_purification.py
---

# Front-1 Research Note: Inhomogeneous Real-Space State and Non-Averaging Local Extraction (Task T05)

## 1. Context and Problem Statement

In previous Phase 1 milestones, real-space cluster solvers (`RealSpaceTorusCluster`) evaluated ground-state Majorana covariance matrices $\mathbf{\Omega} \in \mathbb{R}^{M \times M}$ on finite tori of $N = 2 L^2$ sites ($M = 8 L^2$). However, the extraction routine `extract_mean_field_state` performed an explicit spatial average across all $L^2$ unit cells:
$$U_{\text{avg}} = \frac{1}{L^2} \sum_{m_1, m_2} \mathbf{\Omega}|_{\text{bond}(m_1, m_2)}, \quad V_{\text{avg}} = \frac{1}{L^2} \sum_{m_1, m_2} \mathbf{\Omega}|_{\text{site}(m_1, m_2)}.$$
While appropriate for verifying bulk thermodynamic invariants ($U \to 0.524865$, $E_{\text{site}} \to -0.787$), this spatial averaging completely washes out localized topological excitations (such as isolated vison pairs, vortex branch cuts, and local boundary screening clouds) and site/bond defects.

To study the true physics of topological defects, non-Abelian Majorana zero modes (MZMs), and disorder relaxation without artificial translational symmetry enforcement, we must introduce a native inhomogeneous mean-field state representation `InhomogeneousMeanFieldState` and corresponding assembly/extraction services that retain full spatial resolution.

## 2. Mathematical Specification and Geometry

### 2.1 Lattice Indexing and Directed Bond Geometry
Consider an $L \times L$ torus of $2 L^2$ sites with periodic boundary conditions. Each unit cell $\mathbf{r} = (m_1, m_2)$ where $m_1, m_2 \in \{0, \dots, L-1\}$ contains two sublattices:
- Sublattice $A$ (index $0$): physical site index $s_A = 2(m_1 + L m_2) + 0$.
- Sublattice $B$ (index $1$): physical site index $s_B = 2(m_1 + L m_2) + 1$.

In our canonical orientation, all physical bonds are directed from $A$ to $B$:
1. $z$-bond ($\gamma = 3$, index $0$): $A(m_1, m_2) \to B(m_1, m_2)$, displacement $(0, 0)$.
2. $x$-bond ($\gamma = 1$, index $1$): $A(m_1, m_2) \to B(m_1, m_2 - 1 \pmod L)$, displacement $(0, -1)$.
3. $y$-bond ($\gamma = 2$, index $2$): $A(m_1, m_2) \to B(m_1 - 1 \pmod L, m_2)$, displacement $(-1, 0)$.

Thus, each unit cell $c = m_1 + L m_2 \in \{0, \dots, L^2 - 1\}$ hosts exactly $3$ directed bonds emanating from its $A$ site.

### 2.2 Tensor Architecture of `InhomogeneousMeanFieldState`
To maximize cache locality, prevent stride indirection, and enable seamless vectorization:
- **Bond Correlator Tensor $U$:**
  Shape `(L^2, 3, 4, 4)`, float64.
  $U[c, b, \mu, \nu] = \langle \mathrm{i} c_{A(c)}^\mu c_{B(c, b)}^\nu \rangle = \mathbf{\Omega}[4 s_A + \mu, 4 s_B(c, b) + \nu]$.
- **On-Site Skew-Symmetric Correlator Tensor $V$:**
  Shape `(L^2, 2, 4, 4)`, float64.
  $V[c, \text{sub}, \mu, \nu] = \langle \mathrm{i} c_{s}^\mu c_{s}^\nu \rangle = \mathbf{\Omega}[4 s + \mu, 4 s + \nu]$.
  Subject to strict skew-symmetry: $V[c, \text{sub}]^T = -V[c, \text{sub}]$.
- **Lagrange Multipliers $\boldsymbol{\lambda}$:**
  Shape `(L^2, 2, 3)`, float64.
  $\boldsymbol{\lambda}[c, \text{sub}, \gamma]$ enforces local gauge neutrality $\langle c_s^T G^\gamma c_s \rangle = 0$ on each physical site.

### 2.3 Variational Manifold Dimension
The unconstrained variational parameter dimension is:
$$D = L^2 \times (3 \times 16 + 2 \times 16 + 2 \times 3) = L^2 \times (48 + 32 + 6) = 86 L^2.$$
On an $L = 6$ cluster ($N = 72$ sites): $D = 86 \times 36 = 3096$.
On an $L = 12$ cluster ($N = 288$ sites): $D = 86 \times 144 = 12384$.

## 3. Hamiltonian Assembly and Extraction Invariants

1. **Assembly:**
   Given $U, V, \boldsymbol{\lambda}$, the real skew-symmetric Hamiltonian $A \in \mathbb{R}^{M \times M}$ where $H = \frac{\mathrm{i}}{4} c^T A c$ is assembled locally:
   - For each site $s$:
     $$B_s = \sum_{j \in \mathcal{N}(s)} \sum_{\alpha, \beta} J_{s j}^{\alpha \beta} N^\alpha \operatorname{tr}(V_j^T N^\beta) + \sum_{\gamma=1}^3 (\lambda_s^\gamma G^\gamma - 2 h^\gamma N^\gamma),$$
     placed in $A[\text{slice}(s), \text{slice}(s)]$.
   - For each directed bond $b = (c, \text{type})$ connecting $s_A$ to $s_B$:
     $$A_b = 2 \sum_{\alpha, \beta} J_b^{\alpha \beta} N^\alpha U_{c, \text{type}} N^\beta,$$
     placed in $A[\text{slice}(s_A), \text{slice}(s_B)] = A_b$ and $A[\text{slice}(s_B), \text{slice}(s_A)] = -A_b^T$.
2. **Extraction:**
   Given ground-state covariance $\mathbf{\Omega}$, local extraction populates:
   $$U[c, b] = \mathbf{\Omega}[\text{slice}(s_A), \text{slice}(s_B)], \quad V[c, \text{sub}] = \frac{1}{2}\left(\mathbf{\Omega}[\text{slice}(s), \text{slice}(s)] - \mathbf{\Omega}[\text{slice}(s), \text{slice}(s)]^T\right).$$
   Zero spatial averaging is performed. Any spatial perturbation (e.g. inverted gauge link, vacancy, or localized Zeeman field) is preserved with point-wise fidelity.
