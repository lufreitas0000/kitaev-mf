---
id: t06_unpinned_vortex_self_consistency_and_relaxation
title: Unpinned Vortex Self-Consistency, Topological Trapping vs Healing, and Inhomogeneous Disorder Relaxation
status: completed
date: 2026-11-07
questions:
  - "Does a two-vortex branch cut in the inhomogeneous state relax stably without an artificial post-seed flipped_bonds pinning mask (topological trapping vs dynamical healing)?"
  - "Under what conditions on Zeeman field h and lattice size L do the local gauge links u_b = <i c_A^gamma c_B^gamma> remain topologically trapped in the pi-flux sector (W_{p1} = W_{p2} = -1) vs healing to the vacuum (W_p = +1)?"
  - "How do local non-magnetic site vacancies and exchange bond disorder relax within the unpinned inhomogeneous manifold R^{86 L^2}?"
depends_on:
  - doc/rnd/phase_1/t05_inhomogeneous_state_and_local_extraction.md
contracts:
  - K-020: Exact local Majorana covariance extraction without spatial averaging
  - K-021: Inhomogeneous Hamiltonian assembly preserving local bond and on-site variations
  - K-022: Topological stability criterion for unpinned two-vison states
artifacts:
  - src/kitaev_mf/services/real_space.py
  - test/unit/test_real_space_vortices.py
---

# Front-1 Research Note: Unpinned Vortex Self-Consistency, Trapping vs. Healing, and Inhomogeneous Disorder Relaxation (Task T06)

## 1. Context and Problem Statement

In previous implementations, topological vison defects were introduced by explicitly passing a branch cut mask `flipped_bonds = cut` into `assemble_skew_hamiltonian` on *every* self-consistent iteration. This external constraint artificially pinned the gauge link signs:
$$\operatorname{sign}(u_b) = -1.0 \quad \forall b \in \text{cut},$$
preventing the physical mean-field variational manifold from testing whether an isolated vortex pair represents a genuine, self-sustaining metastable local minimum of the real Grassmannian $\operatorname{Gr}(4 L^2, 8 L^2)$.

To resolve Defect `[B1]` Part 2 and establish Contract `K-022`, we must answer:
1. When a branch cut is inverted solely in the initial variational seed:
   $$U_{c, b}^{\gamma, \gamma} \to -U_{c, b}^{\gamma, \gamma} \quad \text{for } b \in \text{cut}(p_1, p_2),$$
   and iterated self-consistently with `flipped_bonds = None` (unpinned), does the state remain topologically trapped in the two-vison sector ($W_{p1} = W_{p2} = -1$), or does it dynamically heal back to the zero-flux vacuum ($W_p = +1$)?
2. At finite field $h \parallel [111]$, what is the stability threshold between topological trapping and high-field field-induced healing?
3. How do non-magnetic site vacancies and exchange bond disorder relax locally without translational dilution?

## 2. Mathematical Formalism and Topological Trapping

### 2.1 Pure Kitaev Limit ($h = 0$)
In the pure Kitaev honeycomb model at zero magnetic field, the gauge Majoranas $c_i^\gamma$ and $c_j^\gamma$ on bond $\langle ij \rangle_\gamma$ decouple algebraically from all other lattice sites in the mean-field Hamiltonian $A \in \mathbb{R}^{8 L^2 \times 8 L^2}$. The $2 \times 2$ subspace spanned by $(c_{s_A}^\gamma, c_{s_B}^\gamma)$ forms an isolated two-level skew block:
$$A|_{\text{gauge bond } b} = \begin{pmatrix} 0 & 2 K U_b^{\gamma, \gamma} \\ -2 K U_b^{\gamma, \gamma} & 0 \end{pmatrix}.$$
Because $K < 0$, if $U_b^{\gamma, \gamma} = -1.0$, the eigenvalues of $\mathrm{i} A$ are $\pm 2 |K| U_b^{\gamma, \gamma} = \mp 2 |K|$. The negative-energy ground-state projector therefore yields:
$$\mathbf{\Omega}[\text{slice}(s_A), \text{slice}(s_B)]_{\gamma, \gamma} = \operatorname{sgn}(K U_b^{\gamma, \gamma}) = -1.0.$$
Consequently, in the absence of a magnetic field, the inverted gauge links are exact stationary fixed points of the unpinned mean-field map:
$$\Phi(U_b^{\gamma, \gamma} = -1.0) = -1.0.$$
The two visons at $p_1$ and $p_2$ are rigorously trapped with exact $\pi$-fluxes:
$$W_{p_1} = -1.0, \quad W_{p_2} = -1.0, \quad W_p = +1.0 \quad (\forall p \notin \{p_1, p_2\}).$$

### 2.2 Finite Zeeman Field ($h > 0$): Trapping vs. Healing
When $h \parallel [111] > 0$, the on-site Zeeman term $-2 h^\gamma N^\gamma$ and local Lagrange multipliers $\lambda_s^\gamma G^\gamma$ couple the gauge Majoranas $c^\gamma$ to the itinerant matter Majorana $c^0$. This induces virtual gauge fluctuations:
$$u_b(h) = \langle \mathrm{i} c_{s_A}^\gamma c_{s_B}^\gamma \rangle = -1 + \mathcal{O}(h^2).$$
- **Topological Trapping Regime ($0 \le h < h_{\text{heal}}$):** The energy barrier separating the $\pi$-flux sector from the flux-free vacuum is finite and positive. The unpinned self-consistent iteration converges to a metastable saddle/minimum with $W_{p_1}, W_{p_2} < 0$.
- **Dynamical Healing Regime ($h \ge h_{\text{heal}}$):** As the field approaches the polarized spinodal transition ($h \sim h_{c2}^+ \approx 0.40$), the topological barrier collapses, and the unpinned branch cut dynamically heals: $u_b \to 0$ or $+1$, driving all Wilson loops $W_p \to +1$.
