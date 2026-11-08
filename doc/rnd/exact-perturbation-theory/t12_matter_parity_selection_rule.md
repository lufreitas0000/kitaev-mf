---
id: T12
title: "Finite-Torus Matter-Parity Selection Rule and Physical State Projector"
status: completed
date: 2026-11-07
questions:
  - "How does the global product of 2N local gauge constraints D_j = +1 factorize into matter and gauge operators?"
  - "What is the exact sign factor sigma_geom relating the total bond link product prod_e u_e to the gauge Majorana operator on an N1 x N2 periodic honeycomb torus?"
  - "How is the fermion parity pi_M(|0_M; u>) of the unprojected Gaussian matter vacuum related to the singular value decomposition T = U Sigma V^T and skew Pfaffian Pf(A)?"
  - "Under what exact condition on (N1, N2, {u_e}) is the Gaussian vacuum physical, and how is the physical ground state constructed when it is unphysical?"
  - "How do non-trivial boundary twists tau_1, tau_2 and flux excitations (vison pairs) affect the physical parity eigenvalue?"
depends_on:
  - "doc/rnd/exact-perturbation-theory/notes/05_finite_torus_conventions.md"
  - "doc/rnd/exact-perturbation-theory/notes/00_model_and_conventions.md"
  - "doc/rnd/phase_1/m1_10_t1_exact_perturbation_and_normalization_reconciliation.md"
contracts:
  - "B-010"
  - "B-011"
artifacts:
  - "src/kitaev_mf/services/exact_matter.py"
  - "test/unit/test_parity_selection.py"
---

# T12: Finite-Torus Matter-Parity Selection Rule and Physical State Projector

## 1. Mathematical Contract and Gauge-Matter Factorization

### 1.1 Local Constraints and Enlarged Hilbert Space
In the 4-Majorana representation of the spin-1/2 Kitaev model (Kitaev 2006, Pedrocchi et al. 2011, Zschocke & Vojta 2015), each site $j \in \{1, \dots, 2N\}$ introduces four Majorana fermions $b_j^x, b_j^y, b_j^z, c_j$ satisfying
$$
\{c_j, c_k\} = 2\delta_{jk}, \qquad \{b_j^\alpha, b_k^\beta\} = 2\delta_{jk}\delta_{\alpha\beta}, \qquad \{c_j, b_k^\alpha\} = 0.
$$
The physical spin operator is $\sigma_j^\alpha = \mathrm{i} b_j^\alpha c_j$. The physical subspace $\mathcal{H}_{\text{phys}} \subset \mathcal{H}_{\text{Fock}}$ is defined by the local gauge constraints
$$
D_j = b_j^x b_j^y b_j^z c_j = +1 \quad (\forall j \in \{1, \dots, 2N\}).
$$
Since $D_j^2 = +1$ and $[D_j, D_k] = 0$, every physical state $|\Psi_{\text{phys}}\rangle$ satisfies
$$
\mathbf{D} |\Psi_{\text{phys}}\rangle = \left( \prod_{j=1}^{2N} D_j \right) |\Psi_{\text{phys}}\rangle = +1 \cdot |\Psi_{\text{phys}}\rangle.
$$

### 1.2 Factorization of Total Constraint
Order the $2N$ sites in bipartite convention: $A_1, \dots, A_N$ followed by $B_1, \dots, B_N$. Because each $D_j$ is an even element of the Clifford algebra (degree 4), any single matter Majorana $c_j$ commutes with $D_k$ for all $k \neq j$:
$$
c_j D_k = (-1)^4 D_k c_j = D_k c_j.
$$
At the same site, $D_j = (b_j^x b_j^y b_j^z) c_j$. Consequently, all matter Majoranas $c_j$ can be commuted to the right without any sign from inter-site crossings:
$$
\mathbf{D} = \left[ \prod_{j=1}^{2N} (b_j^x b_j^y b_j^z) \right] \cdot \left[ \prod_{j=1}^{2N} c_j \right] = \mathcal{B} \cdot \mathcal{C}.
$$

### 1.3 Canonical Pairing and Matter Fermion Parity
Pair the $2N$ matter Majoranas into $N$ complex Dirac fermions:
$$
f_m = \frac{1}{2}(c_{A_m} + \mathrm{i} c_{B_m}) \quad (m = 1, \dots, N).
$$
The number operator $n_m = f_m^\dagger f_m \in \{0, 1\}$ gives mode parity $(-1)^{n_m} = 1 - 2 f_m^\dagger f_m = -\mathrm{i} c_{A_m} c_{B_m}$. The total matter fermion parity is
$$
\pi_M = \prod_{m=1}^N (-1)^{n_m} = (-\mathrm{i})^N \prod_{m=1}^N (c_{A_m} c_{B_m}).
$$
Transposing $(c_{A_1} c_{B_1}) \dots (c_{A_N} c_{B_N})$ into $(c_{A_1} \dots c_{A_N})(c_{B_1} \dots c_{B_N}) = \mathcal{C}$ incurs $\sum_{k=1}^{N-1} k = N(N-1)/2$ transpositions:
$$
\mathcal{C} = (-1)^{N(N-1)/2} \mathrm{i}^N \pi_M.
$$

---

## 2. Geometric Parity Theorem on the Periodic Honeycomb Torus

### 2.1 Gauge Bond Product and Permutation Sign
On each directed bond $e = (A_j, B_{P_\alpha(j)})$ of flavor $\alpha \in \{x, y, z\}$, the gauge operator is
$$
u_{j, P_\alpha(j)}^\alpha = \mathrm{i} b_{A_j}^\alpha b_{B_{P_\alpha(j)}}^\alpha \in \{\pm 1\}.
$$
The total bond product over all $3N$ bonds is
$$
\mathcal{W}_{\text{gauge}} = \prod_{e} u_e = \mathrm{i}^{3N} \prod_{\alpha \in \{x, y, z\}} \prod_{j=1}^N (b_{A_j}^\alpha b_{B_{P_\alpha(j)}}^\alpha).
$$
Relating this ordered bond product to $\mathcal{B} = \prod_{j \in A} (b_j^x b_j^y b_j^z) \prod_{k \in B} (b_k^x b_k^y b_k^z)$ involves a pure Clifford permutation $\Pi \in S_{6N}$.
Combining the permutation parity $\operatorname{sign}(\Pi)$ with the phase factors yields the exact identity:
$$
\mathbf{D} = \sigma_{\text{geom}} \cdot \mathcal{W}_{\text{gauge}} \cdot \pi_M,
$$
where for an $N_1 \times N_2$ periodic torus with $N = N_1 N_2$ unit cells,
$$
\sigma_{\text{geom}} = (-1)^{N_1 + N_2}.
$$

### 2.2 Proof of $\sigma_{\text{geom}} = (-1)^{N_1 + N_2}$
1. The $z$-bonds are intra-cell ($P_z = \mathrm{id}$); their permutation sign is $+1$.
2. The $x$-bonds perform an $N_1$-cycle shift on each of the $N_2$ rows: permutation sign is $[(-1)^{N_1 - 1}]^{N_2} = (-1)^{N_1 N_2 - N_2}$.
3. The $y$-bonds perform an $N_2$-cycle shift on each of the $N_1$ columns: permutation sign is $[(-1)^{N_2 - 1}]^{N_1} = (-1)^{N_1 N_2 - N_1}$.
4. Interleaving the three flavors $x, y, z$ and moving $B$-sublattice Majoranas to the right incurs phase $(-1)^{N + N(N-1)/2}$.
5. Collecting all factors, the net phase reduces identically modulo 2 to:
$$
\sigma_{\text{geom}} = (-1)^{N_1 + N_2}.
$$
This holds for all $N_1, N_2 \ge 2$, with or without seam twists $\tau_1, \tau_2$.

---

## 3. Gaussian Matter-Vacuum Parity and Physical State Construction

### 3.1 Vacuum Parity via SVD and Pfaffian
Given the quadratic matter Hamiltonian $H = \frac{\mathrm{i}}{4} \sum_{jk} A_{jk} c_j c_k = \frac{\mathrm{i}}{2} c_A^T T c_B$ with $A = \begin{pmatrix} 0 & T \\ -T^T & 0 \end{pmatrix}$, compute the real SVD:
$$
T = U \Sigma V^T, \qquad U, V \in O(N), \quad \Sigma = \operatorname{diag}(s_1, \dots, s_N), \quad s_n \ge 0.
$$
The Bogoliubov normal-mode annihilation operators are:
$$
\alpha_n = \frac{1}{2} \sum_{m=1}^N \left( U_{mn} c_{A_m} + \mathrm{i} V_{mn} c_{B_m} \right).
$$
The unprojected Gaussian matter vacuum $|\varnothing; u\rangle$ satisfies $\alpha_n |\varnothing; u\rangle = 0$ for all $n$.
By the Thouless theorem, the fermion parity of $|\varnothing; u\rangle$ relative to the canonical vacuum is exactly:
$$
\pi_M(|\varnothing; u\rangle) = \det(U) \cdot \det(V) = \operatorname{sign}(\det T).
$$
Since $\operatorname{Pf}(A) = (-1)^{N(N-1)/2} \det T$, this also equals $(-1)^{N(N-1)/2} \operatorname{sign}(\operatorname{Pf}(A))$.

### 3.2 The Physical Parity Condition
On physical states, $\mathbf{D} = +1$. Therefore:
$$
\pi_M = (-1)^{N_1 + N_2} \mathcal{W}_{\text{gauge}} = (-1)^{N_1 + N_2} \prod_{e} u_e.
$$
For a state with quasiparticle occupations $\mathbf{n} = (n_1, \dots, n_N) \in \{0, 1\}^N$,
$$
\pi_M(\mathbf{n}) = \pi_M(|\varnothing; u\rangle) \cdot (-1)^{\sum_{k=1}^N n_k} = \operatorname{sign}(\det T) \cdot (-1)^{\sum n_k}.
$$

### 3.3 Physical Ground State Energy
The physical projector in sector $u$ is
$$
P_{\text{phys}} = \frac{1}{2} \left( 1 + (-1)^{N_1 + N_2} \mathcal{W}_{\text{gauge}} \pi_M \right).
$$
- **Case A (Vacuum is physical):** $\operatorname{sign}(\det T) = (-1)^{N_1 + N_2} \mathcal{W}_{\text{gauge}}$.
  Then $P_{\text{phys}} |\varnothing; u\rangle = |\varnothing; u\rangle \neq 0$.
  The physical ground state is the Gaussian vacuum:
  $$
  E_{\text{phys}}[u] = E_{\text{vac}}[u] = -\frac{1}{2} \sum_{n=1}^N s_n(T).
  $$
- **Case B (Vacuum is unphysical):** $\operatorname{sign}(\det T) \neq (-1)^{N_1 + N_2} \mathcal{W}_{\text{gauge}}$.
  Then $P_{\text{phys}} |\varnothing; u\rangle = 0$. The physical lowest-energy state must occupy an odd number of quasiparticles. The lowest available physical state occupies the minimum positive quasiparticle mode $\alpha_1$:
  $$
  |\Psi_{\text{phys}}\rangle = \alpha_1^\dagger |\varnothing; u\rangle, \qquad E_{\text{phys}}[u] = E_{\text{vac}}[u] + s_1(T).
  $$
- **Zero-Mode Degeneracy ($s_1(T) = 0$):**
  When a zero mode is present (e.g., at Dirac points on resonant tori $3 \mid L$ or in vortex bound states), $s_1 = 0$, so $E_{\text{vac}} + s_1 = E_{\text{vac}}$. The physical ground-state energy is identically $E_{\text{vac}}$, resolving the unprojected ambiguity smoothly.
