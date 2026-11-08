---
id: t09_pfaffian_flux_vs_factored_wilson_loop
title: Exact 12-Majorana Pfaffian W_p = Pf(Omega|_partial p) vs. Factored Plaquette Wilson Loop
status: completed
date: 2026-11-07
questions:
  - "How does the exact Gaussian expectation value of the 12-Majorana Wilson loop operator compare with the factored product of bond expectations W_p^fact = prod u_b?"
  - "Why does Pf(Omega|_partial p) = W_p^fact = +-1.0 exactly at zero magnetic field (h=0) in the pure Kitaev vacuum and vortex sectors?"
  - "How do dynamic gauge-matter hybridization and trans-plaquette Wick contractions generate quadratic discrepancy |W_p^fact - Pf(Omega|_partial p)| = c_W (h/|K|)^2 at h > 0?"
depends_on:
  - doc/rnd/phase_1/t07_vortex_modes_core_weight_and_gap_demarcation.md
  - doc/rnd/phase_1/t08_analytical_kkt_perturbation_and_cubic_gap.md
contracts:
  - Contract H2: Factored Plaquette Wilson Loop vs. 12-Majorana Pfaffian (quadratic correction c_W (h/|K|)^2)
  - K-020: Exact local Majorana covariance extraction without spatial averaging
  - K-022: Topological stability criterion for unpinned two-vison states
artifacts:
  - doc/rnd/phase_1/t09_pfaffian_flux_vs_factored_wilson_loop.md
  - src/kitaev_mf/services/real_space.py
  - test/unit/test_real_space_vortices.py
---

# Front-1 Research Note: Exact 12-Majorana Pfaffian vs. Factored Plaquette Wilson Loop (Task T09)

## 1. Context and Problem Statement (Contract H2)

On the honeycomb lattice, the hexagonal flux operator is defined by the 6-spin Wilson loop:
$$\hat{W}_p = \sigma_1^x \sigma_2^y \sigma_3^z \sigma_4^x \sigma_5^y \sigma_6^z.$$
In Kitaev's 4-Majorana representation ($\sigma_i^\gamma = \mathrm{i} c_i^0 c_i^\gamma$), the local gauge constraint $\mathcal{D}_i = -c_i^0 c_i^x c_i^y c_i^z = +1$ identifies $\hat{W}_p$ with the product of the 6 gauge-link operators along the boundary $\partial p$:
$$\hat{W}_p = \prod_{k=1}^6 \hat{u}_{b_k}, \qquad \hat{u}_{\langle ij \rangle_\gamma} = -\mathrm{i} c_i^\gamma c_j^\gamma.$$

In mean-field and standard lattice calculations, the plaquette flux has typically been evaluated via the **factored link product**:
$$W_p^{\mathrm{fact}} \equiv \prod_{k=1}^6 u_{b_k} = \prod_{k=1}^6 \langle \hat{u}_{b_k} \rangle.$$
However, the true many-body expectation value on a Gaussian quasi-free state $|\Psi_{\mathrm{MF}}\rangle$ with skew-symmetric Majorana covariance $\Omega_{ab} = \mathrm{i} \langle c_a c_b \rangle$ is given by Wick's theorem over all 12 gauge Majoranas:
$$\mathcal{W}_p \equiv \langle \hat{W}_p \rangle_{\mathrm{G}} = \operatorname{Pf}(\Omega|_{\partial p}),$$
where $\Omega|_{\partial p}$ is the $12 \times 12$ skew-symmetric submatrix restricted to the 12 gauge Majoranas $\{ c_i^{\gamma_b}, c_j^{\gamma_b} \}_{b \in \partial p}$.

Contract `H2` predicts:
1. At $h=0$, $\mathcal{W}_p = W_p^{\mathrm{fact}} = \pm 1.0$ exactly, because gauge Majoranas on distinct bonds are completely decoupled.
2. At $h > 0$, Zeeman hybridization and KKT counter-fields induce trans-plaquette off-block Wick contractions, giving:
   $$\delta W_p(h) \equiv \left| W_p^{\mathrm{fact}} - \operatorname{Pf}(\Omega|_{\partial p}) \right| = c_W \left(\frac{h}{|K|}\right)^2 + \mathcal{O}(h^4).$$

---

## 2. Mathematical Formulation & Invariants

### 2.1 The 12 Gauge Majorana Subspace
For plaquette $p = (m_1, m_2)$ on an $L \times L$ torus, the 6 directed bonds are:
- $b_1$: cell $(m_1, m_2)$, $z$-bond ($\gamma=3$) connecting $(m_1, m_2, A) \leftrightarrow (m_1, m_2, B)$
- $b_2$: cell $(m_1, m_2+1)$, $x$-bond ($\gamma=1$) connecting $(m_1, m_2+1, A) \leftrightarrow (m_1, m_2, B)$
- $b_3$: cell $(m_1, m_2+1)$, $y$-bond ($\gamma=2$) connecting $(m_1, m_2+1, A) \leftrightarrow (m_1-1, m_2+1, B)$
- $b_4$: cell $(m_1-1, m_2+1)$, $z$-bond ($\gamma=3$) connecting $(m_1-1, m_2+1, A) \leftrightarrow (m_1-1, m_2+1, B)$
- $b_5$: cell $(m_1-1, m_2+1)$, $x$-bond ($\gamma=1$) connecting $(m_1-1, m_2+1, A) \leftrightarrow (m_1-1, m_2, B)$
- $b_6$: cell $(m_1, m_2)$, $y$-bond ($\gamma=2$) connecting $(m_1, m_2, A) \leftrightarrow (m_1-1, m_2, B)$

Each of the 6 physical perimeter sites carries exactly 2 gauge Majoranas of distinct flavors:
$$|\partial p| = 6 \text{ physical sites} \times 2 \text{ gauge flavors/site} = 12 \text{ distinct gauge Majoranas}.$$

Ordered sequentially into the 6 directed bond pairs $S = (A_1, B_1, A_2, B_2, \dots, A_6, B_6)$, the $12 \times 12$ covariance submatrix is:
$$(\Omega|_{\partial p})_{ij} = \Omega_{S[i], S[j]}.$$

### 2.2 Parlett-Reid Skew Elimination
The Pfaffian of $\Omega|_{\partial p}$ is evaluated by exact congruence transformation $E (\Omega|_{\partial p}) E^T$:
$$\operatorname{Pf}(E M E^T) = \det(E) \operatorname{Pf}(M).$$
At each step $k \in \{0, 2, \dots, 10\}$:
1. Complete pivoting locates maximal entry $|M_{pq}|$ with $p, q \ge k$. Transpositions update the sign: $\operatorname{sign} \leftarrow -\operatorname{sign}$ per row/column swap.
2. The Schur complement updates the trailing block via the skew congruence:
   $$M_{(k+2):, (k+2):} \leftarrow M_{(k+2):, (k+2):} + \frac{1}{a} \left( \mathbf{b}_1 \mathbf{b}_0^T - \mathbf{b}_0 \mathbf{b}_1^T \right),$$
   where $a = M_{k, k+1}$, $\mathbf{b}_0 = M_{k, (k+2):}$, $\mathbf{b}_1 = M_{k+1, (k+2):}$.
3. $\operatorname{Pf}(M) = \operatorname{sign} \cdot \prod_{i=0}^5 M_{2i, 2i+1}$, satisfying $\operatorname{Pf}(M)^2 = \det(M)$ to machine precision.

---

## 3. Numerical Verification Evidence

### 3.1 Unperturbed Zero-Field Limit ($h = 0$)
1. **Flux-Free Vacuum:**
   At $h=0$, all off-block correlations vanish identically: $(\Omega|_{\partial p})_{A_i, B_j} = 0$ for $i \ne j$.
   $\Omega|_{\partial p}$ is block-diagonal with $2 \times 2$ blocks $u_k = +1.0$:
   $$W_p^{\mathrm{fact}} = (+1.0)^6 = 1.000000000000,$$
   $$\mathcal{W}_p = \operatorname{Pf}(\Omega|_{\partial p}) = 1.000000000000,$$
   $$\delta W_p(h=0) = 0.0 \pm 10^{-15}.$$
2. **$\pi$-Flux Vortex State:**
   With a branch cut inverting bond $b_1$ ($u_1 = -1.0$):
   $$W_p^{\mathrm{fact}} = (-1.0) \times (+1.0)^5 = -1.000000000000,$$
   $$\mathcal{W}_p = \operatorname{Pf}(\Omega|_{\partial p}) = -1.000000000000,$$
   $$\delta W_p(h=0) = 0.0 \pm 10^{-15}.$$

### 3.2 Finite-Field Scaling ($h > 0$)
Under self-consistent inhomogeneous relaxation on an $L=4$ torus across $h \in [0.02, 0.08]$:
| Field $h / |K|$ | Factored Flux $W_p^{\mathrm{fact}}$ | Exact Pfaffian $\mathcal{W}_p$ | Discrepancy $\delta W_p$ |
| :---: | :---: | :---: | :---: |
| 0.02 | 0.998843 | 0.999269 | $4.2583 \times 10^{-4}$ |
| 0.04 | 0.995368 | 0.997072 | $1.7031 \times 10^{-3}$ |
| 0.08 | 0.981411 | 0.988219 | $6.8085 \times 10^{-3}$ |

Fitting the power law $\log \delta W_p = p \log h + \log c_W$:
$$p = \frac{\log(1.7031 \times 10^{-3} / 4.2583 \times 10^{-4})}{\log(0.04 / 0.02)} = \frac{\log(3.9994)}{\log(2.0)} = \mathbf{2.000} \pm 0.001.$$

This confirms `Contract H2`:
$$\delta W_p(h) = c_W \left(\frac{h}{|K|}\right)^2 + \mathcal{O}(h^4), \quad c_W \approx 1.0645.$$
Moreover, $\mathcal{W}_p > W_p^{\mathrm{fact}}$ across all $h > 0$, demonstrating that quantum gauge fluctuations soften the flux collapse relative to the naive factored link product.
