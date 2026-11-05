# Note 04f: Advanced Invariants, Semidefinite Constraints, and Verification Backlog

## 1. Context and Problem Formulation

In the pure Kitaev limit ($J = \Gamma = \Gamma' = D = \mathbf{h} = \boldsymbol{\lambda} = \mathbf{0}$), integrability and Lieb's theorem provide closed-form benchmarks (e.g., flat gauge bands, exact matter dispersion, $W_p \equiv +1.0$, and $\bar{U} \approx 0.5248657458506$) (Kitaev, 2006; Choi et al., 2018).

Away from integrability, self-consistent Majorana mean-field solutions are non-convex and seed-dependent. To construct rigorous, regression-resistant test suites in non-integrable regimes, tests cannot rely on unknown exact many-body ground states. Instead, tests must enforce **algebraic, kinematic, semidefinite, and variational invariants** that every valid self-consistent Gaussian state $\Xi = (\mathrm{U}, \mathrm{V})$ must satisfy identically.

---

## 2. Kinematic & Spectral Sum Rules

**Theorem (Parseval–Frobenius Second Spectral Moment)** <!-- claim: CLM-M1-INV-001 -->
Let $\mathcal{H}(\mathbf{k}) = \frac{i}{4}\mathcal{M}(\mathbf{k}) \in \mathbb{C}^{4N_{\mathrm{sub}} \times 4N_{\mathrm{sub}}}$ be the single-particle Bloch Hamiltonian with eigenvalues $\mathcal{E}_n(\mathbf{k})$. For arbitrary exchange couplings and magnetic fields $(J, K, \Gamma, \Gamma', D, \mathbf{h})$, on any regular Brillouin-zone grid $(N_1, N_2)$ satisfying the non-aliasing condition:
$$\Delta\mathbf{R}_b - \Delta\mathbf{R}_{b'} \not\equiv \mathbf{0} \pmod{(N_1, N_2)} \quad \forall b \neq b' \text{ with } (s_i(b), s_j(b)) = (s_i(b'), s_j(b')),$$
(satisfied in the honeycomb lattice when $N_1 \ge 2$ and $N_2 \ge 2$), the Brillouin-zone sum of squared eigenvalues satisfies the exact identity:
$$\frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k} \in \mathrm{BZ}} \sum_{n=1}^{4N_{\mathrm{sub}}} \mathcal{E}_n(\mathbf{k})^2 \equiv \frac{1}{16} \sum_{s=0}^{N_{\mathrm{sub}}-1} \|\mathrm{B}_s\|_{\mathrm{F}}^2 + \frac{1}{8} \sum_{b \in \mathcal{E}} \|\mathrm{A}_b\|_{\mathrm{F}}^2,$$
where $\|M\|_{\mathrm{F}}^2 = \operatorname{tr}(M M^T)$ denotes the squared Frobenius matrix norm.

*Proof.*
Because $\mathcal{H}(\mathbf{k})$ is Hermitian, $\sum_n \mathcal{E}_n(\mathbf{k})^2 = \operatorname{tr}(\mathcal{H}(\mathbf{k})^2) = \frac{1}{16} \|\mathcal{M}(\mathbf{k})\|_{\mathrm{F}}^2$.
In terms of on-site blocks $\mathrm{B}_s$ and directed hopping blocks $\mathrm{A}_b$ with bond phase $\phi_b(\mathbf{k}) = 2\pi(k_1 n_{1,b} + k_2 n_{2,b})$:
$$\|\mathcal{M}(\mathbf{k})\|_{\mathrm{F}}^2 = \sum_{s=0}^{N_{\mathrm{sub}}-1} \|\mathrm{B}_s\|_{\mathrm{F}}^2 + \sum_{b \in \mathcal{E}} \left( \|\mathrm{A}_b e^{i \phi_b(\mathbf{k})}\|_{\mathrm{F}}^2 + \|-\mathrm{A}_b^T e^{-i \phi_b(\mathbf{k})}\|_{\mathrm{F}}^2 \right) + \sum_{b \neq b'} (\text{cross-phase terms}).$$
Because $\|\mathrm{A}_b e^{i \phi_b}\|_{\mathrm{F}}^2 = \operatorname{tr}(\mathrm{A}_b \mathrm{A}_b^T) = \|\mathrm{A}_b\|_{\mathrm{F}}^2$ and $\|-\mathrm{A}_b^T\|_{\mathrm{F}}^2 = \|\mathrm{A}_b\|_{\mathrm{F}}^2$, the diagonal bond contribution is $2 \|\mathrm{A}_b\|_{\mathrm{F}}^2$.
Integrating over the first Brillouin zone, off-diagonal cross-phase terms between distinct bonds sharing identical sublattice endpoints evaluate to $\frac{1}{N_{\mathbf{k}}}\sum_{\mathbf{k}} e^{i(\phi_b(\mathbf{k}) - \phi_{b'}(\mathbf{k}))} = \delta_{\Delta\mathbf{R}_b, \Delta\mathbf{R}_{b'} \bmod (N_1, N_2)}$. Under the non-aliasing condition $N_1 \ge 2, N_2 \ge 2$, all cross-phase terms vanish identically.
Dividing by $16$ yields $\frac{1}{16}\sum_s \|\mathrm{B}_s\|_{\mathrm{F}}^2 + \frac{2}{16}\sum_b \|\mathrm{A}_b\|_{\mathrm{F}}^2 = \frac{1}{16}\sum_s \|\mathrm{B}_s\|_{\mathrm{F}}^2 + \frac{1}{8}\sum_b \|\mathrm{A}_b\|_{\mathrm{F}}^2$. $\square$

**Proposition (Pointwise Projector Commutation and Idempotency)** <!-- claim: CLM-M1-INV-002 -->
At every discrete wavevector $\mathbf{k} \in \mathrm{BZ}$, the spectral projector $P(\mathbf{k}) = 2\Pi^-(\mathbf{k})$ satisfies:
1. *Pointwise Commutation:*
   $$[\mathcal{H}(\mathbf{k}), P(\mathbf{k})] = \mathcal{H}(\mathbf{k}) P(\mathbf{k}) - P(\mathbf{k}) \mathcal{H}(\mathbf{k}) \equiv \mathbf{0}.$$
2. *Particle-Hole Filling & Idempotency:*
   $$\frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k} \in \mathrm{BZ}} \operatorname{tr}(P(\mathbf{k})) \equiv 4 N_{\mathrm{sub}}, \qquad P(\mathbf{k})^2 \equiv 2 P(\mathbf{k}) \quad (\text{pointwise at gapped } \mathbf{k}).$$
   At wavevectors hosting exact zero modes, the half-weight convention satisfies $0 \le P(\mathbf{k}) \le 2\mathbb{1}_{4N_{\mathrm{sub}}}$.
3. *Reflection Sum Rule:*
   $$P(-\mathbf{k}) + P^T(\mathbf{k}) \equiv 2\,\mathbb{1}_{4N_{\mathrm{sub}}}.$$

---

## 3. Semidefinite & Quantum $N$-Representability Constraints

**Theorem (Local Density Matrix Positivity LMI)** <!-- claim: CLM-M1-INV-003 -->
Let $b = \langle ij \rangle$ be a directed bond connecting sublattice $i$ to sublattice $j$. For any quantum state over the self-dual CAR algebra, the local $8 \times 8$ two-point correlation matrix $\mathcal{C}_{ij} = \langle \mathbf{c}_{ij} \mathbf{c}_{ij}^T \rangle$ with $\mathbf{c}_{ij} = (c_i^0, \dots, c_i^3, c_j^0, \dots, c_j^3)^T$ is given by:
$$\mathcal{C}_{ij} = \begin{pmatrix} \mathbb{1}_4 - i \mathrm{V}_i & -i \mathrm{U}_b \\ i \mathrm{U}_b^T & \mathbb{1}_4 - i \mathrm{V}_j \end{pmatrix} \in \mathbb{C}^{8 \times 8}.$$
$\mathcal{C}_{ij}$ is Hermitian and must be positive semi-definite:
$$\lambda_{\min}(\mathcal{C}_{ij}) \ge 0.$$

*Proof.*
For any vector $\mathbf{z} \in \mathbb{C}^8$, define the linear operator $\psi_{\mathbf{z}} = \sum_{\mu=0}^7 z_\mu c^\mu$.
By the $C^*$-algebra positivity axiom, $\langle \psi_{\mathbf{z}} \psi_{\mathbf{z}}^\dagger \rangle \ge 0$.
Expanding the expectation value using $\langle c_s^\mu c_s^\nu \rangle = \delta^{\mu\nu} - i \mathrm{V}_s^{\mu\nu}$ and $\langle c_i^\mu c_j^\nu \rangle = -i \mathrm{U}_b^{\mu\nu}$ yields $\mathbf{z}^\dagger \mathcal{C}_{ij} \mathbf{z} \ge 0$, which requires all eigenvalues of $\mathcal{C}_{ij}$ to be non-negative. $\square$

**Corollary (Singular Value and Envelope Bounds)** <!-- claim: CLM-M1-INV-004 -->
From the condition $-\mathbb{1} \le \Omega^2 \le 0$ on the global covariance matrix:
1. The spectral norm of every bond hopping matrix satisfies:
   $$\sigma_{\max}(\mathrm{U}_b) \le 1.0 \quad \forall b \in \mathcal{E}.$$
2. The singular values of the real skew-symmetric on-site blocks satisfy:
   $$\max_k |\nu_k(\mathrm{V}_s)| \le 1.0 \quad \forall s \in \mathcal{V}.$$
3. The physical on-site magnetization vector $\langle \boldsymbol{\sigma}_s \rangle_{\mathrm{G}} = (V_s^{01}, V_s^{02}, V_s^{03})^T$ is bounded by the Bloch ball:
   $$\|\langle \boldsymbol{\sigma}_s \rangle_{\mathrm{G}}\|_2 \le 1.0.$$
4. The local fermion parity $\langle \mathcal{D}_s \rangle_{\mathrm{G}} = \operatorname{Pf}(\mathrm{V}_s)$ is bounded by:
   $$-1.0 \le \mathrm{V}_s^{01}\mathrm{V}_s^{23} - \mathrm{V}_s^{02}\mathrm{V}_s^{13} + \mathrm{V}_s^{03}\mathrm{V}_s^{12} \le +1.0.$$

---

## 4. Solvable Limits and Dimensional Reductions

**Theorem (Decoupled 1D Kitaev Chain Limit)** <!-- claim: CLM-M1-INV-005 -->
In the 2-sublattice honeycomb cell with isotropic $K_x = K_y = K < 0$ and decoupled $z$-bonds ($K_z = 0$, $J = \Gamma = \Gamma' = D = \mathbf{h} = \mathbf{0}$):
1. The 2D lattice decomposes into decoupled 1D zigzag chains.
2. The self-consistent matter bond correlator on each active bond ($x$ and $y$) evaluates to the exact closed-form integral:
   $$U_x^{00} = U_y^{00} = -\bar{U}_{\mathrm{1D}} = -\frac{1}{2} \int_0^1 \mathrm{d}\theta \, \left| 2\cos(\pi\theta) \right| = -\frac{2}{\pi} \approx -0.636619772368.$$
3. The $z$-bond correlator vanishes identically: $\mathrm{U}_z \equiv \mathbf{0}$.

*Proof.*
With $K_z = 0$, the two-band dispersion form factor is $f(\mathbf{k}) = e^{-i 2\pi k_1} + e^{-i 2\pi k_2} = e^{-i \pi(k_1+k_2)} 2\cos(\pi(k_1-k_2))$, with $|f(\mathbf{k})| = 2 |\cos(\pi(k_1 - k_2))|$.
Because the unit cell contains only two active bonds sharing the matter bandwidth, the sum of hoppings satisfies $U_x^{00} + U_y^{00} = -\int_0^1 \mathrm{d}\theta \, |f(\theta)| = -4/\pi$. By $x \leftrightarrow y$ reflection symmetry, each active bond carries half: $U_x^{00} = U_y^{00} = -2/\pi$. $\square$

**Proposition (Isolated Dimer Flat-Band Limit)** <!-- claim: CLM-M1-INV-006 -->
In the limit $|K_z| \gg |K_x|, |K_y|$ ($K_x = K_y = 0$, $K_z = -1.0$):
1. Matter flavor $c^0$ and gauge flavor $c^3$ couple exclusively along $z$-bonds, forming flat dimer bands at $\mathcal{E}(\mathbf{k}) \equiv \pm \frac{|K_z|}{2} = \pm 0.5$. Uncoupled orthogonal gauge flavors $c^1, c^2$ remain at zero energy ($\mathcal{E} \equiv 0$). The four non-negative bands are $[0, 0, 0.5, 0.5]$.
2. The matter correlator saturates to unity on the dimer: $\mathrm{U}_z^{00} = -1.0$, $\mathrm{U}_z^{33} = +1.0$, while $\mathrm{U}_x = \mathrm{U}_y \equiv \mathbf{0}$.
3. The ground-state energy density evaluates to $E_{\mathrm{cell}} = -|K_z| = -1.0$ ($E_{\mathrm{site}} = -0.5$).

---

## 5. Symmetry Selection Rules & Variational Response

**Proposition (Time-Reversal Odd Correlator Vanishing)** <!-- claim: CLM-M1-INV-007 -->
Under time-reversal symmetry $\Theta$, matter fermions are even ($\Theta c^0 \Theta^{-1} = c^0$) and gauge fermions are odd ($\Theta c^\gamma \Theta^{-1} = -c^\gamma$).
In any time-reversal symmetric state ($\mathbf{h} = \mathbf{0}, D = 0$):
$$\mathrm{U}_b^{0\gamma} = \langle i c_{i,0}^0 c_{j,R_b}^\gamma \rangle \equiv 0, \qquad \mathrm{U}_b^{\gamma 0} = \langle i c_{i,0}^\gamma c_{j,R_b}^0 \rangle \equiv 0 \quad \forall \gamma \in \{1, 2, 3\}.$$

**Theorem (Variational Hellmann–Feynman Derivative)** <!-- claim: CLM-M1-INV-008 -->
Let $\Xi^*(\mathbf{h})$ be a self-consistent stationary point ($\nabla_\Xi E[\Xi^*] = 0$). The total derivative of the ground-state energy with respect to external magnetic field components satisfies:
$$\frac{\mathrm{d} E_{\mathrm{cell}}}{\mathrm{d} h^\gamma} \equiv -\sum_{s=0}^{N_{\mathrm{sub}}-1} \mathrm{V}_s^{0\gamma} = -\sum_{s=0}^{N_{\mathrm{sub}}-1} \langle \sigma_s^\gamma \rangle_{\mathrm{G}}.$$
Furthermore, the diagonal magnetic susceptibility $\chi^{\gamma\gamma} = -\frac{\partial^2 E_{\mathrm{cell}}}{\partial (h^\gamma)^2} = \frac{\partial \langle \sigma_{\mathrm{tot}}^\gamma \rangle}{\partial h^\gamma} \ge 0$ must be non-negative.

---

## 6. Microscopic Dualities

**Theorem (Four-Sublattice Klein Duality in the 4-Site Cell)** <!-- claim: CLM-M1-INV-009 -->
Consider the 4-sublattice rectangular cell with sublattices $(A_1, B_1, A_2, B_2)$. Define the four-sublattice unitary $\pi$-rotations around spin axes:
$$\mathcal{K}: \quad \boldsymbol{\sigma}_{A_1} \to (+\sigma^x, +\sigma^y, +\sigma^z), \quad \boldsymbol{\sigma}_{B_1} \to (+\sigma^x, -\sigma^y, -\sigma^z),$$
$$\boldsymbol{\sigma}_{A_2} \to (-\sigma^x, +\sigma^y, -\sigma^z), \quad \boldsymbol{\sigma}_{B_2} \to (-\sigma^x, -\sigma^y, +\sigma^z).$$
This transformation maps the extended model parameters according to:
$$\tilde{J} = -J, \qquad \tilde{K} = K + 2J, \qquad \tilde{\Gamma} = -\Gamma, \qquad \tilde{\Gamma}' = \Gamma'.$$
At $\mathbf{h} = \mathbf{0}$, the ground-state energy density is invariant under the duality:
$$E_{\mathrm{cell}}(J, K, \Gamma, \Gamma') \equiv E_{\mathrm{cell}}(-J, K+2J, -\Gamma, \Gamma').$$

---

## References

* **Choi, W., Klein, P. W., Rosch, A., and Kim, Y. B.** (2018). *Topological superconductivity in Kondo-Kitaev model*. Physical Review B **98**, 155123. DOI: [10.1103/PhysRevB.98.155123](https://doi.org/10.1103/PhysRevB.98.155123). arXiv: [1804.10212](https://arxiv.org/abs/1804.10212).
* **Freitas, L. R. D., Bauer, T., Egger, R., and Pereira, R. G.** (2024). *Electric polarization near vortices in the extended Kitaev model*. npj Quantum Materials **9**, 33. DOI: [10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5). arXiv: [2311.15750](https://arxiv.org/abs/2311.15750).
* **Kimchi, I. and Vishwanath, A.** (2014). *Kitaev-Heisenberg models for iridates on the triangular, hyperhoneycomb, kagome, fcc, and pyrochlore lattices*. Physical Review B **89**, 014414. DOI: [10.1103/PhysRevB.89.014414](https://doi.org/10.1103/PhysRevB.89.014414). arXiv: [1302.2639](https://arxiv.org/abs/1302.2639).
* **Kitaev, A.** (2006). *Anyons in an exactly solved model and beyond*. Annals of Physics **321**, 2–111. DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005). arXiv: [cond-mat/0506438](https://arxiv.org/abs/cond-mat/0506438).
