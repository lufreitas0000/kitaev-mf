# Note 01: Hartree–Fock–Bogoliubov Decoupling, Covariance Manifolds, and Variational Mean-Field Symmetries

## 1. Context and Quasi-Free State Manifold

Following Note 00, the physical spin Hamiltonian is mapped onto an exact quartic Majorana Hamiltonian $H_\lambda$ on the CAR algebra $\mathcal{A}_{\mathrm{CAR}}$ over graph $\Lambda = (\mathcal{V}, \mathcal{E})$ with local Lagrange multipliers $\boldsymbol{\lambda}_i \in \mathbb{R}^3$ (Freitas et al., 2024).

Because finding exact ground states of quartic interacting fermionic systems is exponentially hard, mean-field theory restricts the variational search to the manifold of **quasi-free (Gaussian) states** $\mathcal{S}_{\mathrm{QF}}$ over $\mathcal{A}_{\mathrm{CAR}}$. By the generalized Wick theorem, any quasi-free state $\rho \in \mathcal{S}_{\mathrm{QF}}$ is uniquely and completely characterized by its real two-point covariance matrices.

**Definition (Mean-Field Covariance Configuration $\Xi$)** <!-- claim: CLM-M1-MF-001 -->
For any quasi-free state $\rho \in \mathcal{S}_{\mathrm{QF}}$, we define the mean-field covariance data $\Xi = (\mathrm{U}, \mathrm{V})$ consisting of real matrix blocks:
$$\mathrm{U}_{ij}^{\mu\nu} \equiv \langle i c_i^\mu c_j^\nu \rangle_\rho \quad (i \neq j), \qquad \mathrm{V}_i^{\mu\nu} \equiv \langle i c_i^\mu c_i^\nu \rangle_\rho - i \delta^{\mu\nu}.$$
The matrices satisfy the exact algebraic symmetry relations:
$$\mathrm{U}_{ji} = -\mathrm{U}_{ij}^T, \qquad \mathrm{V}_i^T = -\mathrm{V}_i \in \mathfrak{so}(4).$$

*Remark (On-Site Subtraction Convention):*
In Freitas et al. (2024), the on-site correlator was written as $\langle i c_i^\mu c_i^\nu \rangle$, which contains an imaginary diagonal $+i\delta^{\mu\nu}$ originating from the Clifford algebra $\{c_i^\mu, c_i^\nu\} = 2\delta^{\mu\nu}\mathbb{1}$. By explicitly defining $\mathrm{V}_i^{\mu\nu} = \langle i c_i^\mu c_i^\nu \rangle - i\delta^{\mu\nu}$, $\mathrm{V}_i$ is strictly real and skew-symmetric ($\mathrm{V}_i \in \mathfrak{so}(4)$). Because the Kitaev matrices $\mathrm{N}^\gamma$ are traceless,
$$\operatorname{tr}\!\left( \mathrm{N}^{\gamma\,T} (\mathrm{V}_i + i\lambda \mathbb{1}) \right) = \operatorname{tr}\!\left( \mathrm{N}^{\gamma\,T} \mathrm{V}_i \right),$$
this subtraction alters neither the Hartree exchange fields nor the energy shifts, while ensuring all arrays stored in memory are real float64 tensors.

---

## 2. Generalized Wick Contraction and Energy Functional

**Lemma (Wick Contraction for Quadratic Bilinears)** <!-- claim: CLM-M1-MF-002 -->
Let $Q_i^a \equiv i c_i^T \mathrm{N}^a c_i = 2\sigma_i^a$ denote the local spin bilinears. On any Gaussian state $\rho \in \mathcal{S}_{\mathrm{QF}}$, the expectation values evaluate identically to:
$$\langle Q_i^a \rangle_\rho = \operatorname{tr}\!\left( \mathrm{V}_i^T \mathrm{N}^a \right) = 2 \mathrm{V}_i^{0a},$$
$$\langle Q_i^a Q_j^b \rangle_\rho = \langle Q_i^a \rangle_\rho \langle Q_j^b \rangle_\rho + 2 \operatorname{tr}\!\left( \mathrm{U}_{ij}^T \mathrm{N}^a \mathrm{U}_{ij} \mathrm{N}^b \right) \quad (i \neq j).$$

*Proof.*
The single-site expectation value follows from $\langle i c_i^\mu c_i^\nu \rangle = \mathrm{V}_i^{\mu\nu} + i\delta^{\mu\nu}$ and $(\mathrm{N}^a)_{\mu\nu} = \delta_{\mu 0}\delta_{\nu a} - \delta_{\nu 0}\delta_{\mu a}$:
$$\langle Q_i^a \rangle = \sum_{\mu,\nu} (\mathrm{N}^a)_{\mu\nu} \langle i c_i^\mu c_i^\nu \rangle = \operatorname{tr}\!\left( (\mathrm{N}^a)^T \mathrm{V}_i \right) = \operatorname{tr}\!\left( \mathrm{V}_i^T \mathrm{N}^a \right) = 2\mathrm{V}_i^{0a}.$$
For the inter-site two-body term, Wick's theorem for four Majorana operators decomposes the product into direct and exchange pairings:
$$\langle (c_i^T \mathrm{N}^a c_i)(c_j^T \mathrm{N}^b c_j) \rangle = \langle c_i^T \mathrm{N}^a c_i \rangle \langle c_j^T \mathrm{N}^b c_j \rangle - 2 \sum_{\mu,\nu,\rho,\sigma} (\mathrm{N}^a)_{\mu\nu} (\mathrm{N}^b)_{\rho\sigma} \langle c_i^\mu c_j^\rho \rangle \langle c_i^\nu c_j^\sigma \rangle.$$
Substituting $\langle c_i^\mu c_j^\rho \rangle = -i \mathrm{U}_{ij}^{\mu\rho}$ gives $(-i)^2 = -1$, canceling the negative sign and yielding $+2 \operatorname{tr}(\mathrm{U}_{ij}^T \mathrm{N}^a \mathrm{U}_{ij} \mathrm{N}^b)$. $\square$

**Definition (Hartree–Fock Energy Functional)** <!-- claim: CLM-M1-MF-003 -->
Evaluating the expectation value $\langle H_\lambda \rangle_\rho$ of the quartic Hamiltonian from Note 00 over a quasi-free state with covariance data $\Xi = (\mathrm{U}, \mathrm{V})$ yields the variational energy functional $E[\Xi, \boldsymbol{\lambda}]$:
$$E[\Xi, \boldsymbol{\lambda}] \equiv \frac{1}{8} \sum_{i \in \mathcal{V}} \sum_{j \in \mathcal{V}_i} \sum_{a,b=1}^3 \mathbf{J}_{ij}^{ab} \left[ \operatorname{tr}\!\left(\mathrm{V}_i^T \mathrm{N}^a\right) \operatorname{tr}\!\left(\mathrm{V}_j^T \mathrm{N}^b\right) + 2\operatorname{tr}\!\left(\mathrm{U}_{ij}^T \mathrm{N}^a \mathrm{U}_{ij} \mathrm{N}^b\right) \right] - \frac{1}{2} \sum_{i \in \mathcal{V}} \sum_{\gamma=1}^3 \left( h^\gamma \operatorname{tr}\!\left(\mathrm{V}_i^T \mathrm{N}^\gamma\right) - \frac{1}{2}\lambda_i^\gamma \operatorname{tr}\!\left((\mathrm{G}^\gamma)^T \mathrm{V}_i\right) \right).$$

---

## 3. Quadratic Mean-Field Hamiltonian and Tangent Property

**Proposition (Majorana HFB Hamiltonian Matrix)** <!-- claim: CLM-M1-MF-004 -->
For any given covariance data $\Xi = (\mathrm{U}, \mathrm{V})$, the effective quadratic mean-field Hamiltonian is (Freitas et al., 2024):
$$H_{\mathrm{MF}}[\Xi] = \sum_{i \in \mathcal{V}} \sum_{j \in \mathcal{V}_i} \frac{i}{4} c_i^T \mathrm{A}_{ij} c_j + \sum_{i \in \mathcal{V}} \frac{i}{4} c_i^T \mathrm{B}_i c_i - \mathrm{C},$$
where the real inter-site hopping matrices $\mathrm{A}_{ij}$, real skew-symmetric on-site matrices $\mathrm{B}_i$, and scalar energy shift $\mathrm{C}$ are:
$$\mathrm{A}_{ij} = 2 \sum_{a,b=1}^3 \mathbf{J}_{ij}^{ab} \mathrm{N}^a \mathrm{U}_{ij} \mathrm{N}^b,$$
$$\mathrm{B}_i = \sum_{j \in \mathcal{V}_i} \sum_{a,b=1}^3 \mathbf{J}_{ij}^{ab} \mathrm{N}^a \operatorname{tr}\!\left(\mathrm{V}_j^T \mathrm{N}^b\right) + \sum_{\gamma=1}^3 \left( \lambda_i^\gamma \mathrm{G}^\gamma - 2 h^\gamma \mathrm{N}^\gamma \right),$$
$$\mathrm{C} = \frac{1}{8} \sum_{i \in \mathcal{V}} \sum_{j \in \mathcal{V}_i} \sum_{a,b=1}^3 \mathbf{J}_{ij}^{ab} \left[ \operatorname{tr}\!\left(\mathrm{V}_i^T \mathrm{N}^a\right) \operatorname{tr}\!\left(\mathrm{V}_j^T \mathrm{N}^b\right) + 2\operatorname{tr}\!\left(\mathrm{U}_{ij}^T \mathrm{N}^a \mathrm{U}_{ij} \mathrm{N}^b\right) \right].$$

**Theorem (Properties of the Mean-Field Decoupling)** <!-- claim: CLM-M1-MF-005 -->
1. *Hermiticity & Skew Symmetries:* $\mathrm{A}_{ji} = -\mathrm{A}_{ij}^T$ and $\mathrm{B}_i^T = -\mathrm{B}_i$, ensuring $H_{\mathrm{MF}}$ is Hermitian.
2. *Expectation Value Consistency:* For any state $\rho$ whose covariance matches $\Xi$, $\langle H_{\mathrm{MF}}[\Xi] \rangle_\rho = E[\Xi, \boldsymbol{\lambda}]$, proving that $\mathrm{C}$ is the exact double-counting scalar subtraction.
3. *Tangent Property (Variational Gradient):* For any arbitrary trial configuration $\Xi' = (\mathrm{U}', \mathrm{V}')$,
   $$\langle H_{\mathrm{MF}}[\Xi] \rangle_{\Xi'} = E[\Xi, \boldsymbol{\lambda}] + \mathrm{d}E_\Xi[\Xi' - \Xi].$$
   Thus, $H_{\mathrm{MF}}[\Xi]$ is the exact affine first-order Taylor expansion of the many-body energy functional $E$ about $\Xi$.

*Proof.*
(1) Using $(\mathrm{N}^a)^T = -\mathrm{N}^a$ and $\mathbf{J}_{ji} = \mathbf{J}_{ij}^T$:
$$\mathrm{A}_{ji} = 2 \sum_{a,b} \mathbf{J}_{ji}^{ab} \mathrm{N}^a \mathrm{U}_{ji} \mathrm{N}^b = 2 \sum_{a,b} \mathbf{J}_{ij}^{ba} \mathrm{N}^a (-\mathrm{U}_{ij}^T) \mathrm{N}^b = -2 \sum_{a,b} \mathbf{J}_{ij}^{ba} (\mathrm{N}^b \mathrm{U}_{ij} \mathrm{N}^a)^T = -\mathrm{A}_{ij}^T.$$
(2) Taking the expectation value of $H_{\mathrm{MF}}$ on state $\Xi$:
$$\sum_{j \in \mathcal{V}_i} \left\langle \frac{i}{4} c_i^T \mathrm{A}_{ij} c_j \right\rangle = \frac{1}{4} \operatorname{tr}\!\left( \mathrm{A}_{ij}^T \mathrm{U}_{ij} \right) = \frac{1}{2} \sum_{a,b} \mathbf{J}_{ij}^{ab} \operatorname{tr}\!\left( \mathrm{U}_{ij}^T \mathrm{N}^a \mathrm{U}_{ij} \mathrm{N}^b \right).$$
Subtracting $\mathrm{C}$ removes exactly half of the quadratic terms, restoring the coefficient $\frac{1}{8} \times 2 = \frac{1}{4}$ of $E[\Xi]$.
(3) The functional $E[\Xi]$ is quadratic in $(\mathrm{U}, \mathrm{V})$. Differentiating $E[\Xi]$ with respect to $\mathrm{U}_{ij}$ and $\mathrm{V}_i$ yields the coefficient matrices $\mathrm{A}_{ij}$ and $\mathrm{B}_i$, proving that $\langle H_{\mathrm{MF}}[\Xi] \rangle_{\Xi'}$ evaluates to the tangent hyperplane $E[\Xi] + \nabla_\Xi E \cdot (\Xi' - \Xi)$. $\square$

**Corollary (Self-Consistency as Variational Stationarity)** <!-- claim: CLM-M1-MF-006 -->
Let $\rho_0[\Xi]$ denote the ground state obtained by diagonalizing $H_{\mathrm{MF}}[\Xi]$. The self-consistent fixed-point condition $\Xi = \Xi[\rho_0[\Xi]]$ is equivalent to the first-order Karush–Kuhn–Tucker (KKT) stationarity condition:
$$\nabla_\Xi E[\Xi, \boldsymbol{\lambda}] = 0.$$
In the pure Kitaev limit ($J = \Gamma = \Gamma' = D = \mathbf{h} = \boldsymbol{\lambda} = \mathbf{0}$), this mean-field ground state matches the exact Lieb bound.

---

## References

* **Freitas, L. R. D., Bauer, T., Egger, R., and Pereira, R. G.** (2024). *Electric polarization near vortices in the extended Kitaev model*. npj Quantum Materials **9**, 33. DOI: [10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5). arXiv: [2311.15750](https://arxiv.org/abs/2311.15750).
* **Kitaev, A.** (2006). *Anyons in an exactly solved model and beyond*. Annals of Physics **321**, 2–111. DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005). arXiv: [cond-mat/0506438](https://arxiv.org/abs/cond-mat/0506438).
