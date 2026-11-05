---
id: RND-P1-M1.4.5-T3
title: "Quasi-Free Covariance Manifolds, Random Stress-Testing, and Variational Saddle-Point Stability"
phase: 1
milestone: M1.4.5
task: T3
status: completed
claims:
  - CLM-M1-STAB-001
  - CLM-M1-STAB-002
  - CLM-M1-STAB-003
  - CLM-M1-STAB-004
---

# Note 05: Quasi-Free Covariance Manifolds, Random Stress-Testing, and Variational Saddle-Point Stability

## 1. Context, Conventions, and Disambiguation

Let $\Lambda = (\mathcal{V}, \mathcal{E})$ be the bipartite honeycomb lattice graph with $N_{\mathrm{sub}}$ sublattices and directed bonds $b = \langle ij \rangle_\gamma$ oriented from sublattice $A$ to sublattice $B$ (Freitas et al., 2024).

### 1.1 Strict Notation Contracts
To eliminate collisions between state parameters, exchange couplings, and correlation matrices:
- **Mean-Field Covariance Configuration $\Xi$:**
  $$\Xi \equiv (\mathrm{U}, \mathrm{V}),$$
  where $\mathrm{U} \in \mathbb{R}^{N_{\mathrm{bonds}} \times 4 \times 4}$ denotes directed-bond inter-site correlators $\mathrm{U}_b^{\mu\nu} = \langle i c_i^\mu c_j^\nu \rangle_{\mathrm{G}}$, and $\mathrm{V} \in \mathfrak{so}(4)^{N_{\mathrm{sub}}}$ denotes real skew-symmetric on-site correlators $\mathrm{V}_s^{\mu\nu} = \langle i c_s^\mu c_s^\nu \rangle_{\mathrm{G}} - i\delta^{\mu\nu}$. The symbol $\Xi$ is the sole collective identifier for mean-field state parameters.
- **Bilinear Exchange Tensor $\mathbf{J}$ and Coupling $\Gamma$:**
  The symbol $\Gamma$ (and $\Gamma'$) is strictly reserved for the bond-directional symmetric off-diagonal exchange interaction in the extended Kitaev model:
  $$H_{\mathrm{ex}} \supset \Gamma_\gamma \left( \sigma_i^{\gamma^+} \sigma_j^{\gamma^-} + \sigma_i^{\gamma^-} \sigma_j^{\gamma^+} \right).$$
- **Local Two-Site Correlation Matrix $\mathcal{C}_{ij}$:**
  The local $8 \times 8$ two-site Majorana correlation matrix is denoted by $\mathcal{C}_{ij}$ (or $\mathcal{M}_{ij}$), never $\Gamma$:
  $$\mathcal{C}_{ij} \equiv \langle \mathbf{c}_{ij} \mathbf{c}_{ij}^T \rangle = \begin{pmatrix} \mathbb{1}_4 - i \mathrm{V}_i & -i \mathrm{U}_b \\ i \mathrm{U}_b^T & \mathbb{1}_4 - i \mathrm{V}_j \end{pmatrix} \in \mathbb{C}^{8 \times 8}.$$

---

## 2. Geometry of the Quasi-Free State Manifold

**Definition (Real Orthogonal Grassmannian $\mathcal{M}$)** <!-- claim: CLM-M1-STAB-001 -->
Let $V \cong \mathbb{R}^{8N}$ be the Euclidean space spanned by the Majorana generators. A quasi-free state $\omega$ over the self-dual CAR algebra $\mathcal{A}_{\mathrm{CAR}}(V)$ is uniquely specified by its real skew-symmetric covariance matrix $\Omega \in \mathfrak{so}(8N)$:
$$\omega(c_i^\mu c_j^\nu) = \delta_{ij}\delta^{\mu\nu} + i \Omega_{(i\mu),(j\nu)}.$$
Fermionic quantum $N$-representability (positivity of the density matrix $\rho \ge 0$) requires:
$$-\mathbb{1}_{8N} \le \Omega^2 \le 0 \iff \mathbb{1}_{8N} + i \Omega \ge 0.$$
A state is pure if and only if $\Omega^2 = -\mathbb{1}_{8N}$. The manifold of pure quasi-free states is the real orthogonal Grassmannian:
$$\mathcal{M} = \{ \Omega \in \mathfrak{so}(8N) : \Omega^2 = -\mathbb{1}_{8N} \} \cong \frac{\mathcal{O}(8N)}{U(4N)}.$$
The tangent space to $\mathcal{M}$ at point $\Omega$ is:
$$T_\Omega \mathcal{M} = \{ \Delta \in \mathfrak{so}(8N) : \{\Omega, \Delta\} = \Omega \Delta + \Delta \Omega = 0 \}.$$

**Theorem (Exact Affine Tangent Hyperplane on Arbitrary States)** <!-- claim: CLM-M1-STAB-002 -->
Let $\Xi_1 = (\mathrm{U}_1, \mathrm{V}_1)$ and $\Xi_2 = (\mathrm{U}_2, \mathrm{V}_2)$ be two arbitrary quasi-free covariance configurations. The single-particle mean-field Hamiltonian $H_{\mathrm{MF}}[\Xi_1]$ evaluated on state $\Xi_2$ satisfies:
$$\langle H_{\mathrm{MF}}[\Xi_1] \rangle_{\Xi_2} \equiv E[\Xi_1, \boldsymbol{\lambda}] + \mathrm{d}E_{\Xi_1}[\Xi_2 - \Xi_1].$$
Because $E[\Xi, \boldsymbol{\lambda}]$ is strictly quadratic in $\Xi$, the tangent representation is exact to machine precision ($< 10^{-14}$) across arbitrary non-stationary and non-converged configurations.

*Proof.*
Expanding the variational energy functional from Note 01:
$$E[\Xi] = \frac{1}{2} \sum_{b, a, b'} J_b^{ab'} \left[ \frac{1}{4} m_i^a(\mathrm{V}) m_j^{b'}(\mathrm{V}) - \frac{1}{2} \operatorname{tr}(N^a \mathrm{U}_b N^{b'} \mathrm{U}_b^T) \right] + E_h[\mathrm{V}] + E_\lambda[\mathrm{V}].$$
The terms are bilinear in $(\mathrm{V}_i, \mathrm{V}_j)$ and $(\mathrm{U}_b, \mathrm{U}_b)$.
For any quadratic form $f(x) = \frac{1}{2} x^T K x + b^T x$, the affine identity $x_2^T K x_1 + b^T x_2 = f(x_1) + \nabla f(x_1)^T (x_2 - x_1)$ holds identically.
Matching coefficients with the definition of $H_{\mathrm{MF}}[\Xi_1] = \sum \frac{i}{4} c_i^T A_{ij}(\mathrm{U}_1) c_j + \sum \frac{i}{4} c_i^T B_i(\mathrm{V}_1) c_i - C(\Xi_1)$ confirms that taking the expectation value $\langle \cdot \rangle_{\Xi_2}$ reproduces the exact first-order Taylor expansion about $\Xi_1$. $\square$

---

## 3. Sampling and Stress-Testing via Random Quasi-Free States

**Proposition (Lie-Algebra Gaussian Sampling of Pure Quasi-Free States on the Grassmannian)** <!-- claim: CLM-M1-STAB-003 -->
Let $W \in \mathfrak{so}(4N_{\mathrm{sub}})$ be a real skew-symmetric matrix ($W^T = -W$) with Gaussian-distributed independent entries $W_{ab} \sim \mathcal{N}(0, \sigma^2)$ ($a < b$). The matrix exponential:
$$Q = \exp(W) \in \mathcal{SO}(4N_{\mathrm{sub}})$$
is a real special orthogonal matrix. Given a canonical reference complex structure $\Omega_0 = \bigoplus_{k=1}^{2N_{\mathrm{sub}}} \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$, the transformed matrix:
$$\Omega_{\mathrm{rand}} = Q \Omega_0 Q^T$$
is a pure quasi-free covariance matrix satisfying $\Omega_{\mathrm{rand}}^T = -\Omega_{\mathrm{rand}}$ and $\Omega_{\mathrm{rand}}^2 = -\mathbb{1}$.
*Note on Measure:* Pushing forward an isotropic Gaussian measure on the Lie algebra $\mathfrak{so}(n)$ via the exponential map does not produce the uniform Haar measure on $\mathcal{SO}(n)$. However, $\Omega_{\mathrm{rand}}$ constitutes a mathematically exact extremal point of $\mathcal{M}$ suitable for non-stationary stress testing.

**Protocol (Non-Stationary Stress-Testing Suite):**
Feeding $\Omega_{\mathrm{rand}}$ into the MMFT services allows non-iterative mathematical testing:
1. **Kinematic Parseval Conservation:** Verify $\frac{1}{N_{\mathbf{k}}}\sum \mathcal{E}_n^2 = \frac{1}{16}\sum \Vert{}B_s\Vert{}_{\mathrm{F}}^2 + \frac{1}{8}\sum \Vert{}A_b\Vert{}_{\mathrm{F}}^2$ on disordered configurations with non-zero cross-hoppings $\mathrm{U}_b^{0\gamma} \neq 0$.
2. **Spectral Projector Retraction:** Verify that passing a randomized Bloch Hamiltonian $\mathcal{H}_{\mathrm{rand}}(\mathbf{k})$ through `compute_spectral_projector` outputs a projector $P(\mathbf{k})$ that strictly satisfies idempotency $P^2 = 2P$ and trace filling $\operatorname{tr} P = 4N_{\mathrm{sub}}$.
3. **Casimir-Parity Envelope:** Verify that the extracted on-site blocks satisfy $0 \le \langle \mathbf{T}_s^2 \rangle_{\mathrm{G}} = \frac{3}{8}(1 - \operatorname{Pf}(\mathrm{V}_s)) \le \frac{3}{4}$.

---

## 4. Formal Definition of Variational Saddle-Point Stability

The self-consistent mean-field solution is defined as the stationary point of the min-max problem:
$$\min_{\Xi \in \mathcal{M}} \max_{\boldsymbol{\lambda} \in \mathbb{R}^{3N_{\mathrm{sub}}}} \mathcal{L}(\Xi, \boldsymbol{\lambda}), \qquad \mathcal{L}(\Xi, \boldsymbol{\lambda}) = E_{\mathrm{spin}}[\Xi] - \frac{1}{4} \sum_{s,\gamma} \lambda_s^\gamma \operatorname{tr}(\mathrm{G}^\gamma \mathrm{V}_s).$$

**Definition (Four Criteria for Variational Stability)** <!-- claim: CLM-M1-STAB-004 -->
A stationary point $(\Xi^*, \boldsymbol{\lambda}^*)$ satisfying KKT conditions $\nabla_\Xi \mathcal{L} = 0, \nabla_{\boldsymbol{\lambda}} \mathcal{L} = 0$ is stable if and only if:

1. **Primal Tangent Stability (Local Energy Minimum):**
   The Hessian of the energy restricted to the tangent space of the Grassmannian $T_{\Xi^*} \mathcal{M}$ and to the gauge-neutral manifold $\ker(\nabla_\Xi \mathbf{g})$ is positive semi-definite:
   $$\delta \Xi^T \left( \nabla_\Xi^2 \mathcal{E}(\Xi^*) \right) \delta \Xi \ge 0 \quad \forall \delta \Xi \in T_{\Xi^*} \mathcal{M} \cap \ker(\nabla_\Xi \mathbf{g}).$$
   Equivalently, all particle-hole RPA excitation energies are real and non-negative (absence of negative-energy Bogoliubov modes).

2. **Dual Multiplier Stability (Gauge Concavity):**
   The concentrated dual energy functional $E(\boldsymbol{\lambda}) \equiv \min_{\Xi \in \mathcal{M}} \mathcal{L}(\Xi, \boldsymbol{\lambda})$ is strictly concave in the multipliers:
   $$\mathcal{J}_{(s\alpha), (s'\beta)} \equiv \frac{\partial^2 E}{\partial \lambda_s^\alpha \partial \lambda_{s'}^\beta} = -\chi_{s\alpha, s'\beta}^{\mathrm{gauge}}(\omega=0) \le 0.$$
   The static gauge susceptibility matrix $\chi^{\mathrm{gauge}} \ge 0$ is positive semi-definite.

3. **Thermodynamic Stability (Magnetic Concavity):**
   The ground-state energy density is concave with respect to external magnetic fields:
   $$\chi^{\alpha\beta} \equiv -\frac{\partial^2 E_{\mathrm{cell}}}{\partial h^\alpha \partial h^\beta} = \frac{\partial \langle \sigma_{\mathrm{tot}}^\alpha \rangle}{\partial h^\beta} \ge 0.$$

4. **Algorithmic Fixed-Point Contraction:**
   Let $\Phi(X)$ denote the self-consistent map $X^{(n+1)} = \Phi(X^{(n)})$ for state vector $X = (\operatorname{vec}(\mathrm{U}), \operatorname{vec}(\mathrm{V}), \boldsymbol{\lambda})^T$.
   The Jacobian $\mathbf{J}_\Phi = D\Phi(X^*)$ restricted to the DIIS iterative subspace satisfies:
   $$\rho(\mathbf{J}_\Phi) < 1.0,$$
   where $\rho(\cdot)$ denotes the spectral radius.

---

## 5. Numerical Measurement Protocols for Stability Metrics

### Protocol A: Finite-Difference Thermodynamic Susceptibility
For a converged ground state at field $\mathbf{h}_0$:
1. Perturb the field along Cartesian direction $\alpha \in \{x, y, z\}$: $\mathbf{h}_\pm = \mathbf{h}_0 \pm \delta \hat{\mathbf{e}}_\alpha$ with $\delta = 10^{-3} \vert{}K\vert{}$.
2. Re-converge the solver to tolerance $10^{-8}$ yielding total magnetizations $m_\pm^\alpha = \sum_s \langle \sigma_s^\alpha \rangle_{\mathbf{h}_\pm}$.
3. Evaluate the diagonal susceptibility:
   $$\chi^{\alpha\alpha} = \frac{m_+^\alpha - m_-^\alpha}{2\delta} \ge -10^{-6}.$$
   A negative value certifies that the solution converged to an unstable saddle point.

### Protocol B: Directional Tangent Energy Curvature
1. Given converged state $\Xi^* = (\mathrm{U}^*, \mathrm{V}^*)$, generate an infinitesimal tangent perturbation $\Delta = (\delta \mathrm{U}, \delta \mathrm{V}) \in T_{\Xi^*} \mathcal{M}$.
2. Verify gauge-neutrality: $\operatorname{tr}(\mathrm{G}^\gamma \delta \mathrm{V}_s) = 0$.
3. Evaluate the symmetric energy second difference:
   $$\kappa_\Delta = \frac{E[\Xi^* + \epsilon \Delta] - 2E[\Xi^*] + E[\Xi^* - \epsilon \Delta]}{\epsilon^2} \ge -10^{-8}.$$

### Protocol C: Asymptotic DIIS Residual Contraction Rate
During the self-consistency loop, record the residual history $r_n = \max(\Vert{}\Delta \mathrm{U}\Vert{}_\infty, \Vert{}\Delta \mathrm{V}\Vert{}_\infty, \Vert{}\mathbf{g}\Vert{}_\infty)$.
Compute the empirical contraction quotient over the final 5 iterations:
$$\mu \equiv \frac{1}{5} \sum_{k=0}^4 \frac{r_{N-k}}{r_{N-k-1}}.$$
For a stable fixed point, $\mu < 1.0$. If $\mu \ge 1.0$, the solver is trapped in a limit cycle or expanding orbit.

---

## References

* **Freitas, L. R. D.** (2024). *Majorana mean-field theory with $\mathbb{Z}_2$ vortices: electrical tests in the extended Kitaev model*. Ph.D. thesis, Universidade Federal do Rio Grande do Norte (UFRN) / Heinrich-Heine-Universität Düsseldorf.
* **Freitas, L. R. D., Bauer, T., Egger, R., and Pereira, R. G.** (2024). *Electric polarization near vortices in the extended Kitaev model*. npj Quantum Materials **9**, 33. DOI: [10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5). arXiv: [2311.15750](https://arxiv.org/abs/2311.15750).
* **Kitaev, A.** (2006). *Anyons in an exactly solved model and beyond*. Annals of Physics **321**, 2–111. DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005). arXiv: [cond-mat/0506438](https://arxiv.org/abs/cond-mat/0506438).
* **Yilmaz, F., Kampf, A. P., and Yip, S. K.** (2022). *Phase diagrams of Kitaev models for arbitrary magnetic field orientations*. Physical Review Research **4**, 043024. DOI: [10.1103/PhysRevResearch.4.043024](https://doi.org/10.1103/PhysRevResearch.4.043024). arXiv: [2202.06158](https://arxiv.org/abs/2202.06158).
