---
id: RND-P1-M1.3-T1
title: "Reciprocal-Space Spectral Projector, Real-Space Correlator Integration, and Joint Anderson Acceleration"
phase: 1
milestone: M1.3
task: T1
status: completed
claims:
  - CLM-M1-COR-001
  - CLM-M1-COR-002
  - CLM-M1-COR-003
  - CLM-M1-COR-004
  - CLM-M1-COR-005
  - CLM-M1-COR-006
  - CLM-M1-COR-007
  - CLM-M1-COR-008
---

# Note 03: Reciprocal-Space Spectral Projector, Real-Space Correlator Integration, and Joint Anderson Acceleration

## 1. Context, Reciprocal Majorana Operators, and Spectral Projector

Let $\Lambda = (\mathcal{V}, \mathcal{E})$ denote a periodic honeycomb lattice graph where the unit cell contains $N_{\mathrm{sub}}$ sublattices ($s \in \{0, \dots, N_{\mathrm{sub}}-1\}$) and $N_{\mathrm{bonds}}$ directed bonds $b = \langle ij \rangle_\gamma \in \mathcal{E}$ with integer cell offsets $\Delta\mathbf{R}_b = (n_{1,b}, n_{2,b}) \in \mathbb{Z}^2$ and bond types $\gamma \in \{1, 2, 3\} \equiv \{x, y, z\}$ (Kitaev, 2006; Freitas et al., 2024).

As established in Note 02, the Fourier-transformed four-flavor Majorana operators satisfy the self-dual CAR algebra involution and anticommutation relations (Burnell and Nayak, 2011; Seifert et al., 2018):
$$c_{s,\mathbf{k}}^\dagger = c_{s,-\mathbf{k}}, \qquad \{ c_{s,\mathbf{k}}^\mu, (c_{s',\mathbf{k}'}^\nu)^\dagger \} = 2 \delta_{ss'} \delta^{\mu\nu} \delta_{\mathbf{k}\mathbf{k}'} \mathbb{1}.$$

The quadratic mean-field Hamiltonian in reciprocal space is:
$$H_{\mathrm{MF}} = \frac{1}{2} \sum_{\mathbf{k} \in \mathrm{BZ}} \boldsymbol{\Psi}_{\mathbf{k}}^\dagger \mathcal{H}(\mathbf{k}) \boldsymbol{\Psi}_{\mathbf{k}} - \mathrm{C},$$
where $\boldsymbol{\Psi}_{\mathbf{k}} = (c_{0,\mathbf{k}}^T, \dots, c_{N_{\mathrm{sub}}-1,\mathbf{k}}^T)^T$ is the $4N_{\mathrm{sub}}$-component Nambu-Majorana spinor and $\mathcal{H}(\mathbf{k}) = \frac{i}{4}\mathcal{M}(\mathbf{k})$ is Hermitian ($\mathcal{H}^\dagger(\mathbf{k}) = \mathcal{H}(\mathbf{k})$) with particle-hole symmetry $\mathcal{H}(-\mathbf{k}) = -\mathcal{H}^*(\mathbf{k}) = -\mathcal{H}^T(\mathbf{k})$.

**Definition (Reciprocal Two-Point Correlation Matrix)** <!-- claim: CLM-M1-COR-001 -->
The reciprocal Majorana two-point correlation matrix $P(\mathbf{k}) \in \mathbb{C}^{4N_{\mathrm{sub}} \times 4N_{\mathrm{sub}}}$ is defined by its matrix elements:
$$P_{s\mu, s'\nu}(\mathbf{k}) \equiv \langle c_{s,\mathbf{k}}^\mu (c_{s',\mathbf{k}}^\nu)^\dagger \rangle.$$

**Theorem (Zero-Temperature Spectral Projector)** <!-- claim: CLM-M1-COR-002 -->
Let $\mathcal{H}(\mathbf{k}) \mathbf{w}_n(\mathbf{k}) = \mathcal{E}_n(\mathbf{k}) \mathbf{w}_n(\mathbf{k})$ denote the spectral decomposition of $\mathcal{H}(\mathbf{k})$ with orthonormal eigenvectors $\mathbf{w}_n^\dagger(\mathbf{k}) \mathbf{w}_m(\mathbf{k}) = \delta_{nm}$. At zero temperature ($T = 0$), the reciprocal correlation matrix is given by:
$$P(\mathbf{k}) = 2 \sum_{n=1}^{4N_{\mathrm{sub}}} \Theta(-\mathcal{E}_n(\mathbf{k}))\, \mathbf{w}_n(\mathbf{k}) \mathbf{w}_n^\dagger(\mathbf{k}) = 2 \Pi^-(\mathbf{k}),$$
where $\Pi^-(\mathbf{k})$ is the orthogonal projector onto the occupied negative-energy subspace ($\mathcal{E}_n < 0$), with zero-energy modes ($\mathcal{E}_n = 0$) weighted by $\frac{1}{2}$ to preserve exact particle-hole symmetry.

*Proof.*
Expand the Majorana field in terms of canonical complex Bogoliubov quasiparticle operators $\eta_{n,\mathbf{k}}$:
$$c_{s,\mathbf{k}}^\mu = \sqrt{2} \sum_{n=1}^{4N_{\mathrm{sub}}} w_{n, s\mu}(\mathbf{k}) \eta_{n,\mathbf{k}}, \qquad \eta_{n,\mathbf{k}} = \frac{1}{\sqrt{2}} \sum_{s,\mu} w_{n, s\mu}^*(\mathbf{k}) c_{s,\mathbf{k}}^\mu.$$
Using $\{c_{s,\mathbf{k}}^\mu, (c_{s',\mathbf{k}'}^\nu)^\dagger\} = 2\delta_{ss'}\delta^{\mu\nu}\delta_{\mathbf{k}\mathbf{k}'}\mathbb{1}$ and the orthonormality $\sum_{s,\mu} w_{n, s\mu}^*(\mathbf{k}) w_{m, s\mu}(\mathbf{k}) = \delta_{nm}$, direct substitution verifies that the quasiparticle operators satisfy standard fermion anticommutation relations:
$$\{ \eta_{n,\mathbf{k}}, \eta_{m,\mathbf{k}'}^\dagger \} = \frac{1}{2} \sum_{s,\mu} \sum_{s',\nu} w_{n,s\mu}^*(\mathbf{k}) w_{m,s'\nu}(\mathbf{k}') \{ c_{s,\mathbf{k}}^\mu, (c_{s',\mathbf{k}'}^\nu)^\dagger \} = \delta_{nm} \delta_{\mathbf{k}\mathbf{k}'} \mathbb{1}.$$
In the many-body ground state $|0\rangle$, positive-energy modes ($\mathcal{E}_n > 0$) correspond to empty excitations ($\eta_{n,\mathbf{k}}|0\rangle = 0$, $\langle \eta_{n,\mathbf{k}}^\dagger \eta_{n,\mathbf{k}} \rangle = 0$), while negative-energy modes ($\mathcal{E}_n < 0$) are fully occupied ($\eta_{n,\mathbf{k}}^\dagger|0\rangle = 0$, $\langle \eta_{n,\mathbf{k}} \eta_{n,\mathbf{k}}^\dagger \rangle = 1$).
Therefore:
$$\langle c_{s,\mathbf{k}}^\mu (c_{s',\mathbf{k}}^\nu)^\dagger \rangle = 2 \sum_{n,m} w_{n,s\mu}(\mathbf{k}) w_{m,s'\nu}^*(\mathbf{k}) \langle \eta_{n,\mathbf{k}} \eta_{m,\mathbf{k}}^\dagger \rangle = 2 \sum_{\mathcal{E}_n(\mathbf{k}) < 0} w_{n,s\mu}(\mathbf{k}) w_{n,s'\nu}^*(\mathbf{k}),$$
which is identically $2 \Pi^-(\mathbf{k})$. $\square$

**Proposition (Symmetries and CAR Sum Rules)** <!-- claim: CLM-M1-COR-003 -->
1. By construction, $P(\mathbf{k})$ is positive semi-definite and Hermitian: $P^\dagger(\mathbf{k}) = P(\mathbf{k})$.
2. Under particle-hole symmetry $\mathcal{H}(-\mathbf{k}) = -\mathcal{H}^*(\mathbf{k})$, if $\mathbf{w}_n(\mathbf{k})$ has eigenvalue $\mathcal{E}_n(\mathbf{k})$, then $\mathbf{w}_n^*(-\mathbf{k})$ has eigenvalue $-\mathcal{E}_n(-\mathbf{k})$. Consequently:
   $$\Pi^-(-\mathbf{k}) = \mathbb{1} - \left(\Pi^-(\mathbf{k})\right)^* \implies P(-\mathbf{k}) = 2\mathbb{1} - P^*(\mathbf{k}) = 2\mathbb{1} - P^T(\mathbf{k}).$$
3. Taking the real part and averaging over the first Brillouin zone:
   $$\frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k} \in \mathrm{BZ}} \operatorname{Re}\!\left[ P_{s\mu, s'\nu}(\mathbf{k}) \right] = \delta_{ss'} \delta^{\mu\nu} \mathbb{1},$$
   reproducing the local CAR trace normalization $\frac{1}{4} \operatorname{tr}(\langle c_{s,\mathbf{R}} c_{s,\mathbf{R}}^T \rangle) = 1.0$ across every sublattice $s$.

---

## 2. Real-Space Correlator Integration

**Theorem (Brillouin-Zone Correlator Formulae)** <!-- claim: CLM-M1-COR-004 -->
Let $b = \langle ij \rangle_\gamma$ be a directed bond oriented from sublattice $s_i$ to sublattice $s_j$ with Bravais offset $\Delta\mathbf{R}_b = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2$. In terms of fractional crystal momenta $\mathbf{k} = (k_1, k_2) \in [0, 1)^2$, the bond phase is $\phi_b(\mathbf{k}) = 2\pi (k_1 n_1 + k_2 n_2)$. The real-space correlators $\mathrm{U}_b^{\mu\nu} = \langle i c_{s_i,\mathbf{0}}^\mu c_{s_j,\Delta\mathbf{R}_b}^\nu \rangle$ and $\mathrm{V}_s^{\mu\nu} = \langle i c_{s}^\mu c_{s}^\nu \rangle - i\delta^{\mu\nu}$ evaluate to:
$$\mathrm{U}_b^{\mu\nu} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \left( \operatorname{Im}\!\left[ P_{s_i\mu, s_j\nu}(\mathbf{k}) \right] \cos(\phi_b(\mathbf{k})) - \operatorname{Re}\!\left[ P_{s_i\mu, s_j\nu}(\mathbf{k}) \right] \sin(\phi_b(\mathbf{k})) \right),$$
$$\mathrm{V}_s^{\mu\nu} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Im}\!\left[ P_{s\mu, s\nu}(\mathbf{k}) \right].$$

*Proof.*
In reciprocal space, $P_{s_i\mu, s_j\nu}(\mathbf{k}) = \langle c_{s_i,\mathbf{k}}^\mu (c_{s_j,\mathbf{k}}^\nu)^\dagger \rangle = \langle c_{s_i,\mathbf{k}}^\mu c_{s_j,-\mathbf{k}}^\nu \rangle$.
By the CAR anticommutation relation $\{c_{s_i,\mathbf{0}}^\mu, c_{s_j,\Delta\mathbf{R}_b}^\nu\} = 2\delta_{s_i s_j}\delta^{\mu\nu}\delta_{\Delta\mathbf{R}_b,\mathbf{0}}\mathbb{1}$, the real-space operator ordering is related to its reverse by:
$$\langle c_{s_i,\mathbf{0}}^\mu c_{s_j,\Delta\mathbf{R}_b}^\nu \rangle = 2\delta_{s_i s_j}\delta^{\mu\nu}\delta_{\Delta\mathbf{R}_b,\mathbf{0}} - \langle c_{s_j,\Delta\mathbf{R}_b}^\nu c_{s_i,\mathbf{0}}^\mu \rangle.$$
Expanding the second term via the Fourier transform yields:
$$\langle c_{s_j,\Delta\mathbf{R}_b}^\nu c_{s_i,\mathbf{0}}^\mu \rangle = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k},\mathbf{k}'} \langle c_{s_j,\mathbf{k}'}^\nu c_{s_i,\mathbf{k}}^\mu \rangle e^{i \mathbf{k}' \cdot \Delta\mathbf{R}_b} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \langle c_{s_j,-\mathbf{k}}^\nu c_{s_i,\mathbf{k}}^\mu \rangle e^{-i \mathbf{k} \cdot \Delta\mathbf{R}_b} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} P_{s_i\mu, s_j\nu}(\mathbf{k}) e^{-i \phi_b(\mathbf{k})}.$$
Therefore, the general two-point expectation value is:
$$\langle c_{s_i,\mathbf{0}}^\mu c_{s_j,\Delta\mathbf{R}_b}^\nu \rangle = 2\delta_{s_i s_j}\delta^{\mu\nu}\delta_{\Delta\mathbf{R}_b,\mathbf{0}} - \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} P_{s_i\mu, s_j\nu}(\mathbf{k}) e^{-i \phi_b(\mathbf{k})}.$$
For distinct sites on directed bonds ($s_i \neq s_j$), the Kronecker delta vanishes:
$$\langle c_{s_i,\mathbf{0}}^\mu c_{s_j,\Delta\mathbf{R}_b}^\nu \rangle = - \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} P_{s_i\mu, s_j\nu}(\mathbf{k}) e^{-i \phi_b(\mathbf{k})}.$$
Multiplying by $i$ and taking the real part:
$$\mathrm{U}_b^{\mu\nu} = \operatorname{Re}\!\langle i c_{s_i,\mathbf{0}}^\mu c_{s_j,\Delta\mathbf{R}_b}^\nu \rangle = \operatorname{Re}\!\left[ -\frac{i}{N_{\mathbf{k}}} \sum_{\mathbf{k}} P_{s_i\mu, s_j\nu}(\mathbf{k}) e^{-i \phi_b(\mathbf{k})} \right] = \operatorname{Im}\!\left[ \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} P_{s_i\mu, s_j\nu}(\mathbf{k}) e^{-i \phi_b(\mathbf{k})} \right].$$
Expanding $P = \operatorname{Re}[P] + i \operatorname{Im}[P]$ and $e^{-i\phi_b} = \cos\phi_b - i\sin\phi_b$:
$$\operatorname{Im}\left[ (\operatorname{Re}[P] + i\operatorname{Im}[P])(\cos\phi_b - i\sin\phi_b) \right] = \operatorname{Im}[P]\cos\phi_b - \operatorname{Re}[P]\sin\phi_b,$$
which yields the bond correlator formula:
$$\mathrm{U}_b^{\mu\nu} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \left( \operatorname{Im}\!\left[ P_{s_i\mu, s_j\nu}(\mathbf{k}) \right] \cos(\phi_b(\mathbf{k})) - \operatorname{Re}\!\left[ P_{s_i\mu, s_j\nu}(\mathbf{k}) \right] \sin(\phi_b(\mathbf{k})) \right).$$
For on-site correlators ($\Delta\mathbf{R} = \mathbf{0}, s_i = s_j = s$), $\phi(\mathbf{k}) = 0$ and the Kronecker delta is retained:
$$\langle c_s^\mu c_s^\nu \rangle = 2\delta^{\mu\nu} - \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} P_{s\mu, s\nu}(\mathbf{k}).$$
Multiplying by $i$:
$$\langle i c_s^\mu c_s^\nu \rangle = 2i\delta^{\mu\nu} - \frac{i}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \left( \operatorname{Re}[P_{s\mu, s\nu}(\mathbf{k})] + i \operatorname{Im}[P_{s\mu, s\nu}(\mathbf{k})] \right).$$
By Proposition CLM-M1-COR-003, $\frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Re}[P_{s\mu, s\nu}(\mathbf{k})] = \delta^{\mu\nu}$. Substituting this trace identity gives:
$$\langle i c_s^\mu c_s^\nu \rangle = 2i\delta^{\mu\nu} - i\delta^{\mu\nu} + \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Im}[P_{s\mu, s\nu}(\mathbf{k})] = i\delta^{\mu\nu} + \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Im}[P_{s\mu, s\nu}(\mathbf{k})].$$
Subtracting $i\delta^{\mu\nu}$ yields:
$$\mathrm{V}_s^{\mu\nu} = \langle i c_s^\mu c_s^\nu \rangle - i\delta^{\mu\nu} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Im}\!\left[ P_{s\mu, s\nu}(\mathbf{k}) \right].$$
Because $P(\mathbf{k})$ is Hermitian, $P_{s\mu, s\nu}^*(\mathbf{k}) = P_{s\nu, s\mu}(\mathbf{k})$, which implies $\operatorname{Im}[P_{s\mu, s\nu}(\mathbf{k})] = -\operatorname{Im}[P_{s\nu, s\mu}(\mathbf{k})]$. Thus, $\mathrm{V}_s^T = -\mathrm{V}_s$ is strictly real and skew-symmetric. $\square$

**Corollary (Exact Gauge Neutrality for Flat Zero-Modes)** <!-- claim: CLM-M1-COR-005 -->
When $u_{\mathrm{matter}} = 0$ in the pure Kitaev state, all six gauge Majorana bands reside at zero energy ($\mathcal{E}_{c^\gamma}(\mathbf{k}) = 0$). Under the particle-hole symmetric assignment $\Theta(0) = \frac{1}{2}$, the spectral projector satisfies:
$$P_{s\gamma, s\gamma}(\mathbf{k}) = 1.0, \qquad P_{s\gamma, s\nu}(\mathbf{k}) = 0 \quad (\gamma \neq \nu).$$
Consequently, $\operatorname{Im}[P_{s\mu, s\nu}(\mathbf{k})] = 0$ for all on-site blocks, ensuring that $\mathrm{V}_s = \mathbf{0}$ and all gauge residuals $g_s^\gamma = \operatorname{tr}(\mathrm{G}^\gamma \mathrm{V}_s) = 0$ vanish identically and independently of the $\mathbf{k}$-grid discretization.

---

## 3. Pure Kitaev Self-Consistent Ground-State Correlator (K-006)

**Theorem (Matter Correlator Invariant K-006)** <!-- claim: CLM-M1-COR-006 -->
In the pure Kitaev limit ($K < 0$, $J = \Gamma = \Gamma' = D = \mathbf{h} = \boldsymbol{\lambda} = \mathbf{0}$) with locked zero-flux gauge background $\mathrm{U}_b^{\gamma\gamma} = u_{\mathrm{gauge}} = +1.0$:
1. The on-site correlation matrix vanishes identically: $\mathrm{V}_s = \mathbf{0}$, and all gauge residuals satisfy $g_s^\gamma = \operatorname{tr}(\mathrm{G}^\gamma \mathrm{V}_s) = 0$.
2. The self-consistent matter bond correlator satisfies the exact sign structure:
   $$u_{\mathrm{matter}} = \operatorname{sgn}(K) \operatorname{sgn}(u_{\mathrm{gauge}}) \bar{U},$$
   where $\bar{U} = \frac{1}{3}\langle |f(\mathbf{k})| \rangle_{\mathrm{BZ}} > 0$. For $K < 0$ and $u_{\mathrm{gauge}} = +1.0$:
   $$\mathrm{U}_b^{00} = -\bar{U} = -\int_0^1 \mathrm{d}k_1 \int_0^1 \mathrm{d}k_2 \frac{1 + \cos(2\pi k_1) + \cos(2\pi k_2)}{\sqrt{3 + 2\cos(2\pi k_1) + 2\cos(2\pi k_2) + 2\cos(2\pi(k_1 - k_2))}} \approx -0.5248645.$$
3. The gauge bond correlator is preserved: $\mathrm{U}_b^{\gamma\gamma} = \operatorname{sgn}(K u_{\mathrm{matter}}) = +1.0$.

*Proof.*
From Theorem CLM-M1-BLO-006 of Note 02, the $2 \times 2$ matter sector in the 2-sublattice cell is governed by $\mathcal{H}_{A0, B0}(\mathbf{k}) = \frac{i|K|}{2} f(\mathbf{k})$, where $f(\mathbf{k}) = 1 + e^{-i2\pi k_1} + e^{-i2\pi k_2} = |f(\mathbf{k})|e^{i\theta_{\mathbf{k}}}$.
The negative-energy state has eigenvector $\mathbf{w}_-(\mathbf{k}) = \frac{1}{\sqrt{2}}(1, i e^{-i\theta_{\mathbf{k}}})^T$.
The off-diagonal matter element of $P(\mathbf{k}) = 2 \mathbf{w}_-(\mathbf{k}) \mathbf{w}_-^\dagger(\mathbf{k})$ is:
$$P_{A0, B0}(\mathbf{k}) = 2 \left( \frac{1}{\sqrt{2}} \right) \left( -\frac{i}{\sqrt{2}} e^{i\theta_{\mathbf{k}}} \right) = -i e^{i\theta_{\mathbf{k}}} = -i \frac{f(\mathbf{k})}{|f(\mathbf{k})|}.$$
Writing $e^{i\theta_{\mathbf{k}}} = \cos\theta_{\mathbf{k}} + i\sin\theta_{\mathbf{k}}$, we have:
$$\operatorname{Im}\!\left[ P_{A0, B0}(\mathbf{k}) \right] = \operatorname{Im}\!\left[ -i (\cos\theta_{\mathbf{k}} + i\sin\theta_{\mathbf{k}}) \right] = \operatorname{Im}\!\left[ \sin\theta_{\mathbf{k}} - i\cos\theta_{\mathbf{k}} \right] = -\cos\theta_{\mathbf{k}} = -\operatorname{Re}\!\left[\frac{f(\mathbf{k})}{|f(\mathbf{k})|}\right].$$
For the $z$-bond ($\Delta\mathbf{R}_z = \mathbf{0}$):
$$\mathrm{U}_z^{00} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Im}[P_{A0, B0}(\mathbf{k})] = -\frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Re}\!\left[\frac{f(\mathbf{k})}{|f(\mathbf{k})|}\right] = -\frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \frac{1 + \cos(2\pi k_1) + \cos(2\pi k_2)}{|f(\mathbf{k})|}.$$
By $C_3$ symmetry (unimodular relabeling on the torus), the $x$- and $y$-bonds integrate to the exact same value. In the thermodynamic limit $N_{\mathbf{k}} \to \infty$, high-precision quadrature evaluates this integral to $\mathrm{U}_b^{00} \approx -0.5248645$. For the gauge sector, evaluating the flat bands at $\mathcal{E}_{c^\gamma} = - \frac{|K|}{2}|u_{\mathrm{matter}}|$ gives $\mathrm{U}_b^{\gamma\gamma} = \operatorname{sgn}(K u_{\mathrm{matter}}) = +1.0$. $\square$

---

## 4. Local Gauge Residuals and Joint Anderson Acceleration (DIIS)

**Lemma (Hellmann–Feynman Multiplier Gradient)** <!-- claim: CLM-M1-COR-007 -->
Let $g_s^\gamma \equiv \operatorname{tr}(\mathrm{G}^\gamma \mathrm{V}_s)$ denote the local gauge constraint residual. The partial derivative of the ground-state mean-field energy $E_0$ with respect to the Lagrange multiplier $\lambda_s^\gamma$ is:
$$\frac{\partial E_0}{\partial \lambda_s^\gamma} = -\frac{1}{4} g_s^\gamma.$$
Consequently, updating the Lagrange multipliers via $\tilde{\lambda}_s^\gamma = \lambda_s^\gamma - \alpha_\lambda g_s^\gamma$ with $\alpha_\lambda > 0$ corresponds to gradient ascent on the concave dual energy functional, driving $g_s^\gamma \to 0$.

*Proof.*
In the quadratic Hamiltonian $H_{\mathrm{MF}}$, the multiplier term is $\sum_{s,\gamma} \lambda_s^\gamma \frac{i}{4} c_s^T \mathrm{G}^\gamma c_s$. By the Hellmann–Feynman theorem:
$$\frac{\partial E_0}{\partial \lambda_s^\gamma} = \left\langle \frac{i}{4} c_s^T \mathrm{G}^\gamma c_s \right\rangle = \frac{1}{4} \sum_{\mu,\nu} (\mathrm{G}^\gamma)_{\mu\nu} \langle i c_s^\mu c_s^\nu \rangle = \frac{1}{4} \operatorname{tr}(\mathrm{G}^\gamma (\mathrm{V}_s + i\mathbb{1})) = \frac{1}{4} \operatorname{tr}(\mathrm{G}^\gamma \mathrm{V}_s^T) = -\frac{1}{4} \operatorname{tr}(\mathrm{G}^\gamma \mathrm{V}_s) = -\frac{1}{4} g_s^\gamma,$$
where we used $\operatorname{tr}(\mathrm{G}^\gamma) = 0$ and the skew-symmetry $\mathrm{V}_s^T = -\mathrm{V}_s, (\mathrm{G}^\gamma)^T = -\mathrm{G}^\gamma$. Because the saddle-point problem requires maximizing $E_0$ with respect to $\boldsymbol{\lambda}$, gradient ascent steps in the direction $+\frac{\partial E_0}{\partial \boldsymbol{\lambda}}$, giving $\Delta \lambda_s^\gamma \propto -g_s^\gamma$. $\square$

**Proposition (Joint Anderson Acceleration / DIIS KKT System)** <!-- claim: CLM-M1-COR-008 -->
Let the mean-field state be represented as a concatenated parameter vector:
$$X = \begin{pmatrix} \operatorname{vec}(\mathrm{U}) \\ \operatorname{vec}(\mathrm{V}) \\ \operatorname{vec}(\boldsymbol{\lambda}) \end{pmatrix} \in \mathbb{R}^{D}, \qquad D = 16 N_{\mathrm{bonds}} + 16 N_{\mathrm{sub}} + 3 N_{\mathrm{sub}}.$$
Given an iterate $X^{(\ell)}$, the single-particle Hamiltonian $\mathcal{H}(\mathbf{k})$ is assembled and projected to yield target correlators $\tilde{\mathrm{U}}(X^{(\ell)})$ and $\tilde{\mathrm{V}}(X^{(\ell)})$. The gauge residuals $g(X^{(\ell)})$ generate target multipliers $\tilde{\boldsymbol{\lambda}} = \boldsymbol{\lambda} - \alpha_\lambda g$.
Defining the fixed-point target mapping $\Phi(X^{(\ell)}) = (\operatorname{vec}(\tilde{\mathrm{U}}), \operatorname{vec}(\tilde{\mathrm{V}}), \operatorname{vec}(\tilde{\boldsymbol{\lambda}}))^T$, the iteration residual is:
$$R_\ell \equiv \Phi(X^{(\ell)}) - X^{(\ell)}.$$
Over a history subspace of depth $m$, the optimal weights $\boldsymbol{\xi} = (\xi_0, \dots, \xi_{m-1})^T$ minimizing $\|\sum_{j=0}^{m-1} \xi_j R_{\ell-j}\|_2^2$ subject to $\sum_{j=0}^{m-1} \xi_j = 1$ satisfy the regularized $(m+1) \times (m+1)$ linear Karush-Kuhn-Tucker system:
$$\begin{pmatrix} \Omega & -\mathbf{1} \\ -\mathbf{1}^T & 0 \end{pmatrix} \begin{pmatrix} \boldsymbol{\xi} \\ \nu \end{pmatrix} = \begin{pmatrix} \mathbf{0} \\ -1 \end{pmatrix},$$
where the Gram matrix elements are regularized against collinearity via Tikhonov parameter $\epsilon_{\mathrm{reg}}$:
$$\Omega_{pq} = \langle R_{\ell-p}, R_{\ell-q} \rangle + \epsilon_{\mathrm{reg}} \delta_{pq} \langle R_{\ell-p}, R_{\ell-p} \rangle.$$
The accelerated iterate is then constructed as $X^{(\ell+1)} = \sum_{j=0}^{m-1} \xi_j \Phi(X^{(\ell-j)})$.

---

## References

* **Burnell, F. J. and Nayak, C.** (2011). *SU(2) slave fermion approach to the Kitaev honeycomb lattice model*. Physical Review B **84**, 125125. DOI: [10.1103/PhysRevB.84.125125](https://doi.org/10.1103/PhysRevB.84.125125). arXiv: [1104.5485](https://arxiv.org/abs/1104.5485).
* **Freitas, L. R. D.** (2024). *Majorana mean-field theory with $\mathbb{Z}_2$ vortices: electrical tests in the extended Kitaev model*. Ph.D. thesis, Universidade Federal do Rio Grande do Norte (UFRN) / Heinrich-Heine-Universität Düsseldorf.
* **Freitas, L. R. D., Bauer, T., Egger, R., and Pereira, R. G.** (2024). *Electric polarization near vortices in the extended Kitaev model*. npj Quantum Materials **9**, 33. DOI: [10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5). arXiv: [2311.15750](https://arxiv.org/abs/2311.15750).
* **Kitaev, A.** (2006). *Anyons in an exactly solved model and beyond*. Annals of Physics **321**, 2–111. DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005). arXiv: [cond-mat/0506438](https://arxiv.org/abs/cond-mat/0506438).
* **Pulay, P.** (1980). *Convergence acceleration of iterative sequences. the case of scf iteration*. Chemical Physics Letters **73**, 393–398. DOI: [10.1016/0009-2614(80)80396-4](https://doi.org/10.1016/0009-2614(80)80396-4).
* **Seifert, U. F. P., Meng, T., and Vojta, M.** (2018). *Fractionalized Fermi liquids and exotic superconductivity in the Kitaev-Kondo lattice*. Physical Review B **97**, 085118. DOI: [10.1103/PhysRevB.97.085118](https://doi.org/10.1103/PhysRevB.97.085118). arXiv: [1710.05018](https://arxiv.org/abs/1710.05018).
* **Walker, H. F. and Ni, P.** (2011). *Anderson acceleration for fixed-point iterations*. SIAM Journal on Numerical Analysis **49**, 1715–1735. DOI: [10.1137/10078356X](https://doi.org/10.1137/10078356X).
* **Yilmaz, F., Kampf, A. P., and Yip, S. K.** (2022). *Phase diagrams of Kitaev models for arbitrary magnetic field orientations*. Physical Review Research **4**, 043024. DOI: [10.1103/PhysRevResearch.4.043024](https://doi.org/10.1103/PhysRevResearch.4.043024). arXiv: [2202.06158](https://arxiv.org/abs/2202.06158).
