# Note 03: Reciprocal-Space Spectral Projector, Real-Space Correlator Integration, and Joint Anderson Acceleration

## 1. Context, Reciprocal Majorana Operators, and Spectral Projector

Let $\Lambda = (\mathcal{V}, \mathcal{E})$ denote a periodic honeycomb lattice graph where the unit cell contains $N_{\mathrm{sub}}$ sublattices ($s \in \{0, \dots, N_{\mathrm{sub}}-1\}$) and $N_{\mathrm{bonds}}$ directed bonds $b = \langle ij \rangle_\gamma \in \mathcal{E}$ with integer cell offsets $\Delta\mathbf{R}_b = (n_{1,b}, n_{2,b}) \in \mathbb{Z}^2$ and bond types $\gamma \in \{1, 2, 3\} \equiv \{x, y, z\}$ (Kitaev, 2006; Freitas et al., 2024).

As established in Note 02, the Fourier-transformed four-flavor Majorana operators satisfy the self-dual CAR algebra involution and anticommutation relations (Burnell and Nayak, 2011; Seifert et al., 2018):
$$c_{s,\mathbf{k}}^\dagger = c_{s,-\mathbf{k}}, \qquad \{ c_{s,\mathbf{k}}^\mu, (c_{s',\mathbf{k}'}^\nu)^\dagger \} = 2 \delta_{ss'} \delta^{\mu\nu} \delta_{\mathbf{k}\mathbf{k}'} \mathbb{1}.$$

The quadratic mean-field Hamiltonian in reciprocal space is:
$$H_{\mathrm{MF}} = \frac{1}{2} \sum_{\mathbf{k} \in \mathrm{BZ}} \boldsymbol{\Psi}_{\mathbf{k}}^\dagger \mathcal{H}(\mathbf{k}) \boldsymbol{\Psi}_{\mathbf{k}} - \mathrm{C},$$
where $\boldsymbol{\Psi}_{\mathbf{k}} = (c_{0,\mathbf{k}}^T, \dots, c_{N_{\mathrm{sub}}-1,\mathbf{k}}^T)^T$ is the $4N_{\mathrm{sub}}$-component Nambu-Majorana spinor and $\mathcal{H}(\mathbf{k}) = \frac{i}{4}\mathcal{M}(\mathbf{k})$ is Hermitian ($\mathcal{H}^\dagger(\mathbf{k}) = \mathcal{H}(\mathbf{k})$) with particle-hole symmetry $\mathcal{H}(-\mathbf{k}) = -\mathcal{H}^*(\mathbf{k}) = -\mathcal{H}^T(\mathbf{k})$.

**Definition 1.1 (Reciprocal Two-Point Correlation Matrix):**
We define the reciprocal Majorana two-point correlation matrix $P(\mathbf{k}) \in \mathbb{C}^{4N_{\mathrm{sub}} \times 4N_{\mathrm{sub}}}$ by its matrix elements:
$$P_{s\mu, s'\nu}(\mathbf{k}) \equiv \langle c_{s,\mathbf{k}}^\mu (c_{s',\mathbf{k}}^\nu)^\dagger \rangle.$$

**Theorem 1.2 (Zero-Temperature Spectral Projector):**
Let $\mathcal{H}(\mathbf{k}) \mathbf{w}_n(\mathbf{k}) = E_n(\mathbf{k}) \mathbf{w}_n(\mathbf{k})$ denote the spectral decomposition of $\mathcal{H}(\mathbf{k})$ with orthonormal eigenvectors $\mathbf{w}_n^\dagger(\mathbf{k}) \mathbf{w}_m(\mathbf{k}) = \delta_{nm}$. At zero temperature ($T = 0$), the reciprocal correlation matrix is given by:
$$P(\mathbf{k}) = 2 \sum_{n=1}^{4N_{\mathrm{sub}}} \Theta(-E_n(\mathbf{k}))\, \mathbf{w}_n(\mathbf{k}) \mathbf{w}_n^\dagger(\mathbf{k}) = 2 \Pi^-(\mathbf{k}),$$
where $\Pi^-(\mathbf{k})$ is the orthogonal projector onto the occupied negative-energy subspace ($E_n < 0$), with zero-energy modes ($E_n = 0$) weighted by $\frac{1}{2}$ to preserve exact particle-hole symmetry.

*Proof:*
Expand the Majorana field in terms of canonical complex Bogoliubov quasiparticle operators $\gamma_{n,\mathbf{k}}$:
$$c_{s,\mathbf{k}}^\mu = \sqrt{2} \sum_{n=1}^{4N_{\mathrm{sub}}} w_{n, s\mu}(\mathbf{k}) \gamma_{n,\mathbf{k}}, \qquad \gamma_{n,\mathbf{k}} = \frac{1}{\sqrt{2}} \sum_{s,\mu} w_{n, s\mu}^*(\mathbf{k}) c_{s,\mathbf{k}}^\mu.$$
Using $\{c_{s,\mathbf{k}}^\mu, (c_{s',\mathbf{k}'}^\nu)^\dagger\} = 2\delta_{ss'}\delta^{\mu\nu}\delta_{\mathbf{k}\mathbf{k}'}\mathbb{1}$ and the orthonormality $\sum_{s,\mu} w_{n, s\mu}^*(\mathbf{k}) w_{m, s\mu}(\mathbf{k}) = \delta_{nm}$, direct substitution verifies that the quasiparticle operators satisfy standard fermion anticommutation relations:
$$\{ \gamma_{n,\mathbf{k}}, \gamma_{m,\mathbf{k}'}^\dagger \} = \frac{1}{2} \sum_{s,\mu} \sum_{s',\nu} w_{n,s\mu}^*(\mathbf{k}) w_{m,s'\nu}(\mathbf{k}') \{ c_{s,\mathbf{k}}^\mu, (c_{s',\mathbf{k}'}^\nu)^\dagger \} = \delta_{nm} \delta_{\mathbf{k}\mathbf{k}'} \mathbb{1}.$$
In the many-body ground state $|0\rangle$, positive-energy modes ($E_n > 0$) correspond to empty excitations ($\gamma_{n,\mathbf{k}}|0\rangle = 0$, $\langle \gamma_{n,\mathbf{k}}^\dagger \gamma_{n,\mathbf{k}} \rangle = 0$), while negative-energy modes ($E_n < 0$) are fully occupied ($\gamma_{n,\mathbf{k}}^\dagger|0\rangle = 0$, $\langle \gamma_{n,\mathbf{k}} \gamma_{n,\mathbf{k}}^\dagger \rangle = 1$).
Therefore:
$$\langle c_{s,\mathbf{k}}^\mu (c_{s',\mathbf{k}}^\nu)^\dagger \rangle = 2 \sum_{n,m} w_{n,s\mu}(\mathbf{k}) w_{m,s'\nu}^*(\mathbf{k}) \langle \gamma_{n,\mathbf{k}} \gamma_{m,\mathbf{k}}^\dagger \rangle = 2 \sum_{E_n(\mathbf{k}) < 0} w_{n,s\mu}(\mathbf{k}) w_{n,s'\nu}^*(\mathbf{k}),$$
which is identically $2 \Pi^-(\mathbf{k})$. $\square$

**Proposition 1.3 (Symmetries and CAR Sum Rules):**
1. By construction, $P(\mathbf{k})$ is positive semi-definite and Hermitian: $P^\dagger(\mathbf{k}) = P(\mathbf{k})$.
2. Under particle-hole symmetry $\mathcal{H}(-\mathbf{k}) = -\mathcal{H}^*(\mathbf{k})$, if $\mathbf{w}_n(\mathbf{k})$ has eigenvalue $E_n(\mathbf{k})$, then $\mathbf{w}_n^*(-\mathbf{k})$ has eigenvalue $-E_n(-\mathbf{k})$. Consequently:
   $$\Pi^-(-\mathbf{k}) = \mathbb{1} - \left(\Pi^-(\mathbf{k})\right)^* \implies P(-\mathbf{k}) = 2\mathbb{1} - P^*(\mathbf{k}) = 2\mathbb{1} - P^T(\mathbf{k}).$$
3. Taking the real part and averaging over the first Brillouin zone:
   $$\frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k} \in \mathrm{BZ}} \operatorname{Re}\!\left[ P_{s\mu, s'\nu}(\mathbf{k}) \right] = \delta_{ss'} \delta^{\mu\nu} \mathbb{1},$$
   reproducing the local CAR trace normalization $\frac{1}{4} \operatorname{tr}(\langle c_{s,\mathbf{R}} c_{s,\mathbf{R}}^T \rangle) = 1.0$ across every sublattice $s$.

---

## 2. Real-Space Correlator Integration

To close the self-consistent loop, the real-space expectation values $\mathrm{U}_b^{\mu\nu} = \langle i c_{s_i,\mathbf{0}}^\mu c_{s_j,\Delta\mathbf{R}_b}^\nu \rangle$ and $\mathrm{V}_s^{\mu\nu} = \langle i c_{s}^\mu c_{s}^\nu \rangle - i\delta^{\mu\nu}$ must be synthesized from $P(\mathbf{k})$.

**Theorem 2.1 (Brillouin-Zone Correlator Formulae):**
Let $b = \langle ij \rangle_\gamma$ be a directed bond oriented from sublattice $s_i$ to sublattice $s_j$ with Bravais offset $\Delta\mathbf{R}_b = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2$. In terms of fractional crystal momenta $\mathbf{k} = (k_1, k_2) \in [0, 1)^2$, the bond phase is $\phi_b(\mathbf{k}) = 2\pi (k_1 n_1 + k_2 n_2)$. The real-space correlators are:
$$\mathrm{U}_b^{\mu\nu} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \left( \operatorname{Im}\!\left[ P_{s_i\mu, s_j\nu}(\mathbf{k}) \right] \cos(\phi_b(\mathbf{k})) - \operatorname{Re}\!\left[ P_{s_i\mu, s_j\nu}(\mathbf{k}) \right] \sin(\phi_b(\mathbf{k})) \right),$$
$$\mathrm{V}_s^{\mu\nu} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Im}\!\left[ P_{s\mu, s\nu}(\mathbf{k}) \right].$$

*Proof:*
Expanding the real-space operators via inverse Fourier transform:
$$c_{s_i,\mathbf{0}}^\mu = \frac{1}{\sqrt{N_{\mathbf{k}}}} \sum_{\mathbf{k}} c_{s_i,\mathbf{k}}^\mu, \qquad c_{s_j,\Delta\mathbf{R}_b}^\nu = \frac{1}{\sqrt{N_{\mathbf{k}}}} \sum_{\mathbf{k}'} c_{s_j,\mathbf{k}'}^\nu e^{i \mathbf{k}' \cdot \Delta\mathbf{R}_b}.$$
Using the involution $c_{s_j,\mathbf{k}'}^\nu = (c_{s_j,-\mathbf{k}'}^\nu)^\dagger$, the expectation value becomes:
$$\langle c_{s_i,\mathbf{0}}^\mu c_{s_j,\Delta\mathbf{R}_b}^\nu \rangle = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}, \mathbf{k}'} \langle c_{s_i,\mathbf{k}}^\mu (c_{s_j,-\mathbf{k}'}^\nu)^\dagger \rangle e^{i \mathbf{k}' \cdot \Delta\mathbf{R}_b} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \langle c_{s_i,\mathbf{k}}^\mu (c_{s_j,\mathbf{k}}^\nu)^\dagger \rangle e^{-i \mathbf{k} \cdot \Delta\mathbf{R}_b} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} P_{s_i\mu, s_j\nu}(\mathbf{k}) e^{-i \phi_b(\mathbf{k})}.$$
Multiplying by $i$:
$$\langle i c_{s_i,\mathbf{0}}^\mu c_{s_j,\Delta\mathbf{R}_b}^\nu \rangle = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} i \left( \operatorname{Re}[P] + i \operatorname{Im}[P] \right) \left( \cos\phi_b - i \sin\phi_b \right).$$
Taking the real part yields:
$$\operatorname{Re}\!\langle i c_{s_i,\mathbf{0}}^\mu c_{s_j,\Delta\mathbf{R}_b}^\nu \rangle = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \left( \operatorname{Im}[P_{s_i\mu, s_j\nu}(\mathbf{k})] \cos\phi_b(\mathbf{k}) - \operatorname{Re}[P_{s_i\mu, s_j\nu}(\mathbf{k})] \sin\phi_b(\mathbf{k}) \right).$$
For the on-site correlator at $\Delta\mathbf{R} = \mathbf{0}$, $\phi(\mathbf{k}) = 0$, giving:
$$\langle i c_s^\mu c_s^\nu \rangle = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \left( i \operatorname{Re}[P_{s\mu, s\nu}(\mathbf{k})] - \operatorname{Im}[P_{s\mu, s\nu}(\mathbf{k})] \right).$$
By Proposition 1.3, $\frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Re}[P_{s\mu, s\nu}(\mathbf{k})] = \delta^{\mu\nu}$. Subtracting $i\delta^{\mu\nu}$ yields:
$$\mathrm{V}_s^{\mu\nu} = \langle i c_s^\mu c_s^\nu \rangle - i\delta^{\mu\nu} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Im}\!\left[ P_{s\mu, s\nu}(\mathbf{k}) \right].$$
Because $P(\mathbf{k})$ is Hermitian, $P_{s\mu, s\nu}^*(\mathbf{k}) = P_{s\nu, s\mu}(\mathbf{k})$, which implies $\operatorname{Im}[P_{s\mu, s\nu}(\mathbf{k})] = -\operatorname{Im}[P_{s\nu, s\mu}(\mathbf{k})]$. Thus, $\mathrm{V}_s^T = -\mathrm{V}_s$ is strictly real and skew-symmetric. $\square$

---

## 3. Pure Kitaev Self-Consistent Ground-State Correlator (K-006)

**Theorem 3.1 (Matter Correlator Invariant K-006):**
In the pure Kitaev limit ($K < 0$, $J = \Gamma = \Gamma' = D = \mathbf{h} = \boldsymbol{\lambda} = \mathbf{0}$) with locked zero-flux gauge background $\mathrm{U}_b^{\gamma\gamma} = +1.0$ and any positive matter seed $u_{\mathrm{matter}} > 0$:
1. The on-site correlation matrix vanishes identically: $\mathrm{V}_s = \mathbf{0}$, and all gauge residuals satisfy $g_s^\gamma = \operatorname{tr}(\mathrm{G}^\gamma \mathrm{V}_s) = 0$.
2. The gauge bond correlator is preserved: $\mathrm{U}_b^{\gamma\gamma} = +1.0$.
3. The self-consistent matter bond correlator evaluates to the exact two-dimensional Brillouin-zone integral (Freitas et al., 2024):
   $$\mathrm{U}_b^{00} = \int_0^1 \mathrm{d}k_1 \int_0^1 \mathrm{d}k_2 \frac{1 + \cos(2\pi k_1) + \cos(2\pi k_2)}{\sqrt{3 + 2\cos(2\pi k_1) + 2\cos(2\pi k_2) + 2\cos(2\pi(k_1 - k_2))}} \approx 0.5248645.$$

*Proof:*
From Theorem 4.1 of Note 02, the $2 \times 2$ matter sector in the 2-sublattice cell is governed by $\mathcal{H}_{A0, B0}(\mathbf{k}) = \frac{i|K|}{2} f(\mathbf{k})$, where $f(\mathbf{k}) = 1 + e^{-i2\pi k_1} + e^{-i2\pi k_2} = |f(\mathbf{k})|e^{i\theta_{\mathbf{k}}}$.
The occupied negative-energy state has eigenvector $\mathbf{w}_-(\mathbf{k}) = \frac{1}{\sqrt{2}}(1, i e^{-i\theta_{\mathbf{k}}})^T$.
The off-diagonal matter element of $P(\mathbf{k}) = 2 \mathbf{w}_-(\mathbf{k}) \mathbf{w}_-^\dagger(\mathbf{k})$ is:
$$P_{A0, B0}(\mathbf{k}) = 2 \left( \frac{1}{\sqrt{2}} \right) \left( -\frac{i}{\sqrt{2}} e^{i\theta_{\mathbf{k}}} \right) = -i e^{i\theta_{\mathbf{k}}} = -i \frac{f(\mathbf{k})}{|f(\mathbf{k})|}.$$
Multiplying by $i$ for the correlator: $i P_{A0, B0}(\mathbf{k}) = \frac{f(\mathbf{k})}{|f(\mathbf{k})|}$.
For the $z$-bond ($\Delta\mathbf{R}_z = \mathbf{0}$):
$$\mathrm{U}_z^{00} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Im}[P_{A0, B0}(\mathbf{k})] = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Re}\!\left[\frac{f(\mathbf{k})}{|f(\mathbf{k})|}\right] = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \frac{1 + \cos(2\pi k_1) + \cos(2\pi k_2)}{|f(\mathbf{k})|}.$$
For the $x$-bond ($\Delta\mathbf{R}_x = (0, -1)$, $\phi_x = -2\pi k_2$):
$$\mathrm{U}_x^{00} = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \operatorname{Re}\!\left[\frac{f(\mathbf{k}) e^{i 2\pi k_2}}{|f(\mathbf{k})|}\right] = \frac{1}{N_{\mathbf{k}}} \sum_{\mathbf{k}} \frac{\cos(2\pi k_2) + \cos(2\pi(k_1 - k_2)) + 1}{|f(\mathbf{k})|},$$
which evaluates to the identical scalar value by $C_3$ rotational symmetry. High-precision numerical quadrature over $[0, 1)^2$ yields $\mathrm{U}_b^{00} = 0.52486453\dots$, matching the zero-flux ground state in Freitas et al. (2024). $\square$

---

## 4. Local Gauge Residuals and Joint Anderson Acceleration (DIIS)

### Gauge Invariance Constraints
Local physical spin representations require the state to reside in the gauge-singlet sector $\langle c_s^T \mathrm{G}^\gamma c_s \rangle = 0$ for all sublattices $s$ and $\gamma \in \{1, 2, 3\}$.
In terms of the on-site correlation matrix $\mathrm{V}_s^{\mu\nu} = \langle i c_s^\mu c_s^\nu \rangle - i\delta^{\mu\nu}$, the gauge residual vector is:
$$g_s^\gamma \equiv \operatorname{tr}(\mathrm{G}^\gamma \mathrm{V}_s) = \sum_{\mu,\nu=0}^3 (\mathrm{G}^\gamma)_{\nu\mu} \mathrm{V}_s^{\mu\nu} = 0.$$
Because $\mathrm{G}^\gamma$ is skew-symmetric and traceless, the $-i\delta^{\mu\nu}$ offset drops out, ensuring $g_s^\gamma$ is strictly real.

### Unified State Representation and Fixed-Point Mapping
We formulate the self-consistent saddle point as a joint root-finding problem in the concatenated parameter space:
$$X = \begin{pmatrix} \operatorname{vec}(\mathrm{U}) \\ \operatorname{vec}(\mathrm{V}) \\ \operatorname{vec}(\boldsymbol{\lambda}) \end{pmatrix} \in \mathbb{R}^{D}, \qquad D = 16 N_{\mathrm{bonds}} + 16 N_{\mathrm{sub}} + 3 N_{\mathrm{sub}}.$$

Given an iterate $X^{(k)}$, the single-particle Hamiltonian $\mathcal{H}(\mathbf{k})$ is assembled and projected onto negative energies to yield updated correlators $\tilde{\mathrm{U}}(X^{(k)})$ and $\tilde{\mathrm{V}}(X^{(k)})$. The gauge residuals $g(X^{(k)})$ generate the updated multiplier targets $\tilde{\boldsymbol{\lambda}} = \boldsymbol{\lambda} - \alpha_\lambda g$.
Defining the target mapping $\Phi(X^{(k)}) = (\operatorname{vec}(\tilde{\mathrm{U}}), \operatorname{vec}(\tilde{\mathrm{V}}), \operatorname{vec}(\tilde{\boldsymbol{\lambda}}))^T$, the iteration residual is:
$$R_k \equiv \Phi(X^{(k)}) - X^{(k)}.$$

### KKT Subspace Inversion
Direct Inversion in the Iterative Subspace (DIIS; Pulay, 1980; Walker and Ni, 2011) retains a history of the last $m$ iterates $\{X^{(k-j)}\}_{j=0}^{m-1}$ and residuals $\{R_{k-j}\}_{j=0}^{m-1}$.
An optimal extrapolated state $X_{\mathrm{opt}} = \sum_{j=0}^{m-1} c_j X^{(k-j)}$ is sought that minimizes the residual norm $\|\sum_{j=0}^{m-1} c_j R_{k-j}\|_2^2$ subject to $\sum_{j=0}^{m-1} c_j = 1$.
The stationary solution satisfies the $(m+1) \times (m+1)$ linear Karush-Kuhn-Tucker system:
$$\begin{pmatrix} B & -\mathbf{1} \\ -\mathbf{1}^T & 0 \end{pmatrix} \begin{pmatrix} \mathbf{c} \\ \mu \end{pmatrix} = \begin{pmatrix} \mathbf{0} \\ -1 \end{pmatrix},$$
where the overlap matrix elements are regularized against collinearity:
$$B_{pq} = \langle R_{k-p}, R_{k-q} \rangle + \epsilon_{\mathrm{reg}} \delta_{pq} \langle R_{k-p}, R_{k-p} \rangle.$$

Once the weights $\mathbf{c} = (c_0, \dots, c_{m-1})^T$ are obtained, the accelerated next iterate is:
$$X^{(k+1)} = \sum_{j=0}^{m-1} c_j \Phi(X^{(k-j)}) = \sum_{j=0}^{m-1} c_j \left( X^{(k-j)} + R_{k-j} \right).$$

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
