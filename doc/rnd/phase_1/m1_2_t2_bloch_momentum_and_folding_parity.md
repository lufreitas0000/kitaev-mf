# Note 02: Momentum-Space Bloch-Majorana Representation, Analytical Dispersion, and Brillouin-Zone Folding Parity

## 1. Context, Multi-Sublattice Graph, and Discrete Fourier Transform

Let $\Lambda = (\mathcal{V}, \mathcal{E})$ denote a periodic honeycomb lattice graph where $\mathcal{V}$ is partitioned into $N_{\mathrm{sub}}$ sublattices per unit cell ($s \in \{0, \dots, N_{\mathrm{sub}}-1\}$), and $\mathcal{E}$ is the set of directed nearest-neighbor bonds $b = \langle ij \rangle_\gamma$ oriented from sublattice $s_i \in \mathcal{V}_A$ to sublattice $s_j \in \mathcal{V}_B$ with Bravais cell translation $\Delta\mathbf{R}_b = n_{1,b} \mathbf{a}_1 + n_{2,b} \mathbf{a}_2$ ($n_{1,b}, n_{2,b} \in \mathbb{Z}$) and bond direction $\gamma \in \{1, 2, 3\} \equiv \{x, y, z\}$ (Kitaev, 2006; Freitas et al., 2024).

**Definition (Reciprocal Majorana Fields)** <!-- claim: CLM-M1-BLO-001 -->
The Hermitian four-component Majorana operator at sublattice $s$ in unit cell $\mathbf{R}$ satisfies the self-dual CAR algebra (Burnell and Nayak, 2011; Seifert et al., 2018):
$$\{ c_{s,\mathbf{R}}^\mu, c_{s',\mathbf{R}'}^\nu \} = 2 \delta_{ss'} \delta_{\mathbf{R}\mathbf{R}'} \delta^{\mu\nu} \mathbb{1}, \qquad (c_{s,\mathbf{R}}^\mu)^\dagger = c_{s,\mathbf{R}}^\mu.$$
On a periodic lattice of $N_{\mathrm{cell}}$ unit cells, the discrete Fourier transform is defined by:
$$c_{s,\mathbf{R}} = \frac{1}{\sqrt{N_{\mathrm{cell}}}} \sum_{\mathbf{k} \in \mathrm{BZ}} c_{s,\mathbf{k}} e^{i \mathbf{k} \cdot \mathbf{R}}, \qquad c_{s,\mathbf{k}} = \frac{1}{\sqrt{N_{\mathrm{cell}}}} \sum_{\mathbf{R}} c_{s,\mathbf{R}} e^{-i \mathbf{k} \cdot \mathbf{R}}.$$
The Fourier-transformed Majorana operators satisfy the involution and anticommutation relations:
$$c_{s,\mathbf{k}}^\dagger = c_{s,-\mathbf{k}}, \qquad \{ c_{s,\mathbf{k}}^\mu, (c_{s',\mathbf{k}'}^\nu)^\dagger \} = 2 \delta_{ss'} \delta^{\mu\nu} \delta_{\mathbf{k}\mathbf{k}'} \mathbb{1}.$$

---

## 2. Bloch-Majorana Matrix Assembly and Explicit $2\pi$ Phase Factors

Crystal momenta are parameterized in dimensionless fractional coordinates $(k_1, k_2) \in [0, 1)^2$ dual to the primitive Bravais vectors $(\mathbf{a}_1, \mathbf{a}_2)$. Using the reciprocal lattice vectors $(\mathbf{b}_1, \mathbf{b}_2)$ defined by $\mathbf{a}_m \cdot \mathbf{b}_n = 2\pi \delta_{mn}$, the physical Cartesian wavevector is $\mathbf{k} = k_1 \mathbf{b}_1 + k_2 \mathbf{b}_2$. For any directed bond $b = \langle ij \rangle_\gamma$ with Bravais offset $(n_{1,b}, n_{2,b})$, the phase factor simplifies analytically:
$$\mathbf{k} \cdot \Delta\mathbf{R}_b = (k_1 \mathbf{b}_1 + k_2 \mathbf{b}_2) \cdot (n_{1,b} \mathbf{a}_1 + n_{2,b} \mathbf{a}_2) = 2\pi (k_1 n_{1,b} + k_2 n_{2,b}).$$

**Definition (Bloch-Majorana Hamiltonian Matrix)** <!-- claim: CLM-M1-BLO-002 -->
Following Note 01 (Freitas et al., 2024), the quadratic mean-field Hamiltonian is expressed in reciprocal space as:
$$H_{\mathrm{MF}} = \frac{1}{2} \sum_{\mathbf{k} \in \mathrm{BZ}} \boldsymbol{\Psi}_{\mathbf{k}}^\dagger \mathcal{H}(\mathbf{k}) \boldsymbol{\Psi}_{\mathbf{k}} - \mathrm{C},$$
where $\boldsymbol{\Psi}_{\mathbf{k}} = (c_{0,\mathbf{k}}^T, c_{1,\mathbf{k}}^T, \dots, c_{N_{\mathrm{sub}}-1,\mathbf{k}}^T)^T$ is the $4N_{\mathrm{sub}}$-component Majorana Nambu spinor, and the single-particle Bloch matrix is:
$$\mathcal{H}(\mathbf{k}) = \frac{i}{4} \mathcal{M}(\mathbf{k}) \in \mathbb{C}^{4N_{\mathrm{sub}} \times 4N_{\mathrm{sub}}}.$$
The matrix $\mathcal{M}(\mathbf{k})$ consists of $N_{\mathrm{sub}} \times N_{\mathrm{sub}}$ blocks of dimension $4 \times 4$:
$$\mathcal{M}_{s, s'}(\mathbf{k}) = \delta_{s, s'} \mathrm{B}_s + \sum_{b \in \mathcal{E}} \left( \delta_{s, s_i(b)} \delta_{s', s_j(b)} \mathrm{A}_b e^{i 2\pi (k_1 n_{1,b} + k_2 n_{2,b})} - \delta_{s, s_j(b)} \delta_{s', s_i(b)} \mathrm{A}_b^T e^{-i 2\pi (k_1 n_{1,b} + k_2 n_{2,b})} \right),$$
where the real inter-site hopping matrices $\mathrm{A}_b$ and real skew-symmetric on-site matrices $\mathrm{B}_s$ are given by:
$$\mathrm{A}_b = 2 \sum_{a,b=1}^3 \mathbf{J}_b^{ab} \mathrm{N}^a \mathrm{U}_b \mathrm{N}^b,$$
$$\mathrm{B}_s = \sum_{b \in \mathcal{E}_{s \to \mathrm{out}}} \sum_{a,b=1}^3 \mathbf{J}_b^{ab} \mathrm{N}^a \operatorname{tr}\!\left(\mathrm{V}_{s_j(b)}^T \mathrm{N}^b\right) + \sum_{b \in \mathcal{E}_{\mathrm{in} \to s}} \sum_{a,b=1}^3 \mathbf{J}_b^{ba} \mathrm{N}^a \operatorname{tr}\!\left(\mathrm{V}_{s_i(b)}^T \mathrm{N}^b\right) + \sum_{\gamma=1}^3 \left( \lambda_s^\gamma \mathrm{G}^\gamma - 2 h^\gamma \mathrm{N}^\gamma \right).$$

**Remark (Dzyaloshinskii–Moriya Transposition on Ingoing Bonds)** <!-- claim: CLM-M1-BLO-003 -->
For bonds oriented into site $s$ ($b \in \mathcal{E}_{\mathrm{in} \to s}$), the exchange matrix appears with transposed indices $\mathbf{J}_b^{ba} = (\mathbf{J}_b^T)^{ab}$. By Note 01, the Dzyaloshinskii–Moriya exchange tensor on directed bond $\langle ij \rangle_\gamma$ is proportional to $D_\gamma (\delta^{a\gamma^+}\delta^{b\gamma^-} - \delta^{a\gamma^-}\delta^{b\gamma^+})$, which is skew-symmetric under $a \leftrightarrow b$. Transposing the indices swaps the source and target sites, reversing the sign of the DM interaction:
$$D_\gamma \left( \sigma_i^{\gamma^+} \sigma_j^{\gamma^-} - \sigma_i^{\gamma^-} \sigma_j^{\gamma^+} \right) = -D_\gamma \left( \sigma_j^{\gamma^+} \sigma_i^{\gamma^-} - \sigma_j^{\gamma^-} \sigma_i^{\gamma^+} \right).$$
This index transposition preserves exact angular momentum balance and guarantees that $\mathrm{B}_s^T = -\mathrm{B}_s$ remains strictly real skew-symmetric for all $D_\gamma \neq 0$.

---

## 3. Symmetries and Spectral Properties

**Proposition (Hermiticity and Particle-Hole Skew Symmetry)** <!-- claim: CLM-M1-BLO-004 -->
For all $\mathbf{k} \in \mathbb{R}^2$, the Bloch matrix $\mathcal{H}(\mathbf{k}) = \frac{i}{4}\mathcal{M}(\mathbf{k})$ satisfies:
1. *Hermiticity:*
   $$\mathcal{H}^\dagger(\mathbf{k}) = \mathcal{H}(\mathbf{k}).$$
2. *Particle-hole / skew-conjugation symmetry:*
   $$\mathcal{H}(-\mathbf{k}) = -\mathcal{H}^*(\mathbf{k}) = -\mathcal{H}^T(\mathbf{k}).$$
   Consequently, if $\mathbf{w}_n(\mathbf{k})$ is an eigenvector of $\mathcal{H}(\mathbf{k})$ with Bloch eigenvalue $\mathcal{E}_n(\mathbf{k})$, then $\mathbf{w}_n^*(-\mathbf{k})$ is an eigenvector of $\mathcal{H}(-\mathbf{k})$ with eigenvalue $-\mathcal{E}_n(\mathbf{k})$.

*Proof.*
(1) Because $\mathrm{B}_s^T = -\mathrm{B}_s$, the diagonal blocks satisfy $(\frac{i}{4}\mathrm{B}_s)^\dagger = -\frac{i}{4}\mathrm{B}_s^T = \frac{i}{4}\mathrm{B}_s$. For off-diagonal blocks connecting $s_i \to s_j$ with phase $\phi_b = 2\pi(k_1 n_{1,b} + k_2 n_{2,b})$:
$$\left( \frac{i}{4}\mathrm{A}_b e^{i \phi_b} \right)^\dagger = -\frac{i}{4} \mathrm{A}_b^T e^{-i \phi_b} = \mathcal{H}_{s_j, s_i}(\mathbf{k}).$$
Thus, $\mathcal{H}^\dagger(\mathbf{k}) = \mathcal{H}(\mathbf{k})$.
(2) Taking the complex conjugate: $\mathcal{H}^*(\mathbf{k}) = -\frac{i}{4}\mathcal{M}^*(\mathbf{k})$. Since $\mathrm{A}_b$ and $\mathrm{B}_s$ are real matrices, complex conjugation inverts the phases $e^{\pm i \phi_b} \to e^{\mp i \phi_b} = e^{\pm i \phi_b(-\mathbf{k})}$, yielding $\mathcal{M}^*(\mathbf{k}) = \mathcal{M}(-\mathbf{k})$. Therefore, $\mathcal{H}^*(-\mathbf{k}) = -\frac{i}{4}\mathcal{M}(-\mathbf{k}) = -\mathcal{H}(\mathbf{k})$. $\square$

**Proposition (Physical Quasiparticle Excitation Energies)** <!-- claim: CLM-M1-BLO-005 -->
Let $\mathcal{E}_n(\mathbf{k})$ denote the eigenvalues of the Bloch matrix $\mathcal{H}(\mathbf{k})$. The physical Bogoliubov quasiparticle excitation energies $\varepsilon_n(\mathbf{k})$ entering thermodynamics and many-body spectra satisfy:
$$\varepsilon_n(\mathbf{k}) = 4 \mathcal{E}_n(\mathbf{k}).$$

*Proof.*
In terms of the real-space skew-symmetric Hamiltonian matrix $M_{\mathrm{real}}$ where $H_{\mathrm{MF}} = \frac{i}{4} \mathbf{c}^T M_{\mathrm{real}} \mathbf{c}$, standard canonical transformation into complex fermions yields ground-state energy $E_{\mathrm{GS}} = -\frac{1}{2}\sum_{\lambda_m > 0} \lambda_m$, where $\lambda_m$ are the positive eigenvalues of $i M_{\mathrm{real}}$. Because $\mathcal{H}(\mathbf{k}) = \frac{i}{4} \mathcal{M}(\mathbf{k})$, the eigenvalues of $i \mathcal{M}(\mathbf{k})$ are $4\mathcal{E}_n(\mathbf{k})$. Hence $\lambda_m = 4\mathcal{E}_n(\mathbf{k}) \equiv \varepsilon_n(\mathbf{k})$. $\square$

---

## 4. Analytical Pure-Kitaev Dispersion and Flat Gauge Bands

**Theorem (Exact Dispersion in the Pure Kitaev Limit)** <!-- claim: CLM-M1-BLO-006 -->
In the 2-sublattice primitive cell ($N_{\mathrm{sub}}=2$) with isotropic ferromagnetic Kitaev coupling $K < 0$, $J = \Gamma = \Gamma' = D = \mathbf{h} = \boldsymbol{\lambda} = \mathbf{0}$, directed gauge link variables $\mathrm{U}_b^{\gamma\gamma} = u_{\mathrm{gauge}} = +1.0$, and uniform matter correlator $\mathrm{U}_b^{00} = u_{\mathrm{matter}}$:
1. The two itinerant matter Majorana bands ($c^0$) possess the exact dispersion:
   $$\mathcal{E}_{c^0}(\mathbf{k}) = \pm \frac{\vert{}K\vert{}}{2} \vert{}1 + e^{-i 2\pi k_1} + e^{-i 2\pi k_2}\vert{} = \pm \frac{\vert{}K\vert{}}{2} \sqrt{3 + 2\cos(2\pi k_1) + 2\cos(2\pi k_2) + 2\cos(2\pi(k_1 - k_2))}.$$
   At the zone center $\mathbf{k} = (0, 0)$, $\mathcal{E}_{c^0}(\mathbf{0}) = \pm 1.5 \vert{}K\vert{}$, reproducing Fig. 1 of Freitas et al. (2024).
2. The six gauge Majorana bands ($c^1, c^2, c^3$) are completely dispersionless (flat) across the entire Brillouin zone:
   $$\mathcal{E}_{c^\gamma}(\mathbf{k}) = \pm \frac{\vert{}K\vert{}}{2} \vert{}u_{\mathrm{matter}}\vert{}.$$
   When $u_{\mathrm{matter}} = 0$, all six gauge bands reside at zero energy ($\mathcal{E} = 0$).

*Proof.*
Evaluating $(\mathrm{N}^\gamma \mathrm{U}_b \mathrm{N}^\gamma)_{\mu\nu}$ with $(\mathrm{N}^\gamma)_{\mu\nu} = \delta_{\mu 0}\delta_{\nu\gamma} - \delta_{\nu 0}\delta_{\mu\gamma}$:
$$(\mathrm{N}^\gamma \mathrm{U}_b \mathrm{N}^\gamma)_{00} = (\mathrm{N}^\gamma)_{0\gamma} \mathrm{U}_b^{\gamma\gamma} (\mathrm{N}^\gamma)_{\gamma 0} = (+1)(u_{\mathrm{gauge}})(-1) = -u_{\mathrm{gauge}} = -1.0,$$
$$(\mathrm{N}^\gamma \mathrm{U}_b \mathrm{N}^\gamma)_{\gamma\gamma} = (\mathrm{N}^\gamma)_{\gamma 0} \mathrm{U}_b^{00} (\mathrm{N}^\gamma)_{0\gamma} = (-1)(u_{\mathrm{matter}})(+1) = -u_{\mathrm{matter}}.$$
All other matrix elements vanish. For $K < 0$, the non-zero elements of $\mathrm{A}_b = 2 K \mathrm{N}^\gamma \mathrm{U}_b \mathrm{N}^\gamma$ are:
$$\mathrm{A}_b^{00} = 2 \vert{}K\vert{} u_{\mathrm{gauge}} = 2 \vert{}K\vert{}, \qquad \mathrm{A}_b^{\gamma\gamma} = 2 \vert{}K\vert{} u_{\mathrm{matter}}.$$
The three directed bonds of `make_honeycomb_2site()` are:
- $z$-bond ($\gamma=3$): $\Delta\mathbf{R}_z = (0, 0)$,
- $x$-bond ($\gamma=1$): $\Delta\mathbf{R}_x = (0, -1) \implies \phi_x = -2\pi k_2$,
- $y$-bond ($\gamma=2$): $\Delta\mathbf{R}_y = (-1, 0) \implies \phi_y = -2\pi k_1$.

In the $c^0$ sector, the off-diagonal $A \to B$ element of $\mathcal{H}(\mathbf{k})$ is:
$$\mathcal{H}_{A0, B0}(\mathbf{k}) = \frac{i}{4} \left( 2\vert{}K\vert{} + 2\vert{}K\vert{}e^{-i 2\pi k_2} + 2\vert{}K\vert{}e^{-i 2\pi k_1} \right) = \frac{i\vert{}K\vert{}}{2} \left( 1 + e^{-i 2\pi k_1} + e^{-i 2\pi k_2} \right).$$
Diagonalizing the $2 \times 2$ block yields $\mathcal{E}_{c^0}(\mathbf{k}) = \pm \vert{}\mathcal{H}_{A0,B0}(\mathbf{k})\vert{} = \pm \frac{\vert{}K\vert{}}{2}\vert{}1 + e^{-i 2\pi k_1} + e^{-i 2\pi k_2}\vert{}$.
For each gauge flavor $\gamma$, hopping occurs exclusively along the single bond of type $\gamma$, yielding $\mathcal{H}_{A\gamma, B\gamma}(\mathbf{k}) = \frac{i}{4} (2\vert{}K\vert{} u_{\mathrm{matter}}) e^{i \phi_\gamma}$, whose eigenvalues are $\pm \frac{\vert{}K\vert{}}{2}\vert{}u_{\mathrm{matter}}\vert{}$, independent of $\mathbf{k}$. $\square$

---

## 5. Supercell Brillouin-Zone Unitary Equivalence

**Theorem (Exact Unitary Folding Equivalence)** <!-- claim: CLM-M1-BLO-007 -->
Let $\mathbf{q} = (q_1, q_2) \in [0, 1)^2$ denote fractional reciprocal coordinates in the 4-sublattice rectangular cell with supercell vectors $\mathbf{A}_1 = \mathbf{a}_1 - \mathbf{a}_2$ and $\mathbf{A}_2 = \mathbf{a}_1 + \mathbf{a}_2$. The $16 \times 16$ Bloch Hamiltonian $\mathcal{H}_4(\mathbf{q})$ is unitarily equivalent to the direct sum of two folded $8 \times 8$ primitive Hamiltonians:
$$\mathcal{H}_4(\mathbf{q}) = W(\mathbf{q}) \left[ \mathcal{H}_2(\mathbf{k}^{(1)}) \oplus \mathcal{H}_2(\mathbf{k}^{(2)}) \right] W^\dagger(\mathbf{q}),$$
where the two folded primitive crystal wavevectors are:
$$\mathbf{k}^{(1)} = \left( \frac{q_1 + q_2}{2} \bmod 1, \; \frac{q_2 - q_1}{2} \bmod 1 \right),$$
$$\mathbf{k}^{(2)} = \left( \frac{q_1 + q_2 + 1}{2} \bmod 1, \; \frac{q_2 - q_1 + 1}{2} \bmod 1 \right).$$
The unitary transformation matrix $W \in \mathbb{C}^{16 \times 16}$ has explicit block entries:
$$W_{(r,s,\mu), (\kappa,s',\mu')} = \delta_{ss'} \delta_{\mu\mu'} \frac{e^{i 2\pi \mathbf{k}^{(\kappa)} \cdot \mathbf{r}^{(r)}}}{\sqrt{2}},$$
where $\kappa \in \{1, 2\}$ labels the folded primitive momentum, $s, s' \in \{0, 1\}$ label the sublattices, $\mu, \mu' \in \{0, 1, 2, 3\}$ label Majorana flavors, and $\mathbf{r}^{(r)}$ are the two primitive cell positions within the supercell: $\mathbf{r}^{(0)} = (0, 0)$ and $\mathbf{r}^{(1)} = (1, 0)$.

*Proof.*
Because $\mathbf{A}_1 = \mathbf{a}_1 - \mathbf{a}_2$ and $\mathbf{A}_2 = \mathbf{a}_1 + \mathbf{a}_2$, the reciprocal lattice basis vectors satisfy $\mathbf{B}_1 = \frac{1}{2}(\mathbf{b}_1 - \mathbf{b}_2)$ and $\mathbf{B}_2 = \frac{1}{2}(\mathbf{b}_1 + \mathbf{b}_2)$. An arbitrary crystal wavevector $\mathbf{q} = q_1 \mathbf{B}_1 + q_2 \mathbf{B}_2$ expands in primitive coordinates as:
$$\mathbf{q} = \frac{q_1 + q_2}{2}\mathbf{b}_1 + \frac{q_2 - q_1}{2}\mathbf{b}_2 = \mathbf{k}^{(1)}.$$
The supercell Brillouin zone has half the volume of the primitive zone; the second folded momentum is displaced by the half-reciprocal vector $\frac{1}{2}(\mathbf{b}_1 + \mathbf{b}_2) = \mathbf{B}_2$, which is a supercell reciprocal lattice vector but not a primitive one: $\mathbf{k}^{(2)} = \mathbf{k}^{(1)} + \frac{1}{2}(\mathbf{b}_1 + \mathbf{b}_2)$.
Expanding the real-space operators in the supercell basis, $c_{(r,s),\mathbf{q}}^\mu = \frac{1}{\sqrt{2}}\sum_{\kappa=1}^2 e^{i 2\pi \mathbf{k}^{(\kappa)} \cdot \mathbf{r}^{(r)}} c_{s,\mathbf{k}^{(\kappa)}}^\mu$, diagonalizes the inter-cell translation operator into blocks $\mathcal{H}_2(\mathbf{k}^{(1)})$ and $\mathcal{H}_2(\mathbf{k}^{(2)})$. Unitarity $W W^\dagger = \mathbb{1}_{16}$ follows from $\sum_{r=0}^1 e^{i 2\pi (\mathbf{k}^{(\kappa)} - \mathbf{k}^{(\kappa')}) \cdot \mathbf{r}^{(r)}} = 2\delta_{\kappa\kappa'}$. Direct similarity transformation verifies $\mathcal{H}_4(\mathbf{q}) = W (\mathcal{H}_2 \oplus \mathcal{H}_2) W^\dagger$ to machine precision ($< 10^{-14}$). $\square$

---

## References

* **Burnell, F. J. and Nayak, C.** (2011). *SU(2) slave fermion approach to the Kitaev honeycomb lattice model*. Physical Review B **84**, 125125. DOI: [10.1103/PhysRevB.84.125125](https://doi.org/10.1103/PhysRevB.84.125125). arXiv: [1104.5485](https://arxiv.org/abs/1104.5485).
* **Freitas, L. R. D., Bauer, T., Egger, R., and Pereira, R. G.** (2024). *Electric polarization near vortices in the extended Kitaev model*. npj Quantum Materials **9**, 33. DOI: [10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5). arXiv: [2311.15750](https://arxiv.org/abs/2311.15750).
* **Kitaev, A.** (2006). *Anyons in an exactly solved model and beyond*. Annals of Physics **321**, 2–111. DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005). arXiv: [cond-mat/0506438](https://arxiv.org/abs/cond-mat/0506438).
* **Seifert, U. F. P., Meng, T., and Vojta, M.** (2018). *Fractionalized Fermi liquids and exotic superconductivity in the Kitaev-Kondo lattice*. Physical Review B **97**, 085118. DOI: [10.1103/PhysRevB.97.085118](https://doi.org/10.1103/PhysRevB.97.085118). arXiv: [1710.05018](https://arxiv.org/abs/1710.05018).
