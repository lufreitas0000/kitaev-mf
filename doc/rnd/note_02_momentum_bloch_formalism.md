# Note 02: Momentum-Space Bloch-Majorana Representation, Analytical Dispersion, and Brillouin-Zone Folding Parity

## 1. Context, Multi-Sublattice Graph, and Discrete Fourier Transform

Let $\Lambda = (\mathcal{V}, \mathcal{E})$ denote a periodic honeycomb lattice graph where $\mathcal{V}$ is partitioned into $N_{\mathrm{sub}}$ sublattices per unit cell, and $\mathcal{E}$ is the set of directed nearest-neighbor bonds $b = \langle ij \rangle_\gamma$ oriented from sublattice $s_i \in \mathcal{V}$ to sublattice $s_j \in \mathcal{V}$ with Bravais cell translation $\Delta\mathbf{R}_b = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2$ ($n_1, n_2 \in \mathbb{Z}$) and bond direction $\gamma \in \{1, 2, 3\} \equiv \{x, y, z\}$ (Kitaev, 2006; Freitas et al., 2024).

The real four-component Majorana operator at sublattice $s \in \{0, \dots, N_{\mathrm{sub}}-1\}$ in unit cell $\mathbf{R}$ satisfies the self-dual CAR algebra (Burnell and Nayak, 2011; Seifert et al., 2018):
$$\{ c_{s,\mathbf{R}}^\mu, c_{s',\mathbf{R}'}^\nu \} = 2 \delta_{ss'} \delta_{\mathbf{R}\mathbf{R}'} \delta^{\mu\nu} \mathbb{1}, \qquad (c_{s,\mathbf{R}}^\mu)^\dagger = c_{s,\mathbf{R}}^\mu.$$

Defining the discrete Fourier transform on a lattice of $N_{\mathrm{cell}}$ unit cells with periodic boundary conditions:
$$c_{s,\mathbf{R}} = \frac{1}{\sqrt{N_{\mathrm{cell}}}} \sum_{\mathbf{k} \in \mathrm{BZ}} c_{s,\mathbf{k}} e^{i \mathbf{k} \cdot \mathbf{R}}, \qquad c_{s,\mathbf{k}} = \frac{1}{\sqrt{N_{\mathrm{cell}}}} \sum_{\mathbf{R}} c_{s,\mathbf{R}} e^{-i \mathbf{k} \cdot \mathbf{R}},$$
the Majorana field in reciprocal space satisfies the involution relation:
$$c_{s,\mathbf{k}}^\dagger = c_{s,-\mathbf{k}}, \qquad \{ c_{s,\mathbf{k}}^\mu, (c_{s',\mathbf{k}'}^\nu)^\dagger \} = 2 \delta_{ss'} \delta^{\mu\nu} \delta_{\mathbf{k}\mathbf{k}'} \mathbb{1}.$$

---

## 2. Bloch-Majorana Matrix Assembly and Explicit $2\pi$ Phase Factors

Crystal momenta are parameterized in dimensionless fractional coordinates $(k_1, k_2) \in [0, 1)^2$ dual to the primitive Bravais lattice basis vectors $(\mathbf{a}_1, \mathbf{a}_2)$. Using the reciprocal lattice basis vectors $(\mathbf{b}_1, \mathbf{b}_2)$ defined by $\mathbf{a}_m \cdot \mathbf{b}_n = 2\pi \delta_{mn}$, the physical Cartesian wavevector is $\mathbf{k} = k_1 \mathbf{b}_1 + k_2 \mathbf{b}_2$.

For any directed bond $b = \langle ij \rangle_\gamma$ with Bravais integer offset $(n_1, n_2)$, the inner product simplifies analytically:
$$\mathbf{k} \cdot \Delta\mathbf{R}_b = (k_1 \mathbf{b}_1 + k_2 \mathbf{b}_2) \cdot (n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2) = 2\pi (k_1 n_1 + k_2 n_2).$$

**Definition 2.1 (Bloch-Majorana Hamiltonian Matrix):**
Following Proposition 5.1 of Note 01 (Freitas et al., 2024), the quadratic mean-field Hamiltonian is expressed in reciprocal space as:
$$H_{\mathrm{MF}} = \sum_{\mathbf{k} \in \mathrm{BZ}} \boldsymbol{\Psi}_{\mathbf{k}}^\dagger \mathcal{H}(\mathbf{k}) \boldsymbol{\Psi}_{\mathbf{k}} - \mathrm{C},$$
where $\boldsymbol{\Psi}_{\mathbf{k}} = (c_{0,\mathbf{k}}^T, c_{1,\mathbf{k}}^T, \dots, c_{N_{\mathrm{sub}}-1,\mathbf{k}}^T)^T$ is the $4N_{\mathrm{sub}}$-component Majorana Nambu spinor, and the single-particle Bloch matrix is:
$$\mathcal{H}(\mathbf{k}) = \frac{i}{4} \mathcal{M}(\mathbf{k}) \in \mathbb{C}^{4N_{\mathrm{sub}} \times 4N_{\mathrm{sub}}}.$$

The matrix $\mathcal{M}(\mathbf{k})$ consists of $N_{\mathrm{sub}} \times N_{\mathrm{sub}}$ blocks of dimension $4 \times 4$:
$$\mathcal{M}_{s, s'}(\mathbf{k}) = \delta_{s, s'} \mathrm{B}_s + \sum_{b \in \mathcal{E}} \left( \delta_{s, s_i(b)} \delta_{s', s_j(b)} \mathrm{A}_b e^{i 2\pi (k_1 n_{1,b} + k_2 n_{2,b})} - \delta_{s, s_j(b)} \delta_{s', s_i(b)} \mathrm{A}_b^T e^{-i 2\pi (k_1 n_{1,b} + k_2 n_{2,b})} \right),$$
where the real inter-site hopping matrices $\mathrm{A}_b$ and real skew-symmetric on-site matrices $\mathrm{B}_s$ are given by:
$$\mathrm{A}_b = 2 \sum_{\alpha,\beta=1}^3 \mathbf{J}_b^{\alpha\beta} \mathrm{N}^\alpha \mathrm{U}_b \mathrm{N}^\beta,$$
$$\mathrm{B}_s = \sum_{b \in \mathcal{E}_{s \to \mathrm{out}}} \sum_{\alpha,\beta=1}^3 \mathbf{J}_b^{\alpha\beta} \mathrm{N}^\alpha \mathrm{tr}\!\left(\mathrm{V}_{s_j(b)}^T \mathrm{N}^\beta\right) + \sum_{b \in \mathcal{E}_{\mathrm{in} \to s}} \sum_{\alpha,\beta=1}^3 \mathbf{J}_b^{\beta\alpha} \mathrm{N}^\alpha \mathrm{tr}\!\left(\mathrm{V}_{s_i(b)}^T \mathrm{N}^\beta\right) + \sum_{\gamma=1}^3 \left( \lambda_s^\gamma \mathrm{G}^\gamma - 2 h^\gamma \mathrm{N}^\gamma \right).$$

Notice that for bonds oriented into site $s$, the exchange matrix is transposed ($\mathbf{J}_{ji} = \mathbf{J}_{ij}^T$), preserving exact angular momentum balance under antisymmetric Dzyaloshinskii-Moriya interactions ($D_\gamma \neq 0$).

---

## 3. Symmetries and Spectral Properties

**Proposition 3.1 (Hermiticity):** The Bloch matrix $\mathcal{H}(\mathbf{k}) = \frac{i}{4}\mathcal{M}(\mathbf{k})$ is strictly Hermitian for all $\mathbf{k} \in \mathbb{R}^2$:
$$\mathcal{H}^\dagger(\mathbf{k}) = \mathcal{H}(\mathbf{k}).$$

*Proof:*
Because $\mathrm{B}_s^T = -\mathrm{B}_s$ is real skew-symmetric, the diagonal blocks satisfy:
$$\left(\frac{i}{4}\mathrm{B}_s\right)^\dagger = -\frac{i}{4}\mathrm{B}_s^T = -\frac{i}{4}(-\mathrm{B}_s) = \frac{i}{4}\mathrm{B}_s.$$
For off-diagonal blocks connecting $s_i \to s_j$ with phase $\phi_b = 2\pi(k_1 n_1 + k_2 n_2)$:
$$\left(\frac{i}{4}\mathrm{A}_b e^{i \phi_b}\right)^\dagger = -\frac{i}{4} \mathrm{A}_b^T e^{-i \phi_b} = \frac{i}{4} \left(-\mathrm{A}_b^T\right) e^{-i \phi_b} = \mathcal{H}_{s_j, s_i}(\mathbf{k}).$$
Thus, $\mathcal{H}^\dagger(\mathbf{k}) = \mathcal{H}(\mathbf{k})$. $\square$

**Proposition 3.2 (Particle-Hole / Skew-Conjugation Symmetry):** In the Majorana basis, the Bloch matrix satisfies:
$$\mathcal{H}(-\mathbf{k}) = -\mathcal{H}^*(\mathbf{k}) = -\mathcal{H}^T(\mathbf{k}).$$
Consequently, the single-particle excitation spectrum occurs in symmetric pairs $\pm E_n(\mathbf{k})$.

*Proof:*
Taking the complex conjugate of $\mathcal{H}(\mathbf{k})$:
$$\mathcal{H}^*(\mathbf{k}) = -\frac{i}{4} \mathcal{M}^*(\mathbf{k}).$$
Because $\mathrm{B}_s$ and $\mathrm{A}_b$ are purely real matrices, $\mathcal{M}^*(\mathbf{k})$ conjugates only the phase factors $e^{\pm i \phi_b} \to e^{\mp i \phi_b}$, which evaluates identically to $\mathcal{M}(-\mathbf{k})$. Therefore, $\mathcal{H}^*(\mathbf{k}) = -\frac{i}{4}\mathcal{M}(-\mathbf{k}) = -\mathcal{H}(-\mathbf{k})$. $\square$

---

## 4. Analytical Pure-Kitaev Dispersion and Flat Gauge Bands

**Theorem 4.1 (Exact Dispersion in the Pure Kitaev Limit):**
In the 2-sublattice primitive cell ($N_{\mathrm{sub}}=2$) with isotropic ferromagnetic Kitaev coupling $K < 0$, $J = \Gamma = \Gamma' = D = \mathbf{h} = \boldsymbol{\lambda} = \mathbf{0}$, directed gauge link variables $\mathrm{U}_b^{\gamma\gamma} = u_{\mathrm{gauge}} = +1.0$, and uniform matter correlator $\mathrm{U}_b^{00} = u_{\mathrm{matter}} \ge 0$:
1. The two itinerant matter Majorana bands ($c^0$) possess the exact dispersion:
   $$E_{c^0}(\mathbf{k}) = \pm \frac{\vert{}K\vert{}}{2} \vert{}1 + e^{-i 2\pi k_1} + e^{-i 2\pi k_2}\vert{} = \pm \frac{\vert{}K\vert{}}{2} \sqrt{3 + 2\cos(2\pi k_1) + 2\cos(2\pi k_2) + 2\cos(2\pi(k_1 - k_2))}.$$
   At the zone center $\mathbf{k} = (0, 0)$, $E_{c^0}(\mathbf{0}) = 1.5 \vert{}K\vert{}$, reproducing Fig. 1 of Freitas et al. (2024).
2. The six gauge Majorana bands ($c^1, c^2, c^3$) are completely dispersionless (flat) across the entire Brillouin zone:
   $$E_{c^\gamma}(\mathbf{k}) = \pm \frac{\vert{}K\vert{}}{2} \vert{}u_{\mathrm{matter}}\vert{}.$$
   When $u_{\mathrm{matter}} = 0$, all six gauge bands reside at zero energy ($E = 0$).

*Proof:*
Evaluating $(\mathrm{N}^\gamma \mathrm{U}_b \mathrm{N}^\gamma)_{\mu\nu}$ with $(\mathrm{N}^\gamma)_{\mu\nu} = \delta_{\mu 0}\delta_{\nu\gamma} - \delta_{\nu 0}\delta_{\mu\gamma}$:
$$(\mathrm{N}^\gamma \mathrm{U}_b \mathrm{N}^\gamma)_{00} = (\mathrm{N}^\gamma)_{0\gamma} \mathrm{U}_b^{\gamma\gamma} (\mathrm{N}^\gamma)_{\gamma 0} = (+1)(u_{\mathrm{gauge}})(-1) = -u_{\mathrm{gauge}} = -1.0,$$
$$(\mathrm{N}^\gamma \mathrm{U}_b \mathrm{N}^\gamma)_{\gamma\gamma} = (\mathrm{N}^\gamma)_{\gamma 0} \mathrm{U}_b^{00} (\mathrm{N}^\gamma)_{0\gamma} = (-1)(u_{\mathrm{matter}})(+1) = -u_{\mathrm{matter}}.$$
All other matrix elements vanish. For $K < 0$, the non-zero elements of $\mathrm{A}_b = 2 K \mathrm{N}^\gamma \mathrm{U}_b \mathrm{N}^\gamma$ are:
$$\mathrm{A}_b^{00} = 2 \vert{}K\vert{} u_{\mathrm{gauge}} = 2 \vert{}K\vert{}, \qquad \mathrm{A}_b^{\gamma\gamma} = 2 \vert{}K\vert{} u_{\mathrm{matter}}.$$
The three directed bonds of `make_honeycomb_2site()` are:
- $z$-bond ($\gamma=3$): $\Delta\mathbf{R}_z = (0, 0)$
- $x$-bond ($\gamma=1$): $\Delta\mathbf{R}_x = (0, -1) \implies \phi_x = -2\pi k_2$
- $y$-bond ($\gamma=2$): $\Delta\mathbf{R}_y = (-1, 0) \implies \phi_y = -2\pi k_1$

In the $c^0$ sector, the off-diagonal $A \to B$ element of $\mathcal{H}(\mathbf{k})$ is:
$$\mathcal{H}_{A0, B0}(\mathbf{k}) = \frac{i}{4} \left( 2\vert{}K\vert{} + 2\vert{}K\vert{}e^{-i 2\pi k_2} + 2\vert{}K\vert{}e^{-i 2\pi k_1} \right) = \frac{i\vert{}K\vert{}}{2} \left( 1 + e^{-i 2\pi k_1} + e^{-i 2\pi k_2} \right).$$
Diagonalizing the $2 \times 2$ block $\begin{pmatrix} 0 & \mathcal{H}_{A0,B0} \\ \mathcal{H}_{A0,B0}^* & 0 \end{pmatrix}$ yields $E_{c^0}(\mathbf{k}) = \pm \vert{}\mathcal{H}_{A0,B0}(\mathbf{k})\vert{} = \pm \frac{\vert{}K\vert{}}{2}\vert{}1 + e^{-i 2\pi k_1} + e^{-i 2\pi k_2}\vert{}$.
For each gauge flavor $\gamma$, hopping occurs exclusively on the single bond of direction $\gamma$, which carries no interference from other bonds. Thus, $\mathcal{H}_{A\gamma, B\gamma}(\mathbf{k}) = \frac{i}{4} (2\vert{}K\vert{} u_{\mathrm{matter}}) e^{i \phi_\gamma}$, whose eigenvalues are $\pm \frac{\vert{}K\vert{}}{2}\vert{}u_{\mathrm{matter}}\vert{}$, independent of $\mathbf{k}$. $\square$

---

## 5. Brillouin-Zone Folding Parity Theorem

**Theorem 5.1 (Supercell Spectral Invariance):**
Let $(K_1, K_2) \in [0, 1)^2$ denote fractional reciprocal coordinates in the 4-sublattice rectangular cell with lattice vectors $\mathbf{A}_1 = \mathbf{a}_1 - \mathbf{a}_2$ and $\mathbf{A}_2 = \mathbf{a}_1 + \mathbf{a}_2$. The 16 eigenvalues of $\mathcal{H}_4(K_1, K_2)$ match the union of the 8 eigenvalues of $\mathcal{H}_2(\mathbf{k}^{(1)})$ and the 8 eigenvalues of $\mathcal{H}_2(\mathbf{k}^{(2)})$ to machine precision ($< 10^{-14}$), where the two folded primitive momenta are:
$$\mathbf{k}^{(1)} = \left( \frac{K_1 + K_2}{2} \bmod 1, \; \frac{K_2 - K_1}{2} \bmod 1 \right),$$
$$\mathbf{k}^{(2)} = \left( \frac{K_1 + K_2 + 1}{2} \bmod 1, \; \frac{K_2 - K_1 + 1}{2} \bmod 1 \right).$$

---

## References

* **Burnell, F. J. and Nayak, C.** (2011). *SU(2) slave fermion approach to the Kitaev honeycomb lattice model*. Physical Review B **84**, 125125. DOI: [10.1103/PhysRevB.84.125125](https://doi.org/10.1103/PhysRevB.84.125125). arXiv: [1104.5485](https://arxiv.org/abs/1104.5485).
* **Freitas, L. R. D., Bauer, T., Egger, R., and Pereira, R. G.** (2024). *Electric polarization near vortices in the extended Kitaev model*. npj Quantum Materials **9**, 33. DOI: [10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5). arXiv: [2311.15750](https://arxiv.org/abs/2311.15750).
* **Kitaev, A.** (2006). *Anyons in an exactly solved model and beyond*. Annals of Physics **321**, 2–111. DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005). arXiv: [cond-mat/0506438](https://arxiv.org/abs/cond-mat/0506438).
* **Seifert, U. F. P., Meng, T., and Vojta, M.** (2018). *Fractionalized Fermi liquids and exotic superconductivity in the Kitaev-Kondo lattice*. Physical Review B **97**, 085118. DOI: [10.1103/PhysRevB.97.085118](https://doi.org/10.1103/PhysRevB.97.085118). arXiv: [1710.05018](https://arxiv.org/abs/1710.05018).
