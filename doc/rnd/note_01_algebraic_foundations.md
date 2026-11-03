# Note 01: Algebraic Foundations of the Majorana Parton Representation, $\mathfrak{so}(4)$ Gauge Structure, and Mean-Field Decoupling

## 1. Context, Lattice Graph, and CAR $C^*$-Algebra

Let $\Lambda = (\mathcal{V}, \mathcal{E})$ denote a tricoordinated honeycomb lattice graph with vertex set $\mathcal{V}$ of cardinality $N_{\mathrm{sites}} = |\mathcal{V}|$ and directed edge set $\mathcal{E}$ (Kitaev, 2006; Freitas et al., 2024). Vertices are indexed by Latin letters $i, j, k, l \in \mathcal{V}$, and $\mathcal{V}_i = \{j \in \mathcal{V} : \langle ij \rangle_\gamma \in \mathcal{E} \text{ or } \langle ji \rangle_\gamma \in \mathcal{E}\}$ denotes the three nearest neighbors of vertex $i$. Each directed nearest-neighbor bond $\langle ij \rangle_\gamma \in \mathcal{E}$ is oriented from sublattice $A$ ($i$) to sublattice $B$ ($j$) along bond direction $\gamma \in \{1, 2, 3\} \equiv \{x, y, z\}$.

**Definition 1.1 (Self-Dual CAR Algebra):** To each vertex $i \in \mathcal{V}$, we assign four real Majorana fermion operators $c_i = (c_i^0, c_i^1, c_i^2, c_i^3)^T$ generating the unital self-dual Canonical Anticommutation Relations (CAR) $C^*$-algebra $\mathcal{A}_{\mathrm{CAR}}$ satisfying (Kitaev, 2006; Fu et al., 2018; Freitas et al., 2024):
$$(c_i^\mu)^\dagger = c_i^\mu, \qquad \{ c_i^\mu, c_j^\nu \} = 2 \delta_{ij} \delta^{\mu\nu} \mathbb{1},$$
for all $i, j \in \mathcal{V}$ and flavor indices $\mu, \nu \in \{0, 1, 2, 3\}$ (contrasted with Cartesian indices $\alpha, \beta, \gamma \in \{1, 2, 3\}$).

Defining two complex Dirac fermion annihilation operators per vertex, $d_{i, 1} = \frac{1}{2}(c_i^0 + i c_i^1)$ and $d_{i, 2} = \frac{1}{2}(c_i^2 + i c_i^3)$ (Burnell and Nayak, 2011; You et al., 2012), the local unprojected Fock space $\mathcal{H}_i \cong \mathbb{C}^4$ has dimension $4$, doubling the dimension of the physical pseudospin-$1/2$ Hilbert space $\mathcal{H}_{\mathrm{phys}, i} \cong \mathbb{C}^2$.

---

## 2. Lie Algebra $\mathfrak{so}(4) \cong \mathfrak{su}(2)_{\mathrm{spin}} \oplus \mathfrak{su}(2)_{\mathrm{gauge}}$

**Lemma 2.1 (Quadratic Bilinear Commutator):** Let $\mathcal{O}_A = \frac{i}{4} c_i^T A c_i$ and $\mathcal{O}_B = \frac{i}{4} c_i^T B c_i$ be on-site quadratic forms with real skew-symmetric matrices $A, B \in \mathfrak{so}(4)$ ($A^T = -A$). By the CAR anticommutation relations, the operator commutator is a Lie algebra homomorphism (Seifert et al., 2018; Freitas, 2024):
$$[\mathcal{O}_A, \mathcal{O}_B] = \frac{i}{4} c_i^T [A, B] c_i.$$

**Definition 2.2 (Kronecker Basis of $\mathfrak{so}(4)$):** Using the real $2 \times 2$ matrices
$$\epsilon = i \sigma_y = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}, \quad I_2 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}, \quad \sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad \sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}, \quad \Pi_+ = \frac{1}{2}(I_2 + \sigma_z) = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix},$$
we define the anti-self-dual gauge generators $\mathrm{G}^\gamma$ and self-dual chiral spin generators $\mathrm{M}^\gamma$ ('t Hooft, 1976; Fu et al., 2018; Seifert et al., 2018; Freitas, 2024):
$$\mathrm{G}^1 = I_2 \otimes \epsilon, \qquad \mathrm{G}^2 = \epsilon \otimes \sigma_z, \qquad \mathrm{G}^3 = \epsilon \otimes \sigma_x,$$
$$\mathrm{M}^1 = \sigma_z \otimes \epsilon, \qquad \mathrm{M}^2 = \epsilon \otimes I_2, \qquad \mathrm{M}^3 = \sigma_x \otimes \epsilon.$$

Explicitly, in the Majorana flavor basis $(c_i^0, c_i^1, c_i^2, c_i^3)^T$, the matrix elements of $\mathrm{G}^\gamma$ are $(\mathrm{G}^\gamma)_{\mu\nu} = \delta_{\mu 0}\delta_{\nu\gamma} - \delta_{\nu 0}\delta_{\mu\gamma} + \sum_{\alpha\beta} \epsilon^{\alpha\beta\gamma}\delta_{\mu\alpha}\delta_{\nu\beta}$, yielding the bilinear identity (Freitas et al., 2024):
$$c_i^T \mathrm{G}^\gamma c_i = c_i^0 c_i^\gamma - c_i^\gamma c_i^0 + \sum_{\alpha,\beta=1}^3 \epsilon^{\alpha\beta\gamma} c_i^\alpha c_i^\beta.$$

**Proposition 2.3 (Orthogonal $\mathfrak{su}(2)$ Subalgebras):** The matrices $\mathrm{G}^\gamma$ and $\mathrm{M}^\gamma$ satisfy:
$$[\mathrm{G}^\alpha, \mathrm{G}^\beta] = -2 \sum_{\gamma=1}^3 \epsilon^{\alpha\beta\gamma} \mathrm{G}^\gamma, \qquad [\mathrm{M}^\alpha, \mathrm{M}^\beta] = +2 \sum_{\gamma=1}^3 \epsilon^{\alpha\beta\gamma} \mathrm{M}^\gamma, \qquad [\mathrm{M}^\alpha, \mathrm{G}^\beta] = 0,$$
and trace orthogonality $\mathrm{tr}((\mathrm{M}^\alpha)^T \mathrm{G}^\beta) = 0$, $\mathrm{tr}((\mathrm{G}^\alpha)^T \mathrm{G}^\beta) = \mathrm{tr}((\mathrm{M}^\alpha)^T \mathrm{M}^\beta) = 4\delta^{\alpha\beta}$. Consequently, the local gauge isospin operators $T_i^\gamma = \frac{i}{4} c_i^T \mathrm{G}^\gamma c_i$ and chiral spin operators $S_i^\gamma = \frac{i}{4} c_i^T \mathrm{M}^\gamma c_i$ generate two mutually commuting $\mathfrak{su}(2)$ algebras:
$$[T_i^\alpha, T_i^\beta] = i \sum_{\gamma=1}^3 \epsilon^{\alpha\beta\gamma} T_i^\gamma, \qquad [S_i^\alpha, S_i^\beta] = i \sum_{\gamma=1}^3 \epsilon^{\alpha\beta\gamma} S_i^\gamma, \qquad [S_i^\alpha, T_j^\beta] = 0.$$

---

## 3. Local Fermion Parity and the Exact Casimir Theorem

**Theorem 3.1 (Casimir-Parity Identity):** Let $D_i = c_i^0 c_i^1 c_i^2 c_i^3$ denote the local fermion parity operator ($D_i^\dagger = D_i$, $D_i^2 = \mathbb{1}$). The quadratic Casimir operator of the gauge isospin algebra $\mathbf{T}_i^2 = \sum_{\gamma=1}^3 (T_i^\gamma)^2$ with $T_i^\gamma = \frac{i}{4} c_i^T \mathrm{G}^\gamma c_i$ (normalized to $j=1/2$ doublet eigenvalues) satisfies (Seifert et al., 2018; Freitas, 2024):
$$\mathbf{T}_i^2 = \frac{3}{8}(\mathbb{1} - D_i), \qquad \mathbf{S}_i^2 = \frac{3}{8}(\mathbb{1} + D_i).$$

*Proof:* Expanding $2 T_i^1 = \frac{i}{2}(c_i^0 c_i^1 + c_i^2 c_i^3)$ and squaring gives $(2 T_i^1)^2 = -\frac{1}{4}[(c_i^0 c_i^1)^2 + (c_i^2 c_i^3)^2 + 2 c_i^0 c_i^1 c_i^2 c_i^3] = \frac{1}{2}(\mathbb{1} - D_i)$. Identical contractions hold for $\gamma = 2, 3$. Summing over $\gamma \in \{1, 2, 3\}$ with the $\frac{1}{2}$ spin-$1/2$ normalization yields $\mathbf{T}_i^2 = \frac{3}{8}(\mathbb{1} - D_i)$. $\square$

**Corollary 3.2 (Equivalence of Physical Projection and Gauge Neutrality):** Because $\mathbf{T}_i^2$ is positive semi-definite, a state $|\psi\rangle \in \mathcal{H}_{\mathrm{Fock}}$ lies in the physical subspace $D_i |\psi\rangle = +|\psi\rangle$ if and only if it is a gauge singlet annihilated by all three gauge generators (Seifert et al., 2018; Ralko and Merino, 2020; Yilmaz et al., 2022; Freitas et al., 2024):
$$D_i |\psi\rangle = +|\psi\rangle \iff \mathbf{T}_i^2 |\psi\rangle = 0 \iff c_i^T \mathrm{G}^\gamma c_i |\psi\rangle = 0, \quad \forall \gamma \in \{1, 2, 3\}.$$

---

## 4. The Kitaev Representation and Extended Exchange Hamiltonian

**Definition 4.1 (Kitaev Projection Matrices $\mathrm{N}^\gamma$):** In the Kitaev representation (Kitaev, 2006; Freitas et al., 2024), the Pauli operators are represented as:
$$\sigma_i^\gamma = i c_i^0 c_i^\gamma = \frac{i}{2} c_i^T \mathrm{N}^\gamma c_i, \qquad \mathrm{N}^\gamma = \frac{1}{2}\left(\mathrm{M}^\gamma + \mathrm{G}^\gamma\right),$$
where explicitly:
$$\mathrm{N}^1 = \Pi_+ \otimes \epsilon, \qquad \mathrm{N}^2 = \epsilon \otimes \Pi_+, \qquad \mathrm{N}^3 = \frac{1}{2}\left(\sigma_x \otimes \epsilon + \epsilon \otimes \sigma_x\right),$$
with matrix elements $(\mathrm{N}^\gamma)_{\mu\nu} = \delta_{\mu 0}\delta_{\nu\gamma} - \delta_{\nu 0}\delta_{\mu\gamma}$ and trace norm $\mathrm{tr}((\mathrm{N}^\alpha)^T \mathrm{N}^\beta) = 2\delta^{\alpha\beta}$. Unlike arbitrary linear combinations of $\mathrm{M}^\gamma$ and $\mathrm{G}^\gamma$, the Kitaev operators $\sigma_i^\gamma = i c_i^0 c_i^\gamma$ satisfy $(\sigma_i^\gamma)^2 = \mathbb{1}$ identically on the entire unprojected Fock space and preserve the exact integrability of the pure Kitaev model via conserved $\mathbb{Z}_2$ bond operators $\hat{u}_{ij}^\gamma = i c_i^\gamma c_j^\gamma$.

**Definition 4.2 (Extended $J\text{-}K\text{-}\Gamma\text{-}\Gamma'\text{-}D\text{-}\mathbf{h}$ Spin Hamiltonian):** On the honeycomb lattice $\Lambda = (\mathcal{V}, \mathcal{E})$, the general bilinear spin Hamiltonian in an external Zeeman field $\mathbf{h} = (h^1, h^2, h^3)^T$ is (Rau et al., 2014; Winter et al., 2016; Freitas et al., 2024):
$$H = \frac{1}{2} \sum_{i \in \mathcal{V}} \sum_{j \in \mathcal{V}_i} \sum_{\alpha,\beta=1}^3 \sigma_i^\alpha \mathbf{J}_{ij}^{\alpha\beta} \sigma_j^\beta - \sum_{i \in \mathcal{V}} \mathbf{h} \cdot \boldsymbol{\sigma}_i,$$
where for a directed bond $\langle ij \rangle_\gamma \in \mathcal{E}$ of type $\gamma \in \{1, 2, 3\}$ with $(\alpha, \beta, \gamma)$ a cyclic permutation of $(1, 2, 3)$, the exchange matrix $\mathbf{J}_{\langle ij \rangle_\gamma}$ is parameterized by bond-dependent couplings $(J_\gamma, K_\gamma, \Gamma_\gamma, \Gamma'_\gamma, D_\gamma)$:
$$\mathbf{J}_{\langle ij \rangle_\gamma}^{\alpha'\beta'} = J_\gamma \delta^{\alpha'\beta'} + K_\gamma \delta^{\alpha'\gamma}\delta^{\beta'\gamma} + \Gamma_\gamma \left(\delta^{\alpha'\alpha}\delta^{\beta'\beta} + \delta^{\alpha'\beta}\delta^{\beta'\alpha}\right) + \Gamma'_\gamma \left(\delta^{\alpha'\alpha}\delta^{\beta'\gamma} + \delta^{\alpha'\gamma}\delta^{\beta'\alpha} + \delta^{\alpha'\beta}\delta^{\beta'\gamma} + \delta^{\alpha'\gamma}\delta^{\beta'\beta}\right) + D_\gamma \left(\delta^{\alpha'\alpha}\delta^{\beta'\beta} - \delta^{\alpha'\beta}\delta^{\beta'\alpha}\right),$$
and $\mathbf{J}_{ji} = \mathbf{J}_{ij}^T$ for the reverse orientation $j \to i$.

---

## 5. Exact Mean-Field Decoupling and Ultralocal State Representation

Substituting $\sigma_i^\gamma = \frac{i}{2} c_i^T \mathrm{N}^\gamma c_i$ and introducing site-dependent Lagrange multipliers $\boldsymbol{\lambda}_i = (\lambda_i^1, \lambda_i^2, \lambda_i^3)^T$ to enforce $\langle c_i^T \mathrm{G}^\gamma c_i \rangle = 0$ symmetrically across all three gauge directions yields the quartic Majorana Hamiltonian (Freitas et al., 2024):
$$H = \frac{1}{8} \sum_{i \in \mathcal{V}} \sum_{j \in \mathcal{V}_i} \sum_{\alpha,\beta=1}^3 \left(i c_i^T \mathrm{N}^\alpha c_i\right) \mathbf{J}_{ij}^{\alpha\beta} \left(i c_j^T \mathrm{N}^\beta c_j\right) - \frac{1}{4} \sum_{i \in \mathcal{V}} \sum_{\gamma=1}^3 \left( 2 h^\gamma i c_i^T \mathrm{N}^\gamma c_i - \lambda_i^\gamma i c_i^T \mathrm{G}^\gamma c_i \right).$$

**Proposition 5.1 (Hartree-Fock-Bogoliubov Mean-Field Hamiltonian):** Defining the real ultralocal bond and on-site correlation matrices (Freitas et al., 2024):
$$\mathrm{U}_{ij}^{\mu\nu} = \langle i c_i^\mu c_j^\nu \rangle = -\mathrm{U}_{ji}^{\nu\mu} \quad (i \neq j), \qquad \mathrm{V}_i^{\mu\nu} = \langle i c_i^\mu c_i^\nu \rangle - i\delta^{\mu\nu} = -\mathrm{V}_i^{\nu\mu},$$
the Wick-decoupled quadratic mean-field Hamiltonian is:
$$H_{\mathrm{MF}} = \sum_{i \in \mathcal{V}} \sum_{j \in \mathcal{V}_i} \frac{i}{4} c_i^T \mathrm{A}_{ij} c_j + \sum_{i \in \mathcal{V}} \frac{i}{4} c_i^T \mathrm{B}_i c_i - \mathrm{C},$$
where the $4 \times 4$ real inter-site hopping matrix $\mathrm{A}_{ij}$ ($\mathrm{A}_{ji} = -\mathrm{A}_{ij}^T$), on-site skew-symmetric matrix $\mathrm{B}_i$ ($\mathrm{B}_i^T = -\mathrm{B}_i$), and scalar energy shift $\mathrm{C}$ are given by Eqs. (25)–(27) of Freitas et al. (2024):
$$\mathrm{A}_{ij} = 2 \sum_{\alpha,\beta=1}^3 \mathbf{J}_{ij}^{\alpha\beta} \mathrm{N}^\alpha \mathrm{U}_{ij} \mathrm{N}^\beta,$$
$$\mathrm{B}_i = \sum_{j \in \mathcal{V}_i} \sum_{\alpha,\beta=1}^3 \mathbf{J}_{ij}^{\alpha\beta} \mathrm{N}^\alpha \mathrm{tr}\!\left(\mathrm{V}_j^T \mathrm{N}^\beta\right) + \sum_{\gamma=1}^3 \left(\lambda_i^\gamma \mathrm{G}^\gamma - 2 h^\gamma \mathrm{N}^\gamma\right),$$
$$\mathrm{C} = \frac{1}{8} \sum_{i \in \mathcal{V}} \sum_{j \in \mathcal{V}_i} \sum_{\alpha,\beta=1}^3 \mathbf{J}_{ij}^{\alpha\beta} \left[ \mathrm{tr}\!\left(\mathrm{V}_i^T \mathrm{N}^\alpha\right) \mathrm{tr}\!\left(\mathrm{V}_j^T \mathrm{N}^\beta\right) + 2\,\mathrm{tr}\!\left(\mathrm{U}_{ij}^T \mathrm{N}^\alpha \mathrm{U}_{ij} \mathrm{N}^\beta\right) \right].$$

---

## References

* **Burnell, F. J. and Nayak, C.** (2011). *SU(2) slave fermion approach to the Kitaev honeycomb lattice model*. Physical Review B **84**, 125125. DOI: [10.1103/PhysRevB.84.125125](https://doi.org/10.1103/PhysRevB.84.125125). arXiv: [1104.5485](https://arxiv.org/abs/1104.5485).
* **Freitas, L. R. D.** (2024). *Majorana mean-field theory with $\mathbb{Z}_2$ vortices: electrical tests in the extended Kitaev model*. Ph.D. thesis, Universidade Federal do Rio Grande do Norte (UFRN) / Heinrich-Heine-Universität Düsseldorf.
* **Freitas, L. R. D., Bauer, T., Egger, R., and Pereira, R. G.** (2024). *Electric polarization near vortices in the extended Kitaev model*. npj Quantum Materials **9**, 33. DOI: [10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5). arXiv: [2311.15750](https://arxiv.org/abs/2311.15750).
* **Fu, J., Knolle, J., and Perkins, N. B.** (2018). *Three types of representation of spin in terms of Majorana fermions and an alternative solution of the Kitaev honeycomb model*. Physical Review B **97**, 115142. DOI: [10.1103/PhysRevB.97.115142](https://doi.org/10.1103/PhysRevB.97.115142). arXiv: [1801.09713](https://arxiv.org/abs/1801.09713).
* **Kitaev, A.** (2006). *Anyons in an exactly solved model and beyond*. Annals of Physics **321**, 2–111. DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005). arXiv: [cond-mat/0506438](https://arxiv.org/abs/cond-mat/0506438).
* **Knolle, J., Bhattacharjee, S., and Moessner, R.** (2018). *Dynamics of a quantum spin liquid beyond integrability: The Kitaev-Heisenberg-$\Gamma$ model in an augmented parton mean-field theory*. Physical Review B **97**, 134432. DOI: [10.1103/PhysRevB.97.134432](https://doi.org/10.1103/PhysRevB.97.134432). arXiv: [1801.03774](https://arxiv.org/abs/1801.03774).
* **Ralko, A. and Merino, J.** (2020). *Novel chiral quantum spin liquids in Kitaev magnets*. Physical Review Letters **124**, 217203. DOI: [10.1103/PhysRevLett.124.217203](https://doi.org/10.1103/PhysRevLett.124.217203). arXiv: [2001.08505](https://arxiv.org/abs/2001.08505).
* **Rau, J. G., Lee, E. K.-H., and Kee, H.-Y.** (2014). *Generic spin model for the honeycomb iridates beyond the Kitaev limit*. Physical Review Letters **112**, 077204. DOI: [10.1103/PhysRevLett.112.077204](https://doi.org/10.1103/PhysRevLett.112.077204). arXiv: [1401.5477](https://arxiv.org/abs/1401.5477).
* **Seifert, U. F. P., Meng, T., and Vojta, M.** (2018). *Fractionalized Fermi liquids and exotic superconductivity in the Kitaev-Kondo lattice*. Physical Review B **97**, 085118. DOI: [10.1103/PhysRevB.97.085118](https://doi.org/10.1103/PhysRevB.97.085118). arXiv: [1710.05018](https://arxiv.org/abs/1710.05018).
* **'t Hooft, G.** (1976). *Computation of the quantum effects due to a four-dimensional pseudoparticle*. Physical Review D **14**, 3432–3450. DOI: [10.1103/PhysRevD.14.3432](https://doi.org/10.1103/PhysRevD.14.3432).
* **Winter, S. M., Li, Y., Jeschke, H. O., and Valentí, R.** (2016). *Challenges in design of Kitaev materials: Magnetic interactions from competing energy scales*. Physical Review B **93**, 214431. DOI: [10.1103/PhysRevB.93.214431](https://doi.org/10.1103/PhysRevB.93.214431). arXiv: [1602.05459](https://arxiv.org/abs/1602.05459).
* **Yilmaz, F., Kampf, A. P., and Yip, S. K.** (2022). *Phase diagrams of Kitaev models for arbitrary magnetic field orientations*. Physical Review Research **4**, 043024. DOI: [10.1103/PhysRevResearch.4.043024](https://doi.org/10.1103/PhysRevResearch.4.043024). arXiv: [2202.06158](https://arxiv.org/abs/2202.06158).
* **You, Y.-Z., Kimchi, I., and Vishwanath, A.** (2012). *Doping a spin-orbit Mott insulator: Topological superconductivity from the Kitaev-Heisenberg model*. Physical Review B **86**, 085145. DOI: [10.1103/PhysRevB.86.085145](https://doi.org/10.1103/PhysRevB.86.085145). arXiv: [1109.4155](https://arxiv.org/abs/1109.4155).
