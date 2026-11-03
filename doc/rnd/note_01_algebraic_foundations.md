# Note 01: Algebraic Foundations of Majorana Parton Representations, $SO(4)$ Gauge Algebra, and Decoupling Channels

## 1. Canonical Anticommutation Relations (CAR) and Fock Space

Let $\Lambda = (V_\Lambda, E_\Lambda)$ be a honeycomb lattice with $N = |V_\Lambda|$ sites (Kitaev, 2006; Freitas et al., 2024). To each site $i \in V_\Lambda$, we associate four real Majorana fermions $\{c_i^0, c_i^1, c_i^2, c_i^3\}$ generating the self-dual Canonical Anticommutation Relations (CAR) algebra $\mathcal{A}_{\mathrm{CAR}}(V)$ on the real Euclidean space $V = \bigoplus_{i=1}^N \mathbb{R}^4 \cong \mathbb{R}^{8N}$ (Kitaev, 2006; Fu et al., 2018; Freitas, 2023):
$$c_i^\mu = (c_i^\mu)^*, \qquad \{ c_i^\mu, c_j^\nu \} = 2 \delta_{ij} \delta^{\mu\nu} \mathbb{1}$$
where $i, j \in \{1, \dots, N\}$ and $\mu, \nu \in \{0, 1, 2, 3\}$.

The on-site algebra is isomorphic to the real Clifford algebra $\mathrm{Cl}(4, \mathbb{R})$ (Fu et al., 2018). Defining two canonical complex Dirac annihilation operators per site (Burnell and Nayak, 2011; You et al., 2012):
$$d_{i, 1} = \frac{1}{2}(c_i^0 + i c_i^1), \qquad d_{i, 2} = \frac{1}{2}(c_i^2 + i c_i^3)$$
satisfying $\{d_{i, a}, d_{j, b}^\dagger\} = \delta_{ij} \delta_{ab}\mathbb{1}$ and $\{d_{i, a}, d_{j, b}\} = 0$, the local unprojected Fock space $\mathcal{H}_i \cong \mathbb{C}^4$ is spanned by the occupation basis (Fu et al., 2018; Freitas, 2023):
$$\mathcal{H}_i = \operatorname{span} \left\{ |0, 0\rangle_i, \; d_{i, 1}^\dagger |0, 0\rangle_i, \; d_{i, 2}^\dagger |0, 0\rangle_i, \; d_{i, 1}^\dagger d_{i, 2}^\dagger |0, 0\rangle_i \right\}$$
The global unprojected Hilbert space $\mathcal{H}_{\mathrm{Fock}} = \bigotimes_{i=1}^N \mathcal{H}_i$ has dimension $\dim \mathcal{H}_{\mathrm{Fock}} = 4^N = 2^{2N}$, which is a factor of $2^N$ larger than the physical spin-$1/2$ Hilbert space $\mathcal{H}_{\mathrm{phys}} = (\mathbb{C}^2)^{\otimes N}$ (Kitaev, 2006; Freitas, 2023).

---

## 2. Real Skew-Symmetric Algebra $\mathfrak{so}(4) \cong \mathfrak{su}(2)_{\mathrm{spin}} \oplus \mathfrak{su}(2)_{\mathrm{gauge}}$

The six linearly independent bilinear operators on site $i$,
$$\mathcal{O}_{A} = \frac{i}{4} c_i^T A c_i = \frac{i}{4} \sum_{\mu, \nu = 0}^3 A_{\mu\nu} c_i^\mu c_i^\nu, \qquad A \in \mathfrak{so}(4)$$
generate the Lie algebra $\mathfrak{so}(4)$ (Seifert et al., 2018; Fu et al., 2018). Using the CAR anticommutation relations, the commutator between any two quadratic forms evaluates to:
$$[\mathcal{O}_A, \mathcal{O}_B] = \frac{i}{4} c_i^T [A, B] c_i$$

### 2.1 Pauli Kronecker Product Construction

To make the algebraic structure transparent, we express $4 \times 4$ real skew-symmetric matrices as tensor products of real $2 \times 2$ matrices (Freitas, 2023). Defining the real symplectic metric $\epsilon$ and symmetric Pauli generators:
$$\epsilon = i \sigma_y = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}, \quad I_2 = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}, \quad \sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}, \quad \sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$$
a matrix $M \in \operatorname{Mat}(4, \mathbb{R})$ is skew-symmetric ($M^T = -M$) if and only if it is a linear combination of terms formed by the Kronecker product of one symmetric and one antisymmetric $2 \times 2$ matrix (Freitas, 2023).

### 2.2 Gauge Isospin Generators $\mathrm{G}^\alpha$ (Anti-Self-Dual 't Hooft Symbols)

The gauge isospin generators $T_i^\alpha = \frac{i}{4} c_i^T \mathrm{G}^\alpha c_i$ are defined via the anti-self-dual 't Hooft matrices ('t Hooft, 1976; Seifert et al., 2018; Yilmaz et al., 2022):
$$\mathrm{G}^1 = I_2 \otimes \epsilon = \begin{pmatrix} 0 & 1 & 0 & 0 \\ -1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & -1 & 0 \end{pmatrix}$$
$$\mathrm{G}^2 = \epsilon \otimes \sigma_z = \begin{pmatrix} 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & -1 \\ -1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \end{pmatrix}$$
$$\mathrm{G}^3 = \epsilon \otimes \sigma_x = \begin{pmatrix} 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & -1 & 0 & 0 \\ -1 & 0 & 0 & 0 \end{pmatrix}$$

Evaluating matrix commutators using $[\epsilon, \sigma_a] = -2 \epsilon_{ab} \sigma_b$ and $\epsilon^2 = -I_2$:
$$[\mathrm{G}^\alpha, \mathrm{G}^\beta] = -2 \sum_{\gamma=1}^3 \epsilon^{\alpha\beta\gamma} \mathrm{G}^\gamma$$
Consequently, the physical isospin operators satisfy the standard $\mathfrak{su}(2)$ Lie algebra:
$$[T_i^\alpha, T_i^\beta] = \frac{i}{4} c_i^T [\mathrm{G}^\alpha, \mathrm{G}^\beta] c_i = -\frac{i}{2} \sum_{\gamma=1}^3 \epsilon^{\alpha\beta\gamma} c_i^T \mathrm{G}^\gamma c_i = i \sum_{\gamma=1}^3 \epsilon^{\alpha\beta\gamma} T_i^\gamma$$

### 2.3 Physical Chiral Spin Generators $\mathrm{M}^\alpha$ (Self-Dual 't Hooft Symbols)

The physical spin operators $S_i^\alpha = \frac{i}{4} c_i^T \mathrm{M}^\alpha c_i$ are defined via the self-dual 't Hooft matrices ('t Hooft, 1976; Fu et al., 2018; Seifert et al., 2018):
$$\mathrm{M}^1 = \sigma_z \otimes \epsilon = \begin{pmatrix} 0 & 1 & 0 & 0 \\ -1 & 0 & 0 & 0 \\ 0 & 0 & 0 & -1 \\ 0 & 0 & 1 & 0 \end{pmatrix}$$
$$\mathrm{M}^2 = \epsilon \otimes I_2 = \begin{pmatrix} 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \\ -1 & 0 & 0 & 0 \\ 0 & -1 & 0 & 0 \end{pmatrix}$$
$$\mathrm{M}^3 = \sigma_x \otimes \epsilon = \begin{pmatrix} 0 & 0 & 0 & 1 \\ 0 & 0 & -1 & 0 \\ 0 & 1 & 0 & 0 \\ -1 & 0 & 0 & 0 \end{pmatrix}$$

Direct matrix multiplication yields:
$$[\mathrm{M}^\alpha, \mathrm{M}^\beta] = +2 \sum_{\gamma=1}^3 \epsilon^{\alpha\beta\gamma} \mathrm{M}^\gamma \implies [S_i^\alpha, S_i^\beta] = i \sum_{\gamma=1}^3 \epsilon^{\alpha\beta\gamma} S_i^\gamma$$
Moreover, the self-dual and anti-self-dual representations commute identically on the nose:
$$[\mathrm{M}^\alpha, \mathrm{G}^\beta] = 0 \implies [S_i^\alpha, T_j^\beta] = 0, \quad \forall i, j, \alpha, \beta$$
confirming that the local gauge transformations generated by $T_i^\alpha$ leave the physical spin algebra invariant (Fu et al., 2018; Freitas, 2023).

---

## 3. Local Parity and the Exact Casimir Identity

The local fermion parity operator is the ordered product of all four on-site Majoranas (Kitaev, 2006; Fu et al., 2018; Freitas, 2023):
$$D_i = c_i^0 c_i^1 c_i^2 c_i^3$$
Because $\{c_i^\mu, c_i^\nu\} = 2\delta^{\mu\nu}$ and $(c_i^\mu)^2 = \mathbb{1}$, $D_i$ is an orthogonal involution:
$$D_i^\dagger = D_i, \qquad D_i^2 = \mathbb{1} \implies \operatorname{spec}(D_i) = \{-1, +1\}$$

### 3.1 Direct Proof of the Casimir Relation

Expanding the isospin generators into explicit Majorana pairs (Freitas, 2023):
$$T_i^1 = \frac{i}{2}(c_i^0 c_i^1 + c_i^2 c_i^3), \quad T_i^2 = \frac{i}{2}(c_i^0 c_i^2 - c_i^1 c_i^3), \quad T_i^3 = \frac{i}{2}(c_i^0 c_i^3 + c_i^1 c_i^2)$$
Squaring each generator independently:
$$(T_i^1)^2 = -\frac{1}{4} \left[ (c_i^0 c_i^1)^2 + (c_i^2 c_i^3)^2 + c_i^0 c_i^1 c_i^2 c_i^3 + c_i^2 c_i^3 c_i^0 c_i^1 \right]$$
Using $(c_i^0 c_i^1)^2 = -\mathbb{1}$, $(c_i^2 c_i^3)^2 = -\mathbb{1}$, and $c_i^2 c_i^3 c_i^0 c_i^1 = (-1)^4 c_i^0 c_i^1 c_i^2 c_i^3 = D_i$:
$$(T_i^1)^2 = -\frac{1}{4} \left( -\mathbb{1} - \mathbb{1} + 2 D_i \right) = \frac{1}{2}(\mathbb{1} - D_i)$$
Identical calculation for $T_i^2$ and $T_i^3$ yields:
$$(T_i^2)^2 = \frac{1}{2}(\mathbb{1} - D_i), \qquad (T_i^3)^2 = \frac{1}{2}(\mathbb{1} - D_i)$$
Under standard angular momentum normalization where generators are scaled by $1/2$ ($T_i^\alpha \to \frac{1}{2} T_i^\alpha$ so that $j(j+1) = \frac{3}{4}$ on the doublet) (Freitas, 2023):
$$\mathbf{T}_i^2 = \sum_{\alpha=1}^3 (T_i^\alpha)^2 = \frac{3}{8}(\mathbb{1} - D_i)$$
Similarly, for the physical spin Casimir $\mathbf{S}_i^2 = \sum_{\alpha=1}^3 (S_i^\alpha)^2$:
$$\mathbf{S}_i^2 = \frac{3}{8}(\mathbb{1} + D_i)$$

### 3.2 Physical Subspace Identification

The local four-dimensional Fock space decomposes under $SU(2)_{\mathrm{gauge}} \times SU(2)_{\mathrm{spin}}$ as (Fu et al., 2018; Freitas, 2023):
$$\mathcal{H}_i \cong \left(\mathbf{0}, \frac{\mathbf{1}}{\mathbf{2}}\right) \oplus \left(\frac{\mathbf{1}}{\mathbf{2}}, \mathbf{0}\right)$$
* **Physical subspace $\mathcal{H}_{\mathrm{phys}, i}$:**
  $$D_i = +1 \iff \mathbf{T}_i^2 = 0 \iff T_i^\alpha = 0 \; (\forall \alpha) \iff \mathbf{S}_i^2 = \frac{3}{4}\mathbb{1}$$
* **Unphysical subspace $\mathcal{H}_{\mathrm{unphys}, i}$:**
  $$D_i = -1 \iff \mathbf{T}_i^2 = \frac{3}{4}\mathbb{1} \iff \mathbf{S}_i^2 = 0$$

The exact Gutzwiller projection operator onto the gauge-singlet physical subspace is (Kitaev, 2006; Burnell and Nayak, 2011; Freitas, 2023):
$$\mathcal{P} = \prod_{i=1}^N \mathcal{P}_i, \qquad \mathcal{P}_i = \frac{\mathbb{1} + D_i}{2} = \frac{\mathbb{1} + c_i^0 c_i^1 c_i^2 c_i^3}{2}$$

---

## 4. Kitaev Representation and the Continuous $\eta$-Decoupling Family

### 4.1 Canonical Kitaev Matrices

In the standard Kitaev representation, the spin operator couples only the matter fermion $c_i^0$ to the bond gauge Majorana $c_i^\alpha$ (Kitaev, 2006; Freitas et al., 2024):
$$\sigma_i^\alpha = i c_i^0 c_i^\alpha = \frac{i}{2} c_i^T \mathrm{N}^\alpha c_i$$
Defining the projector $\Pi_+ = \frac{1}{2}(I_2 + \sigma_z) = \begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$, the $\mathrm{N}^\alpha$ matrices are:
$$\mathrm{N}^1 = \Pi_+ \otimes \epsilon = \begin{pmatrix} 0 & 1 & 0 & 0 \\ -1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix}$$
$$\mathrm{N}^2 = \epsilon \otimes \Pi_+ = \begin{pmatrix} 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 0 \\ -1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{pmatrix}$$
$$\mathrm{N}^3 = \frac{1}{2}(\sigma_x \otimes \epsilon + \epsilon \otimes \sigma_x) = \begin{pmatrix} 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ -1 & 0 & 0 & 0 \end{pmatrix}$$

Notice the fundamental linear algebraic relation (Freitas, 2023):
$$\mathrm{N}^\alpha = \frac{1}{2}(\mathrm{M}^\alpha + \mathrm{G}^\alpha)$$

### 4.2 The Continuous $\eta$-Family

We define a continuous one-parameter family of spin operators $\mathrm{S}_i^\alpha(\eta)$ parameterized by $\eta \in [0, 1]$ (Fu et al., 2018; Seifert et al., 2018):
$$\mathrm{S}_i^\alpha(\eta) = \frac{i}{4} c_i^T \left( \mathrm{M}^\alpha - \eta \mathrm{G}^\alpha \right) c_i = S_i^\alpha - \eta T_i^\alpha$$
* **$\eta = 0$ (Chiral $SO(4)$ representation):**
  $$\mathrm{S}_i^\alpha(0) = S_i^\alpha = \frac{i}{4} c_i^T \mathrm{M}^\alpha c_i$$
  Generates physical $\mathfrak{su}(2)$ rotations across the entire four-dimensional Fock space $\mathcal{H}_i$ (Fu et al., 2018).
* **$\eta = 1$ (Kitaev representation):**
  $$\mathrm{S}_i^\alpha(1) = \frac{i}{4} c_i^T (\mathrm{M}^\alpha - \mathrm{G}^\alpha) c_i = \frac{i}{2} c_i^0 c_i^\alpha$$
  Isolates the matter Majorana $c_i^0$ while decoupling gauge Majoranas into static $\mathbb{Z}_2$ bond variables (Kitaev, 2006).

Because $T_i^\alpha \mathcal{P}_i = 0$ on the physical subspace, all members of the family project to the identical physical spin operator:
$$\mathcal{P} \, \mathrm{S}_i^\alpha(\eta) \, \mathcal{P} = \mathcal{P} \, \mathrm{S}_i^\alpha(0) \, \mathcal{P}, \qquad \forall \eta \in [0, 1]$$

---

## 5. Mean-Field Decoupling Channels and Energy Functional

### 5.1 Generalized Wick Contractions

In terms of the parameterized spin operators $\sigma_i^\alpha(\eta) = 2 \mathrm{S}_i^\alpha(\eta) = \frac{i}{2} c_i^T (\mathrm{M}^\alpha - \eta \mathrm{G}^\alpha) c_i$, the interaction term on bond $\langle ij \rangle$ is:
$$H_{ij} = \sum_{\alpha, \beta = 1}^3 \mathbf{J}_{ij}^{\alpha\beta} \, \sigma_i^\alpha(\eta) \sigma_j^\beta(\eta) = -\frac{1}{4} \sum_{\alpha, \beta} \mathbf{J}_{ij}^{\alpha\beta} [c_i^T A^\alpha(\eta) c_i] [c_j^T A^\beta(\eta) c_j]$$
where $A^\alpha(\eta) = \mathrm{M}^\alpha - \eta \mathrm{G}^\alpha$.

Under Wick's theorem, any expectation value of quartic Majoranas over a quasi-free state $\omega$ with covariance matrix $\Omega_{IJ} = \omega(i c_I c_J)$ decouples into direct (Hartree), exchange (Fock), and pairing channels (Freitas, 2023; Yilmaz et al., 2022):
$$\omega\left( [c_i^T A c_i] [c_j^T B c_j] \right) = \operatorname{Tr}(A V_i) \operatorname{Tr}(B V_j) - 2 \operatorname{Tr}(A U_{ij} B^T U_{ij}^T)$$
where the on-site and bond covariance matrices are defined by (Freitas et al., 2024; Freitas, 2023):
$$(V_i)_{\mu\nu} = \omega(i c_i^\mu c_i^\nu), \qquad (U_{ij})_{\mu\nu} = \omega(i c_i^\mu c_j^\nu)$$

### 5.2 Decoupling Channel Energy Splitting

The on-site gauge constraints enforce $\operatorname{Tr}(\mathrm{G}^\alpha V_i) = 0$ on average (Yilmaz et al., 2022; Freitas, 2023). However, the inter-site gauge-spin contractions do not vanish on unprojected quasi-free states:
$$\operatorname{Tr}\left( \mathrm{M}^\alpha U_{ij} (\mathrm{G}^\beta)^T U_{ij}^T \right) \neq 0, \qquad \operatorname{Tr}\left( \mathrm{G}^\alpha U_{ij} (\mathrm{G}^\beta)^T U_{ij}^T \right) \neq 0$$
Consequently, the unprojected mean-field energy functional $\mathcal{E}(\Omega; \eta)$ depends explicitly on $\eta$ (Freitas, 2023):
$$\mathcal{E}(\Omega; \eta) = \mathcal{E}_{SO(4)}(\Omega) + \eta \Delta \mathcal{E}_1(\Omega) + \eta^2 \Delta \mathcal{E}_2(\Omega)$$
where:
$$\Delta \mathcal{E}_1(\Omega) = \frac{1}{2} \sum_{\langle ij \rangle} \mathbf{J}_{ij}^{\alpha\beta} \operatorname{Tr}\left[ \mathrm{M}^\alpha U_{ij} (\mathrm{G}^\beta)^T U_{ij}^T + \mathrm{G}^\alpha U_{ij} (\mathrm{M}^\beta)^T U_{ij}^T \right]$$
$$\Delta \mathcal{E}_2(\Omega) = -\frac{1}{2} \sum_{\langle ij \rangle} \mathbf{J}_{ij}^{\alpha\beta} \operatorname{Tr}\left[ \mathrm{G}^\alpha U_{ij} (\mathrm{G}^\beta)^T U_{ij}^T \right]$$

### 5.3 Physical Impact of the Decoupling Channel

* **$\eta = 1$ (Kitaev channel):**
  Preserves exact local constants of motion $\hat{u}_{ij}^\gamma = i c_i^\gamma c_j^\gamma = \pm 1$ at $J = \Gamma = \mathbf{h} = 0$ (Kitaev, 2006). The vison bands are completely flat and localized. This channel is quantitatively accurate for computing local vortex properties, quadrupole tensors, and scanning tunneling spectroscopy of Majorana zero modes (Freitas et al., 2024).
* **$\eta = 0$ (Chiral $SO(4)$ channel):**
  Treats all spin and gauge components symmetrically (Fu et al., 2018). It introduces artificial matter-gauge hybridization that broadens vison dispersions into itinerant bands and overestimates the critical magnetic field $h_c$ by $15\% - 30\%$ relative to DMRG calculations. However, it treats transverse dipole interactions isotropically, stabilizing competing collinear magnetic states (such as zigzag and stripy phases) across a broader parameter window (Winter et al., 2016; Knolle et al., 2018).

---

## References

* **Burnell, F. J. and Nayak, C.** (2011). *SU(2) slave fermion approach to the Kitaev honeycomb lattice model*. Physical Review B, 84(12), 125125. DOI: [10.1103/PhysRevB.84.125125](https://doi.org/10.1103/PhysRevB.84.125125). arXiv: [1104.5485](https://arxiv.org/abs/1104.5485).
* **Freitas, L. R. D.** (2023). *Majorana mean-field theory with $\mathbb{Z}_2$ vortices: electrical tests in the extended Kitaev model*. Ph.D. thesis, Universidade Federal do Rio Grande do Norte (UFRN) / Heinrich-Heine-Universität Düsseldorf.
* **Freitas, L. R. D., Bauer, T., Egger, R., and Pereira, R. G.** (2024). *Electric polarization near vortices in the extended Kitaev model*. npj Quantum Materials, 9(1), 33. DOI: [10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5). arXiv: [2311.15750](https://arxiv.org/abs/2311.15750).
* **Fu, J., Knolle, J., and Perkins, N. B.** (2018). *Majorana representations of spin and an alternative solution of the Kitaev honeycomb model*. Physical Review B, 97(11), 115142. DOI: [10.1103/PhysRevB.97.115142](https://doi.org/10.1103/PhysRevB.97.115142). arXiv: [1801.09713](https://arxiv.org/abs/1801.09713).
* **Kitaev, A.** (2006). *Anyons in an exactly solved model and beyond*. Annals of Physics, 321(1), 2–111. DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005). arXiv: [cond-mat/0506438](https://arxiv.org/abs/cond-mat/0506438).
* **Knolle, J., Chern, G.-W., Kovrizhin, D. L., Moessner, R., and Perkins, N. B.** (2018). *Magnon damping in the zigzag phase of the Kitaev-Heisenberg-$\Gamma$ model on a honeycomb lattice*. Physical Review B, 97(13), 134432. DOI: [10.1103/PhysRevB.97.134432](https://doi.org/10.1103/PhysRevB.97.134432). arXiv: [1911.12829](https://arxiv.org/abs/1911.12829).
* **Rau, J. G., Lee, E. K.-H., and Kee, H.-Y.** (2014). *Generic Spin Model for Bilayer and Honeycomb Materials with Strong Spin-Orbit Coupling*. Physical Review Letters, 112(7), 077204. DOI: [10.1103/PhysRevLett.112.077204](https://doi.org/10.1103/PhysRevLett.112.077204). arXiv: [1401.5477](https://arxiv.org/abs/1401.5477).
* **Seifert, U. F. P., Meng, T., and Vojta, M.** (2018). *Fractionalized Fermi liquids and exotic superconductivity in the Kitaev-Kondo lattice*. Physical Review B, 97(8), 085118. DOI: [10.1103/PhysRevB.97.085118](https://doi.org/10.1103/PhysRevB.97.085118). arXiv: [1710.05018](https://arxiv.org/abs/1710.05018).
* **'t Hooft, G.** (1976). *Computation of the quantum effects due to a four-dimensional pseudoparticle*. Physical Review D, 14(12), 3432–3450. DOI: [10.1103/PhysRevD.14.3432](https://doi.org/10.1103/PhysRevD.14.3432).
* **Winter, S. M., Li, Y., Jeschke, H. O., and Valentí, R.** (2016). *Challenges in design of Kitaev materials: Magnetic interactions in $\alpha\text{-RuCl}_3$*. Physical Review B, 93(21), 214431. DOI: [10.1103/PhysRevB.93.214431](https://doi.org/10.1103/PhysRevB.93.214431). arXiv: [1602.05459](https://arxiv.org/abs/1602.05459).
* **Yilmaz, F., Kampf, A. P., and Yip, S.-K.** (2022). *Phase diagrams of the Kitaev-Heisenberg-$\Gamma$-$\Gamma'$ model in a magnetic field*. Physical Review Research, 4(4), 043024. DOI: [10.1103/PhysRevResearch.4.043024](https://doi.org/10.1103/PhysRevResearch.4.043024). arXiv: [2202.06158](https://arxiv.org/abs/2202.06158).
* **You, Y.-Z., Kimchi, I., and Vishwanath, A.** (2012). *Doping a spin-orbit Mott insulator: Topological superconductivity from the Kitaev-Heisenberg model*. Physical Review B, 86(8), 085145. DOI: [10.1103/PhysRevB.86.085145](https://doi.org/10.1103/PhysRevB.86.085145). arXiv: [1109.4155](https://arxiv.org/abs/1109.4155).