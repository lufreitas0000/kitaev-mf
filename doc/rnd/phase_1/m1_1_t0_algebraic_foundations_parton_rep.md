---
id: RND-P1-M1.1-T0
title: "Algebraic Foundations, SO(4) Gauge Group, and Exact Parton Representation"
phase: 1
milestone: M1.1
task: T0
status: completed
claims:
  - CLM-M1-ALG-001
  - CLM-M1-ALG-002
  - CLM-M1-ALG-003
  - CLM-M1-ALG-004
  - CLM-M1-ALG-005
  - CLM-M1-ALG-006
  - CLM-M1-ALG-007
  - CLM-M1-ALG-008
  - CLM-M1-ALG-009
  - CLM-M1-ALG-010
  - CLM-M1-ALG-011
---

# Note 00: Algebraic Foundations, $\mathfrak{so}(4)$ Gauge Group, and Exact Parton Representation

## 0. Conventions and Disambiguated Symbols

- Greek indices $\mu,\nu,\rho,\sigma \in \{0,1,2,3\}$ label Majorana flavors.
- Greek indices $\alpha,\beta,\gamma \in \{1,2,3\} \equiv \{x,y,z\}$ label Cartesian spin and bond directions.
- For bond type $\gamma$, $(\gamma^+,\gamma^-) := (\gamma+1,\gamma+2) \bmod 3$ denotes the cyclic pair in $\{1,2,3\}$.
- Local fermion parity is denoted by $\mathcal{D}_i \equiv c_i^0 c_i^1 c_i^2 c_i^3$, resolving the collision with the Dzyaloshinskii–Moriya coupling $D_\gamma$.
- The $2 \times 2$ projector is $\sigma_{+} := \operatorname{diag}(1, 0) = \frac{1}{2}(I_2 + \sigma_z)$. Replacing $\Pi_0$.
- Quasiparticle Bogoliubov excitations are denoted by $\eta_{n,\mathbf{k}}$.
- Exchange matrices $\mathbf{J}$ are in units where the Kitaev interaction is $H_K = K \sum \sigma_i^\gamma \sigma_j^\gamma$, with $K < 0$ ferromagnetic.

---

## 1. Lattice Graph and Hilbert Space Embedding

**Definition (Lattice Graph)** <!-- claim: CLM-M1-ALG-001 -->
Let $\Lambda = (\mathcal{V}, \mathcal{E})$ be a tricoordinated periodic honeycomb graph with vertex set $\mathcal{V}$ ($|\mathcal{V}| = N_{\mathrm{sites}}$). The graph is bipartite: $\mathcal{V} = \mathcal{V}_A \sqcup \mathcal{V}_B$. Directed nearest-neighbor edges $\langle ij \rangle_\gamma \in \mathcal{E}$ are oriented from $i \in \mathcal{V}_A$ to $j \in \mathcal{V}_B$, with directional type $\gamma \in \{1, 2, 3\} \equiv \{x, y, z\}$. Every vertex belongs to exactly one bond of each type. Nearest neighbors of vertex $i$ are denoted by $\mathcal{V}_i$. Sums $\sum_i \sum_{j \in \mathcal{V}_i}$ run over ordered pairs $(i, j)$, traversing every physical edge twice.

**Definition (Physical Spin Hilbert Space)** <!-- claim: CLM-M1-ALG-002 -->
To each site $i \in \mathcal{V}$, we assign a physical spin-$1/2$ degree of freedom described by the local Hilbert space $\mathcal{H}_{\mathrm{phys}, i} \cong \mathbb{C}^2$. The total physical Hilbert space is:
$$\mathcal{H}_{\mathrm{phys}} = \bigotimes_{i \in \mathcal{V}} \mathcal{H}_{\mathrm{phys}, i} \cong (\mathbb{C}^2)^{\otimes N_{\mathrm{sites}}}, \qquad \dim \mathcal{H}_{\mathrm{phys}} = 2^{N_{\mathrm{sites}}}.$$

**Definition (Self-Dual CAR Algebra and Extended Fock Space)** <!-- claim: CLM-M1-ALG-003 -->
To each vertex $i \in \mathcal{V}$, we assign four Hermitian Majorana fermion operators $c_i = (c_i^0, c_i^1, c_i^2, c_i^3)^T$ generating the self-dual Canonical Anticommutation Relations (CAR) $C^*$-algebra $\mathcal{A}_{\mathrm{CAR}}$:
$$(c_i^\mu)^\dagger = c_i^\mu, \qquad \{ c_i^\mu, c_j^\nu \} = 2 \delta_{ij} \delta^{\mu\nu} \mathbb{1}.$$
Defining two canonical Dirac fermion annihilation operators per site,
$$d_{i, 1} = \frac{1}{2}(c_i^0 + i c_i^1), \qquad d_{i, 2} = \frac{1}{2}(c_i^2 + i c_i^3),$$
satisfying $\{d_{i, a}, d_{j, b}^\dagger\} = \delta_{ij}\delta_{ab}$, the unprojected fermionic Fock space is:
$$\mathcal{H}_{\mathrm{Fock}} = \bigotimes_{i \in \mathcal{V}} \mathcal{F}_i \cong (\mathbb{C}^4)^{\otimes N_{\mathrm{sites}}}, \qquad \dim \mathcal{H}_{\mathrm{Fock}} = 4^{N_{\mathrm{sites}}} = 2^{2N_{\mathrm{sites}}}.$$
Fractionalization doubles the local degrees of freedom from $2$ to $4$, introducing unphysical states that must be projected out by local gauge constraints.

---

## 2. Lie Algebra $\mathfrak{so}(4) \cong \mathfrak{su}(2)_{\mathrm{spin}} \oplus \mathfrak{su}(2)_{\mathrm{gauge}}$

For any real skew-symmetric matrix $A \in \mathfrak{so}(4)$ ($A^T = -A$) and vertex $i$, define the Hermitian quadratic operator:
$$\mathcal{O}_i[A] \equiv \frac{i}{4} c_i^T A c_i = \frac{i}{4} \sum_{\mu,\nu=0}^3 A_{\mu\nu} c_i^\mu c_i^\nu.$$

**Lemma (Commutator Lie Homomorphism)** <!-- claim: CLM-M1-ALG-004 -->
For all $A, B \in \mathfrak{so}(4)$,
$$[\mathcal{O}_i[A], \mathcal{O}_i[B]] = i\, \mathcal{O}_i\big[[A, B]\big].$$

*Proof.*
Expanding the product of quadratic forms using the CAR relation $\{c^\mu, c^\nu\} = 2\delta^{\mu\nu}$ yields:
$$[c^a c^b, c^c c^d] = 2\big(\delta^{bc} c^a c^d - \delta^{ac} c^b c^d + \delta^{bd} c^c c^a - \delta^{ad} c^c c^b\big).$$
Contracting with $A_{ab} B_{cd}$ and using $A^T = -A, B^T = -B$ gives $[c^T A c, c^T B c] = 4 c^T [A, B] c$. Thus:
$$[\mathcal{O}[A], \mathcal{O}[B]] = \left(\frac{i}{4}\right)^2 \cdot 4 c^T [A, B] c = -\frac{1}{4} c^T [A, B] c = i \cdot \frac{i}{4} c^T [A, B] c = i\, \mathcal{O}\big[[A, B]\big]. \quad \square$$

**Definition (Basis Generators of $\mathfrak{so}(4)$)** <!-- claim: CLM-M1-ALG-005 -->
Let $I_2$ be the $2 \times 2$ identity, $\sigma_x, \sigma_z$ the standard Pauli matrices, and $\epsilon = i\sigma_y = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$. Define the gauge generators $\mathrm{G}^\gamma$ (anti-self-dual 't Hooft matrices) and chiral-spin generators $\mathrm{M}^\gamma$ (self-dual 't Hooft matrices) by:
$$\mathrm{G}^1 = I_2 \otimes \epsilon, \qquad \mathrm{G}^2 = \epsilon \otimes \sigma_z, \qquad \mathrm{G}^3 = \epsilon \otimes \sigma_x,$$
$$\mathrm{M}^1 = \sigma_z \otimes \epsilon, \qquad \mathrm{M}^2 = \epsilon \otimes I_2, \qquad \mathrm{M}^3 = \sigma_x \otimes \epsilon.$$
Explicitly, the matrix entries in the basis $(c^0, c^1, c^2, c^3)^T$ satisfy:
$$c_i^T \mathrm{G}^\gamma c_i = 2\left( c_i^0 c_i^\gamma + c_i^{\gamma^+} c_i^{\gamma^-} \right), \qquad c_i^T \mathrm{M}^\gamma c_i = 2\left( c_i^0 c_i^\gamma - c_i^{\gamma^+} c_i^{\gamma^-} \right).$$

**Proposition (Commutation and Trace Orthogonality)** <!-- claim: CLM-M1-ALG-006 -->
1. *Matrix level:*
   $$[\mathrm{G}^\alpha, \mathrm{G}^\beta] = -2 \sum_\gamma \epsilon^{\alpha\beta\gamma} \mathrm{G}^\gamma, \qquad [\mathrm{M}^\alpha, \mathrm{M}^\beta] = +2 \sum_\gamma \epsilon^{\alpha\beta\gamma} \mathrm{M}^\gamma, \qquad [\mathrm{M}^\alpha, \mathrm{G}^\beta] = 0,$$
   $$\operatorname{tr}\!\left( (\mathrm{G}^\alpha)^T \mathrm{G}^\beta \right) = \operatorname{tr}\!\left( (\mathrm{M}^\alpha)^T \mathrm{M}^\beta \right) = 4\delta^{\alpha\beta}, \qquad \operatorname{tr}\!\left( (\mathrm{M}^\alpha)^T \mathrm{G}^\beta \right) = 0.$$
2. *Operator level:* With the spin-$1/2$ normalization:
   $$S_i^\gamma \equiv \frac{\mathcal{O}_i[\mathrm{M}^\gamma]}{2} = \frac{i}{8} c_i^T \mathrm{M}^\gamma c_i, \qquad T_i^\gamma \equiv \frac{\mathcal{O}_i[\mathrm{G}^\gamma]}{2} = \frac{i}{8} c_i^T \mathrm{G}^\gamma c_i,$$
   the operators satisfy two mutually commuting $\mathfrak{su}(2)$ Lie algebras:
   $$[S_i^\alpha, S_i^\beta] = +i \sum_\gamma \epsilon^{\alpha\beta\gamma} S_i^\gamma, \qquad [T_i^\alpha, T_i^\beta] = -i \sum_\gamma \epsilon^{\alpha\beta\gamma} T_i^\gamma, \qquad [S_i^\alpha, T_j^\beta] = 0.$$

---

## 3. Local Fermion Parity, Casimir Identity, and Gutzwiller Projection

**Theorem (Casimir–Parity Identity)** <!-- claim: CLM-M1-ALG-007 -->
Let $\mathcal{D}_i \equiv c_i^0 c_i^1 c_i^2 c_i^3$ denote the local fermion parity operator ($\mathcal{D}_i^\dagger = \mathcal{D}_i, \mathcal{D}_i^2 = \mathbb{1}$). The quadratic Casimir operators of the gauge isospin and spin algebras satisfy:
$$\mathbf{T}_i^2 \equiv \sum_{\gamma=1}^3 (T_i^\gamma)^2 = \frac{3}{8}\left(\mathbb{1} - \mathcal{D}_i\right), \qquad \mathbf{S}_i^2 \equiv \sum_{\gamma=1}^3 (S_i^\gamma)^2 = \frac{3}{8}\left(\mathbb{1} + \mathcal{D}_i\right).$$

*Proof.*
Using $T_i^1 = \frac{i}{4}(c_i^0 c_i^1 + c_i^2 c_i^3)$, we square the generator:
$$(T_i^1)^2 = -\frac{1}{16}\left( (c_i^0 c_i^1)^2 + (c_i^2 c_i^3)^2 + 2 c_i^0 c_i^1 c_i^2 c_i^3 \right) = -\frac{1}{16}\left(-\mathbb{1} - \mathbb{1} + 2\mathcal{D}_i\right) = \frac{1}{8}\left(\mathbb{1} - \mathcal{D}_i\right).$$
Summing over cyclic permutations $\gamma \in \{1, 2, 3\}$ yields $\mathbf{T}_i^2 = \frac{3}{8}(\mathbb{1} - \mathcal{D}_i)$. An identical calculation for $S_i^1 = \frac{i}{4}(c_i^0 c_i^1 - c_i^2 c_i^3)$ flips the cross-term sign, giving $\mathbf{S}_i^2 = \frac{3}{8}(\mathbb{1} + \mathcal{D}_i)$. $\square$

**Corollary (Exact Gauge-Invariant Physical Subspace)** <!-- claim: CLM-M1-ALG-008 -->
The unprojected 4-dimensional local Fock space decomposes under $\mathfrak{su}(2)_{\mathrm{gauge}} \oplus \mathfrak{su}(2)_{\mathrm{spin}}$ as:
$$\mathcal{F}_i \cong \left(\mathbf{0}, \frac{\mathbf{1}}{\mathbf{2}}\right) \oplus \left(\frac{\mathbf{1}}{\mathbf{2}}, \mathbf{0}\right).$$
A state $|\psi\rangle \in \mathcal{H}_{\mathrm{Fock}}$ belongs to the physical spin Hilbert space $\mathcal{H}_{\mathrm{phys}}$ if and only if it is a gauge singlet annihilated by all three gauge generators:
$$\mathcal{D}_i |\psi\rangle = +|\psi\rangle \iff \mathbf{T}_i^2 |\psi\rangle = 0 \iff T_i^\gamma |\psi\rangle = 0, \quad \forall \gamma \in \{1, 2, 3\}.$$
The exact orthogonal Gutzwiller projector onto the physical Hilbert space is:
$$\mathcal{P}_G = \prod_{i \in \mathcal{V}} \mathcal{P}_{G, i}, \qquad \mathcal{P}_{G, i} = \frac{\mathbb{1} + \mathcal{D}_i}{2} = \frac{\mathbb{1} + c_i^0 c_i^1 c_i^2 c_i^3}{2}.$$

**Proposition (Pauli Algebra of the Kitaev Representation)** <!-- claim: CLM-M1-ALG-009 -->
In the Kitaev representation, the spin-$1/2$ Pauli operators are represented as:
$$\sigma_i^\gamma \equiv i c_i^0 c_i^\gamma = \frac{i}{2} c_i^T \mathrm{N}^\gamma c_i, \qquad \mathrm{N}^\gamma \equiv \frac{1}{2}\left(\mathrm{M}^\gamma + \mathrm{G}^\gamma\right).$$
On the full unprojected Fock space $\mathcal{H}_{\mathrm{Fock}}$, the operators satisfy:
$$(\sigma_i^\gamma)^2 = \mathbb{1}, \qquad \{\sigma_i^\alpha, \sigma_i^\beta\} = 2\delta^{\alpha\beta}\mathbb{1}, \qquad \sigma_i^\alpha \sigma_i^\beta = i \sigma_i^\gamma - 4i T_i^\gamma \quad (\text{cyclic } \alpha, \beta, \gamma).$$
Therefore, on the physical subspace ($\mathcal{D}_i = +1, T_i^\gamma = 0$), the operators $\sigma_i^\gamma$ obey the exact $\mathfrak{su}(2)$ Pauli algebra:
$$\mathcal{P}_G \left( \sigma_i^\alpha \sigma_i^\beta - i\epsilon^{\alpha\beta\gamma} \sigma_i^\gamma \right) \mathcal{P}_G = 0.$$

---

## 4. The Extended Spin Model and Exact Quartic Majorana Hamiltonian

**Definition (Extended $J$-$K$-$\Gamma$-$\Gamma'$-$D$-$\mathbf{h}$ Spin Hamiltonian)** <!-- claim: CLM-M1-ALG-010 -->
On the honeycomb lattice $\Lambda = (\mathcal{V}, \mathcal{E})$, the general bilinear spin Hamiltonian in an external magnetic field $\mathbf{h} = (h^1, h^2, h^3)^T$ is:
$$H = \frac{1}{2} \sum_{i \in \mathcal{V}} \sum_{j \in \mathcal{V}_i} \sum_{a,b=1}^3 \sigma_i^a \mathbf{J}_{ij}^{ab} \sigma_j^b - \sum_{i \in \mathcal{V}} \mathbf{h} \cdot \boldsymbol{\sigma}_i, \qquad \mathbf{J}_{ji} = \mathbf{J}_{ij}^T.$$
For a directed bond $\langle ij \rangle_\gamma \in \mathcal{E}$ ($i \in \mathcal{V}_A \to j \in \mathcal{V}_B$), the $3 \times 3$ exchange tensor $\mathbf{J}_{\langle ij \rangle_\gamma}$ is:
$$\mathbf{J}_{\langle ij \rangle_\gamma}^{ab} = J_\gamma \delta^{ab} + K_\gamma \delta^{a\gamma}\delta^{b\gamma} + \Gamma_\gamma \left( \delta^{a\gamma^+}\delta^{b\gamma^-} + \delta^{a\gamma^-}\delta^{b\gamma^+} \right) + \Gamma'_\gamma \left( \delta^{a\gamma^+}\delta^{b\gamma} + \delta^{a\gamma}\delta^{b\gamma^+} + \delta^{a\gamma^-}\delta^{b\gamma} + \delta^{a\gamma}\delta^{b\gamma^-} \right) + D_\gamma \left( \delta^{a\gamma^+}\delta^{b\gamma^-} - \delta^{a\gamma^-}\delta^{b\gamma^+} \right).$$
Explicitly, for a $z$-bond ($\gamma=3, \gamma^+=1, \gamma^-=2$):
$$\mathbf{J}_{\langle ij \rangle_z} = \begin{pmatrix} J & \Gamma + D & \Gamma' \\ \Gamma - D & J & \Gamma' \\ \Gamma' & \Gamma' & J + K \end{pmatrix}.$$

**Definition (Exact Quartic Majorana Hamiltonian with Lagrange Multipliers)** <!-- claim: CLM-M1-ALG-011 -->
Substituting $\sigma_i^\gamma = \frac{i}{2} c_i^T \mathrm{N}^\gamma c_i$ and introducing local Lagrange multipliers $\boldsymbol{\lambda}_i = (\lambda_i^1, \lambda_i^2, \lambda_i^3)^T \in \mathbb{R}^3$ to enforce $\langle T_i^\gamma \rangle = 0$, the exact Hamiltonian on the unprojected Fock space is:
$$H_\lambda = \frac{1}{8} \sum_{i \in \mathcal{V}} \sum_{j \in \mathcal{V}_i} \sum_{a,b=1}^3 \left( i c_i^T \mathrm{N}^a c_i \right) \mathbf{J}_{ij}^{ab} \left( i c_j^T \mathrm{N}^b c_j \right) - \frac{1}{4} \sum_{i \in \mathcal{V}} \sum_{\gamma=1}^3 \left( 2h^\gamma i c_i^T \mathrm{N}^\gamma c_i - \lambda_i^\gamma i c_i^T \mathrm{G}^\gamma c_i \right).$$
Because $\frac{i}{4} c_i^T \mathrm{G}^\gamma c_i = 2 T_i^\gamma$, the multiplier term vanishes identically on the physical subspace $\mathcal{H}_{\mathrm{phys}}$:
$$\mathcal{P}_G H_\lambda \mathcal{P}_G = \mathcal{P}_G H \mathcal{P}_G.$$
This establishes the exact microscopic operator equivalence between the spin-$1/2$ Hamiltonian and the quartic Majorana system before any mean-field approximation.

---

## References

* **Burnell, F. J. and Nayak, C.** (2011). *SU(2) slave fermion approach to the Kitaev honeycomb lattice model*. Physical Review B **84**, 125125. DOI: [10.1103/PhysRevB.84.125125](https://doi.org/10.1103/PhysRevB.84.125125). arXiv: [1104.5485](https://arxiv.org/abs/1104.5485).
* **Freitas, L. R. D.** (2024). *Majorana mean-field theory with $\mathbb{Z}_2$ vortices: electrical tests in the extended Kitaev model*. Ph.D. thesis, Universidade Federal do Rio Grande do Norte (UFRN) / Heinrich-Heine-Universität Düsseldorf.
* **Freitas, L. R. D., Bauer, T., Egger, R., and Pereira, R. G.** (2024). *Electric polarization near vortices in the extended Kitaev model*. npj Quantum Materials **9**, 33. DOI: [10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5). arXiv: [2311.15750](https://arxiv.org/abs/2311.15750).
* **Fu, J., Knolle, J., and Perkins, N. B.** (2018). *Three types of representation of spin in terms of Majorana fermions and an alternative solution of the Kitaev honeycomb model*. Physical Review B **97**, 115142. DOI: [10.1103/PhysRevB.97.115142](https://doi.org/10.1103/PhysRevB.97.115142). arXiv: [1801.09713](https://arxiv.org/abs/1801.09713).
* **Kitaev, A.** (2006). *Anyons in an exactly solved model and beyond*. Annals of Physics **321**, 2–111. DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005). arXiv: [cond-mat/0506438](https://arxiv.org/abs/cond-mat/0506438).
* **Rau, J. G., Lee, E. K.-H., and Kee, H.-Y.** (2014). *Generic spin model for the honeycomb iridates beyond the Kitaev limit*. Physical Review Letters **112**, 077204. DOI: [10.1103/PhysRevLett.112.077204](https://doi.org/10.1103/PhysRevLett.112.077204). arXiv: [1401.5477](https://arxiv.org/abs/1401.5477).
* **'t Hooft, G.** (1976). *Computation of the quantum effects due to a four-dimensional pseudoparticle*. Physical Review D **14**, 3432–3450. DOI: [10.1103/PhysRevD.14.3432](https://doi.org/10.1103/PhysRevD.14.3432).
* **Winter, S. M., Li, Y., Jeschke, H. O., and Valentí, R.** (2016). *Challenges in design of Kitaev materials: Magnetic interactions from competing energy scales*. Physical Review B **93**, 214431. DOI: [10.1103/PhysRevB.93.214431](https://doi.org/10.1103/PhysRevB.93.214431). arXiv: [1602.05459](https://arxiv.org/abs/1602.05459).
* **You, Y.-Z., Kimchi, I., and Vishwanath, A.** (2012). *Doping a spin-orbit Mott insulator: Topological superconductivity from the Kitaev-Heisenberg model*. Physical Review B **86**, 085145. DOI: [10.1103/PhysRevB.86.085145](https://doi.org/10.1103/PhysRevB.86.085145). arXiv: [1109.4155](https://arxiv.org/abs/1109.4155).
