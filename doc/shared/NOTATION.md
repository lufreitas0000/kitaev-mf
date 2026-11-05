# Mathematical & Physical Notation Registry: kitaev-mf

This registry records the canonical symbol conventions, mathematical domains, definitions, and collision-avoidance rules enforced across the Living Lab Book (`doc/labbook/`), Front 1 R&D notes (`doc/rnd/`), test specifications, and the Front 3 manuscript (`doc/manuscript/`).

---

## 1. Primary Notation and Disambiguation Index

| Symbol | Mathematical Domain / Type | Physical Definition | First Used | Anti-Collision / Disambiguation Directives |
| :--- | :--- | :--- | :---: | :--- |
| $\Omega$ | Matrix in $\mathfrak{so}(2n)$ | Global real skew-symmetric covariance matrix of a quasi-free state: $\omega(c_a c_b) = \delta_{ab}\Id + i\Omega_{ab}$. | Ch 1b | **Never** used for DIIS Gram matrix or frequency. Always skew-symmetric with $\|\Omega\|_2 \le 1$. |
| $\Xi$ | Tuple $(\mathrm{U}, \mathrm{V})$ | Collective parameter tuple of bond and on-site mean-field covariances. | Ch 2 | **Sole** collective identifier for mean-field covariances; never overloaded with couplings. |
| $\mathrm{U}_b$ | Matrix in $\mathbb{R}^{4 \times 4}$ | Directed bond Majorana correlation matrix $\mathrm{U}_b^{\mu\nu} \equiv -i\langle c_i^\mu c_j^\nu \rangle$ along directed bond $b = \langle ij \rangle$. | Ch 2 | Lower index $b$ indicates directed bond; upper indices $\mu, \nu \in \{0, 1, 2, 3\}$. |
| $\mathrm{V}_s$ | Matrix in $\mathfrak{so}(4)$ | On-site Majorana correlation matrix $\mathrm{V}_s^{\mu\nu} \equiv -i\langle c_{s,i}^\mu c_{s,i}^\nu \rangle$ on sublattice $s$. | Ch 2 | Skew-symmetric: $\mathrm{V}_s^T = -\mathrm{V}_s$. Transposition introduces an explicit sign. |
| $\Gamma, \Gamma'$ | Scalar in $\mathbb{R}$ | Bond-directional symmetric off-diagonal Kitaev exchange interactions ($H_{\mathrm{ex}} \supset \Gamma_\gamma [S_i^\alpha S_j^\beta + S_i^\beta S_j^\alpha]$). | Ch 1 | **Strictly reserved** for Hamiltonian exchange coupling constants. Never used as a covariance or correlation matrix. |
| $\mathcal{C}_{ij}$ | Matrix in $\mathbb{C}^{8 \times 8}$ | Local two-site Majorana correlation matrix $\begin{pmatrix} \Id_4 - i\mathrm{V}_i & -i\mathrm{U}_b \\ i\mathrm{U}_b^T & \Id_4 - i\mathrm{V}_j \end{pmatrix}$. | Ch 5 | Used for semidefinite positivity and two-site LMI bounds. |
| $\mathbf{B}_{\mathrm{DIIS}}$ | Matrix in $\mathbb{R}^{(m+1) \times (m+1)}$ | Augmented Pulay/DIIS Karush-Kuhn-Tucker Gram matrix for subspace acceleration. | Ch 3 | Denoted $\mathbf{B}_{\mathrm{DIIS}}$ or $\mathbf{G}_{\mathrm{DIIS}}$, never $\Omega$. |
| $b \in \mathcal{E}$ | Lattice element | Directed bond identifier $b = (i, j, \gamma_b)$ on the honeycomb lattice. | Ch 2 | **Never** overloaded as a spin index. Spin indices strictly use Greek letters $\alpha, \beta, \gamma \in \{x, y, z\}$ or primed Latin $a, a' \in \{1, 2, 3\}$. |
| $\alpha, \beta, \gamma$ | Indices in $\{x, y, z\}$ | Cartesian spatial and spin polarization components. | Ch 1 | Reserved for physical spin operators and bond types ($x$-, $y$-, $z$-bonds). |
| $\mu, \nu$ | Indices in $\{0, 1, 2, 3\}$ | Majorana flavor indices: $\mu=0$ is matter fermion, $\mu \in \{1, 2, 3\}$ are gauge fermions. | Ch 1 | Greek indices $\mu, \nu$ always range over four Majorana modes per site. |
| $\mathcal{E}_n(\mathbf{k})$ | Scalar $\ge 0$ | Eigenvalues of positive sector of quadratic Bloch Hamiltonian $\mathcal{H}_2(\mathbf{k})$ or $\mathcal{H}_4(\mathbf{q})$. | Ch 2 | Denoted with calligraphic $\mathcal{E}$. Related to excitation energies by factor 4. |
| $\varepsilon_n(\mathbf{k})$ | Scalar $\ge 0$ | Physical Bogoliubov quasiparticle excitation energy: $\varepsilon_n(\mathbf{k}) = 4\mathcal{E}_n(\mathbf{k})$. | Ch 2 | Denoted with varepsilon $\varepsilon$. Always four times $\mathcal{E}_n$. |
| $\omega$ | Functional $\mathcal{A} \to \mathbb{C}$ | Normalized, positive linear expectation functional (state) on the $C^*$-algebra $\mathcal{A}$. | Ch 1b | Distinct from angular frequency; represents quantum state. |
| $\mathcal{D}_i$ | Operator in $\Cl(4)$ | Local fermion parity / Klein operator $\mathcal{D}_i \equiv -c_i^0 c_i^1 c_i^2 c_i^3$. | Ch 1 | Has eigenvalues $\pm 1$. $\mathcal{D}_i = +1$ identifies the physical spin subspace. |
| $T_i^\gamma$ | Operator in $\Cl(4)$ | Anti-self-dual gauge generator $T_i^\gamma \equiv -\frac{i}{4}\epsilon^{\gamma\alpha\beta} c_i^\alpha c_i^\beta$. | Ch 1 | Related to $G_i^\gamma$ by $G_i^\gamma = 2 T_i^\gamma$. Satisfies $[T^\alpha, T^\beta] = -i\epsilon^{\alpha\beta\gamma}T^\gamma$. |
| $W_p$ | Scalar $\in [-1, 1]$ | Factored plaquette Wilson loop product $W_p \equiv \prod_{b \in \partial p} u_b$. | Ch 4 | Approximates full multi-link operator expectation $\langle \hat{W}_p \rangle$. |

---

## 2. Quantum Information & Algebraic QFT Nomenclature

- **$\mathcal{A}_{\Lambda}$:** Quasi-local $C^*$-algebra associated to spatial region $\Lambda \subset \mathcal{L}$.
- **$\mathcal{S}(\mathcal{A}_{\Lambda})$:** Convex state space of positive normalized linear functionals $\omega: \mathcal{A}_{\Lambda} \to \mathbb{C}$, with $\omega(\Id) = 1$.
- **$\mathcal{M}_{\mathrm{QF}}$:** Manifold of quasi-free (Gaussian) states on $\mathcal{A}_{\mathrm{CAR}}(V)$.
- **$\Gr(n, 2n)$:** Real orthogonal Grassmannian manifold of pure quasi-free states satisfying $\Omega^2 = -\Id$.
- **$D_{\mathrm{tr}}(\rho, \sigma)$:** Trace distance $\frac{1}{2}\Tr|\rho - \sigma|$ between density matrices $\rho, \sigma$.
- **$F(\rho, \sigma)$:** Uhlmann fidelity $(\Tr\sqrt{\sqrt{\rho}\sigma\sqrt{\rho}})^2$.
- **$S(\rho \| \sigma)$:** Quantum relative entropy $\Tr[\rho(\log\rho - \log\sigma)]$.
- **1-RDM / 2-RDM:** One-particle and two-particle reduced density matrices in the fermionic particle-hole basis.
