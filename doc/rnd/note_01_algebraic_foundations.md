# Note 01: Algebraic Foundations of the Majorana Parton Representation, $\mathfrak{so}(4)$ Gauge Structure, and Mean-Field Decoupling

> **Status of this revision.** All statements below were re-derived and cross-checked numerically with `verify_notes_01_03.py` (exact diagonalization on one and two sites, Jordan–Wigner representation). Corrections with respect to the previous version are collected in the *Revision log* at the end. Notation introduced in §0 is shared by Notes 02 and 03.

---

## 0. Conventions and Notation (shared by Notes 01–03)

**Index sets.** Greek indices from the middle of the alphabet, $\mu,\nu,\rho,\sigma \in \{0,1,2,3\}$, label *Majorana flavors*; Greek indices from the start, $\alpha,\beta,\gamma \in \{1,2,3\} \equiv \{x,y,z\}$, label *Cartesian components* and, for $\gamma$, also the *type of a nearest-neighbor bond*. Dummy indices of $3\times 3$ exchange matrices are written $a,b \in \{1,2,3\}$. For a bond type $\gamma$ we write $(\gamma^+,\gamma^-) := (\gamma+1,\gamma+2) \bmod 3$ (values in $\{1,2,3\}$), so that $(\gamma^+,\gamma^-,\gamma)$ is a cyclic permutation of $(1,2,3)$. $\epsilon^{\alpha\beta\gamma}$ is the Levi-Civita symbol with $\epsilon^{123}=+1$.

**Symbols that were overloaded in the previous version and are now disambiguated.**

| Quantity | Symbol (this revision) | Replaces / clashes with |
|---|---|---|
| local fermion parity $c^0c^1c^2c^3$ | $\mathcal{Z}_i$ (called $D_i$ in Freitas et al., 2024) | $D_\gamma$ = Dzyaloshinskii–Moriya coupling |
| $\mathrm{diag}(1,0)$ | $e_{11}$ | $\Pi_+$ (clashes with spectral projectors $\Pi^\pm(\mathbf{k})$) |
| Bogoliubov quasiparticle operators | $\eta_{n,\mathbf{k}}$ | $\gamma_{n,\mathbf{k}}$ (clashes with bond type $\gamma$) |
| supercell crystal momentum | $\mathbf{q}=(q_1,q_2)$ | $(K_1,K_2)$ (clashes with the Kitaev coupling $K$) |
| iteration counter (Note 03) | $\ell$ | $k$ (clashes with momentum $\mathbf{k}$) |
| DIIS weights / Lagrange multiplier / Gram matrix | $\xi_j$ / $\nu$ / $\Omega$ | $c_j$ / $\mu$ / $B$ (clash with Majorana $c$, flavor $\mu$, on-site matrix $\mathrm{B}$) |
| eigenvalues of the Bloch matrix $\mathcal{H}(\mathbf{k})$ | $\mathcal{E}_n(\mathbf{k})$ ("Bloch energies") | $E_n(\mathbf{k})$ |
| physical excitation energies | $\varepsilon_n(\mathbf{k}) = 4\,\mathcal{E}_n(\mathbf{k})$ (Note 02, Prop. 3.5) | — |
| Pauli-matrix (Kronecker) factors | $I_2,\sigma_x,\sigma_z,\epsilon=\mathrm{i}\sigma_y$ | — |

**Normalization of spins.** All spin operators are Pauli operators, $\sigma^\gamma_i$ with $(\sigma^\gamma_i)^2=\mathbb{1}$ (not $S=\sigma/2$). Consequently the exchange matrices $\mathbf{J}$ are in units where the pure Kitaev Hamiltonian reads $H_K=K\sum_{\langle ij\rangle_\gamma}\sigma_i^\gamma\sigma_j^\gamma$, as in Kitaev (2006) and Freitas et al. (2024). The physical regime of $\alpha$-RuCl$_3$ is $K<0$ (ferromagnetic).

**Real/imaginary conventions.** $c_i^\mu$ are Hermitian operators. A *real skew-symmetric* matrix $A\in\mathfrak{so}(4)$ yields the *anti-Hermitian* operator $c^TAc$ (because $(c^TAc)^\dagger=c^TA^Tc=-c^TAc$) and the *Hermitian* operator $\mathrm{i}\,c^TAc$. All mean-field parameters are expectation values of Hermitian operators and are therefore real.

---

## 1. Lattice Graph and the Self-Dual CAR Algebra

**Definition 1.1 (Lattice graph).** Let $\Lambda=(\mathcal{V},\mathcal{E})$ be a finite or periodic honeycomb graph with vertex set $\mathcal{V}$, $N_{\mathrm{sites}}=\vert{}\mathcal{V}\vert{}$. $\mathcal{V}=\mathcal{V}_A\sqcup\mathcal{V}_B$ is bipartite. The edge set $\mathcal{E}$ consists of *directed* nearest-neighbor bonds $\langle ij\rangle_\gamma$ with $i\in\mathcal{V}_A$, $j\in\mathcal{V}_B$, of type $\gamma\in\{1,2,3\}$; every vertex belongs to exactly one bond of each type (periodic case). $\mathcal{V}_i:=\{j:\langle ij\rangle_\gamma\in\mathcal{E}\ \text{or}\ \langle ji\rangle_\gamma\in\mathcal{E}\}$ is the set of the three nearest neighbors of $i$. Sums $\sum_{i}\sum_{j\in\mathcal{V}_i}$ run over **ordered** pairs, i.e. every bond is visited twice (once in each orientation).

**Definition 1.2 (Self-dual CAR algebra).** To each vertex $i$ we assign four Hermitian generators $c_i=(c_i^0,c_i^1,c_i^2,c_i^3)^T$ of the unital algebra $\mathcal{A}_{\mathrm{CAR}}$ with
$$(c_i^\mu)^\dagger=c_i^\mu,\qquad \{c_i^\mu,c_j^\nu\}=2\,\delta_{ij}\delta^{\mu\nu}\,\mathbb{1}.$$
$\mathcal{A}_{\mathrm{CAR}}$ is the complex Clifford algebra on $4N_{\mathrm{sites}}$ generators, $\mathcal{A}_{\mathrm{CAR}}\cong \mathrm{Mat}(2^{2N_{\mathrm{sites}}},\mathbb{C})$. It is a finite-dimensional $C^*$-algebra, so every representation is unitarily equivalent to the Jordan–Wigner (Fock) representation.

**Remark 1.3 (Local Fock space).** The two operators $d_{i,1}=\tfrac12(c_i^0+\mathrm{i}c_i^1)$ and $d_{i,2}=\tfrac12(c_i^2+\mathrm{i}c_i^3)$ satisfy the canonical (complex) anticommutation relations, e.g. $\{d_{i,1},d_{i,1}^\dagger\}=\tfrac14(2+2)=1$ (Burnell and Nayak, 2011; You et al., 2012). The local Fock space is $\mathcal{F}_i\cong\mathbb{C}^4$, twice as large as the physical pseudospin space $\mathbb{C}^2$. The physical subspace is selected in §3.

---

## 2. Lie Algebra $\mathfrak{so}(4)\cong\mathfrak{su}(2)_{\mathrm{spin}}\oplus\mathfrak{su}(2)_{\mathrm{gauge}}$

For $A\in\mathfrak{so}(4)$ (real, $A^T=-A$) and a vertex $i$ define the Hermitian on-site operator
$$\mathcal{O}_i[A]:=\frac{\mathrm{i}}{4}\,c_i^TAc_i=\frac{\mathrm{i}}{4}\sum_{\mu\nu}A_{\mu\nu}\,c_i^\mu c_i^\nu.$$

**Lemma 2.1 (Commutator of quadratic forms).** For $A,B\in\mathfrak{so}(4)$,
$$\boxed{\;[\mathcal{O}_i[A],\mathcal{O}_i[B]]=\mathrm{i}\,\mathcal{O}_i\big[[A,B]\big].\;}$$
*(The previous version omitted the factor $\mathrm{i}$. Without it the identity cannot hold: the left side is anti-Hermitian, the right side $\tfrac{\mathrm{i}}{4}c^T[A,B]c$ is Hermitian.)*

*Proof.* From $\{c^\mu,c^\nu\}=2\delta^{\mu\nu}$ one finds the elementary identity
$$[c^ac^b,c^cc^d]=2\big(\delta^{bc}c^ac^d-\delta^{ac}c^bc^d+\delta^{bd}c^cc^a-\delta^{ad}c^cc^b\big).$$
Contracting with $A_{ab}B_{cd}$ and using $A^T=-A$, $B^T=-B$, the four terms give $2c^TABc,\ 2c^TABc,\ -2c^TBAc,\ -2c^TBAc$, i.e. $[c^TAc,c^TBc]=4\,c^T[A,B]c$. Hence $[\mathcal{O}[A],\mathcal{O}[B]]=(\tfrac{\mathrm{i}}{4})^2\cdot 4\,c^T[A,B]c=-\tfrac14c^T[A,B]c=\mathrm{i}\cdot\tfrac{\mathrm{i}}{4}c^T[A,B]c$. $\square$

**Definition 2.2 (Basis of $\mathfrak{so}(4)$).** With the real $2\times 2$ matrices $I_2$, $\sigma_x$, $\sigma_z$, $\epsilon=\mathrm{i}\sigma_y=\begin{psmallmatrix}0&1\\-1&0\end{psmallmatrix}$ and $e_{11}=\mathrm{diag}(1,0)$, and the flavor ordering $\mu=2a+b$ with $a,b\in\{0,1\}$ (so $\mathbb{R}^4=\mathbb{R}^2\otimes\mathbb{R}^2$), define the *gauge* generators $\mathrm{G}^\gamma$ (anti-self-dual) and *chiral-spin* generators $\mathrm{M}^\gamma$ (self-dual) ('t Hooft, 1976; Fu et al., 2018):
$$\mathrm{G}^1=I_2\otimes\epsilon,\quad \mathrm{G}^2=\epsilon\otimes\sigma_z,\quad \mathrm{G}^3=\epsilon\otimes\sigma_x,$$
$$\mathrm{M}^1=\sigma_z\otimes\epsilon,\quad \mathrm{M}^2=\epsilon\otimes I_2,\quad \mathrm{M}^3=\sigma_x\otimes\epsilon.$$
Their entries are, in the flavor basis $(c^0,c^1,c^2,c^3)$,
$$(\mathrm{G}^\gamma)_{\mu\nu}=\delta_{\mu0}\delta_{\nu\gamma}-\delta_{\nu0}\delta_{\mu\gamma}+\sum_{\alpha\beta}\epsilon^{\alpha\beta\gamma}\delta_{\mu\alpha}\delta_{\nu\beta},$$
$$(\mathrm{M}^\gamma)_{\mu\nu}=\delta_{\mu0}\delta_{\nu\gamma}-\delta_{\nu0}\delta_{\mu\gamma}-\sum_{\alpha\beta}\epsilon^{\alpha\beta\gamma}\delta_{\mu\alpha}\delta_{\nu\beta}.$$
For example $\mathrm{G}^1$ has non-zero entries $(0,1)=(2,3)=+1$, $(1,0)=(3,2)=-1$. Consequently,
$$c^T\mathrm{G}^\gamma c=c^0c^\gamma-c^\gamma c^0+\sum_{\alpha\beta}\epsilon^{\alpha\beta\gamma}c^\alpha c^\beta=2\Big(c^0c^\gamma+c^{\gamma^+}c^{\gamma^-}\Big).\tag{2.1}$$

**Proposition 2.3 (Orthogonal $\mathfrak{su}(2)$ subalgebras).**
(i) *Matrix level:*
$$[\mathrm{G}^\alpha,\mathrm{G}^\beta]=-2\sum_\gamma\epsilon^{\alpha\beta\gamma}\mathrm{G}^\gamma,\quad [\mathrm{M}^\alpha,\mathrm{M}^\beta]=+2\sum_\gamma\epsilon^{\alpha\beta\gamma}\mathrm{M}^\gamma,\quad [\mathrm{M}^\alpha,\mathrm{G}^\beta]=0,$$
$$\mathrm{tr}\big((\mathrm{G}^\alpha)^T\mathrm{G}^\beta\big)=\mathrm{tr}\big((\mathrm{M}^\alpha)^T\mathrm{M}^\beta\big)=4\delta^{\alpha\beta},\qquad \mathrm{tr}\big((\mathrm{M}^\alpha)^T\mathrm{G}^\beta\big)=0.$$
(ii) *Operator level.* With the **spin-$\tfrac12$ normalization**
$$\boxed{\,S_i^\gamma:=\frac{\mathcal{O}_i[\mathrm{M}^\gamma]}{2}=\frac{\mathrm{i}}{8}c_i^T\mathrm{M}^\gamma c_i,\qquad T_i^\gamma:=\frac{\mathcal{O}_i[\mathrm{G}^\gamma]}{2}=\frac{\mathrm{i}}{8}c_i^T\mathrm{G}^\gamma c_i,\,}$$
the following relations hold:
$$[S_i^\alpha,S_i^\beta]=+\mathrm{i}\sum_\gamma\epsilon^{\alpha\beta\gamma}S_i^\gamma,\qquad [T_i^\alpha,T_i^\beta]=-\mathrm{i}\sum_\gamma\epsilon^{\alpha\beta\gamma}T_i^\gamma,\qquad [S_i^\alpha,T_j^\beta]=0.$$

*Remarks.* (a) The previous version defined $T^\gamma=\tfrac{\mathrm{i}}{4}c^T\mathrm{G}^\gamma c$ and still claimed the commutators $\pm\mathrm{i}\epsilon T$ and the Casimir $\tfrac38(1-\mathcal{Z})$; those statements are correct only for the $\tfrac{\mathrm{i}}{8}$ normalization used here (with $\tfrac{\mathrm{i}}{4}$ one gets structure constants $\mp 2\mathrm{i}$ and $\mathbf{T}^2=\tfrac32(1-\mathcal{Z})$). (b) The gauge algebra has the *opposite* sign of the structure constants; $\tilde T^\gamma:=-T^\gamma$ obeys the standard relation $[\tilde T^\alpha,\tilde T^\beta]=+\mathrm{i}\epsilon^{\alpha\beta\gamma}\tilde T^\gamma$. Casimirs are unaffected.

*Proof.* (i) is a finite matrix computation (checked in `test_so4_algebra_and_traces`). (ii) By Lemma 2.1, $[S^\alpha,S^\beta]=\tfrac14\,\mathrm{i}\,\mathcal{O}[[\mathrm{M}^\alpha,\mathrm{M}^\beta]]=\tfrac14\mathrm{i}\cdot2\epsilon^{\alpha\beta\gamma}\mathcal{O}[\mathrm{M}^\gamma]=\mathrm{i}\epsilon^{\alpha\beta\gamma}S^\gamma$; the gauge case is identical with $-2$. $[S^\alpha_i,T^\beta_j]=0$ follows from $[\mathrm{M},\mathrm{G}]=0$ for $i=j$ and from CAR (even operators on different sites commute) for $i\ne j$. $\square$

---

## 3. Local Fermion Parity, Casimir Identity, and Physical Subspace

**Theorem 3.1 (Casimir–parity identity).** Let $\mathcal{Z}_i:=c_i^0c_i^1c_i^2c_i^3$ ($\mathcal{Z}_i^\dagger=\mathcal{Z}_i$, $\mathcal{Z}_i^2=\mathbb{1}$). Then
$$\mathbf{T}_i^2=\sum_{\gamma}(T_i^\gamma)^2=\frac{3}{8}\big(\mathbb{1}-\mathcal{Z}_i\big),\qquad \mathbf{S}_i^2=\frac{3}{8}\big(\mathbb{1}+\mathcal{Z}_i\big).$$

*Proof.* By (2.1), $T^1=\tfrac{\mathrm{i}}{4}(c^0c^1+c^2c^3)$. The two bilinears commute, $(c^0c^1)^2=(c^2c^3)^2=-\mathbb{1}$ and $(c^0c^1)(c^2c^3)=\mathcal{Z}$. Hence,
$$(T^1)^2=-\frac{1}{16}\big[(c^0c^1)^2+(c^2c^3)^2+2\,c^0c^1c^2c^3\big]=\frac{1}{8}(\mathbb{1}-\mathcal{Z}).$$
The cases $\gamma=2,3$ follow by relabeling $(1,2,3)\to(2,3,1)\to(3,1,2)$ [use (2.1) with $c^TG^\gamma c=2(c^0c^\gamma+c^{\gamma^+}c^{\gamma^-})$], so summing gives $\tfrac38(\mathbb{1}-\mathcal{Z})$. For $\mathbf{S}^2$ the same computation with $c^TM^\gamma c=2(c^0c^\gamma-c^{\gamma^+}c^{\gamma^-})$ flips the sign of the cross term. $\square$

*Consequence.* On $\mathcal{Z}_i=+1$ one has $\mathbf{S}^2=\tfrac34$ (a spin-$\tfrac12$ doublet) and $\mathbf{T}^2=0$; on $\mathcal{Z}_i=-1$ the roles are exchanged. The physical pseudospin is the $\mathcal{Z}_i=+1$ sector.

**Corollary 3.2 (Physical projection $\Leftrightarrow$ gauge neutrality).** For $\vert{}\psi\rangle\in\mathcal{F}$:
$$\mathcal{Z}_i\vert{}\psi\rangle=+\vert{}\psi\rangle\iff \mathbf{T}_i^2\vert{}\psi\rangle=0\iff T_i^\gamma\vert{}\psi\rangle=0\ \ \forall\gamma\iff c_i^T\mathrm{G}^\gamma c_i\vert{}\psi\rangle=0\ \ \forall\gamma.$$
Moreover each *single* condition $T_i^\gamma\vert{}\psi\rangle=0$ already implies $\mathcal{Z}_i\vert{}\psi\rangle=\vert{}\psi\rangle$, because $(T^\gamma)^2=\tfrac18(\mathbb{1}-\mathcal{Z})$ for every $\gamma$ separately.

*Proof.* The first equivalence is Theorem 3.1. For the second, $\langle\psi\vert{}\mathbf{T}^2\vert{}\psi\rangle=\sum_\gamma\Vert{}T^\gamma\psi\Vert{}^2$ because the $T^\gamma$ are Hermitian. The last equivalence is $T^\gamma=\tfrac{\mathrm{i}}{8}c^T\mathrm{G}^\gamma c$. The single-$\gamma$ statement uses the identity $(T^\gamma)^2=\tfrac18(\mathbb{1}-\mathcal{Z})$ proved above. $\square$

**Remark 3.3 (Exact versus averaged constraints).** Corollary 3.2 is a statement about *states*. Mean-field theory only controls *expectation values*, $\langle T^\gamma_i\rangle=0$, which is strictly weaker (e.g. a Gaussian state can have $\langle T^\gamma\rangle=0$ and $\langle\mathbf{T}^2\rangle>0$). Since the three conditions are exactly redundant but only approximately redundant on the mean-field level, they are imposed *symmetrically* through three multipliers $\lambda_i^\gamma$ (Freitas et al., 2024, Eq. (3) and Methods). The violation is monitored by the residual $g_i^\gamma$ of Note 03, §4.

**Proposition 3.4 (Pauli algebra holds modulo gauge generators).** For the Kitaev operators $\sigma_i^\gamma:=\mathrm{i}c_i^0c_i^\gamma$ one has, on the **full** Fock space and for $(\alpha,\beta,\gamma)$ cyclic,
$$\sigma_i^\alpha\sigma_i^\beta=\mathrm{i}\,\sigma_i^\gamma-4\mathrm{i}\,T_i^\gamma,\qquad \sigma_i^\gamma=2\big(S_i^\gamma+T_i^\gamma\big),\qquad (\sigma_i^\gamma)^2=\mathbb{1},\quad \{\sigma_i^\alpha,\sigma_i^\beta\}=2\delta^{\alpha\beta}.$$
Hence on the gauge-neutral subspace ($T_i^\gamma=0$) the $\sigma_i^\gamma$ obey the Pauli algebra $\sigma^\alpha\sigma^\beta=\mathrm{i}\epsilon^{\alpha\beta\gamma}\sigma^\gamma$ and $\sigma_i^\gamma=2S_i^\gamma$, while on $\mathcal{Z}_i=-1$ one finds $\sigma^\alpha\sigma^\beta=-\mathrm{i}\epsilon^{\alpha\beta\gamma}\sigma^\gamma$.

*Proof.* $\sigma^\alpha\sigma^\beta=(\mathrm{i}c^0c^\alpha)(\mathrm{i}c^0c^\beta)=c^\alpha c^\beta$ for $\alpha\ne\beta$. By (2.1), $c^\alpha c^\beta=\tfrac12c^T\mathrm{G}^\gamma c-c^0c^\gamma=-4\mathrm{i}T^\gamma+\mathrm{i}\sigma^\gamma$. The decomposition $\sigma^\gamma=\tfrac{\mathrm{i}}{4}c^T(\mathrm{M}^\gamma+\mathrm{G}^\gamma)c=2(S^\gamma+T^\gamma)$ follows from $\mathrm{N}^\gamma=\tfrac12(\mathrm{M}^\gamma+\mathrm{G}^\gamma)$ (Definition 4.1). $\square$

*(The statement "$\sigma_i^\gamma$ satisfy the spin algebra on the entire unprojected Fock space" that could be read into the previous version is false: $[\sigma^x,\sigma^y]=2\mathrm{i}\sigma^z$ fails off the gauge-neutral subspace.)*

---

## 4. The Kitaev Representation and the Extended Spin Hamiltonian

**Definition 4.1 (Kitaev matrices $\mathrm{N}^\gamma$).**
$$\sigma_i^\gamma=\mathrm{i}c_i^0c_i^\gamma=\frac{\mathrm{i}}{2}c_i^T\mathrm{N}^\gamma c_i,\qquad \mathrm{N}^\gamma=\frac{1}{2}(\mathrm{M}^\gamma+\mathrm{G}^\gamma),\qquad (\mathrm{N}^\gamma)_{\mu\nu}=\delta_{\mu0}\delta_{\nu\gamma}-\delta_{\nu0}\delta_{\mu\gamma}.$$
In Kronecker form $\mathrm{N}^1=e_{11}\otimes\epsilon$, $\mathrm{N}^2=\epsilon\otimes e_{11}$, $\mathrm{N}^3=\tfrac12(\sigma_x\otimes\epsilon+\epsilon\otimes\sigma_x)$, and $\mathrm{tr}((\mathrm{N}^\alpha)^T\mathrm{N}^\beta)=2\delta^{\alpha\beta}$. Two properties will be used repeatedly (both follow from the explicit form):
$$(\mathrm{N}^\gamma)^T=-\mathrm{N}^\gamma,\qquad \mathrm{tr}\,\mathrm{N}^\gamma=0\ \implies\ \mathrm{tr}\big(\mathrm{N}^{\gamma\,T}(V+\mathrm{i}\lambda\mathbb{1})\big)=\mathrm{tr}(\mathrm{N}^{\gamma\,T}V)\ \ \forall\lambda.\tag{4.1}$$
The Kitaev operators $\sigma_i^\gamma$ are special among linear combinations of $\mathrm{M}^\gamma,\mathrm{G}^\gamma$: each squares to $\mathbb{1}$ on the whole Fock space, and in the pure Kitaev limit the bond operators $\hat u_{ij}=\mathrm{i}c_i^\gamma c_j^\gamma$ commute with $H$, which gives exact solvability.

**Definition 4.2 (Extended $J$-$K$-$\Gamma$-$\Gamma'$-$D$-$\mathbf{h}$ Hamiltonian).** On $\Lambda$, with Zeeman field $\mathbf{h}=(h^1,h^2,h^3)^T$,
$$H=\frac12\sum_{i\in\mathcal{V}}\sum_{j\in\mathcal{V}_i}\sum_{a,b=1}^3\sigma_i^a\,\mathbf{J}_{ij}^{ab}\,\sigma_j^b-\sum_{i\in\mathcal{V}}\mathbf{h}\cdot\boldsymbol{\sigma}_i,\qquad \mathbf{J}_{ji}=\mathbf{J}_{ij}^T.\tag{4.2}$$
For a bond $\langle ij\rangle_\gamma\in\mathcal{E}$ (oriented $A\to B$) the exchange matrix is
$$\mathbf{J}^{ab}_{\langle ij\rangle_\gamma}=J_\gamma\delta^{ab}+K_\gamma\delta^{a\gamma}\delta^{b\gamma}+\Gamma_\gamma\big(\delta^{a\gamma^+}\delta^{b\gamma^-}+\delta^{a\gamma^-}\delta^{b\gamma^+}\big)+\Gamma'_\gamma\big(\delta^{a\gamma^+}\delta^{b\gamma}+\delta^{a\gamma}\delta^{b\gamma^+}+\delta^{a\gamma^-}\delta^{b\gamma}+\delta^{a\gamma}\delta^{b\gamma^-}\big)+D_\gamma\big(\delta^{a\gamma^+}\delta^{b\gamma^-}-\delta^{a\gamma^-}\delta^{b\gamma^+}\big).$$
Written out for a $z$-bond ($\gamma=3$, $\gamma^+=1$, $\gamma^-=2$), in the basis $(x,y,z)$:
$$\mathbf{J}_{\langle ij\rangle_z}=\begin{pmatrix}J&\Gamma+D&\Gamma'\\ \Gamma-D&J&\Gamma'\\ \Gamma'&\Gamma'&J+K\end{pmatrix};$$
$x$- and $y$-bonds follow by cyclic permutation of the indices. For $\Gamma'=D=0$ this is the $J$-$K$-$\Gamma$ matrix of Freitas et al. (2024), Eq. (1) (Rau et al., 2014; Winter et al., 2016). The $D_\gamma$ term reproduces $D_{ij}(\sigma^\alpha_i\sigma^\beta_j-\sigma^\beta_i\sigma^\alpha_j)$ of Eq. (13) in that paper. The prefactor $\tfrac12$ in (4.2) compensates the double counting of ordered pairs; the rule $\mathbf{J}_{ji}=\mathbf{J}_{ij}^T$ makes $H$ symmetric under relabeling of the two ends of a bond and is essential for the antisymmetric ($D$) part.

---

## 5. Exact Mean-Field Decoupling

**Definition 5.1 (Quartic Majorana Hamiltonian with gauge multipliers).** Insert $\sigma_i^\gamma=\tfrac{\mathrm{i}}{2}c_i^T\mathrm{N}^\gamma c_i$ and add site-dependent multipliers $\boldsymbol{\lambda}_i=(\lambda_i^1,\lambda_i^2,\lambda_i^3)^T$:
$$H_\lambda=\frac18\sum_{i}\sum_{j\in\mathcal{V}_i}\sum_{a,b}\big(\mathrm{i}c_i^T\mathrm{N}^ac_i\big)\mathbf{J}_{ij}^{ab}\big(\mathrm{i}c_j^T\mathrm{N}^bc_j\big)-\frac14\sum_i\sum_\gamma\Big(2h^\gamma\,\mathrm{i}c_i^T\mathrm{N}^\gamma c_i-\lambda_i^\gamma\,\mathrm{i}c_i^T\mathrm{G}^\gamma c_i\Big).$$
Since $\tfrac{\mathrm{i}}{4}c^T\mathrm{G}^\gamma c=2T^\gamma$, the multiplier term is $+2\sum_{i,\gamma}\lambda_i^\gamma T_i^\gamma$: it vanishes identically on the physical (gauge-neutral) subspace and therefore does not change physical matrix elements.

**Definition 5.2 (Mean-field data).** For a Gaussian (quasi-free) state $\rho$ define the real matrices
$$\mathrm{U}_{ij}^{\mu\nu}:=\langle\mathrm{i}c_i^\mu c_j^\nu\rangle_\rho\ \ (i\ne j),\qquad \mathrm{V}_i^{\mu\nu}:=\langle\mathrm{i}c_i^\mu c_i^\nu\rangle_\rho-\mathrm{i}\,\delta^{\mu\nu}.$$
They satisfy $\mathrm{U}_{ji}=-\mathrm{U}_{ij}^T$ and $\mathrm{V}_i^T=-\mathrm{V}_i$ (on-site: $\{c^\mu,c^\nu\}=2\delta^{\mu\nu}$ gives $\langle\mathrm{i}c^\mu c^\nu\rangle+\langle\mathrm{i}c^\nu c^\mu\rangle=2\mathrm{i}\delta^{\mu\nu}$). Write $\Gamma=(\mathrm{U},\mathrm{V})$ for the full set.

> **Convention difference to the paper.** Freitas et al. (2024), Eq. (5), define $\mathrm{V}_i^{\mu\nu}=\langle\mathrm{i}c_i^\mu c_i^\nu\rangle$ *without* subtracting $\mathrm{i}\delta^{\mu\nu}$, so that $\mathrm{V}+\mathrm{V}^T=2\mathrm{i}\delta$. Our $\mathrm{V}$ is the real, skew-symmetric part. By (4.1) this changes neither $\mathrm{B}_i$ nor $\mathrm{C}$, which depend on $\mathrm{V}$ only through $\mathrm{tr}(\mathrm{V}^T\mathrm{N}^\beta)$. The physical content is $\langle\sigma_i^\gamma\rangle=\mathrm{V}_i^{0\gamma}$ and $\mathrm{V}_i^{\alpha\beta}=-\mathrm{V}_i^{0\gamma}$ on the gauge-neutral subspace.

**Lemma 5.3 (Wick contraction).** Let $Q_i^a:=\mathrm{i}c_i^T\mathrm{N}^ac_i=2\sigma_i^a$. In a Gaussian state, for $i\ne j$,
$$\langle Q_i^a\rangle=\mathrm{tr}\big((\mathrm{V}_i)^T\mathrm{N}^a\big)=2\mathrm{V}_i^{0a},\qquad
\langle Q_i^aQ_j^b\rangle=\langle Q_i^a\rangle\langle Q_j^b\rangle+2\,\mathrm{tr}\big(\mathrm{U}_{ij}^T\mathrm{N}^a\mathrm{U}_{ij}\mathrm{N}^b\big).$$

*Proof.* The first relation is (4.1). The second is Wick's theorem for four Majorana operators: the fully disconnected term gives the product; the two connected pairings $\langle c_i^\mu c_j^\rho\rangle\langle c_i^\nu c_j^\sigma\rangle$ and its mirror contribute equal amounts (use $\mathrm{N}^T=-\mathrm{N}$), and the two ordering signs combine with the factor $\mathrm{i}^2$ to give $+2\,\mathrm{tr}(\mathrm{U}^T\mathrm{N}^a\mathrm{U}\mathrm{N}^b)$. $\square$ (Verified to $10^{-15}$ in `test_wick_and_hfb_two_site`.)

**Definition 5.4 (Hartree–Fock energy functional).** For a bond-translation-invariant or inhomogeneous configuration $\Gamma$,
$$E[\Gamma]:=\frac18\sum_i\sum_{j\in\mathcal{V}_i}\sum_{a,b}\mathbf{J}_{ij}^{ab}\Big[\mathrm{tr}(\mathrm{V}_i^T\mathrm{N}^a)\mathrm{tr}(\mathrm{V}_j^T\mathrm{N}^b)+2\,\mathrm{tr}(\mathrm{U}_{ij}^T\mathrm{N}^a\mathrm{U}_{ij}\mathrm{N}^b)\Big]-\frac12\sum_{i,\gamma}\Big(h^\gamma\,\mathrm{tr}(\mathrm{V}_i^T\mathrm{N}^\gamma)-\frac12\lambda_i^\gamma\,\mathrm{tr}((\mathrm{G}^\gamma)^T\mathrm{V}_i)\Big).$$
By Lemma 5.3, $E[\Gamma]=\langle H_\lambda\rangle_\rho$ for every Gaussian state $\rho$ with data $\Gamma$.

**Proposition 5.5 (Majorana HFB Hamiltonian).** For fixed $\Gamma$ define the $4\times4$ real matrices (Freitas et al., 2024, Eqs. (25)–(27))
$$\mathrm{A}_{ij}=2\sum_{a,b}\mathbf{J}_{ij}^{ab}\,\mathrm{N}^a\mathrm{U}_{ij}\mathrm{N}^b,$$
$$\mathrm{B}_i=\sum_{j\in\mathcal{V}_i}\sum_{a,b}\mathbf{J}_{ij}^{ab}\,\mathrm{N}^a\,\mathrm{tr}(\mathrm{V}_j^T\mathrm{N}^b)+\sum_\gamma\big(\lambda_i^\gamma\mathrm{G}^\gamma-2h^\gamma\mathrm{N}^\gamma\big),$$
$$\mathrm{C}=\frac18\sum_i\sum_{j\in\mathcal{V}_i}\sum_{a,b}\mathbf{J}_{ij}^{ab}\Big[\mathrm{tr}(\mathrm{V}_i^T\mathrm{N}^a)\mathrm{tr}(\mathrm{V}_j^T\mathrm{N}^b)+2\,\mathrm{tr}(\mathrm{U}_{ij}^T\mathrm{N}^a\mathrm{U}_{ij}\mathrm{N}^b)\Big],$$
and
$$H_{\mathrm{MF}}[\Gamma]=\sum_i\sum_{j\in\mathcal{V}_i}\frac{\mathrm{i}}{4}c_i^T\mathrm{A}_{ij}c_j+\sum_i\frac{\mathrm{i}}{4}c_i^T\mathrm{B}_ic_i-\mathrm{C}.$$
Then:

1. $\mathrm{A}_{ji}=-\mathrm{A}_{ij}^T$ and $\mathrm{B}_i^T=-\mathrm{B}_i$. Hence $H_{\mathrm{MF}}$ is Hermitian and each bond contributes $\tfrac{\mathrm{i}}{2}c_i^T\mathrm{A}_{ij}c_j$ (both orientations are equal).
2. $\langle H_{\mathrm{MF}}[\Gamma]\rangle_\rho=E[\Gamma]=\langle H_\lambda\rangle_\rho$ for every Gaussian $\rho$ with data $\Gamma$ (so $\mathrm{C}$ is the double-counting correction).
3. (*Tangent property.*) For *any* data $\Gamma'=(\mathrm{U}',\mathrm{V}')$,
$$\langle H_{\mathrm{MF}}[\Gamma]\rangle_{\Gamma'}=E[\Gamma]+\mathrm{d}E_\Gamma[\Gamma'-\Gamma],$$
i.e. $H_{\mathrm{MF}}[\Gamma]$ is the first-order Taylor expansion of $E$ about $\Gamma$ evaluated on $\Gamma'$ (here $\langle\tfrac{\mathrm{i}}{4}c_i^T\mathrm{A}c_j\rangle_{\Gamma'}=\tfrac14\mathrm{tr}(\mathrm{A}^T\mathrm{U}'_{ij})$, $\langle\tfrac{\mathrm{i}}{4}c_i^T\mathrm{B}c_i\rangle_{\Gamma'}=\tfrac14\mathrm{tr}(\mathrm{B}^T\mathrm{V}'_i)$).

*Proof.* (1) $\mathrm{A}_{ji}=2\sum\mathbf{J}_{ij}^{ba}\mathrm{N}^a(-\mathrm{U}_{ij}^T)\mathrm{N}^b=-2\sum\mathbf{J}_{ij}^{ba}(\mathrm{N}^b\mathrm{U}_{ij}\mathrm{N}^a)^T=-\mathrm{A}_{ij}^T$, using $\mathrm{N}^T=-\mathrm{N}$ twice. $\mathrm{B}_i^T=-\mathrm{B}_i$ because $\mathrm{N}^\gamma,\mathrm{G}^\gamma$ are skew. For two sites, $c_j^T(-\mathrm{A}^T)c_i=+c_i^T\mathrm{A}c_j$ (distinct-site Majoranas anticommute), so the two orientations add. (2) and (3): $E$ is a quadratic polynomial in $\Gamma$; $\mathrm{A}_{ij}$ and $\mathrm{B}_i$ are, up to the factors $\tfrac14\cdot 2$ and $\tfrac14$ from the pairing of orientations, the partial derivatives $\partial E/\partial\mathrm{U}_{ij}$ and $\partial E/\partial\mathrm{V}_i$ obtained by differentiating Definition 5.4 term by term, and $\mathrm{C}$ is fixed by requiring equality at $\Gamma'=\Gamma$. The multiplier and Zeeman terms are linear in $\mathrm{V}$ and drop out of $\mathrm{C}$. This was verified numerically for a Gaussian two-site state and for random (not necessarily physical) $\Gamma'$ in `test_wick_and_hfb_two_site`. $\square$

**Corollary 5.6 (Self-consistency = stationarity).** Let $\rho_0[\Gamma]$ denote a ground state of $H_{\mathrm{MF}}[\Gamma]$. Because $\rho\mapsto\langle H_{\mathrm{MF}}[\Gamma]\rangle_\rho$ is a *linear* functional of $\Gamma'$ minimized at $\rho_0[\Gamma]$, the self-consistency condition $\Gamma=\Gamma[\rho_0[\Gamma]]$ is exactly the first-order stationarity condition of $E$ on the manifold of Gaussian states (a necessary condition for a variational minimum, not a sufficient one). In the pure Kitaev limit the mean-field energy is exact (Note 03, Prop. 3.3).

---

## Revision log (Note 01)

1. **Lemma 2.1:** missing factor $\mathrm{i}$ restored: $[\mathcal{O}[A],\mathcal{O}[B]]=\mathrm{i}\,\mathcal{O}[[A,B]]$. A full proof was added.
2. **Prop. 2.3 / Thm. 3.1:** $S^\gamma,T^\gamma$ now defined with prefactor $\tfrac{\mathrm{i}}{8}$ (spin-$\tfrac12$), which is the normalization for which the stated commutators and the Casimir $\tfrac38(1\mp\mathcal{Z})$ hold; the *chirality* of the gauge algebra ($-\mathrm{i}\epsilon$) is now stated. The old proof's intermediate formula $2T^1=\tfrac{\mathrm{i}}{2}(\cdots)$ was inconsistent with its own definition.
3. **Cor. 3.2:** equivalence proved (sum of squares of Hermitian operators), single-$\gamma$ redundancy shown, and the exact-versus-averaged distinction (Remark 3.3) added.
4. **New Prop. 3.4:** the Pauli algebra of $\sigma^\gamma=\mathrm{i}c^0c^\gamma$ holds only on the gauge-neutral subspace; exact operator identity given.
5. **Def. 4.2:** dummy indices renamed ($a,b$ instead of primed $\alpha',\beta'$ clashing with the cyclic triple), explicit $3\times3$ matrix for the $z$-bond added.
6. **§5:** convention difference of $\mathrm{V}$ to the paper made explicit; Wick lemma, energy functional, tangent property of $H_{\mathrm{MF}}$ and the stationarity corollary added (the old note only quoted Eqs. (25)–(27)).
7. Symbol clashes resolved (§0 table).

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
