# Note 01: Algebraic Foundations of the Majorana Parton Representation, $\mathfrak{so}(4)$ Gauge Structure, and Mean-Field Decoupling

> **Status of this revision.** All statements below were re-derived and cross-checked numerically with exact diagonalization on one and two sites (Jordan–Wigner representation). Notation introduced in §0 is shared across Notes 01–03.

---

## 0. Conventions and Notation (shared by Notes 01–03)

**Index sets.** Greek indices from the middle of the alphabet, $\mu,\nu,\rho,\sigma \in \{0,1,2,3\}$, label *Majorana flavors*; Greek indices from the start, $\alpha,\beta,\gamma \in \{1,2,3\} \equiv \{x,y,z\}$, label *Cartesian components* and, for $\gamma$, also the *type of a nearest-neighbor bond*. Dummy indices of $3\times 3$ exchange matrices are written $a,b \in \{1,2,3\}$. For a bond type $\gamma$ we write $(\gamma^+,\gamma^-) := (\gamma+1,\gamma+2) \bmod 3$ (values in $\{1,2,3\}$), so that $(\gamma^+,\gamma^-,\gamma)$ is a cyclic permutation of $(1,2,3)$. $\epsilon^{\alpha\beta\gamma}$ is the Levi-Civita symbol with $\epsilon^{123}=+1$.

**Disambiguated Symbols Table.**

| Quantity | Symbol (this revision) | Replaces / clashes with |
|---|---|---|
| local fermion parity $c^0c^1c^2c^3$ | $\mathcal{Z}_i$ (called $D_i$ in Freitas et al., 2024) | $D_\gamma$ = Dzyaloshinskii–Moriya coupling |
| $\mathrm{diag}(1,0)$ | $e_{11}$ | $\Pi_+$ (clashes with spectral projectors $\Pi^\pm(\mathbf{k})$) |
| Bogoliubov quasiparticle operators | $\eta_{n,\mathbf{k}}$ | $\gamma_{n,\mathbf{k}}$ (clashes with bond type $\gamma$) |
| supercell crystal momentum | $\mathbf{q}=(q_1,q_2)$ | $(K_1,K_2)$ (clashes with the Kitaev coupling $K$) |
| iteration counter (Note 03) | $\ell$ | $k$ (clashes with momentum $\mathbf{k}$) |
| DIIS weights / Lagrange multiplier / Gram matrix | $\xi_j$ / $\nu$ / $\Omega$ | $c_j$ / $\mu$ / $B$ (clash with Majorana $c$, flavor $\mu$, on-site matrix $\mathrm{B}$) |
| eigenvalues of the Bloch matrix $\mathcal{H}(\mathbf{k})$ | $\mathcal{E}_n(\mathbf{k})$ ("Bloch energies") | $E_n(\mathbf{k})$ |
| physical excitation energies | $\varepsilon_n(\mathbf{k}) = 4\,\mathcal{E}_n(\mathbf{k})$ (Note 02) | — |
| Pauli-matrix (Kronecker) factors | $I_2,\sigma_x,\sigma_z,\epsilon=\mathrm{i}\sigma_y$ | — |

**Normalization of spins.** All spin operators are Pauli operators, $\sigma^\gamma_i$ with $(\sigma^\gamma_i)^2=\mathbb{1}$ (not $S=\sigma/2$). Consequently the exchange matrices $\mathbf{J}$ are in units where the pure Kitaev Hamiltonian reads $H_K=K\sum_{\langle ij\rangle_\gamma}\sigma_i^\gamma\sigma_j^\gamma$, as in Kitaev (2006) and Freitas et al. (2024). The physical regime of $\alpha$-RuCl$_3$ is $K<0$ (ferromagnetic).

**Real/imaginary conventions.** $c_i^\mu$ are Hermitian operators. A *real skew-symmetric* matrix $A\in\mathfrak{so}(4)$ yields the *anti-Hermitian* operator $c^TAc$ (because $(c^TAc)^\dagger=c^TA^Tc=-c^TAc$) and the *Hermitian* operator $\mathrm{i}\,c^TAc$. All mean-field parameters are expectation values of Hermitian operators and are therefore real.

---

## 1. Lattice Graph and the Self-Dual CAR Algebra

**Definition (Lattice Graph)** <!-- claim: CLM-M1-ALG-001 -->
Let $\Lambda=(\mathcal{V},\mathcal{E})$ be a finite or periodic honeycomb graph with vertex set $\mathcal{V}$, $N_{\mathrm{sites}}=\vert{}\mathcal{V}\vert{}$. $\mathcal{V}=\mathcal{V}_A\sqcup\mathcal{V}_B$ is bipartite. The edge set $\mathcal{E}$ consists of *directed* nearest-neighbor bonds $\langle ij\rangle_\gamma$ with $i\in\mathcal{V}_A$, $j\in\mathcal{V}_B$, of type $\gamma\in\{1,2,3\}$; every vertex belongs to exactly one bond of each type (periodic case)[cite: 2]. $\mathcal{V}_i:=\{j:\langle ij\rangle_\gamma\in\mathcal{E}\ \text{or}\ \langle ji\rangle_\gamma\in\mathcal{E}\}$ is the set of the three nearest neighbors of $i$[cite: 2]. Sums $\sum_{i}\sum_{j\in\mathcal{V}_i}$ run over **ordered** pairs, i.e. every bond is visited twice (once in each orientation)[cite: 4].

**Definition (Self-Dual CAR Algebra)** <!-- claim: CLM-M1-ALG-002 -->
To each vertex $i$ we assign four Hermitian generators $c_i=(c_i^0,c_i^1,c_i^2,c_i^3)^T$ of the unital algebra $\mathcal{A}_{\mathrm{CAR}}$ with:
$$(c_i^\mu)^\dagger=c_i^\mu,\qquad \{c_i^\mu,c_j^\nu\}=2\,\delta_{ij}\delta^{\mu\nu}\,\mathbb{1}.$$
$\mathcal{A}_{\mathrm{CAR}}$ is the complex Clifford algebra on $4N_{\mathrm{sites}}$ generators, $\mathcal{A}_{\mathrm{CAR}}\cong \mathrm{Mat}(2^{2N_{\mathrm{sites}}},\mathbb{C})$.

**Remark (Local Fock Space)**
The two operators $d_{i,1}=\tfrac12(c_i^0+\mathrm{i}c_i^1)$ and $d_{i,2}=\tfrac12(c_i^2+\mathrm{i}c_i^3)$ satisfy canonical complex anticommutation relations: $\{d_{i,1},d_{i,1}^\dagger\}=1$ (Burnell and Nayak, 2011; You et al., 2012)[cite: 2, 4]. The local Fock space is $\mathcal{F}_i\cong\mathbb{C}^4$, twice as large as the physical pseudospin space $\mathbb{C}^2$[cite: 4].

---

## 2. Lie Algebra $\mathfrak{so}(4)\cong\mathfrak{su}(2)_{\mathrm{spin}}\oplus\mathfrak{su}(2)_{\mathrm{gauge}}$

For $A\in\mathfrak{so}(4)$ (real, $A^T=-A$) and a vertex $i$ define the Hermitian on-site operator:
$$\mathcal{O}_i[A]:=\frac{\mathrm{i}}{4}\,c_i^TAc_i=\frac{\mathrm{i}}{4}\sum_{\mu\nu}A_{\mu\nu}\,c_i^\mu c_i^\nu.$$

**Lemma (Commutator of Quadratic Forms)** <!-- claim: CLM-M1-ALG-003 -->
For $A,B\in\mathfrak{so}(4)$,
$$[\mathcal{O}_i[A],\mathcal{O}_i[B]]=\mathrm{i}\,\mathcal{O}_i\big[[A,B]\big].$$

*Proof.*
From $\{c^\mu,c^\nu\}=2\delta^{\mu\nu}$, the product of bilinears expands as:
$$[c^a c^b, c^c c^d] = 2\big(\delta^{bc}c^a c^d - \delta^{ac}c^b c^d + \delta^{bd}c^c c^a - \delta^{ad}c^c c^b\big),$$
which contracts to $[c^TAc,c^TBc]=4\,c^T[A,B]c$. Hence:
$$[\mathcal{O}[A],\mathcal{O}[B]]=\left(\frac{\mathrm{i}}{4}\right)^2 \cdot 4\,c^T[A,B]c = -\frac{1}{4}c^T[A,B]c = \mathrm{i}\cdot\frac{\mathrm{i}}{4}c^T[A,B]c = \mathrm{i}\,\mathcal{O}\big[[A,B]\big].\quad\square$$

**Definition (Basis of $\mathfrak{so}(4)$)** <!-- claim: CLM-M1-ALG-004 -->
With $I_2$, $\sigma_x$, $\sigma_z$, $\epsilon=\mathrm{i}\sigma_y=\begin{psmallmatrix}0&1\\-1&0\end{psmallmatrix}$ and $e_{11}=\mathrm{diag}(1,0)$, define the gauge generators $\mathrm{G}^\gamma$ and chiral-spin generators $\mathrm{M}^\gamma$ ('t Hooft, 1976; Fu et al., 2018)[cite: 2, 4]:
$$\mathrm{G}^1=I_2\otimes\epsilon,\quad \mathrm{G}^2=\epsilon\otimes\sigma_z,\quad \mathrm{G}^3=\epsilon\otimes\sigma_x,$$
$$\mathrm{M}^1=\sigma_z\otimes\epsilon,\quad \mathrm{M}^2=\epsilon\otimes I_2,\quad \mathrm{M}^3=\sigma_x\otimes\epsilon.$$
Consequently,
$$c^T\mathrm{G}^\gamma c=c^0c^\gamma-c^\gamma c^0+\sum_{\alpha\beta}\epsilon^{\alpha\beta\gamma}c^\alpha c^\beta=2\Big(c^0c^\gamma+c^{\gamma^+}c^{\gamma^-}\Big).\tag{2.1}$$

**Proposition (Orthogonal $\mathfrak{su}(2)$ Subalgebras)** <!-- claim: CLM-M1-ALG-005 -->
1. *Matrix level:*
   $$[\mathrm{G}^\alpha,\mathrm{G}^\beta]=-2\sum_\gamma\epsilon^{\alpha\beta\gamma}\mathrm{G}^\gamma,\quad [\mathrm{M}^\alpha,\mathrm{M}^\beta]=+2\sum_\gamma\epsilon^{\alpha\beta\gamma}\mathrm{M}^\gamma,\quad [\mathrm{M}^\alpha,\mathrm{G}^\beta]=0,$$
   $$\operatorname{tr}\big((\mathrm{G}^\alpha)^T\mathrm{G}^\beta\big)=\operatorname{tr}\big((\mathrm{M}^\alpha)^T\mathrm{M}^\beta\big)=4\delta^{\alpha\beta},\qquad \operatorname{tr}\big((\mathrm{M}^\alpha)^T\mathrm{G}^\beta\big)=0.$$
2. *Operator level:* With the spin-$\tfrac12$ normalization
   $$S_i^\gamma:=\frac{\mathcal{O}_i[\mathrm{M}^\gamma]}{2}=\frac{\mathrm{i}}{8}c_i^T\mathrm{M}^\gamma c_i,\qquad T_i^\gamma:=\frac{\mathcal{O}_i[\mathrm{G}^\gamma]}{2}=\frac{\mathrm{i}}{8}c_i^T\mathrm{G}^\gamma c_i,$$
   the following commutation relations hold:
   $$[S_i^\alpha,S_i^\beta]=+\mathrm{i}\sum_\gamma\epsilon^{\alpha\beta\gamma}S_i^\gamma,\qquad [T_i^\alpha,T_i^\beta]=-\mathrm{i}\sum_\gamma\epsilon^{\alpha\beta\gamma}T_i^\gamma,\qquad [S_i^\alpha,T_j^\beta]=0.$$

---

## 3. Local Fermion Parity, Casimir Identity, and Physical Subspace

**Theorem (Casimir–Parity Identity)** <!-- claim: CLM-M1-ALG-006 -->
Let $\mathcal{Z}_i:=c_i^0c_i^1c_i^2c_i^3$ ($\mathcal{Z}_i^\dagger=\mathcal{Z}_i$, $\mathcal{Z}_i^2=\mathbb{1}$). Then:
$$\mathbf{T}_i^2=\sum_{\gamma}(T_i^\gamma)^2=\frac{3}{8}\big(\mathbb{1}-\mathcal{Z}_i\big),\qquad \mathbf{S}_i^2=\frac{3}{8}\big(\mathbb{1}+\mathcal{Z}_i\big).$$

*Proof.*
By (2.1), $T^1=\tfrac{\mathrm{i}}{4}(c^0c^1+c^2c^3)$[cite: 4]. The bilinears commute with $(c^0c^1)^2=(c^2c^3)^2=-\mathbb{1}$ and $(c^0c^1)(c^2c^3)=\mathcal{Z}$[cite: 4]. Hence:
$$(T^1)^2=-\frac{1}{16}\big[(c^0c^1)^2+(c^2c^3)^2+2\,c^0c^1c^2c^3\big]=\frac{1}{8}(\mathbb{1}-\mathcal{Z}).$$
Summing over cyclic permutations yields $\mathbf{T}^2 = \frac{3}{8}(\mathbb{1}-\mathcal{Z})$. $\square$

**Corollary (Physical Projection $\Leftrightarrow$ Gauge Neutrality)** <!-- claim: CLM-M1-ALG-007 -->
For $\vert{}\psi\rangle\in\mathcal{F}$:
$$\mathcal{Z}_i\vert{}\psi\rangle=+\vert{}\psi\rangle\iff \mathbf{T}_i^2\vert{}\psi\rangle=0\iff T_i^\gamma\vert{}\psi\rangle=0\ \ \forall\gamma\iff c_i^T\mathrm{G}^\gamma c_i\vert{}\psi\rangle=0\ \ \forall\gamma.$$
Each single condition $T_i^\gamma\vert{}\psi\rangle=0$ separately implies $\mathcal{Z}_i\vert{}\psi\rangle=+\vert{}\psi\rangle$.

**Proposition (Pauli Algebra Holds Modulo Gauge Generators)** <!-- claim: CLM-M1-ALG-008 -->
For $\sigma_i^\gamma:=\mathrm{i}c_i^0c_i^\gamma$ and $(\alpha,\beta,\gamma)$ cyclic:
$$\sigma_i^\alpha\sigma_i^\beta=\mathrm{i}\,\sigma_i^\gamma-4\mathrm{i}\,T_i^\gamma,\qquad \sigma_i^\gamma=2\big(S_i^\gamma+T_i^\gamma\big),\qquad (\sigma_i^\gamma)^2=\mathbb{1},\quad \{\sigma_i^\alpha,\sigma_i^\beta\}=2\delta^{\alpha\beta}.$$
On the gauge-neutral subspace ($T_i^\gamma=0$), the $\sigma_i^\gamma$ satisfy the standard $\mathfrak{su}(2)$ algebra.

---

## 4. The Kitaev Representation and the Extended Spin Hamiltonian

**Definition (Kitaev Matrices $\mathrm{N}^\gamma$)** <!-- claim: CLM-M1-ALG-009 -->
$$\sigma_i^\gamma=\mathrm{i}c_i^0c_i^\gamma=\frac{\mathrm{i}}{2}c_i^T\mathrm{N}^\gamma c_i,\qquad \mathrm{N}^\gamma=\frac{1}{2}(\mathrm{M}^\gamma+\mathrm{G}^\gamma),\qquad (\mathrm{N}^\gamma)_{\mu\nu}=\delta_{\mu0}\delta_{\nu\gamma}-\delta_{\nu0}\delta_{\mu\gamma}.$$
In Kronecker form $\mathrm{N}^1=e_{11}\otimes\epsilon$, $\mathrm{N}^2=\epsilon\otimes e_{11}$, $\mathrm{N}^3=\tfrac12(\sigma_x\otimes\epsilon+\epsilon\otimes\sigma_x)$, with:
$$(\mathrm{N}^\gamma)^T=-\mathrm{N}^\gamma,\qquad \operatorname{tr}\!\left(\mathrm{N}^\gamma\right)=0\implies \operatorname{tr}\!\left(\mathrm{N}^{\gamma\,T}(V+\mathrm{i}\lambda\mathbb{1})\right)=\operatorname{tr}\!\left(\mathrm{N}^{\gamma\,T}V\right)\quad \forall\lambda.\tag{4.1}$$

**Definition (Extended $J$-$K$-$\Gamma$-$\Gamma'$-$D$-$\mathbf{h}$ Hamiltonian)** <!-- claim: CLM-M1-ALG-010 -->
On $\Lambda$, with Zeeman field $\mathbf{h}=(h^1,h^2,h^3)^T$:
$$H=\frac12\sum_{i\in\mathcal{V}}\sum_{j\in\mathcal{V}_i}\sum_{a,b=1}^3\sigma_i^a\,\mathbf{J}_{ij}^{ab}\,\sigma_j^b-\sum_{i\in\mathcal{V}}\mathbf{h}\cdot\boldsymbol{\sigma}_i,\qquad \mathbf{J}_{ji}=\mathbf{J}_{ij}^T.\tag{4.2}$$
On a directed bond $\langle ij\rangle_\gamma\in\mathcal{E}$ (oriented $A\to B$):
$$\mathbf{J}^{ab}_{\langle ij\rangle_\gamma}=J_\gamma\delta^{ab}+K_\gamma\delta^{a\gamma}\delta^{b\gamma}+\Gamma_\gamma\big(\delta^{a\gamma^+}\delta^{b\gamma^-}+\delta^{a\gamma^-}\delta^{b\gamma^+}\big)+\Gamma'_\gamma\big(\delta^{a\gamma^+}\delta^{b\gamma}+\delta^{a\gamma}\delta^{b\gamma^+}+\delta^{a\gamma^-}\delta^{b\gamma}+\delta^{a\gamma}\delta^{b\gamma^-}\big)+D_\gamma\big(\delta^{a\gamma^+}\delta^{b\gamma^-}-\delta^{a\gamma^-}\delta^{b\gamma^+}\big).$$

---

## 5. Exact Mean-Field Decoupling

**Definition (Quartic Majorana Hamiltonian with Multipliers)** <!-- claim: CLM-M1-ALG-011 -->
$$H_\lambda=\frac18\sum_{i}\sum_{j\in\mathcal{V}_i}\sum_{a,b}\big(\mathrm{i}c_i^T\mathrm{N}^ac_i\big)\mathbf{J}_{ij}^{ab}\big(\mathrm{i}c_j^T\mathrm{N}^bc_j\big)-\frac14\sum_i\sum_\gamma\Big(2h^\gamma\,\mathrm{i}c_i^T\mathrm{N}^\gamma c_i-\lambda_i^\gamma\,\mathrm{i}c_i^T\mathrm{G}^\gamma c_i\Big).$$

**Definition (Mean-Field Correlators)** <!-- claim: CLM-M1-ALG-012 -->
For a Gaussian quasi-free state $\rho$:
$$\mathrm{U}_{ij}^{\mu\nu}:=\langle\mathrm{i}c_i^\mu c_j^\nu\rangle_\rho\quad (i\ne j),\qquad \mathrm{V}_i^{\mu\nu}:=\langle\mathrm{i}c_i^\mu c_i^\nu\rangle_\rho-\mathrm{i}\,\delta^{\mu\nu}.$$
They satisfy $\mathrm{U}_{ji}=-\mathrm{U}_{ij}^T$ and $\mathrm{V}_i^T=-\mathrm{V}_i$[cite: 4].

**Lemma (Wick Contraction Identity)** <!-- claim: CLM-M1-ALG-013 -->
Let $Q_i^a:=\mathrm{i}c_i^T\mathrm{N}^ac_i=2\sigma_i^a$. In a Gaussian state, for $i\ne j$[cite: 4]:
$$\langle Q_i^a\rangle=\operatorname{tr}\!\left(\mathrm{V}_i^T\mathrm{N}^a\right)=2\mathrm{V}_i^{0a},\qquad \langle Q_i^aQ_j^b\rangle=\langle Q_i^a\rangle\langle Q_j^b\rangle+2\,\operatorname{tr}\!\left(\mathrm{U}_{ij}^T\mathrm{N}^a\mathrm{U}_{ij}\mathrm{N}^b\right).$$

**Definition (Hartree–Fock Energy Functional)** <!-- claim: CLM-M1-ALG-014 -->
$$E[\Gamma]:=\frac18\sum_i\sum_{j\in\mathcal{V}_i}\sum_{a,b}\mathbf{J}_{ij}^{ab}\Big[\operatorname{tr}\!\left(\mathrm{V}_i^T\mathrm{N}^a\right)\operatorname{tr}\!\left(\mathrm{V}_j^T\mathrm{N}^b\right)+2\,\operatorname{tr}\!\left(\mathrm{U}_{ij}^T\mathrm{N}^a\mathrm{U}_{ij}\mathrm{N}^b\right)\Big]-\frac12\sum_{i,\gamma}\Big(h^\gamma\,\operatorname{tr}\!\left(\mathrm{V}_i^T\mathrm{N}^\gamma\right)-\frac12\lambda_i^\gamma\,\operatorname{tr}\!\left((\mathrm{G}^\gamma)^T\mathrm{V}_i\right)\Big).$$

**Proposition (Majorana HFB Hamiltonian)** <!-- claim: CLM-M1-ALG-015 -->
For fixed $\Gamma = (\mathrm{U}, \mathrm{V})$, the quadratic mean-field Hamiltonian is:
$$H_{\mathrm{MF}}[\Gamma]=\sum_i\sum_{j\in\mathcal{V}_i}\frac{\mathrm{i}}{4}c_i^T\mathrm{A}_{ij}c_j+\sum_i\frac{\mathrm{i}}{4}c_i^T\mathrm{B}_ic_i-\mathrm{C},$$
where:
$$\mathrm{A}_{ij}=2\sum_{a,b}\mathbf{J}_{ij}^{ab}\,\mathrm{N}^a\mathrm{U}_{ij}\mathrm{N}^b,$$
$$\mathrm{B}_i=\sum_{j\in\mathcal{V}_i}\sum_{a,b}\mathbf{J}_{ij}^{ab}\,\mathrm{N}^a\,\operatorname{tr}\!\left(\mathrm{V}_j^T\mathrm{N}^b\right)+\sum_\gamma\big(\lambda_i^\gamma\mathrm{G}^\gamma-2h^\gamma\mathrm{N}^\gamma\big),$$
$$\mathrm{C}=\frac18\sum_i\sum_{j\in\mathcal{V}_i}\sum_{a,b}\mathbf{J}_{ij}^{ab}\Big[\operatorname{tr}\!\left(\mathrm{V}_i^T\mathrm{N}^a\right)\operatorname{tr}\!\left(\mathrm{V}_j^T\mathrm{N}^b\right)+2\,\operatorname{tr}\!\left(\mathrm{U}_{ij}^T\mathrm{N}^a\mathrm{U}_{ij}\mathrm{N}^b\right)\Big].$$
Then:
1. $\mathrm{A}_{ji}=-\mathrm{A}_{ij}^T$ and $\mathrm{B}_i^T=-\mathrm{B}_i$, ensuring $H_{\mathrm{MF}}$ is Hermitian[cite: 4].
2. $\langle H_{\mathrm{MF}}[\Gamma]\rangle_\rho=E[\Gamma]$ for every Gaussian state $\rho$ with data $\Gamma$[cite: 4].
3. (*Tangent property:*) For any trial configuration $\Gamma'$,
   $$\langle H_{\mathrm{MF}}[\Gamma]\rangle_{\Gamma'}=E[\Gamma]+\mathrm{d}E_\Gamma[\Gamma'-\Gamma].$$

**Corollary (Self-Consistency as Variational Stationarity)** <!-- claim: CLM-M1-ALG-016 -->
The self-consistency condition $\Gamma = \Gamma[\rho_0[\Gamma]]$ is the exact first-order stationarity condition of $E[\Gamma]$ on the manifold of Gaussian states[cite: 4].

---

## References

* **Burnell, F. J. and Nayak, C.** (2011). *SU(2) slave fermion approach to the Kitaev honeycomb lattice model*. Physical Review B **84**, 125125. DOI: [10.1103/PhysRevB.84.125125](https://doi.org/10.1103/PhysRevB.84.125125)[cite: 2]. arXiv: [1104.5485](https://arxiv.org/abs/1104.5485)[cite: 2].
* **Freitas, L. R. D.** (2024). *Majorana mean-field theory with $\mathbb{Z}_2$ vortices: electrical tests in the extended Kitaev model*. Ph.D. thesis, Universidade Federal do Rio Grande do Norte (UFRN) / Heinrich-Heine-Universität Düsseldorf[cite: 2].
* **Freitas, L. R. D., Bauer, T., Egger, R., and Pereira, R. G.** (2024). *Electric polarization near vortices in the extended Kitaev model*. npj Quantum Materials **9**, 33. DOI: [10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5)[cite: 2]. arXiv: [2311.15750](https://arxiv.org/abs/2311.15750)[cite: 2].
* **Fu, J., Knolle, J., and Perkins, N. B.** (2018). *Three types of representation of spin in terms of Majorana fermions and an alternative solution of the Kitaev honeycomb model*. Physical Review B **97**, 115142. DOI: [10.1103/PhysRevB.97.115142](https://doi.org/10.1103/PhysRevB.97.115142)[cite: 2]. arXiv: [1801.09713](https://arxiv.org/abs/1801.09713)[cite: 2].
* **Kitaev, A.** (2006). *Anyons in an exactly solved model and beyond*. Annals of Physics **321**, 2–111. DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005)[cite: 2]. arXiv: [cond-mat/0506438](https://arxiv.org/abs/cond-mat/0506438)[cite: 2].
* **Rau, J. G., Lee, E. K.-H., and Kee, H.-Y.** (2014). *Generic spin model for the honeycomb iridates beyond the Kitaev limit*. Physical Review Letters **112**, 077204. DOI: [10.1103/PhysRevLett.112.077204](https://doi.org/10.1103/PhysRevLett.112.077204)[cite: 2]. arXiv: [1401.5477](https://arxiv.org/abs/1401.5477)[cite: 2].
* **'t Hooft, G.** (1976). *Computation of the quantum effects due to a four-dimensional pseudoparticle*. Physical Review D **14**, 3432–3450. DOI: [10.1103/PhysRevD.14.3432](https://doi.org/10.1103/PhysRevD.14.3432)[cite: 2].
* **Winter, S. M., Li, Y., Jeschke, H. O., and Valentí, R.** (2016). *Challenges in design of Kitaev materials: Magnetic interactions from competing energy scales*. Physical Review B **93**, 214431. DOI: [10.1103/PhysRevB.93.214431](https://doi.org/10.1103/PhysRevB.93.214431)[cite: 2]. arXiv: [1602.05459](https://arxiv.org/abs/1602.05459)[cite: 2].
* **You, Y.-Z., Kimchi, I., and Vishwanath, A.** (2012). *Doping a spin-orbit Mott insulator: Topological superconductivity from the Kitaev-Heisenberg model*. Physical Review B **86**, 085145. DOI: [10.1103/PhysRevB.86.085145](https://doi.org/10.1103/PhysRevB.86.085145)[cite: 2]. arXiv: [1109.4155](https://arxiv.org/abs/1109.4155)[cite: 2].
