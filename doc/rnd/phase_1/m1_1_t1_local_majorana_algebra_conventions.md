---
id: RND-P1-M1.1-T1
title: "Finite-Lattice Conventions and Local Majorana Algebra"
phase: 1
milestone: M1.1
task: T1
status: completed
contracts:
  - A01
  - A02
  - A03
  - A04
  - A05
  - A06
  - A07
claims:
  - CLM-M1-ALG-003
  - CLM-M1-ALG-006
---

# Note 00a: Finite-Lattice Conventions and Local Majorana Algebra

**Status:** M1.1a R&D specification
**Parent note:** [Note 00](note_00_algebraic_qft_parton_representation.md)
**Scope:** Clarify conventions and specify independently testable local-algebra contracts. This note supplements Note 00; it does not replace or amend it.

## 1. Scope and conventions

This batch concerns only the finite-dimensional local algebra used by the M1 construction. Representation decomposition, the parity constraint, physical projection, exchange tensors, and the lattice Hamiltonian are deferred to later batches.

### 1.1 Finite lattice

Throughout this R&D, the honeycomb lattice is a finite periodic graph
\[
\Lambda=(\mathcal V,\mathcal E),\qquad |\mathcal V|=N_{\mathrm{sites}}<\infty.
\]
All tensor-product Hilbert spaces and operator algebras considered here are therefore finite-dimensional. No infinite-volume CAR algebra or thermodynamic-limit construction is assumed in this batch.

The graph is bipartite, \(\mathcal V=\mathcal V_A\sqcup\mathcal V_B\), and each directed nearest-neighbor edge is oriented from A to B. A vertex has one incident bond of each type \(\gamma\in\{1,2,3\}\equiv\{x,y,z\}\). A sum over ordered neighboring pairs counts each undirected physical edge twice.

### 1.2 Indices and cyclic ordering

Majorana-flavor indices are \(\mu,\nu,\rho,\sigma\in\{0,1,2,3\}\). Spin and bond indices are \(\alpha,\beta,\gamma\in\{1,2,3\}\).

For each bond type, the cyclic pair used in the Majorana bilinears is fixed by
\[
\begin{array}{c|cc}
\gamma & \gamma^+ & \gamma^-\\ \hline
1&2&3\\
2&3&1\\
3&1&2
\end{array}
\]
Equivalently, with labels represented by the integers \(1,2,3\),
\[
\gamma^+=1+(\gamma\bmod 3),\qquad
\gamma^-=1+((\gamma+1)\bmod 3).
\]
This table is the operative definition. It makes explicit the intended cyclic ordering without changing the contractions or generator definitions in Note 00.

The symbol \(\epsilon\) is retained for both the two-index antisymmetric tensor/matrix and the three-index Levi-Civita tensor. Their index structure and surrounding expression distinguish them; no separate symbol is introduced in this batch.

Local fermion parity remains \(\mathcal D_i\). The Dzyaloshinskii–Moriya coupling remains \(D_\gamma\).

The rank-one matrix appearing in the \(N^\gamma\) constructions is denoted
\[
e_{11}:=|1\rangle\langle1|
=\begin{pmatrix}1&0\\0&0\end{pmatrix}
=\frac12(I_2+\sigma_z).
\]
The corresponding Python constant is to be named \`E_11\`. This is a naming contract only in this batch; the existing implementation is not changed here.

### 1.3 Local Hilbert spaces and Majorana normalization

At each site, the physical spin space is \(\mathcal H_{\mathrm{phys},i}\cong\mathbb C^2\). The four Hermitian Majorana operators obey
\[
(c_i^\mu)^\dagger=c_i^\mu,\qquad
\{c_i^\mu,c_j^\nu\}=2\delta_{ij}\delta^{\mu\nu}\mathbb 1.
\]
The local Dirac operators are
\[
d_{i,1}=\frac12(c_i^0+i c_i^1),\qquad
d_{i,2}=\frac12(c_i^2+i c_i^3),
\]
with \(\{d_{i,a},d_{j,b}^\dagger\}=\delta_{ij}\delta_{ab}\). Thus the local unprojected Fock space has dimension four. The relation of this four-dimensional space to the physical spin representation is not derived in this note; it is the subject of a later batch.

## 2. Local bilinear map

For a real skew-symmetric \(4\times4\) matrix \(A\), define
\[
\mathcal O_i[A]=\frac{i}{4}c_i^T A c_i.
\]
Because \(A^T=-A\) and the Majoranas are Hermitian, \(\mathcal O_i[A]\) is Hermitian.

The local CAR imply
\[
[c^a c^b,c^c c^d]
=2\left(
\delta^{bc}c^a c^d-\delta^{ac}c^b c^d
+\delta^{bd}c^c c^a-\delta^{ad}c^c c^b
\right).
\]
Contracting with real skew matrices \(A,B\) gives
\[
[c^T A c,c^T B c]=4c^T[A,B]c,
\]
and hence
\[
[\mathcal O_i[A],\mathcal O_i[B]]
=i\,\mathcal O_i[[A,B]].
\]
The factor and sign in this identity are part of the normalization contract.

For distinct sites, even local bilinears commute:
\[
[\mathcal O_i[A],\mathcal O_j[B]]=0,\qquad i\ne j.
\]
This follows because each bilinear contains an even number of Majorana operators and the site-local Majorana sets anticommute.

## 3. Generator matrices and local operator identities

Use the real \(2\times2\) matrices
\[
I_2=\begin{pmatrix}1&0\\0&1\end{pmatrix},\quad
\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
\sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\quad
\epsilon=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\]
The six real \(4\times4\) matrices are
\[
\begin{aligned}
G^1&=I_2\otimes\epsilon,&
G^2&=\epsilon\otimes\sigma_z,&
G^3&=\epsilon\otimes\sigma_x,\\
M^1&=\sigma_z\otimes\epsilon,&
M^2&=\epsilon\otimes I_2,&
M^3&=\sigma_x\otimes\epsilon.
\end{aligned}
\]
They are real and skew-symmetric. With the cyclic ordering in Section 1.2, their contractions are
\[
c_i^T G^\gamma c_i
=2(c_i^0c_i^\gamma+c_i^{\gamma^+}c_i^{\gamma^-}),
\]
\[
c_i^T M^\gamma c_i
=2(c_i^0c_i^\gamma-c_i^{\gamma^+}c_i^{\gamma^-}).
\]

Define the normalized Hermitian bilinears
\[
S_i^\gamma=\frac12\mathcal O_i[M^\gamma]
=\frac{i}{8}c_i^TM^\gamma c_i,\qquad
T_i^\gamma=\frac12\mathcal O_i[G^\gamma]
=\frac{i}{8}c_i^TG^\gamma c_i.
\]
The matrix commutators are
\[
[G^\alpha,G^\beta]
=-2\epsilon^{\alpha\beta\gamma}G^\gamma,\qquad
[M^\alpha,M^\beta]
=+2\epsilon^{\alpha\beta\gamma}M^\gamma,\qquad
[M^\alpha,G^\beta]=0,
\]
where the repeated \(\gamma\) is summed. The corresponding operator relations are
\[
[S_i^\alpha,S_i^\beta]
=i\epsilon^{\alpha\beta\gamma}S_i^\gamma,\qquad
[T_i^\alpha,T_i^\beta]
=-i\epsilon^{\alpha\beta\gamma}T_i^\gamma,\qquad
[S_i^\alpha,T_i^\beta]=0.
\]
The opposite sign in the \(T\)-commutator is intentional and must not be normalized away by changing the definition of \(G^\gamma\) or \(T^\gamma\).

The trace contractions are
\[
\operatorname{tr}((G^\alpha)^TG^\beta)=4\delta^{\alpha\beta},\qquad
\operatorname{tr}((M^\alpha)^TM^\beta)=4\delta^{\alpha\beta},\qquad
\operatorname{tr}((M^\alpha)^TG^\beta)=0.
\]

## 4. M1.1a test contracts

The following tests are required to validate the local algebra independently of the later representation-theory and projection batches. They specify mathematical behavior, not a particular test-framework implementation.

### T00a-01 — Majorana CAR

Construct a concrete four-dimensional local Fock representation of the four Majoranas and verify
\[
(c^\mu)^\dagger=c^\mu,\qquad
c^\mu c^\nu+c^\nu c^\mu=2\delta^{\mu\nu}I_4.
\]
All sixteen ordered pairs \((\mu,\nu)\) must be checked. This test must not use the production \(G,M,N\) matrices to construct its expected CAR result.

### T00a-02 — Dirac CAR and inverse relations

From the represented Majoranas, construct
\[
d_1=(c^0+i c^1)/2,\qquad d_2=(c^2+i c^3)/2.
\]
Verify the two-mode canonical anticommutation relations, including vanishing same-type and cross-mode anticommutators. Also verify the inverse relations \(c^0=d_1+d_1^\dagger\), \(c^1=-i(d_1-d_1^\dagger)\), and analogously for \(c^2,c^3\).

### T00a-03 — Bilinear Hermiticity and commutator identity

For a spanning set of real skew \(4\times4\) matrices (the six elementary matrices \(E_{ab}-E_{ba}\), \(a<b\)), verify Hermiticity of \(\mathcal O[A]\) and the exact finite-matrix identity
\[
[\mathcal O[A],\mathcal O[B]]=i\mathcal O[[A,B]].
\]
Testing all pairs in the six-element basis is sufficient by bilinearity. Include at least one pair with nonzero matrix commutator so the test is not vacuous.

### T00a-04 — Explicit \(G,M\) contractions

For every \(\gamma\in\{1,2,3\}\), verify both bilinear contraction formulas against the represented Majorana operators. This catches basis-ordering and cyclic-index errors that matrix commutator tests alone do not detect.

### T00a-05 — Matrix Lie algebra and mutual commutation

Verify skew-symmetry of all six matrices and all nine \((\alpha,\beta)\) cases of
\[
[G^\alpha,G^\beta]=-2\epsilon^{\alpha\beta\gamma}G^\gamma,\quad
[M^\alpha,M^\beta]=+2\epsilon^{\alpha\beta\gamma}M^\gamma,\quad
[M^\alpha,G^\beta]=0.
\]
Use the same three-index orientation as the note; do not silently switch to a different Levi-Civita convention.

### T00a-06 — Trace orthogonality and normalization

Verify the complete \(3\times3\) Gram matrices:
\[
\operatorname{tr}((G^\alpha)^TG^\beta)=4I_3,\quad
\operatorname{tr}((M^\alpha)^TM^\beta)=4I_3,\quad
\operatorname{tr}((M^\alpha)^TG^\beta)=0.
\]
The full Gram matrices, rather than only their diagonal entries, must be checked.

### T00a-07 — Cross-site even-bilinear commutation

For a two-site finite Fock representation, verify that even bilinears on distinct sites commute. This is a separate CAR consequence and guards against an implementation that accidentally treats site-local Majorana matrices as commuting odd operators.

## 5. Numerical and exactness policy

These identities are finite-dimensional algebraic statements. Where the test implementation uses floating-point NumPy arrays, tolerances are a numerical verification device, not part of the mathematical claims. Prefer exact integer-valued matrix constructions where practical; otherwise use a tight, explicitly stated tolerance and test Hermiticity/antisymmetry separately.

For operator-level tests, use a concrete Majorana representation constructed independently of the production generator matrices. For matrix-level tests, compare against the explicit Kronecker definitions and the stated commutator/trace identities. Avoid tests whose expected values are computed from the same production expression being tested.

## 6. Explicit exclusions from this batch

This note does not establish or test:
- the decomposition of the four-dimensional local Fock space into gauge and spin irreducible representations;
- the Casimir-parity identities or the local parity operator;
- the equivalence between the gauge-singlet condition and the physical constraint;
- the Gutzwiller projector or projected Pauli algebra;
- the exchange tensor, quartic Hamiltonian, multiplier terms, or bond-counting prefactors;
- any mean-field decoupling, self-consistency equation, or numerical phase result.

Those subjects require separate derivations and test contracts in subsequent M1.1 batches.
