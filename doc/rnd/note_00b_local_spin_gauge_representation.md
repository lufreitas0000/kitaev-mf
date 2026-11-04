# Note 00b: Local Parity Sectors and the Spin–Gauge Representation

**Status:** M1.1b R&D specification
**Parent notes:** [Note 00](note_00_algebraic_qft_parton_representation.md), [Note 00a](note_00a_local_majorana_algebra_and_conventions.md)
**Scope:** Derive the decomposition of the four-dimensional local Fock representation under the commuting spin and gauge algebras. This note supplements the parent notes; it does not replace or amend them.

## 1. Algebraic setting and sign convention

At one site, the four Hermitian Majoranas satisfy
\[
\{c^\mu,c^\nu\}=2\delta^{\mu\nu}\mathbb 1,\qquad \mu,\nu\in\{0,1,2,3\}.
\]
The local Fock space \(\mathcal F\) has dimension four. The six Hermitian bilinears are
\[
S^\alpha=\frac{i}{8}c^T M^\alpha c,\qquad
T^\alpha=\frac{i}{8}c^T G^\alpha c.
\]
Their commutators, fixed in Note 00a, are
\[
[S^\alpha,S^\beta]=i\epsilon^{\alpha\beta\gamma}S^\gamma,\qquad
[T^\alpha,T^\beta]=-i\epsilon^{\alpha\beta\gamma}T^\gamma,\qquad
[S^\alpha,T^\beta]=0.
\]
The opposite sign in the gauge commutator matters when assigning conventional representation labels. Define
\[
J_{\mathrm g}^\alpha:=-T^\alpha.
\]
Then
\[
[J_{\mathrm g}^\alpha,J_{\mathrm g}^\beta]
=i\epsilon^{\alpha\beta\gamma}J_{\mathrm g}^\gamma,\qquad
[S^\alpha,J_{\mathrm g}^\beta]=0.
\]
Thus \(S^\alpha\) and \(J_{\mathrm g}^\alpha\) are two commuting sets of Hermitian generators with the standard \(\mathfrak{su}(2)\) convention. The sign change does not alter the gauge quadratic Casimir:
\[
\sum_\alpha (J_{\mathrm g}^\alpha)^2=\sum_\alpha(T^\alpha)^2.
\]

Throughout this note, representation labels are ordered as
\[
(j_{\mathrm{spin}},j_{\mathrm{gauge}}).
\]

## 2. Local parity and its invariant sectors

Define the local Majorana parity operator
\[
\mathcal D:=c^0c^1c^2c^3.
\]
The CAR imply
\[
\mathcal D^\dagger=\mathcal D,\qquad \mathcal D^2=\mathbb 1.
\]
Every \(S^\alpha\) and \(T^\alpha\) is an even Majorana polynomial. Consequently,
\[
[\mathcal D,S^\alpha]=[\mathcal D,T^\alpha]=0.
\]
The local Fock space therefore splits into invariant parity eigenspaces
\[
\mathcal F=\mathcal F_+\oplus\mathcal F_-,\qquad
\mathcal F_\pm:=\ker(\mathcal D\mp\mathbb 1).
\]

The two Dirac modes of Note 00a are
\[
d_1=\frac12(c^0+ic^1),\qquad d_2=\frac12(c^2+ic^3),
\]
with number operators \(n_a=d_a^\dagger d_a\). Directly from these definitions and the CAR,
\[
c^0c^1=i(1-2n_1),\qquad
c^2c^3=i(1-2n_2),
\]
and hence
\[
\mathcal D=-(1-2n_1)(1-2n_2)
=-(-1)^{n_1+n_2}.
\]
In the occupation basis \(\{|n_1n_2\rangle:n_1,n_2\in\{0,1\}\}\),
\[
\mathcal F_+=\operatorname{span}\{|10\rangle,|01\rangle\},\qquad
\mathcal F_-=\operatorname{span}\{|00\rangle,|11\rangle\}.
\]
In particular, each parity eigenspace has dimension two. The minus sign in \(\mathcal D=-(-1)^{n_1+n_2}\) depends on the stated Majorana/Dirac convention and should not be silently replaced by the usual Dirac fermion-parity convention.

## 3. Casimir identities and the representation decomposition

The exact local Casimir–parity identities are
\[
\mathbf T^2:=\sum_{\alpha=1}^3(T^\alpha)^2
=\frac38(\mathbb 1-\mathcal D),\qquad
\mathbf S^2:=\sum_{\alpha=1}^3(S^\alpha)^2
=\frac38(\mathbb 1+\mathcal D).
\]
Restricting these identities to the two parity sectors gives
\[
\begin{array}{c|cc}
&\mathbf S^2&\mathbf J_{\mathrm g}^2\\ \hline
\mathcal F_+ \;(\mathcal D=+1)&\frac34\mathbb 1&0\\
\mathcal F_- \;(\mathcal D=-1)&0&\frac34\mathbb 1
\end{array}
\]
where \(\mathbf J_{\mathrm g}^2=\mathbf T^2\).

For Hermitian generators in a finite-dimensional unitary representation of \(\mathfrak{su}(2)\), the quadratic Casimir on an irreducible spin-\(j\) representation is \(j(j+1)\mathbb 1\). Since \(j(j+1)=3/4\) gives \(j=1/2\), and each parity sector is two-dimensional, the restrictions are irreducible doublets:
- on \(\mathcal F_+\), the spin algebra acts as spin \(1/2\), while the gauge algebra acts trivially;
- on \(\mathcal F_-\), the gauge algebra acts as spin \(1/2\), while the spin algebra acts trivially.

The local representation therefore decomposes as
\[
\boxed{\mathcal F\cong
\left(\frac12,0\right)_{\mathrm{spin,gauge}}
\oplus
\left(0,\frac12\right)_{\mathrm{spin,gauge}}.}
\]
Equivalently, with the order \((j_{\mathrm{gauge}},j_{\mathrm{spin}})\), this is
\[
\mathcal F\cong
\left(0,\frac12\right)_{\mathrm{gauge,spin}}
\oplus
\left(\frac12,0\right)_{\mathrm{gauge,spin}}.
\]
This is a reducible four-dimensional representation of
\(\mathfrak{su}(2)_{\mathrm{spin}}\oplus\mathfrak{su}(2)_{\mathrm{gauge}}\), not a four-dimensional irreducible representation. The physical spin-\(1/2\) space is precisely \(\mathcal F_+\), the gauge-singlet summand.

## 4. Gauge-singlet characterization

On \(\mathcal F_+\), \(\mathbf T^2=0\). Since each \(T^\alpha\) is Hermitian,
\[
0=\langle\psi|\mathbf T^2|\psi\rangle
=\sum_{\alpha=1}^3\|T^\alpha|\psi\rangle\|^2
\quad\Longrightarrow\quad
T^\alpha|\psi\rangle=0
\quad\text{for every }\alpha.
\]
Conversely, if all three \(T^\alpha\) annihilate \(|\psi\rangle\), then \(\mathbf T^2|\psi\rangle=0\), and the Casimir–parity identity gives \(\mathcal D|\psi\rangle=|\psi\rangle\). Therefore
\[
\mathcal F_+
=\bigcap_{\alpha=1}^3\ker T^\alpha
=\ker(\mathcal D-\mathbb 1).
\]
This is a local statement. The global physical subspace is the tensor product of these local sectors, with the local constraints imposed at every site.

## 5. Test contracts

These tests concern the local four-dimensional representation and are independent of the production \(G,M\) matrices when constructing the Majorana operators. They may use the production matrices only to construct and check the stated bilinears.

### T00b-01 — Parity operator

In an independently constructed local Majorana representation, verify
\[
\mathcal D^\dagger=\mathcal D,\qquad
\mathcal D^2=\mathbb 1,\qquad
\operatorname{tr}\mathcal D=0.
\]
The trace condition confirms that the two parity eigenspaces have equal dimension.

### T00b-02 — Invariance of parity sectors

Construct \(S^\alpha,T^\alpha\) from the bilinear definitions and verify
\[
[\mathcal D,S^\alpha]=[\mathcal D,T^\alpha]=0
\]
for all three components.

### T00b-03 — Casimir–parity identities

Verify both operator identities
\[
\sum_\alpha(T^\alpha)^2=\frac38(\mathbb 1-\mathcal D),\qquad
\sum_\alpha(S^\alpha)^2=\frac38(\mathbb 1+\mathcal D).
\]
The two sides must be constructed independently enough that the test does not merely restate the same formula.

### T00b-04 — Parity-sector representation content

Form \(P_\pm=(\mathbb 1\pm\mathcal D)/2\). Verify projector idempotency, orthogonality, ranks two, and
\[
\mathbf S^2P_+=\frac34P_+,\quad \mathbf T^2P_+=0,\qquad
\mathbf S^2P_-=0,\quad \mathbf T^2P_-=\frac34P_-.
\]
Also verify \(T^\alpha P_+=0\) and \(S^\alpha P_-=0\) componentwise.

### T00b-05 — Conventional gauge-generator sign and occupation parity

Set \(J_{\mathrm g}^\alpha=-T^\alpha\). Verify its standard-sign commutator, its commutation with \(S^\beta\), and
\[
\mathcal D=-(-1)^{n_1+n_2}
\]
using the Dirac operators and number operators constructed from the represented Majoranas. This test protects both the gauge sign convention and the parity convention.

## 6. Exactness and scope

All claims in this note are exact finite-dimensional operator statements. Floating-point residuals in the tests are numerical checks of these identities, not proofs.

This batch does not derive the Gutzwiller projector on the full lattice, the projected Pauli algebra, exchange tensors, the quartic Hamiltonian, or any mean-field approximation. Those remain separate M1.1 batches.
