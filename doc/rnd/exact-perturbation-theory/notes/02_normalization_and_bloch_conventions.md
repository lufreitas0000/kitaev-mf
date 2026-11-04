# B-001/B-002 — Operator normalization and Bloch-convention audit

## Scope and status

This note audits the normalization and primitive-cell conventions shared by the exact
Kitaev research track and the production Majorana mean-field implementation on `main`.
It is an operator-level comparison, not a code port. The research convention is

\[
H_K^{\mathrm{R}}=-K_{\mathrm R}\sum_{\langle ij\rangle_\gamma}
\sigma_i^\gamma\sigma_j^\gamma,\qquad K_{\mathrm R}>0,
\]

with each physical bond counted once. The production exchange convention has
\(K_{\mathrm P}<0\) and the pure-Kitaev Hamiltonian
\(H_K^{\mathrm P}=K_{\mathrm P}\sum_{\langle ij\rangle_\gamma}
\sigma_i^\gamma\sigma_j^\gamma\). Therefore the spin-Hamiltonian coupling map is

\[
K_{\mathrm P}=-K_{\mathrm R}.
\]

This relation is exact and does not involve a cell-volume or Fourier normalization.

## B-001 — Operator-level normalization

Both implementations use four Hermitian Majoranas per site with
\(\{c_i^\mu,c_j^\nu\}=2\delta_{ij}\delta^{\mu\nu}\), and the same Kitaev
representation
\[
\sigma_i^\gamma=i c_i^0c_i^\gamma
=\frac{i}{2}c_i^T N^\gamma c_i,
\qquad
(N^\gamma)_{0\gamma}=+1,\quad
(N^\gamma)_{\gamma0}=-1.
\]

For an A-to-B bond, define the oriented gauge bilinear
\[
u_{ij}=i c_i^\gamma c_j^\gamma.
\]
Using anticommutation between Majoranas on distinct sites gives
\[
\begin{aligned}
K_{\mathrm P}\sigma_i^\gamma\sigma_j^\gamma
&=-K_{\mathrm P}c_i^0c_i^\gamma c_j^0c_j^\gamma\\
&=K_{\mathrm P}c_i^0c_j^0c_i^\gamma c_j^\gamma\\
&=-iK_{\mathrm P}u_{ij}c_i^0c_j^0
=iK_{\mathrm R}u_{ij}c_i^0c_j^0.
\end{aligned}
\]
Thus, in the fixed-gauge sector (u_{ij}=+1), the exact matter hopping on
each directed bond is (iK_{\mathrm R}c_i^0c_j^0).

For a real antisymmetric Majorana matrix (A), the convention
\[
H=\frac{i}{4}c^T A c
\]
assigns (A_{ij}=2K_{\mathrm R}u_{ij}) on this bond: the two ordered
matrix entries combine to (iA_{ij}c_i c_j/2). Consequently the bipartite
off-diagonal block is (T_{ij}=2K_{\mathrm R}u_{ij}), exactly the
coefficient used in the research note `01_flux_sectors_and_spectrum.md`.

The production mean-field hopping builder uses
\[
A_b=2\sum_{\alpha,\beta}J_b^{\alpha\beta}
N^\alpha U_b N^\beta.
\]
In the pure Kitaev limit, (J_b^{\gamma\gamma}=K_{\mathrm P}),
(U_b^{\gamma\gamma}=u_{ij}), and
\((N^\gamma U_bN^\gamma)_{00}=-u_{ij}). Hence
\[
(A_b)_{00}=-2K_{\mathrm P}u_{ij}
=2K_{\mathrm R}u_{ij}.
\]
At the level of the real-space antisymmetric Majorana matrix, the production
hopping coefficient and the research (T)-matrix coefficient therefore agree.
There is no residual factor of four in the microscopic coupling or in the
bond-level (A_b/T) coefficient.

### Spin-operator convention

If another source writes (S_i^\gamma=\sigma_i^\gamma/2) and
\(H=-J\sum_{\langle ij\rangle_\gamma}S_i^\gamma S_j^\gamma\), then
\[
K_{\mathrm R}=J/4,\qquad K_{\mathrm P}=-J/4.
\]
This is a separate Pauli-versus-spin normalization. It must not be conflated
with the Bloch-matrix normalization below.

## B-002 — Primitive cell, directed bonds, and Fourier phases

The production primitive cell has
\[
\mathbf a_1=(1/2,\sqrt3/2),\qquad
\mathbf a_2=(-1/2,\sqrt3/2),
\]
with basis positions
\[
\mathbf r_A=(0,0),\qquad
\mathbf r_B=(0,1/\sqrt3).
\]
The three directed A-to-B bonds are:

| Bond | Type | Sublattice pair | Cell offset |
|---|---:|---|---|
| (z) | 3 | (A\to B) | ((0,0)) |
| (x) | 1 | (A\to B) | ((0,-1)) |
| (y) | 2 | (A\to B) | ((-1,0)) |

For fractional reciprocal coordinates ((k_1,k_2)), the production phase
is (e^{i2\pi(k_1n_1+k_2n_2)}). The matter-sector structure factor is
therefore
\[
f_{\mathrm P}(\mathbf k)
=1+e^{-i2\pi k_2}+e^{-i2\pi k_1}.
\]
The research convention writes
\[
f_{\mathrm R}(\boldsymbol\theta)
=1+e^{i\theta_1}+e^{i\theta_2},
\qquad \theta_a=\mathbf q\cdot\mathbf n_a.
\]
Identifying ((\theta_1,\theta_2)=(-2\pi k_1,-2\pi k_2))
makes the two structure factors identical. Reordering the (x)- and
(y)-bond terms is immaterial. In particular,
\[
|f_{\mathrm P}(\mathbf k)|=|f_{\mathrm R}(\boldsymbol\theta)|.
\]
The primitive-cell geometry, directed bond count, and Fourier phases do
not explain the reported factor-four spectral discrepancy.

## B-001/B-002 — Bloch matrix scaling: discrepancy identified

The remaining factor four comes from the matrix used by the production Bloch
builder, not from the coupling or lattice conventions. The real-space
quadratic operator is
\[
H=\frac{i}{4}c^T A c,
\qquad A^T=-A.
\]
Its canonical Majorana equations of motion are
\[
\frac{d c}{dt}=A c,
\]
so the Hermitian generator governing the one-particle excitation frequencies
is (iA), not (iA/4). Equivalently, canonical mode diagonalization gives
\[
H=\sum_{n=1}^{N}\epsilon_n
\left(d_n^\dagger d_n-\frac12\right),
\qquad \epsilon_n=s_n(T),
\]
for the bipartite block matrix
\[
A=\begin{pmatrix}0&T\\-T^T&0\end{pmatrix}.
\]
The positive frequencies are the singular values of (T).

The two-Majorana calculation makes the normalization explicit. For
\(A_{12}=t\), (A_{21}=-t\), and
(d=(c_1+ic_2)/2),
\[
H=\frac{i}{4}c^TAc
=\frac{it}{2}c_1c_2
=t\left(d^\dagger d-\frac12\right).
\]
The excitation energy is (t). The eigenvalues of (iA) are
(pm t), whereas those of (iA/4) are (pm t/4). The latter are not
the canonical excitation energies.

The production builder constructs its matrix using
\[
H_{\mathrm{Bloch}}^{\mathrm P}(\mathbf k)
=\frac{i}{4}A(\mathbf k).
\]
Its test therefore reports
\[
E_{\mathrm{test}}^{\mathrm P}(\mathbf k)
=\pm\frac{|K_{\mathrm P}|}{2}|f_{\mathrm P}(\mathbf k)|.
\]
But the real-space antisymmetric matrix has
(T(\mathbf k)=2K_{\mathrm R}f(\mathbf k)), so its canonical positive
excitation frequency is
\[
\epsilon_{\mathrm{phys}}(\mathbf k)
=s(T(\mathbf k))
=2K_{\mathrm R}|f(\mathbf k)|.
\]
Using (K_{\mathrm P}=-K_{\mathrm R}), the production test eigenvalue is
exactly one quarter of the physical excitation energy. This reproduces the
factor-four discrepancy without changing the coupling map or the lattice
convention.

This identifies a mismatch between the production Bloch matrix's reported
eigenvalues and the exact quasiparticle excitation energies. It does not,
by itself, prove that the production mean-field operator or self-consistent
solver is wrong. The current solver passes the matrix to a negative-energy
spectral-projector routine; multiplying a Hermitian matrix by a positive
constant leaves its occupied eigenspace unchanged, apart from the numerical
zero-mode threshold. Thus the scale cancels from the ideal self-consistency
map. The builder's matrix can serve as a quadratic-form coefficient for
that purpose, but its eigenvalues must not be labelled physical excitation
energies without the corresponding canonical normalization. Its public
contract and downstream use must be audited before changing it.
For an API that returns physical excitation energies, the matrix to
diagonalize must reproduce the spectrum of (iA) (or, equivalently, the
singular values of (T)); the current (iA/4) eigenvalues do not.

**Decision:** B-001's coupling and bond-matrix normalization, and B-002's
primitive-cell geometry and Fourier phases, are reconciled. The source of
the factor-four discrepancy is identified as the production Bloch matrix
scaling. Do not modify production code from the research branch. First
clarify the production builder's intended mathematical contract and verify
the relation between its returned matrix, the many-body quadratic operator,
and its canonical excitation frequencies.

## Required next checks

1. Derive the production Fourier-space quadratic operator from its real-space
   (i/4,c^TAc) expression, retaining the full (mathbf k,-\mathbf k)
   Majorana redundancy.
2. Derive the canonical complex-fermion mode energies from that operator and
   compare them with both the eigenvalues of (iA(\mathbf k)/4) and the
   singular values of (T(\mathbf k)).
3. Add a minimal two-site/two-Majorana normalization test before changing
   `test_pure_kitaev_analytical_dispersion` or the production builder.
4. Only after this is settled, proceed to finite-torus momentum grids and
   the 4-site folding map.

