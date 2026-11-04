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

## B-002 — Rectangular four-site cell as a two-fold supercell

The rectangular Bravais vectors are
[
mathbf A_1=mathbf a_1-mathbf a_2,qquad
mathbf A_2=mathbf a_1+mathbf a_2.
]
In primitive-cell integer coordinates, a rectangular translation
((m_1,m_2)) is therefore
[
(n_1,n_2)=(m_1+m_2,-m_1+m_2).
]
The rectangular cell has twice the primitive area. Its four basis sites are
the two primitive-cell sublattices in each of the two cosets
[
mathbf t_0=(0,0),qquad mathbf t_1=(1,0),
]
where (mathbf t_1) denotes the primitive translation (mathbf a_1).
In Cartesian coordinates the four sites are
[
mathbf r_0=0,quad mathbf r_1=mathbf r_B,quad
mathbf r_2=mathbf a_1,quad
mathbf r_3=mathbf a_1+mathbf r_B.
]
Thus the production ordering ((0,1,2,3)) is precisely
((A_0,B_0,A_{mathbf t_1},B_{mathbf t_1})), with the first two
four-Majorana blocks belonging to coset 0 and the last two to coset 1.

Let (oldsymbolkappa=(kappa_1,kappa_2)) be rectangular fractional
reciprocal coordinates and (mathbf k=(k_1,k_2)) primitive fractional
coordinates. Equality of the Bloch phases under the change of direct
lattice basis gives
[
kappa_1=k_1-k_2,qquad
kappa_2=k_1+k_2.
]
For fixed (oldsymbolkappa), the two primitive momenta in the folded
fiber are
[
mathbf k^{(a)}
=left(rac{kappa_1+kappa_2}{2},
       rac{kappa_2-kappa_1}{2}ight),qquad
mathbf k^{(b)}=mathbf k^{(a)}+(1/2,1/2)
quad(mathrm{mod} 1).
]
Both map to the same rectangular momentum because the shift
((1/2,1/2)) changes (kappa_2) by one and leaves (kappa_1)
unchanged. This is the two-fold Brillouin-zone folding used in the
existing production test.

The folding also has an explicit unitary basis transformation. Group the
four-site Bloch basis by coset, so each coset block has dimension eight
(two sublattices times four Majoranas). Define
[
mathcal U(oldsymbolkappa)
=rac{1}{sqrt2}
egin{pmatrix}
I_8&I_8\
e^{i2pi k^{(a)}_1}I_8&e^{i2pi k^{(b)}_1}I_8
end{pmatrix}
=rac{1}{sqrt2}
egin{pmatrix}
I_8&I_8\
e^{i2pi k^{(a)}_1}I_8&-e^{i2pi k^{(a)}_1}I_8
end{pmatrix}.
]
The columns are orthonormal because the two coset phases differ by
(e^{ipi}=-1); hence (mathcal U^daggermathcal U=I_{16}).
With the same cell-periodic Fourier convention as the production builder,
the matrix identity is
[
H_4(oldsymbolkappa)
=mathcal U(oldsymbolkappa)
egin{pmatrix}
H_2(mathbf k^{(a)})&0\
0&H_2(mathbf k^{(b)})
end{pmatrix}
mathcal U(oldsymbolkappa)^dagger.
]
The phase in the lower coset row is the Fourier phase of the internal
primitive translation (mathbf t_1=mathbf a_1). More generally, for a
bond from coset (s_i) to coset (s_j), its primitive translation is
(mathbf n=mathbf t_{s_j}+m_1(mathbf a_1-mathbf a_2)
+m_2(mathbf a_1+mathbf a_2)-mathbf t_{s_i}). Its primitive Bloch
phase factors into the rectangular cell phase and the two coset phases;
this gives the displayed unitary block decomposition bond by bond.

This identity requires the four-site state and all couplings to be
translations of the same two-site primitive-cell data. It applies to the
uniform pure-Kitaev seed used in the existing tests (including anisotropic
bond couplings and uniform fields), but not to a genuinely enlarged-cell
mean-field state that breaks primitive translations. In that latter case,
the four-site cell describes a physically distinct modulation and folding
into two primitive-cell Hamiltonians is not valid.

The existing tests
`test_brillouin_zone_folding_parity_pure_kitaev` and
`test_brillouin_zone_folding_parity_general_couplings` compare sorted
eigenvalue sets at the two folded momenta. That is useful numerical
evidence for spectral folding, but it does not directly test the stronger
matrix identity above. A direct test should construct (mathcal U) and
check the operator-norm (or Frobenius-norm) residual
[
left|H_4-mathcal U,
mathrm{diag}(H_2(mathbf k^{(a)}),H_2(mathbf k^{(b)}))
,mathcal U^daggeright|<arepsilon
]
at generic momenta, including points on the rectangular Brillouin-zone
boundary. This direct matrix test remains the numerical acceptance check
for B-002; the basis and folding map are now derived analytically.

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

