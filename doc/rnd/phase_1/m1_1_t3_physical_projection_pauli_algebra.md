# Note 00c: Physical Projection and Recovery of the Pauli Algebra

**Status:** M1.1c R&D specification
**Parent notes:** [Note 00](note_00_algebraic_qft_parton_representation.md), [Note 00a](note_00a_local_majorana_algebra_and_conventions.md), [Note 00b](note_00b_local_spin_gauge_representation.md)
**Scope:** Construct the exact local and global physical projectors and establish the Pauli algebra on the projected Hilbert space. This note supplements the parent notes; it does not replace or amend them.

## 1. Local physical projector

At each site, define
\[
\mathcal D_i=c_i^0c_i^1c_i^2c_i^3,\qquad
P_{G,i}:=\frac{\mathbb 1+\mathcal D_i}{2}.
\]
By the local CAR, \(\mathcal D_i^\dagger=\mathcal D_i\) and \(\mathcal D_i^2=\mathbb 1\). Therefore
\[
P_{G,i}^\dagger=P_{G,i},\qquad
P_{G,i}^2=P_{G,i}.
\]
Its range is exactly the positive-parity sector:
\[
\operatorname{Ran}P_{G,i}
=\ker(\mathcal D_i-\mathbb 1)
=\bigcap_{\gamma=1}^3\ker T_i^\gamma.
\]
The final equality is the gauge-singlet characterization established in Note 00b. Since the local Fock space has dimension four and each parity sector has dimension two,
\[
\operatorname{rank}P_{G,i}=2,\qquad
\operatorname{tr}P_{G,i}=2.
\]
Thus \(P_{G,i}\) is the orthogonal projector from the four-dimensional local Fock space onto the two-dimensional physical spin space.

## 2. Global physical projector on a finite lattice

For the finite site set \(\mathcal V\), define
\[
P_G:=\prod_{i\in\mathcal V}P_{G,i}
=\prod_{i\in\mathcal V}\frac{\mathbb 1+\mathcal D_i}{2}.
\]
Each \(\mathcal D_i\) is even in the Majoranas at site \(i\). Consequently, for \(i\ne j\),
\[
[\mathcal D_i,\mathcal D_j]=0,\qquad
[P_{G,i},P_{G,j}]=0.
\]
The product is therefore independent of the ordering of sites and is itself an orthogonal projector:
\[
P_G^\dagger=P_G,\qquad P_G^2=P_G.
\]
Its range is the tensor product of the local physical sectors,
\[
\operatorname{Ran}P_G
=\bigotimes_{i\in\mathcal V}\operatorname{Ran}P_{G,i}
\cong(\mathbb C^2)^{\otimes N_{\mathrm{sites}}}.
\]
In particular,
\[
\operatorname{rank}P_G
=\operatorname{tr}P_G
=2^{N_{\mathrm{sites}}}.
\]
The unprojected Fock space has dimension \(4^{N_{\mathrm{sites}}}\), so the physical projector retains exactly a fraction \(2^{-N_{\mathrm{sites}}}\) of its dimension.

This is a product of local constraints, not a single global-parity constraint. Projecting only onto a prescribed eigenvalue of \(\prod_i\mathcal D_i\) would retain states with compensating local parity violations and would not define the physical Hilbert space.

## 3. Majorana representation of the spin operators

For \(\gamma\in\{1,2,3\}\), define
\[
\sigma_i^\gamma:=i c_i^0c_i^\gamma.
\]
With
\[
N^\gamma:=\frac12(M^\gamma+G^\gamma),
\]
the explicit contractions in Note 00a give
\[
c_i^T N^\gamma c_i=2c_i^0c_i^\gamma,
\qquad
\sigma_i^\gamma=\frac{i}{2}c_i^T N^\gamma c_i.
\]
These are Hermitian operators. Direct use of the Majorana CAR yields, on the full unprojected local Fock space,
\[
(\sigma_i^\gamma)^2=\mathbb 1,
\qquad
\{\sigma_i^\alpha,\sigma_i^\beta\}
=2\delta^{\alpha\beta}\mathbb 1.
\]
Thus the three operators already satisfy the Pauli Clifford relations on the enlarged space. What fails on that space is the correct oriented multiplication law: for distinct components, the product contains a gauge-generator contribution.

## 4. Exact unprojected multiplication identity

Using the definitions
\[
S_i^\gamma=\frac{i}{8}c_i^T M^\gamma c_i,\qquad
T_i^\gamma=\frac{i}{8}c_i^T G^\gamma c_i,
\]
the relation between the bilinears is
\[
S_i^\gamma+T_i^\gamma
=\frac{i}{8}c_i^T(M^\gamma+G^\gamma)c_i
=\frac12\sigma_i^\gamma.
\]
For cyclically ordered, pairwise distinct \((\alpha,\beta,\gamma)\), the local Majorana CAR give
\[
\sigma_i^\alpha\sigma_i^\beta
=i\sigma_i^\gamma-4iT_i^\gamma.
\]
Including the reversed order and equal-component cases, this is equivalently
\[
\boxed{
\sigma_i^\alpha\sigma_i^\beta
=\delta^{\alpha\beta}\mathbb 1
+i\sum_{\gamma=1}^3\epsilon^{\alpha\beta\gamma}
\bigl(\sigma_i^\gamma-4T_i^\gamma\bigr).
}
\]
The gauge term is essential: omitting it would assert the physical Pauli multiplication law on the entire four-dimensional Fock space, which is false. The identity above is an exact operator identity before projection.

## 5. Projection recovers the physical Pauli algebra

Every \(\sigma_i^\gamma\) is an even local Majorana bilinear, so
\[
[\mathcal D_i,\sigma_i^\gamma]=0.
\]
It also commutes with \(\mathcal D_j\) for \(j\ne i\), because even operators on distinct sites commute. Hence
\[
[P_G,\sigma_i^\gamma]=0,\qquad
P_G\sigma_i^\gamma P_G=\sigma_i^\gamma P_G.
\]
From Note 00b, the physical sector is annihilated by every gauge generator:
\[
T_i^\gamma P_G=P_GT_i^\gamma=0.
\]
Sandwiching the exact unprojected multiplication identity by \(P_G\) therefore gives
\[
\boxed{
P_G\sigma_i^\alpha\sigma_i^\beta P_G
=\delta^{\alpha\beta}P_G
+i\sum_{\gamma=1}^3\epsilon^{\alpha\beta\gamma}
P_G\sigma_i^\gamma P_G.
}
\]
Equivalently, define the physical operator on \(\operatorname{Ran}P_G\) by
\[
\sigma_{i,\mathrm{phys}}^\gamma
:=\left.P_G\sigma_i^\gamma P_G\right|_{\operatorname{Ran}P_G}.
\]
Then
\[
\sigma_{i,\mathrm{phys}}^\alpha
\sigma_{i,\mathrm{phys}}^\beta
=\delta^{\alpha\beta}\mathbb 1_{\mathrm{phys}}
+i\sum_\gamma\epsilon^{\alpha\beta\gamma}
\sigma_{i,\mathrm{phys}}^\gamma.
\]
The physical spin operators are \(S_{i,\mathrm{phys}}^\gamma=\sigma_{i,\mathrm{phys}}^\gamma/2\), with the conventional spin-\(1/2\) commutator and Casimir.

For distinct sites, the unprojected \(\sigma_i^\alpha\) and \(\sigma_j^\beta\) commute because they are even local operators. Their projected restrictions therefore commute as well:
\[
[\sigma_{i,\mathrm{phys}}^\alpha,
 \sigma_{j,\mathrm{phys}}^\beta]=0,\qquad i\ne j.
\]
Thus the product physical Hilbert space carries one independent spin-\(1/2\) Pauli algebra at each site.

## 6. Test contracts

The tests use an independently constructed Jordan–Wigner representation of the Majoranas. The production \(G,M,N\) matrices are used only where the contract explicitly compares their bilinear definitions with direct Majorana products.

### T00c-01 — Local orthogonal projector

Verify for one site that \(P_{G,i}^\dagger=P_{G,i}\), \(P_{G,i}^2=P_{G,i}\), and \(\operatorname{tr}P_{G,i}=2\). Also verify \(P_{G,i}=(\mathbb 1+\mathcal D_i)/2\).

### T00c-02 — Global product projector and rank

On a two-site finite Fock space, construct the two local parity projectors independently and verify that they commute. Verify that their product is Hermitian and idempotent, and has trace/rank four. Check that it equals the projector onto the simultaneous positive-parity eigenspace, not merely the positive eigenspace of total parity.

### T00c-03 — Direct Kitaev bilinear and unprojected Clifford relations

For every \(\gamma\), compare \(i c^0c^\gamma\) with \((i/2)c^TN^\gamma c\). Verify Hermiticity, squares equal to identity, and all anticommutators on the full local Fock space.

### T00c-04 — Full-space product identity

For all nine pairs \((\alpha,\beta)\), verify
\[
\sigma^\alpha\sigma^\beta
=\delta^{\alpha\beta}\mathbb 1
+i\epsilon^{\alpha\beta\gamma}
(\sigma^\gamma-4T^\gamma).
\]
The test must include both cyclic and reversed component orderings and equal-component cases.

### T00c-05 — Projected local Pauli multiplication

For all nine component pairs, verify the sandwiched Pauli multiplication identity with the local physical projector. Also verify that the projector commutes with every \(\sigma^\gamma\).

### T00c-06 — Two-site physical tensor-product algebra

On a two-site Fock space, verify the global projector rank and show that each projected site operator obeys the Pauli multiplication law while operators on distinct sites commute. This checks that local projection composes into the intended tensor-product physical spin space.

## 7. Exactness and scope

The projector identities, rank statements, and Pauli multiplication law are exact finite-dimensional operator results. Floating-point matrix residuals in the tests are numerical checks, not proofs.

This batch does not derive exchange tensors, the quartic Majorana Hamiltonian, multiplier prefactors, bond counting, or mean-field decouplings. Those are deferred to M1.1d and later batches.
