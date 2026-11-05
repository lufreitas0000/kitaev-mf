---
id: RND-P1-M1.1-T4
title: "Exchange Tensors and the Exact Quartic Majorana Hamiltonian"
phase: 1
milestone: M1.1
task: T4
status: completed
contracts:
  - A16
  - A17
  - A18
  - A19
  - A20
claims:
  - CLM-M1-ALG-010
  - CLM-M1-ALG-011
---

# Note 00d: Exchange Tensors and the Exact Quartic Majorana Hamiltonian

**Status:** M1.1d R&D specification
**Parent notes:** [Note 00](note_00_algebraic_qft_parton_representation.md), [Note 00a](note_00a_local_majorana_algebra_and_conventions.md), [Note 00b](note_00b_local_spin_gauge_representation.md), [Note 00c](note_00c_physical_projection_and_pauli_algebra.md)
**Scope:** Fix the directed-bond exchange convention, derive the quartic Majorana Hamiltonian and constraint-multiplier normalization, and establish projected operator equivalence. This note supplements the parent notes; it does not replace or amend them.

## 1. Finite directed honeycomb graph and edge counting

Let \(\Lambda=(\mathcal V,\mathcal E)\) be the finite periodic bipartite honeycomb graph fixed in Note 00a. Each undirected nearest-neighbor edge has a unique orientation from its A-sublattice endpoint to its B-sublattice endpoint. Write this set of oriented representatives as \(\overrightarrow{\mathcal E}_{A\to B}\). Each physical edge appears once in this set.

For an oriented edge \(e=(i,j)\) of type \(\gamma\), let \(J_e\) be a real \(3\times3\) exchange tensor. The reverse orientation is not a second independent coupling:
\[
J_{ji}=J_{ij}^{T}.
\]
The ordered-neighbor sum therefore obeys
\[
\frac12\sum_{i\in\mathcal V}\sum_{j\in\mathcal V_i}
\sum_{a,b}\sigma_i^a J_{ij}^{ab}\sigma_j^b
=
\sum_{(i,j)\in\overrightarrow{\mathcal E}_{A\to B}}
\sum_{a,b}\sigma_i^a J_{ij}^{ab}\sigma_j^b.
\]
Indeed, the reverse contribution is
\[
\sum_{a,b}\sigma_j^a J_{ji}^{ab}\sigma_i^b
=\sum_{a,b}\sigma_j^a J_{ij}^{ba}\sigma_i^b
=\sum_{a,b}\sigma_i^a J_{ij}^{ab}\sigma_j^b,
\]
where operators on distinct sites commute. Thus the factor \(1/2\) removes exactly the double counting of each undirected edge. A sum over \(\overrightarrow{\mathcal E}_{A\to B}\) must not carry this extra factor.

The directed-edge convention is important when the exchange tensor has an antisymmetric component, such as a Dzyaloshinskii–Moriya term. Reversing the edge transposes the tensor; it does not leave its components unchanged.

## 2. Exchange tensor and its symmetry content

For a bond of type \(\gamma\in\{1,2,3\}\), use the cyclic labels \((\gamma,\gamma^+,\gamma^-)\) from Note 00a. The exchange tensor is
\[
\begin{aligned}
J_\gamma^{ab}={}&
J_\gamma\delta^{ab}
+K_\gamma\delta^{a\gamma}\delta^{b\gamma}\\
&+\Gamma_\gamma\left(
\delta^{a\gamma^+}\delta^{b\gamma^-}
+\delta^{a\gamma^-}\delta^{b\gamma^+}\right)\\
&+\Gamma'_\gamma\left(
\delta^{a\gamma^+}\delta^{b\gamma}
+\delta^{a\gamma}\delta^{b\gamma^+}
+\delta^{a\gamma^-}\delta^{b\gamma}
+\delta^{a\gamma}\delta^{b\gamma^-}\right)\\
&+D_\gamma\left(
\delta^{a\gamma^+}\delta^{b\gamma^-}
-\delta^{a\gamma^-}\delta^{b\gamma^+}\right).
\end{aligned}
\]
Its symmetric and antisymmetric parts are
\[
J_\gamma^{(s)}=\frac12(J_\gamma+J_\gamma^T),\qquad
J_\gamma^{(a)}=\frac12(J_\gamma-J_\gamma^T).
\]
The \(J,K,\Gamma,\Gamma'\) terms are symmetric under transpose, while the \(D_\gamma\) term is antisymmetric. In the cyclic basis \((\gamma^+,\gamma^-,\gamma)\), the tensor is
\[
J_{\langle ij\rangle_\gamma}=
\begin{pmatrix}
J_\gamma & \Gamma_\gamma+D_\gamma & \Gamma'_\gamma\\
\Gamma_\gamma-D_\gamma & J_\gamma & \Gamma'_\gamma\\
\Gamma'_\gamma & \Gamma'_\gamma & J_\gamma+K_\gamma
\end{pmatrix}.
\]
For \(\gamma=3\), the ordered basis is \((1,2,3)\), reproducing the \(z\)-bond matrix in Note 00.

No symmetry of \(J_{ij}\) under transpose is assumed when \(D_\gamma\ne0\). Hermiticity of the spin Hamiltonian instead follows from real couplings and the relation \(J_{ji}=J_{ij}^T\) on the reversed edge.

## 3. Exact spin Hamiltonian and Majorana substitution

The finite-lattice spin Hamiltonian is
\[
H=
\frac12\sum_{i\in\mathcal V}\sum_{j\in\mathcal V_i}
\sum_{a,b=1}^3
\sigma_i^aJ_{ij}^{ab}\sigma_j^b
-\sum_{i\in\mathcal V}\sum_{a=1}^3h_i^a\sigma_i^a.
\]
The field is written site-dependently here; a uniform field is the special case \(h_i^a=h^a\).

Use the exact local representation
\[
\sigma_i^a=i c_i^0c_i^a
=\frac{i}{2}c_i^TN^ac_i.
\]
It follows that
\[
i c_i^TN^ac_i=2\sigma_i^a.
\]
Substitution into the exchange term gives
\[
\begin{aligned}
H_{\mathrm{ex}}
&=\frac12\sum_{i,j\ {\rm ordered}}\sum_{a,b}
\left(\frac{i}{2}c_i^TN^ac_i\right)
J_{ij}^{ab}
\left(\frac{i}{2}c_j^TN^bc_j\right)\\
&=\frac18\sum_{i,j\ {\rm ordered}}\sum_{a,b}
(i c_i^TN^ac_i)J_{ij}^{ab}(i c_j^TN^bc_j).
\end{aligned}
\]
The factor \(1/8\) is the product of the original ordered-edge factor \(1/2\) and the two factors \(1/2\) from the spin representation. For the field term,
\[
H_h=-\frac12\sum_{i,a}h_i^a\,i c_i^TN^ac_i.
\]
These are exact substitutions, not a mean-field decoupling. The quartic character refers to the product of two local Majorana bilinears on distinct sites.

## 4. Constraint multipliers and normalization

The physical subspace is characterized by
\[
T_i^\gamma=0,\qquad
T_i^\gamma=\frac{i}{8}c_i^TG^\gamma c_i.
\]
Introduce real multipliers \(\lambda_i^\gamma\) through the Hermitian constraint term
\[
H_\lambda^{\mathrm{constr}}
:=2\sum_{i,\gamma}\lambda_i^\gamma T_i^\gamma
=\frac14\sum_{i,\gamma}\lambda_i^\gamma\,i c_i^TG^\gamma c_i.
\]
The factor \(2\) is a convention for the multiplier: it can be absorbed into a redefinition of \(\lambda_i^\gamma\). With this convention the full extended operator is
\[
\boxed{
\begin{aligned}
H_\lambda={}&
\frac18\sum_{i}\sum_{j\in\mathcal V_i}\sum_{a,b}
(i c_i^TN^ac_i)J_{ij}^{ab}(i c_j^TN^bc_j)\\
&-\frac14\sum_{i,a}2h_i^a\,i c_i^TN^ac_i
+\frac14\sum_{i,\gamma}\lambda_i^\gamma\,i c_i^TG^\gamma c_i.
\end{aligned}}
\]
Equivalently, the last two terms can be combined as
\[
-\frac14\sum_{i,\gamma}
\left(2h_i^\gamma i c_i^TN^\gamma c_i
-\lambda_i^\gamma i c_i^TG^\gamma c_i\right)
\]
when the field and constraint component indices are both summed over \(1,2,3\).

The multipliers are not needed to define the exact projected spin Hamiltonian: the local projector already imposes the constraints. They are included to prepare a later constrained mean-field formulation. At this stage they are arbitrary real parameters, not self-consistent fields and not an approximation.

## 5. Exact projected operator equivalence

Let \(P_G=\prod_i(1+\mathcal D_i)/2\). Since
\[
T_i^\gamma P_G=P_GT_i^\gamma=0,
\]
the multiplier term vanishes after projection:
\[
P_G H_\lambda^{\mathrm{constr}}P_G=0.
\]
The exchange and field terms are exact Majorana substitutions of the spin operators, and all these local spin bilinears commute with the local parity projectors. Consequently,
\[
\boxed{P_GH_\lambda P_G=P_GHP_G.}
\]
This is an operator identity on the enlarged Fock space after sandwiching with the physical projector. Equivalently, the restrictions of \(H_\lambda\) and \(H\) to \(\operatorname{Ran}P_G\) coincide.

This statement does not assert \(H_\lambda=H\) on the full unprojected Fock space. Nor does it imply that an unprojected mean-field Hamiltonian obtained by decoupling the quartic term is exact. Mean-field factorization, saddle-point equations, and the treatment of gauge fluctuations are outside M1.1d.

## 6. Test contracts

Use independently constructed Jordan–Wigner Majorana matrices for a two-site finite Fock space. The test must build the spin operators directly as \(i c_i^0c_i^a\); production \(N^\gamma\) matrices may be used only in the explicit substitution comparison. All tolerances must be explicit and use zero relative tolerance.

### T00d-01 — Exchange tensor structure

For all three bond types, verify the tensor entries in the cyclic basis, the symmetric/antisymmetric decomposition, and
\[
J_{ji}=J_{ij}^T.
\]
Use nonzero independent \(J,K,\Gamma,\Gamma',D\) values so that no term is hidden by a zero coupling. In particular, verify that the DM contribution changes sign under transpose while the other four contributions do not.

### T00d-02 — Directed-edge counting

For a single physical edge with arbitrary real \(J_{ij}\), verify that the half-weighted two-orientation ordered sum equals the single A-to-B edge contribution. Include a nonsymmetric tensor to ensure the transpose convention is actually exercised.

### T00d-03 — Bilinear substitution and exchange prefactor

Verify \(i c_i^TN^ac_i=2\sigma_i^a\) for every component and verify that the \(1/8\) ordered-pair quartic expression equals the original \(1/2\)-weighted spin exchange operator before projection.

### T00d-04 — Field prefactor

For arbitrary site-dependent fields, verify that
\[
-\frac14\sum_{i,a}2h_i^a\,i c_i^TN^ac_i
=-\sum_{i,a}h_i^a\sigma_i^a.
\]

### T00d-05 — Constraint multiplier normalization

Verify
\[
\frac14\lambda_i^\gamma i c_i^TG^\gamma c_i
=2\lambda_i^\gamma T_i^\gamma
\]
and that its projection by \(P_G\) vanishes for arbitrary real multipliers.

### T00d-06 — Two-site projected Hamiltonian equivalence

For a two-site edge with generic real exchange tensor (including an antisymmetric component), arbitrary fields, and arbitrary real multipliers, construct \(H\) and \(H_\lambda\) independently and verify
\[
P_GH_\lambda P_G=P_GHP_G.
\]
Also verify that the unprojected operators need not be equal when multipliers are nonzero.

## 7. Exactness and scope

The edge-counting identity, exchange-tensor decomposition, Majorana substitution, multiplier normalization, and projected operator equivalence are exact finite-dimensional statements. Matrix residuals are numerical checks of these identities, not proofs.

This batch stops before any Hubbard–Stratonovich transformation, mean-field channel decomposition, Wick factorization, saddle-point approximation, self-consistency condition, or phase calculation. Those require separate mathematical and approximation contracts.
