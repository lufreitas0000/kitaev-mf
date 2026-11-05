---
id: RND-P1-M1.4-T4
title: "Pure-Kitaev Energy Normalization Audit"
phase: 1
milestone: M1.4
task: T4
status: completed
claims:
  - CLM-M1-ENERGY-NORM-001
  - CLM-M1-ENERGY-NORM-002
  - CLM-M1-ENERGY-NORM-003
---

# Note 04c: Pure-Kitaev Energy Normalization Audit

**Status:** M1.4.3a analytical audit
**Parent notes:** [Note 00d](note_00d_exchange_tensor_and_quartic_hamiltonian.md), [Note 02](note_02_momentum_bloch_formalism.md), [Note 03](note_03_correlators_saddle_point.md), [Note 04](note_04_observables_energetics.md), [Note 04b](note_04b_energetics_test_specification.md)
**Scope:** Reconcile the pure-Kitaev energy with the exact-solution literature, fixing operator and normalization conventions. No production code or existing note is changed.

## 1. Result

The apparent factor-of-two discrepancy identified in the previous branch review is not a demonstrated implementation error. It arose from comparing different spin normalizations and different energy denominators.

The branch Hamiltonian uses Pauli operators,
\[
H_{\mathrm P}=K_{\mathrm P}\sum_{\langle ij\rangle_\gamma}
\sigma_i^\gamma\sigma_j^\gamma,
\qquad K_{\mathrm P}<0,
\]
with each physical edge counted once. The reference exact-solution convention uses physical spin operators,
\[
H_{S}=-K_S\sum_{\langle ij\rangle_\gamma}
S_i^\gamma S_j^\gamma,
\qquad S_i^\gamma=\frac12\sigma_i^\gamma,
\qquad K_S>0.
\]
For the same ferromagnetic model, \(K_{\mathrm P}=-K_S/4\), and hence
\[
H_{\mathrm P}=H_S.
\]
The reference value \(E_0=-0.3936K_S\) is **per two-site unit cell**, not per physical site. Therefore, for \(|K_{\mathrm P}|=1\), the corresponding Pauli-convention values are approximately
\[
\frac{E_0}{N_{\mathrm{cell}}}=-1.5744,
\qquad
\frac{E_0}{N_{\mathrm{site}}}=-0.7872.
\]
These agree, within the precision of the quoted reference decimal, with the branch's values \(-1.5745935\) per cell and \(-0.78729675\) per site.

<!-- claim: CLM-M1-ENERGY-NORM-001 -->
**Conclusion:** retain the Pauli-convention arithmetic in T04b-03. Do not replace its expected values by \(-0.3936\) per site. That number belongs to a different coupling convention and is quoted per unit cell in the cited source.

## 2. Literature convention and correlators

Choi, Klein, Rosch, and Kim define the Kitaev interaction as
\[
H_K=-K_S\sum_{a\text{-links}}S_i^aS_j^a
\]
and represent \(S_i^a=(i/2)\gamma_i^a\gamma_i^0\). They report the exact ground-state energy \(E_0=-0.3936K_S\) **per unit cell**, together with the zero-flux correlators
\[
\langle i\gamma_i^a\gamma_j^a\rangle=1,
\qquad
\langle i\gamma_i^0\gamma_j^0\rangle=-0.524864.
\]
These statements occur together in Sec. III.A, Eqs. (2), (13), and (14) of the cited preprint. The paper explicitly distinguishes this Majorana mean-field representation, which reproduces the exact energy in the Kitaev limit, from a naive Abrikosov-fermion decoupling whose energy is smaller by a factor of four.

Reference:
W. Choi, P. W. Klein, A. Rosch, and Y. B. Kim,
“Topological superconductivity in Kondo-Kitaev model,”
arXiv:1804.10212v1 (2018), Sec. III.A, Eqs. (2), (13), and (14).

In the branch convention, the pure-state covariance is
\[
U_b^{00}=-\bar U,\qquad U_b^{\gamma\gamma}=1,
\qquad V_s=0.
\]
The trace contraction in Note 04 is
\[
C_b^{\gamma\gamma}
=\operatorname{tr}\!\left[(N^\gamma)^TU_bN^\gamma U_b^T\right]
=2U_b^{00}U_b^{\gamma\gamma}=-2\bar U.
\]
Consequently,
\[
\langle\sigma_i^\gamma\sigma_j^\gamma\rangle_G
=-\frac12C_b^{\gamma\gamma}=\bar U,
\qquad
\langle S_i^\gamma S_j^\gamma\rangle_G
=\frac14\bar U.
\]
At \(\bar U\simeq0.5248657\), the latter is approximately \(0.1312164\), consistent with the reference energy: three bonds per two-site cell give
\[
-3K_S\frac{\bar U}{4}
\simeq -0.3936493K_S.
\]
The small difference from \(-0.3936K_S\) is the rounding of the literature value.

<!-- claim: CLM-M1-ENERGY-NORM-002 -->
This comparison is a convention-matched check: the branch's Pauli bond correlator is four times the physical-spin bond correlator, while the coupling multiplying it is correspondingly four times larger in magnitude.

## 3. Independent integral and convergence check

The pure matter dispersion in Note 02 gives
\[
f(\mathbf k)=1+e^{-2\pi i k_1}+e^{-2\pi i k_2},
\qquad
\bar U=\frac13\int_0^1\!dk_1\int_0^1\!dk_2\,
|f(\mathbf k)|.
\]
By symmetry of the three honeycomb bond directions, this is equivalent to the integral written in Note 03,
\[
\bar U=\int_0^1\!dk_1\int_0^1\!dk_2\,
\frac{1+\cos(2\pi k_1)+\cos(2\pi k_2)}
{|f(\mathbf k)|}.
\]
The integrand has isolated Dirac-point singularities in the second representation; the first representation is continuous but non-differentiable at those points. A midpoint Riemann sum of the first expression gives:

| Mesh \(N\times N\) | \(\bar U_N\) |
|---:|---:|
| 576 | 0.524865746616 |
| 1152 | 0.524865745946 |
| 2304 | 0.524865745863 |
| 4096 | 0.524865745850 |

The sequence supports
\[
\bar U\simeq0.5248657459.
\]
This updates the numerical precision of the previously quoted \(0.5248645\): the latter is adequate for the existing loose 36-by-36 regression tolerance, but it is not a high-precision value. The mesh sequence is numerical convergence evidence, not a rigorous error bound. A quadrature routine with an explicit error estimator or analytic treatment of the Dirac neighborhoods remains preferable before assigning certified digits.

<!-- claim: CLM-M1-ENERGY-NORM-003 -->
The integral identity and the converged midpoint values are independent of the energy-per-site convention. They validate the covariance magnitude, not by themselves the mean-field scalar subtraction.

## 4. Spectral and Wick-functional consistency

For the branch's pure-Kitaev Bloch spectrum, the positive eigenvalues are
\[
\mathcal E_m^+(\mathbf k)=\frac{|K_{\mathrm P}|}{2}|f(\mathbf k)|,
\qquad
\mathcal E_g^+(\mathbf k)=\frac{|K_{\mathrm P}|}{2}\bar U
\]
with three gauge flavors. With the full-BZ positive-eigenvalue counting convention of Note 04,
\[
\frac{E_{\mathrm{quad},0}}{N_{\mathrm{cell}}}
=-2\left\langle\mathcal E_m^+
+3\mathcal E_g^+\right\rangle
=-6|K_{\mathrm P}|\bar U.
\]
The scalar from the Gaussian exchange functional is
\[
\frac{\mathrm C}{N_{\mathrm{cell}}}
=3K_{\mathrm P}\bar U=-3|K_{\mathrm P}|\bar U.
\]
Thus
\[
\frac{E_{\mathrm{quad},0}-\mathrm C}{N_{\mathrm{cell}}}
=-3|K_{\mathrm P}|\bar U,
\qquad
\frac{E_0}{N_{\mathrm{site}}}
=-\frac32|K_{\mathrm P}|\bar U.
\]
For \(|K_{\mathrm P}|=1\), this gives approximately \(-1.57459724\) per cell and \(-0.78729862\) per site using the refined integral. The finite-precision test values in Note 04b are consistent with this result at their stated precision.

The spectral expression includes all four positive bands (one matter and three gauge flavors). The scalar subtraction is derived from the same Gaussian Wick functional; it is not an independently adjustable correction. Removing the gauge bands or changing the positive-eigenvalue counting without re-deriving the scalar would invalidate the equality.

## 5. Disposition and remaining work

1. Keep the Pauli-convention energy formulas and expected values in T04b-03. Update the prose in a later editorial batch only if a clearer explicit mapping to \(S=\sigma/2\) is needed.
2. Correct the earlier review statement that treated \(-0.3936\) as a Pauli-convention energy per site. It is a physical-spin-convention energy per two-site cell in the cited reference.
3. Treat \(0.5248645\) as a rounded legacy regression target, not a high-precision benchmark. The midpoint study supports approximately \(0.5248657459\), but does not supply a rigorous error estimate.
4. Keep M1.4.3 awaiting local WSL validation. This audit does not execute or validate the test suite and does not authorize implementation of the production energetics service.

No production service, test, existing note, or project-state entry is changed in this batch.
