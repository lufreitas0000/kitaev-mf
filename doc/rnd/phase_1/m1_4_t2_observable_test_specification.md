---
id: RND-P1-M1.4-T2
title: "Observable Test Specification"
phase: 1
milestone: M1.4
task: T2
status: completed
contracts:
  - T04a-01
  - T04a-02
  - T04a-03
---

# Note 04a: Observable Test Specification

**Status:** M1.4.2 test specification and implementation
**Parent note:** [Note 04](note_04_observables_energetics.md)
**Scope:** Independently validate the Gaussian one-point and intersite spin-observable contractions before adding production observable services.

## 1. Testing strategy

The observable formulas in Note 04 are tested against direct finite-dimensional Fock-space expectations, not against a second implementation of the same covariance formula.

The test fixture constructs eight Jordan–Wigner Majorana operators, corresponding to two sites with four Majoranas per site. A generic real skew-symmetric quadratic form
\[
H_{\mathrm q}=\frac{i}{4}\sum_{\mu,\nu=0}^{7}
A_{\mu\nu}c_\mu c_\nu,\qquad A^T=-A,
\]
is diagonalized exactly. Its nondegenerate lowest eigenvector is a pure fermionic Gaussian state. The deterministic matrix seed fixes reproducibility; the finite-dimensional diagonalization is a numerical oracle for the exact finite-Fock-space expectation values.

The covariance matrices are then measured directly:
\[
V_s^{\mu\nu}=\langle i c_{s}^{\mu}c_s^\nu\rangle,\quad \mu\ne\nu,
\qquad
U^{\mu\nu}=\langle i c_0^\mu c_1^\nu\rangle.
\]
The diagonal of $V_s$ is set to zero, consistent with the subtracted CAR contact term in Note 04. All measured covariance entries are checked to be real within the stated tolerance.

## 2. Test contracts

### T04a-01 — One-point Pauli and spin normalization

For each component $a$, compare the direct Fock expectation of
$\sigma_0^a=i c_0^0c_0^a$ with
\[
\langle\sigma_0^a\rangle=\frac12
\operatorname{tr}(V_0^TN^a).
\]
Then compare the physical-spin expectation of $S_0^a=\sigma_0^a/2$ with
\[
\langle S_0^a\rangle=\frac14
\operatorname{tr}(V_0^TN^a).
\]
This test fixes the Pauli-versus-spin factor independently of a production observable service.

### T04a-02 — Two-site Wick contraction

For all nine component pairs $(a,b)$, compare the direct Fock expectation
$\langle\sigma_0^a\sigma_1^b\rangle$ with
\[
\frac14m_0^am_1^b
-\frac12\operatorname{tr}\!\left[(N^a)^TU N^bU^T\right],
\qquad
m_s^a=\operatorname{tr}(V_s^TN^a).
\]
Also verify the physical-spin correlator is one quarter of the Pauli correlator. The sites are distinct, so no on-site contact-term identity is used.

### T04a-03 — Hartree and exchange terms are both exercised

For a deterministic component pair in the seeded Gaussian state, verify that the Hartree and intersite exchange contractions are individually nonzero and that their sum gives the Wick result. This prevents a test from passing merely because one contribution vanishes for the chosen state.

## 3. Evidence and limitations

These tests validate the stated formulas for a nondegenerate pure Gaussian state on a two-site finite Fock space. The direct Fock-space diagonalization and matrix expectations are numerical evaluations of exact finite-dimensional definitions, with absolute tolerance $2\times10^{-12}$ and zero relative tolerance.

They do not establish equivalence between unprojected Gaussian expectations and expectations in the exactly projected physical spin state. They do not test finite-temperature mixed Gaussian states, arbitrary lattice sizes, or the Brillouin-zone integration routine. They also do not test a production observable API: the observable service is a later M1.4.4 deliverable, and these tests deliberately keep the analytical oracle independent of that future implementation.

The quadratic test Hamiltonian is a fixture, not a physical Kitaev-model benchmark. The tests check operator/covariance normalization and Wick factorization only.
