---
id: W001-T4
wave: W001
task: T4
title: Earlier perturbative contracts and claim-to-evidence review
status: completed
questions:
  - Which original projection and resolvent requirements are actually established by later code and tests?
  - Which finite algebraic results can support the canonical labbook without accepting an exact prefactor claim?
depends_on: [W001-T2, W001-T3-source-assembly]
tags: [conventions, physical-state-space, perturbation, exact-benchmarks]
artifacts:
  - doc/audits/w001_foundation_review.md
  - doc/audits/w001_acceptance.yaml
  - scripts/audits/w001_foundation_checks.py
---

## Contract before diagnostics

Review the original finite physical spin-space projectors, zero-flux block,
external energy window, and energy-dependent Feshbach map. Compare later claims
with primary source spans and implementation/test scope. Preserve the original
notes and all production numerical behavior. Do not repair braiding, change a
physical threshold, or replace the disputed prefactor by another assumed number.

Finite diagnostic objects are 2x2 and 4x4 Hermitian diagonal matrices with commuting
involutions, one-spin Pauli matrices, and the six-spin isolated-hexagon Hamiltonian
already used in the legacy tests. Toy examples test algebra/API properties and
are not claimed to model the thermodynamic honeycomb lattice. The hexagon has
open missing bonds, not periodic bulk geometry. Use standard Pauli matrices and
least-significant-bit site ordering, explicitly distinguished from other tensor
orderings. Only isolated finite computations are authorized by this intake.

## Preregistered evidence and criteria

- Energy-reference shifts must not remove genuine sector eigenvalues. Compare
  fixed-sector compression before and after a constant energy shift.
- A single plaquette sector does not define the joint flux-free sector on a
  multi-plaquette system. Compare an explicit four-state joint projector to a
  single-involution call; report the actual scope instead of a bulk conclusion.
- Check the one-spin Zeeman generator against explicit standard Pauli matrices,
  with exact array equality for the integer/imaginary-unit entries.
- Evaluate second order independently on the finite hexagon using its ground
  eigenspace and complementary restricted spectral sum. Report a nonzero witness
  or its absence, with spectral gaps and roundoff residuals. A numerical witness
  supports only this finite system; the general sign statement has its own proof.
- Read the pinned Kitaev source archive directly. Verify archive/file hashes and
  use TeX labels/line locators rather than invent rendered equation numbers.
- Re-run the existing focused perturbation tests without changing assertions;
  distinguish a passing implementation identity from a validated physical claim.

No Bayesian criterion applies to exact matrix identities or source attribution.
For the small eigensystems, record machine epsilon and residuals; a diagnostic
scale of 1e-10 is used only to group the known finite hexagon's degenerate ground
levels, and the observed gap must exceed that grouping scale decisively. Do not
infer fifth-order constants or thermodynamic convergence from these checks.

## Destination and workflow

Main book: explicitly scoped finite projectors/resolvents and convention glossary.
Wave appendix: full audit, counterexamples, source locators, and rejected claims.
No manuscript prefactor result. T4 stays on the single T3 integration candidate,
per PI steering to foundation review while PR #4 remains draft; this is a bounded
prerequisite review, not a new parallel research stream. Record dispositions in
the acceptance register and do not merge or publish on the basis of this audit.
