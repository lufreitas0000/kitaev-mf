---
id: W001-T3-BUNDLE
wave: W001
task: T3-bundle
status: completed
title: External braiding and exact-perturbation candidate integration
questions:
  - Which supplied components satisfy independent finite contracts without changing the physical model silently?
  - What causes PR 4 CI to fail before the known braiding test runs?
depends_on: [W001-T4]
tags: [defects, perturbation, numerics, physical-state-space]
artifacts:
  - doc/audits/w001_bundle_review.md
---

## Contract and evidence registered before probes

PI requested inspection and implementation of the downloaded W001 bundle on
2026-11-08 UTC. Its SHA-256 is
`e34f26cf74933db5f2a0235121c32b76c64bf9661ec3fbcede1c62afcd16ed41`.
Preserve the original sources and manifest. Ignore bundled bytecode for import;
do not overwrite existing public APIs or the original failing braiding assertion.
Keep historical Phase 4 disapproved. T5 chapter drafts have a separate checkpoint.

GitHub CI at da17ffe fails on two collection imports of missing matplotlib;
lockfile/lint/type stages pass. Align declared runtime dependencies and lockfile
specifications using the existing pip-tools workflow; no unrelated upgrades.

For braiding, distinguish a finite real-skew chiral matter generator with
K>0, fixed kappa, periodic LxL honeycomb graph from the production four-flavor
finite-field mean-field model. Independent checks: local gauge covariance under
arbitrary site signs, flux endpoints and exact final gauge restoration, spectral
separation of the central rank-two subspace, overlap singular values, orthogonal
or unitary closure, and refinement of transport discretization. Reject invalid
paths, zero interpolation count, ambiguous boundaries, singular transport, and
failed closure rather than unconditionally taking a polar factor.

For perturbation, preserve candidate coefficients as unaccepted. Inspect exact
CAR constraint signs; zero-mode projector completeness; compare Gaussian overlap
and insertion formulas with explicit small Fock-space exponentials before any
physical interpretation. A pinned-gauge vacuum-energy difference is not a
physical vison threshold without parity/gauge projection. Fixed quadrature is a
numerical approximation requiring convergence evidence, not exact integration.
No bulk/all-orders coefficient or topology assertion is accepted by import.

Use finite deterministic matrix diagnostics, test original and adapted functions,
and obtain Terra 6.1 adversarial review. Bayesian thresholds are inapplicable to
algebraic identities and package imports. Record failures honestly. Full suite
runs in an ignored disposable copy because existing CLI tests write figures.
Original braiding remains a gate unless the PI explicitly changes its status.

## Outcome

PI authorized retiring the unsupported legacy guarantee from required tests;
its original assertions remain an explicit failing research check. Separate
chiral transport and isolated unprojected quadrature are implemented and reviewed.
No physical SW coefficient, parity map, all-orders closure, or mean-field braid
is accepted. Local full-suite and focused-check evidence is in the verification
receipt; current-head GitHub outcomes must be inspected after push.

## CI follow-up: Phase 3 degenerate-state fixtures

Remote CI after collection repair exposes three inherited failures, while new
bundle tests pass. Cylinder quenched H has an exact null pair, so its vacuum
orientation and det(X) singularity are not fixed by H alone. Torus L=3 has a
twofold-degenerate lowest positive energy; occupied_mode=0 is a chosen state,
not a basis-invariant observable. Do not loosen the original numerical tolerances
or set a backend flag to hide this. Register explicit null-vacuum selection for
the singular-path fixture and replace its arbitrary golden spectral weight by
an independent finite-Fock projection for the same declared initial state.
Keep performance measurements separate from required algebraic correctness.
These fixture contracts do not establish projected-spin or thermodynamic claims.

### Schur block parsing regression contract

A real Schur form can contain separated 1x1 null blocks. Pairing columns in
positions (0,1),(2,3),... is invalid if a null block precedes a 2x2 rotation.
Register an exact 4x4 witness T=diag(0,2J,0), Z=I; it must reconstruct T with
energies [0,2]. Group actual nonzero 2x2 Schur blocks and pair the collected null
columns only afterwards. Reuse one low-level implementation for both existing
entry points, retaining public exports. At a zero mode det(Q) is a chosen basis
orientation, not sign(Pf(A))=0. Verify reconstruction on singular matrices and
reject malformed matrices; do not assume null-state physics is uniquely fixed.
