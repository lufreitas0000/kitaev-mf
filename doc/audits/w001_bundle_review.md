# W001-T3 bundle review and CI repair — 2026-11-08 UTC

## Scope and user decision

PI supplied `kitaev_mf_w001_t3_and_exact_perturbation_bundle.zip` and requested
implementation plus PR #4 CI repair. Archive SHA-256:
`e34f26cf74933db5f2a0235121c32b76c64bf9661ec3fbcede1c62afcd16ed41`.
Raw text members and hashes are preserved under
`doc/rnd/exact-perturbation-theory/imports/w001_bundle/`. Bytecode is not executed.

After the model distinction was explained, PI explicitly chose to retire the
unsupported four-flavor pi-holonomy guarantee from required tests while retaining
its failing research reproduction. This supersedes the earlier instruction to
keep that assertion as an integration gate; it does not make the claim true.
The full assertion body is retained in
`test/research_checks/check_legacy_mean_field_braiding.py`. Run it explicitly
with pytest. No threshold was relaxed and no xfail was added. Public legacy
classes remain available to existing scripts. Historical fitted benchmark rows
are retained as historical data, not independent evidence of a physical braid.

## Actual GitHub failure

At checkpoint `da17ffe3cc8ab03301eab832482a864b3a185078`, CI run 42423051113,
job 128555418191, failed during pytest collection: both Phase 3 plotting test
modules imported unavailable `matplotlib`. Lockfile verification, Ruff, and Mypy
had passed. Blueprint and LaTeX runs at that head succeeded. The lockfile input
omitted a runtime dependency already declared in `pyproject.toml`.

Add `matplotlib>=3.8.0` to `requirements.in` and regenerate both lockfiles using
Python 3.12 pip-tools. The existing locked packages are not upgraded; the diff
adds matplotlib's missing dependency closure and constraint annotations. The
old braiding failure is separate and would become visible after collection.

## Braiding: separate model, validated finite transport

The new service is `services/phase_2/chiral_vortex_transport.py`, not a replacement
of `vortex_braiding.py`. It uses H=(i/4)c^T A c on an LxL periodic honeycomb graph,
NN K>0 and prescribed nonzero NNN kappa. It is not the finite-field four-flavor
self-consistent Hamiltonian, a physically projected spin braid, or a many-body
non-Abelian braid representation.

Useful proposal components retained:

- NNN hopping contains both NN gauge-link factors. Arbitrary independent site
  sign conjugations are verified against direct transformed-link assembly.
- A simple closed dual loop is the bond boundary of a declared site region.
  The final Hamiltonian must close exactly under that site's sign gauge.
- Consecutive central rank-two eigenspaces are parallel-transported with SVD
  polar factors; eigenvalue phases refer to this one-particle subspace.

Additional admission checks reject missing custom data, malformed/disconnected
cycles, path collisions with the stationary vortex, nonmatching site boundaries,
zero interpolation count, absent spectral isolation, insufficient overlap,
failed Hamiltonian closure, and failed subspace closure. The overlap threshold
is a declared discretization criterion (default 0.1), not a physical theorem;
users must refine if it fails. Sampled gap positivity does not prove positivity
between samples or in a thermodynamic limit. Selecting two central modes alone
does not prove their Majorana localization.

For L=10, K=1, kappa=0.12, radius=1, center=(4,4), six and ten interpolation
steps per crossed bond, enclosing fixed plaquette (4,4) yields -I and nonenclosing
fixed plaquette (0,0) yields +I within the registered finite tolerances. Final
Hamiltonian closure is exact for these sign operations; spectral subspace
separation and successive-overlap lower bounds are asserted separately. These
are bounded numerical results for this model and these paths.

## Perturbation: runnable candidate, unresolved interpretation

`doc/rnd/exact-perturbation-theory/code/w001_gaussian_quadrature.py` is an isolated
research implementation adapted from the supplied code. It removes external
absolute-path imports and hardcoded writes, supplies a CLI, uses explicit types,
and rejects invalid orders, reference zero modes, nonpositive vacuum-energy
differences, ill-conditioned matrices, and unresolved determinant phases.
It continues to label all coefficient output as an unaccepted unprojected
candidate. It does not repair an all-orders derivation by relabeling it.

Concrete blockers to accepting the supplied note:

1. It writes D=-i b^x b^y b^z c=+1. With four distinct Hermitian CAR Majoranas,
   the product is Hermitian and squares to +1; multiplying it by -i is
   anti-Hermitian and squares to -1. It cannot be the asserted constraint.
2. A gauge-sector Gaussian vacuum is not necessarily a physical spin state.
   Matter parity, global sectors, and boundary conventions are not supplied by
   an unprojected vacuum-energy difference. A full flux-free projector is not
   automatically a degenerate ground-space projector for canonical SW.
3. On L=3 the candidate has four eigenvalues within 1e-10 of zero. Its proposed
   transpose/complement identity has Frobenius defect about 1.98, so zero-mode
   occupation must be specified; the adapted numerical routine rejects it.
4. Bounded factors do not by themselves prevent singular solves. A principal
   complex square root does not specify a continuous many-body overlap branch.
5. Fixed quadrature orders are approximate. For L=4, orders (12,6) produce total
   13.6815162963 and (24,10) produce 13.5611188552 in the proposal's definitions.
   This is a refinement difference, not a certified error bar or physical
   prefactor. L=8 numbers in the note are not independently certified here.
6. Power counting and a Hartree-Fock closure do not prove that all higher orders
   reduce to two parameters or protect topology at arbitrary coupling. Those
   statements, lifetime laws, and cross-project walk formulas remain unaccepted.

A useful finite positive check is retained: the single-slice overlap and two-point
insertion formulas agree with explicit 4x4 Fock-space matrix exponentials for
a noncommuting four-Majorana example at tau=0,0.3,1.7 (residuals below 7e-16 in
this run). This supports that kernel on those cases, not the multi-slice third
order, physical parity map, or the claimed coefficient interpretation.

## Next scientific work

Before accepting physical perturbative coefficients, derive the correctly
normalized spin/Majorana operator map and the chosen P/Q reduction; establish
physical parity in each intermediate sector; verify multi-slice insertions and
the graded coefficient extraction against independent small Fock/spin matrices;
then perform converged quadrature, size, and boundary studies. The supplied raw
note and runnable prototype make that work inspectable without replacing the
earlier seven notes or endorsing disapproved historical Phase 4 results.

T5 canonical book assembly remains in progress with separately committed chapter
drafts; its include map and PDF are not yet complete. T6 PRL migration and T7
bounded main promotion follow. No main merge or publication is part of this fix.

## Remote numerical follow-up

After the dependency repair, CI run 42437696341 exposed three inherited Phase 3
failures (352 passed, one skipped). The new bundle tests passed. The failed
assertions assumed a unique null-mode vacuum on a cylinder and a basis-specific
truncated weight for one of two degenerate lowest torus modes. Local runs with
1/2/4 BLAS threads did not reproduce the same mismatch, so no thread-only cause
is asserted or used to mask the issue.

An independent exact counterexample also established a shared solver defect:
for a real Schur form diag(0,2J,0), blind consecutive pairing returned [0,0] and
reconstruction error 2. The new common canonical-skew helper collects actual
rotation blocks and null columns before pairing. Both existing entry points use
it, preserving the scalable-module public imports. Rank-deficient spectra are
checked against an independent Hermitian diagonalization. A relative numerical
cutoff is explicit; returned basis orientation at zero is not a physical parity
certificate. The exact-arithmetic reason is that a matrix which is both skew and
real quasi-upper-triangular has only 2x2 skew rotation blocks and 1x1 zeros. Its
null columns span an invariant nullspace and can be paired after collection.

The cylinder SVD fixture now declares its null-mode vacuum and verifies that
swapping the chosen pair acts inside the actual nullspace. It retains the exact
SVD/exhaustive equivalence and overlap thresholds; a regression exercises both
returned null orientations. Environment-dependent wall-clock speed is no longer
a required algebraic assertion. The historical timing measurement is not renewed
as a performance certification.

The torus test now compares the same explicitly chosen one-particle state to an
independent 512-dimensional sparse-CAR construction and projection onto quenched
number 0 or 2. It checks both the returned degenerate basis and a 0.37-radian
rotation inside it. Agreement uses absolute 1e-10, tighter than the old relative
1e-3 golden comparison. This is a truncated finite-Fock check, not proof of a
full physical spin sum rule. A production policy for degenerate physical states
remains research work.
