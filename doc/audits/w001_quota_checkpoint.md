# W001 quota checkpoint and resume handover — 2026-11-07

## Pause trigger and authorization

Lu Freitas requested a smooth pause when the five-hour Codex quota reached 5%
remaining, a three-hour continuation timer, and a Git commit/push checkpoint.
The live app usage-limits tool reported 96% consumed in the 300-minute window
at approximately 22:56 UTC (4% remaining); the trigger was reached. This is an
operational pause, not completion of W001 or acceptance of disputed physics.
The quota observation is dated telemetry, not durable project state.

## Preserved state

Before this checkpoint, the primary checkout was clean at
`e04638833910c07f448da4ec4e16bce92c1c7735`, on
`research/w001-t3-source-integration`, synchronized with its origin branch.
That commit completes the bounded T4 foundation audit. T1 and T2 are complete;
T3 source assembly exists but its integration merge remains gated. PR #4 is draft:
https://github.com/lufreitas0000/kitaev-mf/pull/4, targeting `exploration`.

The PI explicitly selected: continue foundation review and keep braiding
unresolved. Preserve its original failing assertion; no threshold relaxation,
xfail, or main promotion. Historical Phase 4 calculations remain disapproved by
the PI. The earlier exact-perturbation work was ignored, not presumed altered.
Seven original notes were preserved and compared byte for byte during T2.

No T5 source, labbook chapter, mathematical claim, numerical output, or PDF has
been edited in this preflight. The current canonical PDF is an inherited build;
it is not a newly reviewed consolidated book. Main and exploration have not been
advanced by merging this candidate. Original worktrees and archive refs remain.

## Next bounded task: T5 canonical labbook

The requested ten-minute reading interval after T4 has elapsed. Read the live
status, applicable rules, and existing diff before resuming; do not assume the
checkout stayed unchanged during the pause. Keep one active integration stream.
Write the T5 intake in `doc/rnd/waves/` before implementing the book changes.

T5 must produce actual readable chapters and technical appendices, not just an
outline. One root is `doc/labbook/labbook_kitaev_mf.tex`. The main narrative follows
the argument; technical appendices preserve discovery chronology and identify
current claim status. Historical phase/milestone paths remain provenance IDs.
A standalone incremental-note template is a review artifact with base revision,
status, and canonical destination, not a separate official book per branch.

Suggested logical flow, revisable with the science:

1. Physical question, finite FM Kitaev model, Pauli/spin-half conventions, scope.
2. Physical spin space versus enlarged Majorana space and exact projection.
3. Gaussian states, covariance conventions, Wick contractions, and pure-state space.
4. Mean-field stationarity, weak versus exact constraints, numerical reliability.
5. Bulk/defect states and what flux measurements do and do not establish.
6. Perturbative and ED comparison contracts, projectors, resolvents, gap definitions.
7. Thermal ensembles, response normalization, and bounded open questions.

Put complete proofs of newly accepted statements in an explicitly scoped
technical appendix. New foundational statements beyond T4 require their own
proof/referee review; do not promote them merely through an editorial rewrite.
Preserve the eleven existing Phase 1 chapters as historical technical records.
The Phase 3 chapter and five Phase 4 chapters need visible status warnings;
Phase 4 is disapproved, not accepted by inclusion. Preserve equations and source
identity, recording any rendering changes. Map all R&D notes to their canonical
source and incorporation status; do not claim full incorporation of notes that
are only indexed. Select original perturbation and Phase 2 records for faithful
rendering where useful, rather than silently dropping their derivations.

## Preflight findings that affect the build

- Existing root includes only eleven Phase 1 chapters; Phase 3 is commented out.
  Its title still names Phases 1 and 3. It has no active bibliography.
- Seventeen tracked chapter sources total 3,486 lines: eleven Phase 1, one Phase 3,
  five Phase 4. There are 113 Markdown R&D records (about 1 MB); source existence
  is not acceptance. Avoid an indiscriminate main-narrative dump.
- Three Phase 3 figures use ignored PDFs, so the old chapter is not a clean-source
  reproducible include. Their tables are committed under the Phase 3 benchmark
  data area. Render archival displays directly from those tables if included.
- The original Phase 3 captions claim L=10 and localization length 1.42, while
  committed fig3/fig4 headers say L=6 and xi=3.606412. The comparison table uses
  L_exact=6 and L_mf=4. Show the discrepancy explicitly; do not silently revise
  data or present those displays as validated reproduction of a paper.
- Existing pure-state text calls O(2n)/U(n) an ordinary real Grassmannian and
  treats generic momentum blocks too independently. Review conventions before
  reusing this in the logical narrative. Particle-hole constraints link k and -k.
- A Gaussian variational energy is an upper bound only for the declared trial
  state space and Hamiltonian. Physical spin-space bounds require a physical
  trial state or justified projection; average gauge constraints do not supply it.
- For H=(i/4) gamma^T A gamma and Omega_ab=(i/2)<[gamma_a,gamma_b]>, derive the
  single complex-mode block before fixing signs, occupation, thermal covariance,
  Wick factors, or Pfaffian parity conventions. State the Majorana ordering.
- For a fixed physical finite Hamiltonian the Gibbs variance identity applies.
  A temperature-dependent self-consistent A(T) needs the implicit derivatives;
  one cannot copy the fixed-A identity as a complete thermodynamic calculation.
- An independent binary vison toy with cost 0.137 has W=tanh(Delta/(2T)); it does
  not imply the old proposed >90% drop across [0.02,0.15]. This toy is not a
  derivation of the coupled flux/matter ensemble on a torus.

## Authoritative foundation and verification records

Use `doc/audits/w001_foundation_review.md`, `w001_foundation_observations.json`,
`w001_acceptance.yaml`, and `w001_t4_verification.md`. T4 permits the conditional
finite-dimensional Schur/Feshbach theorem and unit substitution, not a computed
bulk prefactor, a full-Q threshold certificate, or an all-orders result.

Existing focused T4 tests: 33 passed. T3 repaired candidate: 349 passed and one
original braiding failure, also reproduced on the untouched Phase 2 source.
Full numerical suite need not be rerun for this checkpoint documentation.
Do not claim latest-head remote CI success without inspecting its actual run.

For T5, build from a disposable ignored source copy to expose reliance on ignored
figures and avoid tests overwriting committed benchmark assets. Use LaTeX build,
PDF text/visual inspection, source/include mapping, bibliography checks, and
independent adversarial review. The existing referee is `/root/w001_policy_review`
(Terra 6.1); continue read-only review through that agent if it is available.
Commit/push coherent changes with Lu Freitas's configured author identity.
If committing a reviewed draft PDF, bind it to the source revision and hash;
publication from main is a separate gate. Do not relabel the old PDF as current.

## References companion and subsequent tasks

References chat `01a1b7f8-1799-7820-af9d-ed9dbbf3e0ee` is user-owned and separately
steered by the PI. Do not create another chat or override its scope. Pilot files:
`/home/lucas/Projects/References/references-transcripts/collections/kitaev-mf/pilot-2026-11-07/`.
The project manifest `doc/reference/literature_manifest.yaml` pins the raw primary
archives and reviewed TeX spans. Metadata reconciliation is not transcript
certification. Preserve originals and keep interpretation in paper cards.

After T5, T6 owns neutral manuscript naming and PRL/REVTeX formatting; T7 owns a
bounded reviewed promotion to main. The goal is not complete at this checkpoint.
At timed continuation, inspect live quota and checkout first. If quota remains at
or below 5%, do not begin a large task or claim it reset. Report the condition and
keep a safe checkpoint; any new scheduling should remain bounded and explicit.

## Timer and checkpoint verification

The app created heartbeat `resume-kitaev-mf-consolidation-after-quota-pause`,
attached to this chat. Its finite schedule has one future occurrence, approximately
2026-11-08 01:57:46 UTC
(2026-11-07 22:57:46 America/Fortaleza).
The scheduler excludes its already-past start occurrence. No recurring
indefinite job was requested. Timed execution must recheck quota and live state.

The checkpoint received an independent read-only referee PASS for accurate task
status and handover scope. No numerical tests were rerun for this documentation
only change. `git diff --check` and staged scope inspection are required before
commit; the task branch is pushed without merging or publishing.
