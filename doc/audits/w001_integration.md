# W001-T3: Provisional source integration

## Intake and contract

- Question: can the preserved research inputs be assembled on one exploration
  line without changing accepted remote main or endorsing disapproved results?
- Prerequisites: W001 T1/T2, archive checkpoints, source inventories, clean primary
  worktree, and current remote refs. T2 reading interval observed before edits.
- Scope: Git integration, evidence-status records, and workflow routing. Preserve
  original physics formulas, archived sources, and other worktrees. Mechanical
  typing/import repairs may adapt the assembled software to existing interfaces. No scientific acceptance or new
  physics calculation. Retain inherited failures visibly rather than changing
  scientific thresholds to make tests pass.
- Evidence: input ancestry, conflict resolutions, source-file accounting, baseline
  and candidate interface tests, state consistency, and workflow trigger checks.
- Destination: chronological consolidation appendix and repository governance.
- Tags: numerics, conventions, perturbation. Integration target: exploration.
- Stop: reviewed provisional integration, pushed research branch/PR and exploration;
  no promotion of the full inherited tree to main and no automatic publication.

## Baseline decision

Observed remote main is `dd7cb44044a6c1fe6b760b92c102001144dc427b`.
Local main `39e733e8e920f24702a1e3413b9ec9f187ddfb09` contains seven additional
inherited commits. Preserve that local line as `archive/w001-local-main`, then
create local main tracking origin/main. This renames rather than resets/discards
local work. Its immutable archive tag already exists remotely. Research and T1/T2
policy changes enter exploration; remote main remains unchanged.

## Pre-integration checks

On T2 commit `a1026d8`: 241 tests passed in 19.22 seconds with single-threaded BLAS;
Ruff reported no errors; Mypy reported no issues in 103 files. Tests establish
only their finite computational contracts, not the inherited physical claims.

## Source selection and conflict resolutions

All three historical tips and the Phase 3 snapshot are ancestors of the candidate.
The source-specific Python, tests, scripts, and compact data were initially
imported byte-for-byte. Ten shared Python paths from the older Phase 3 lineage
retain the newer Phase 2/4 versions, rather than downgrading the shared solver.
The source-accounting manifest records those selections explicitly.

W001 PROJECT_STATE and CHANGELOG were retained over incoming blanket completion
labels. Every unique table row from incoming discovery logs was retained.
The original Phase 3/4 state and changelog remain retrievable at archived refs.
The Phase 3 chapter is present but its canonical include is commented pending
T4/T5 review. The Phase 4 addition to the default Makefile was not adopted; its
new standalone PDF is retained at the archive ref and in ignored local storage,
not among active publication outputs. Original historical notes remain in place.

## Verification findings and bounded repairs

Raw combined tree: 349 tests passed, one braiding test failed; Ruff 165 findings;
Mypy 43 errors across 15 files. Broad imported-code hygiene was repaired without
changing physical formulas or original assertion thresholds: import ordering and
unused imports/locals, type annotations, formatting, missing Sequence import,
and NumPy result typing. Intentional path-bootstrap imports have scoped E402
comments, and historical project imports remain explicit re-exports where needed.
Two benchmark scripts now resolve output relative to their own repository instead
of writing into the old Phase 3 worktree.

The Phase 4 solver now provides tuple histories to the established SolverMetrics
interface instead of mutable lists; this is an explicit container compatibility
change, tested by assertions in the existing regression. The numerical iteration
and stopping criteria are unchanged. The malformed-plaquette zip retains the
original truncation behavior with strict=False; no stricter input contract was
introduced for typing compliance.

The full test suite runs in a disposable ignored copy, because inherited Phase 3
CLI tests regenerate figure/data files. Existing tracked project datasets are
not replaced by test output. Scientific acceptance is not inferred from test
success or static cleanup.

## Unresolved braiding gate and PI steering

The enclosing-loop test expects a phase near pi but returns about 0.5904 radians.
The identical failure reproduces on untouched archive/w001/phase2 in 2.78 seconds.
No tolerance or expectation was changed and no xfail was introduced.
Lu Freitas chose: continue the foundation review, keep braiding unresolved, and
leave this integration PR in draft. Thus T3 source assembly is ready for review,
but final merge into exploration is deliberately gated. T4 foundation review may
proceed on this single candidate line without returning to parallel phase branches.
