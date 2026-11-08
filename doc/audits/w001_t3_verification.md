# W001-T3 verification receipt — 2026-11-07

## Scope

Source assembly and software/interface hygiene are complete on the candidate.
Final integration merge is NOT approved: PI chose foundation review while the
original braiding failure stays unresolved. Main is unchanged remotely.

## Commands and results

- Pre-integration at a1026d8: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
  MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/pytest -q
  -p no:cacheprovider --tb=short test/`: 241 passed in 19.22 seconds.
- Pre-integration `.venv/bin/ruff check --no-cache src test scripts --statistics`:
  PASS; `.venv/bin/mypy --no-incremental src test scripts`: PASS, 103 files.
- Raw combined source suite: 349 passed, 1 failed in 156.93 seconds.
- Repaired combined suite: 349 passed, 1 failed in 170.81 seconds. Run with
  single-threaded BLAS, MPLBACKEND=Agg, PYTHONDONTWRITEBYTECODE=1, PYTHONPATH pointing
  to the disposable copy's src, and the repository venv Python executing
  `-m pytest -q -p no:cacheprovider --tb=short test/ --junitxml=...` from
  `build/w001/integration-check`. Logs/XML stay under ignored build/w001.
- Untouched Phase 2 archive reproducer: the same enclosing-loop test failed in
  2.78 seconds, returning 0.5904 radians against the original near-pi criterion.
- Focused final metrics/path checks: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
  MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/pytest -q
  -p no:cacheprovider --tb=short test/unit/test_phase4_analytical_mf.py
  test/unit/test_phase2_m2_1_atlas.py`: 8 passed in 5.72 seconds. The added tuple
  assertions verify container compatibility, not scientific acceptance.
- Final `.venv/bin/ruff check --no-cache src test scripts` and
  `.venv/bin/ruff format --check src test scripts`: PASS, 167 files formatted.
- Final `.venv/bin/mypy --no-incremental src test scripts`: PASS, 167 files.
- Workflow YAML/routing assertions: PASS for exploration push/PR paths, boolean
  publication default false, and manual-main publication guards.
- `git diff --check`: PASS. Source accounting and preserved-input hashes checked.
- Internal adversarial review: PASS for draft/source-preservation scope after
  restoring re-exports and original zip behavior, documenting tuple adaptation.
  This is explicitly NOT merge approval or independent physics validation.

## Limits and remaining gate

No numerical tolerance was relaxed and no xfail added. The test failure remains
visible. Existing tests have limited scientific scope; their success does not
certify reproduction, non-Abelian braiding, or perturbative prefactors. No TeX
build or document publication was run. No remote CI success is claimed.
The raw source and original mathematical notes remain preserved by ancestry and
archive tags; software normalization is separately enumerated in the source map.
