# W001-T4 verification receipt — 2026-11-07

Scope: bounded mathematical/source audit and finite diagnostics. Production
implementations and test assertions were not edited. T3/PR #4 remains draft-gated;
the original braiding failure is not repaired or reclassified as expected failure.

## Commands and observed outcomes

- `git status --short --branch` and source/working diff inspected before work.
- `PYTHONPATH=orchestrator python3 -m agentctl project kitaev-mf status` from the
  Projects hub: current wave and source tree read. Its workspace canonical-branch
  label remains a stale external manifest field; no hub files were changed.
- `... project kitaev-mf hypothesis --id K-019 --explain`: provenance slice read;
  generated local-model commentary was not used as scientific evidence.
- Original primary archives read without executing TeX: tar-member bytes and
  SHA-256 matched the References pilot for cond-mat/0506438v3/hexagon.tex and
  1105.4573v2/HoneycombPhysical_Arxiv.tex. Selected TeX labels/lines inspected;
  this is not full transcript or rendered-paper certification.
- `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
  PYTHONDONTWRITEBYTECODE=1 .venv/bin/python scripts/audits/w001_foundation_checks.py`:
  observations collected. Legacy mismatches are recorded as mismatches, not PASS.
  Compact observations and source hashes are committed; raw output lives in
  ignored build/rnd/exact-perturbation-theory. Baseline source revision plus
  script/production hashes identify the dirty-tree diagnostic execution honestly.
- `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
  PYTHONDONTWRITEBYTECODE=1 .venv/bin/pytest -q -p no:cacheprovider --tb=short
  test/unit/test_feshbach_bounds.py test/unit/test_exact_perturbation.py
  test/unit/test_parity_selection.py`: 33 passed in 1.41 seconds. This establishes
  only their existing finite computational checks, not their docstring claims.
- `.venv/bin/ruff check --no-cache scripts/audits/w001_foundation_checks.py` and
  `.venv/bin/ruff format --check scripts/audits/w001_foundation_checks.py`: PASS.
- `.venv/bin/mypy --no-incremental scripts/audits/w001_foundation_checks.py`:
  PASS, one file. Broad package checks were not repeated for this new isolated
  diagnostic plus documentation; no production numerical code changed.
- `.venv/bin/python /tmp/kitaev-w001-verify.py`: JSON/YAML, state aliases, source
  inventory blobs, evidence paths, and preserved untracked hashes passed.
- `git diff --check`: PASS. Compact observation hashes and source-manifest
  archive/member hashes checked against actual files.
- Independent adversarial review: PASS for conditional finite algebra and bounded
  audit; reviewer did not execute the diagnostic script. No disputed physical
  prefactor, thermodynamic bound, braiding claim, or integration merge was accepted.

## Remote state at the inspected prior PR head

At `91d59bde6bf3b96e8ac221134d1a64dd659379ae`, GitHub reported CI failure,
blueprint build success, and LaTeX cancellation. These are not current-head
publication receipts and do not override the explicit draft merge gate.

## Handoff

T5 can use the conditional finite theorem, explicit unit map, distinct projector
contracts, and second-order sign identity, with their hypotheses. It must qualify
or exclude other unsupported statements while reorganizing the canonical book.
The required implementation repairs remain explicit backlog, not silently fixed.
The separate References chat is being steered by the PI; its latest metadata
reconciliation does not imply completed conversion or scientific acceptance.
