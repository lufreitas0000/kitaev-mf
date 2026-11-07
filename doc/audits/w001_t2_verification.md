# W001-T2 verification — 2026-11-07

Scope: preservation, corrected provenance, and version-control procedure.
No scientific source was edited; no simulation, TeX build, or remote CI result
is claimed. This task does not resolve the prefactor problem.

- `.venv/bin/python /tmp/kitaev-w001-verify.py`: PASS, parses JSON/YAML, checks
  state aliases/evidence paths, immutable documentation blobs and untracked hashes.
- `.venv/bin/python -` (focused preservation assertions): PASS, all six archive
  refs match recorded commits; snapshot parent equals original Phase 3 tip; exactly
  five added paths match inventory; snapshot content SHA-256 equals source files;
  original Phase 3 HEAD, tracked files, and index unchanged; seven earlier notes
  match the earlier research ref at all three later source tips.
- `git diff --check`: PASS.
- `python3 /home/lucas/Projects/scripts/sync_agents.py --dry-run`: exit 0, preview
  only, no cross-project propagation.
- Internal adversarial review: PASS. Reviewer independently checked snapshot
  parent and exactly five paths, and confirmed source files remain untracked.
- `git push origin refs/tags/archive/w001/policy-base refs/tags/archive/w001/phase2
  refs/tags/archive/w001/phase3 refs/tags/archive/w001/phase4
  refs/tags/archive/w001/earlier-perturbation-ref
  refs/tags/archive/w001/phase3-untracked-inputs`: PASS (one bounded command).
- `git ls-remote --tags origin 'refs/tags/archive/w001/*'`: inspected all six remote
  targets; each matched the inventory. Archive tags are not published releases.

The five unreviewed files are preserved remotely without modifying or cleaning
that worktree. Current assessments now say earlier work was ignored, not presumed
altered. The dated T1 observation remains in history with an explicit correction.
No main-targeted PR was opened: the seven inherited local-main commits still
require T3 baseline selection/reconciliation. Existing CI filters do not yet
support the whole new branch scheme; this is a declared T3 deliverable.

The References pilot has produced five archived source packages and candidate
conversions/cards; two conversions failed parsing and others have warnings.
Fidelity remains pending. The user is directly steering that chat's reference
list, so this task does not override its waiting state or promote its outputs.
