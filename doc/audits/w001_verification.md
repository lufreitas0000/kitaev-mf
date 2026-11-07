# W001-T1 verification receipt — 2026-11-07

Scope: policy and preservation inventory, based on local main
`39e733e` (full SHA recorded as policy_base in w001_source_inventory.json).
This is not a physics acceptance, branch merge, build release, or remote CI receipt.

## Commands and results

- `.venv/bin/python /tmp/kitaev-w001-verify.py`: PASS. The task-local verification
  script parsed PROJECT_STATE and acceptance YAML plus inventory JSON; checked
  matching W001 state aliases, unique review IDs and evidence/target file paths;
  checked all 210 documentation-path blob records against recorded immutable
  commits, existence of source commit objects, and SHA-256 for the five untracked
  Phase 3 preservation inputs. The temporary script is not a production tool.
- `git diff --check`: PASS after removing a trailing blank line.
- `python3 /home/lucas/Projects/scripts/sync_agents.py --dry-run`: exit 0,
  preview only; no workspace policy propagation. The synchronizer has no project
  selector; project-local policy did not authorize changing other repositories.
- Adversarial internal review: PASS after fixing legacy-only naming/state wording
  and adding the earlier perturbation remote ref to inventory. Reviewer did not
  independently rerun the validator. Review applies to policy/inventory only.

## Not run / limits

No numerical suite, Ruff/Mypy, simulation, TeX build, manuscript migration, or
remote CI was run for these documentation-only changes. Existing numerical and
static failures from the previous audit are not resolved by this task. The
standalone References ingestion chat remains separate ongoing work.

The preserved perturbation ref is an observed historical tip, not a claim that
its contents are pristine. T2 must compare original notes and later changes.
T1 does not restore an earlier snapshot over current source. Recent historical
Phase 4 calculations remain disapproved-by-PI. Legacy main-based claims are not
recertified. No main-targeted PR may silently include the seven inherited commits
beyond observed origin/main; baseline reconciliation precedes promotion.
