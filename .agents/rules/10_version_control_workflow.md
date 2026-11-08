# Rule 10: Branch, PR, Data, and Document Version-Control Workflow

Rule 01 owns history safety and data-storage tiers. This rule describes the
operating procedure. It is project-local and model agnostic.

## Meanings

- A commit is a named snapshot of files and its parent history.
- A branch is a movable name for the latest commit in one line of work.
- `origin` is the GitHub remote; `origin/main` is the last observed remote-main
  ref, which must be checked for drift. Local `main` can differ from it.
- A PR compares a source branch with a destination (GitHub's base branch).
  It proposes integration; it does not establish scientific validity.
- A tag is a fixed preservation/release pointer. `archive/wNNN/...` tags preserve
  inputs only; they are neither release tags nor scientific certificates.
- A worktree is another working directory for the same Git repository. It is not
  a separate backup and does not isolate shared refs. Do not switch branches in
  another worker's dirty worktree or overwrite its files.

## Two permanent branches

`main` contains consolidated work. `exploration` contains recorded provisional
research, including negative results and unresolved questions. One task branch
is active at a time in this project. Historical branches remain preservation
inputs during W001; do not delete them or their worktrees automatically.

| Task kind | Start from | Name | PR base |
|---|---|---|---|
| Research experiment or provisional derivation | exploration | research/wNNN-tN-slug | exploration |
| Accepted research promotion | main | promote/wNNN-tN-slug | main |
| Fix, feature, or documentation for consolidated work | main | fix/, feat/, or doc/wNNN-tN-slug | main |

Apply the complete `wNNN-tN-slug` suffix to every prefix. Research implementation
can use feat/ or fix/ when clearer, but its recorded base and PR destination
remain exploration. The PR must state its intended target explicitly.

## Start and finish a task

1. Inspect status/diff, branch/worktree ownership, recorded intake, and relevant
   upstream commits. Preserve unrelated files. Resolve local/remote-main drift
   explicitly; do not assume local main is identical to the reviewed remote base.
2. Create the temporary branch from the recorded intended base. Record dependencies
   on other unmerged work. During W001 only, T2 may stack on T1 before exploration
   exists; record the inherited baseline and do not submit it wholesale to main.
   The PI-authorized T4 prerequisite review remains on the T3 candidate while its
   merge gate is unresolved; it must not open a competing research stream.
3. Work, test, and obtain internal review. Keep scientific status separate from
   software status. Inspect the staged diff and commit with Lu Freitas's requested
   identity. Push the temporary branch as a backup/review candidate.
4. Open a PR to the correct base. Describe the question, actual result, provenance,
   limitations, dependencies, exact checks, and unresolved acceptance decisions.
5. Inspect the final PR diff and remote checks. Re-run affected checks if conflict
   resolution changes the result. No auto-merge merely because CI is green.
6. Merge only within explicit authorization and reviewed acceptance scope. Stop
   work on the temporary branch after integration. Delete branches/worktrees only
   when authorized and after preservation checks.

## Promote exploration to main

Start a promotion branch from the inspected main baseline. Select only accepted
commits or reconstruct a bounded reviewed change when research commits mix
accepted and provisional content. Include required code, tests, data provenance,
and documentation dependencies. Inspect the actual diff; a clean cherry-pick is
not acceptance. Never force a cherry-pick with unsafe conflict resolution or use
an entire exploration merge as a substitute for selection.

An accepted refutation can be promoted. Disapproved calculations remain preserved
with their status and cannot be presented as accepted results. Main acceptance
records may describe rejected work without adopting its scientific conclusions.

After main integration, deliberately synchronize accepted main changes into
exploration, inspect conflicts against provisional work, and verify affected
interfaces. This propagates fixes without promoting exploration back to main.

## Commits, raw data, documents, and releases

- Git stores source, tests, policy, provenance manifests, paper cards, and compact
  reviewed numerical tables. A pushed branch does not prove data are reproducible.
- Heavy arrays and source archives follow Rule 01's external storage/DVC policy;
  retain checksums and versions. Do not migrate or delete originals as cleanup.
- Generated build artifacts live in ignored build areas unless an explicitly
  reviewed deliverable requires preservation. Hashes identify bytes, not validity.
- One canonical labbook source is versioned with the project. Standalone incremental
  task/wave notes must identify their base and eventual canonical destination.
- Release tags `vX.Y.Z`, package version, accepted publication outputs, and build
  receipts are separate from archive checkpoints. Publishing Pages or a release
  requires its own declared gate; do not confuse branch backup with publication.
- PR review and CI do not replace scientific review. Current workflow filters do
  not yet cover every new branch prefix or PR base; W001 T3 must reconcile them.

## Continuous execution and strategy review

A user-requested continuous goal permits sequential tasks after their reading
interval; it does not remove task boundaries, claim review, or publication gates.
Before each task, inspect steering and decide whether its prerequisites still
hold. Ask for PI input when choosing scientific thresholds, changing the physical
scope/observable, abandoning a required deliverable, or resolving an ambiguous
scientific acceptance decision. Do not silently pivot. Independent preservation
or documentation may continue while an answer is pending; dependent work waits.
