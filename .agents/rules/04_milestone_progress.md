# Rule 04: Milestone Progress and Commit-Level Handover

Every commit on a milestone/development branch MUST update both PROJECT_STATE.yaml and doc/MILESTONE_ROADMAP.md in the same Git commit. This applies to code, tests, documentation, refactors, and corrective changes.

## Atomicity
1. Include artifact changes and both progress files in one atomic commit; never defer progress synchronization.
2. With GitHub APIs, create all blobs, one tree, one commit, then advance the branch ref using an expected-head check. Do not use separate Contents API writes for a multi-file change.
3. If atomic multi-file commit is unavailable, stop and report the limitation.

## Roadmap content
For every active milestone, record objective, scope, assumptions, dependencies, definition of done, ordered batches/sub-batches with stable IDs, small intermediate tasks, expected files, validation commands, statuses (planned/in_progress/blocked/awaiting_local_validation/completed), risks, unresolved questions, decisions, non-goals, and the exact next action. Preserve completed history and explain execution-order changes.

## Evidence discipline
1. Update both progress files whenever scope, order, implementation, tests, blockers, or validation changes.
2. Mark work completed only when deliverables exist and acceptance criteria are met.
3. Local WSL validation remains awaiting_local_validation until the user reports exact commands and outcomes. Never infer test success from code inspection or GitHub writes.
4. Distinguish analytical proof, literature agreement, numerical convergence, automated tests, and user-run validation.
5. Record the previous known commit SHA in the next roadmap update when available; never guess the SHA of the commit being created.
6. A fresh agent must read AGENTS.md, this rule, PROJECT_STATE.yaml, and doc/MILESTONE_ROADMAP.md, then inspect branch/head and working-tree state before editing.

## Branch-local ownership
Roadmap and project state describe the current milestone branch. Do not synchronize in-progress state to main, create a PR, merge, or switch branches unless explicitly authorized.
