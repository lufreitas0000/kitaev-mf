# Rule 05: Research Programs, Waves, Sprints, and Tasks

## Planning units

A research program identifies a model family and physical questions. A wave
(`W001`, `W002`, ...) is a bounded deliverable and review gate. A sprint is a
working interval within a wave; it does not imply acceptance. Tasks (`T1`, `T2`,
...) are atomic units within a wave. Execute one atomic task per turn, including
internal revisions and review. Leave a 10-minute reading interval between tasks;
without a continuous-goal instruction, do not advance automatically merely because
the interval has elapsed. An explicit continuous goal authorizes the next task
after that interval and a steering/prerequisite check; new scientific decisions
still require PI input under Rule 10.

New work uses wave IDs, not new numbered phases. Historical phase/milestone IDs,
paths, and citations remain stable aliases. Do not mechanically relabel archives.

## Ownership and compatibility

- `doc/MILESTONE_ROADMAP.md` remains the roadmap entry point; its legacy filename
  is retained for existing links. New sections are waves and future research scopes.
- `PROJECT_STATE.yaml` owns the active wave/tasks and durable blockers. Legacy
  phase/milestone keys remain compatibility aliases during the tool migration.
- New discovery notes use `doc/rnd/waves/wNNN_tN_<slug>.md`, with YAML frontmatter
  and the intake contract in Rule 09. Existing notes keep their locations.
- Subject tags do not dictate exclusive directories. Architecture layers remain
  responsibility based; discovery records remain chronological.

## Completion

A task closes only after its declared evidence is produced, affected checks are
run and reported, and review findings are resolved or explicitly left open.
Scientific rejection/refutation is a valid outcome. Do not change a hypothesis
or threshold to obtain PASS. Human-owned scientific criteria remain authoritative.

A wave closes when its declared acceptance gates are met and unresolved work is
explicitly deferred. Physics waves require relevant numerical verification and
an updated labbook narrative/appendix with a successful build. Documentation-only
waves require documentation checks; do not manufacture numerical acceptance work.

A task's completion, scientific acceptance, PR merge, document publication, and
release are distinct events. Update roadmap and project state together whenever
scope, active tasks, decisions, blockers, or evidence materially change.
