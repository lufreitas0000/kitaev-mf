# Project entry index: reconciliation and scientific foundations

Start here before new work. [Plan and scope](PLAN.md), [foundation register](FOUNDATIONS.yaml),
and [pinned source/Git inventory](INVENTORY.json) describe the reconciliation inputs.
The inventory is a dated observation, not a live branch/CI dashboard.

## Where to find information

| Question | Canonical entry | Role / limitation |
|---|---|---|
| What work is active? | [Project state](../../PROJECT_STATE.yaml), [sprint brief](../../CURRENT_SPRINT.md) | YAML task state; readable scope in the brief |
| Why this order? | [Reconciliation plan](PLAN.md), [roadmap](../MILESTONE_ROADMAP.md) | Current bounded program versus later outcome candidates |
| What ideas are deferred? | [Idea backlog](../../BACKLOG.md) | Incubation, not accepted science |
| What changed our understanding? | [Discovery log](../../LOG_DISCOVERY_AND_PIVOTS.md), linked R&D notes | Dated reasoning, failed approaches and counterexamples |
| What assumptions support a claim? | [Foundation register](FOUNDATIONS.yaml) | Question, assumptions, consequences, source assessments and open gate |
| Where is an older result? | [Catalog](INVENTORY.json) | Every tracked file in the pinned main and W001 inputs; filter by path/kind/ref |
| How does source map to planning? | [Existing planning migration](../PLANNING_RECORDS_MIGRATION.md) | Retain its rationale; do not start a competing idea register |
| What is verified? | [Hypothesis contracts](../HYPOTHESIS_CONTRACTS.md), linked evidence and errata | Historical labels require current qualifications; tests cover only tested domains |
| Where is W001? | [Preserved integration PR #4](https://github.com/lufreitas0000/kitaev-mf/pull/4) | Provisional candidate; its files are indexed at the pinned commit |
| What is the current solver-study design? | [PR #16](https://github.com/lufreitas0000/kitaev-mf/pull/16), [PR #17](https://github.com/lufreitas0000/kitaev-mf/pull/17), [PR #20](https://github.com/lufreitas0000/kitaev-mf/pull/20) | Design/infrastructure/sampling proposal; not a physical sampling certification |

## Information contract for future contributions

Use a stable claim/question ID and one canonical authored source. Record model
family, couplings/field, geometry/BC, representation, normalization, assumptions,
method, epistemic kind (proof/finite calculation/approximation/conjecture),
dependencies, consequences, evidence location and scope, review status, and book
or appendix destination. Tags may overlap. Use source locators and exact Git refs
for material on another branch; do not copy a changing file merely to make a link.

These fields apply proportionally when a source is created or reviewed. The full
historical catalog is not a claim that every old note already satisfies them.
Indexing establishes coverage, not source fidelity or physical correctness.

## Current conflicts to resolve

- Main and W001 carry different task/governance/acceptance states; reconcile them
  explicitly, preserving subsequent main work and earlier perturbation sources.
- Main still advertises broad K-019 reconciliation while its erratum refutes a
  supporting inference. Keep the original statement and link the current assessment.
- The hub canonical-branch/focus pointers and local governance checkout are stale
  relative to inspected remote sources. No hub or governance files are edited here.
- Historical Phase 4 (FM analytical/MF-versus-ED) is disapproved by Lu; newer Wave 4
  (paired solver-study methods) is a distinct workstream.
- Local agent/config files remain untracked preservation inputs until reviewed.
- W001 book drafts are not included in the canonical root; PRL migration remains
  pending. An old compiled PDF is not a current consolidated book.

## Repository boundaries

kitaev-mf owns its approximation and scientific interpretation. References owns
original articles/transcripts; kitaev-numerics owns independent finite-spin oracle
work; kitaev-walk-geometry owns its mathematical kernels. Reuse through explicit
convention/interface maps and pinned evidence. Neighboring repositories' completion
labels are not acceptance of a claim in this model.

Validate the pinned catalog and original untracked inputs with
`python3 -B doc/reconciliation/check_index.py` from the repository root.
A moved input is reported unavailable; a changed input is a dated-observation
difference requiring review, not permission to overwrite it.
