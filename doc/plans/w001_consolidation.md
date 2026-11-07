# W001: All-track consolidation and acceptance recovery

## Task T1 intake: policy and preservation inventory

- Question: how do we converge the existing work without promoting disapproved
  results or losing earlier perturbation derivations?
- Prerequisites: user annotations dated 2026-11-07; current source refs and files;
  repository Rules 00-07 and governance architecture; ADR-002.
- Contract: documentation/governance only. Preserve scientific sources, numerical
  arrays, old branch refs, and Phase 3 untracked files. Record acceptance separately
  from historical completion labels. Do not claim an all-orders derivation exists.
- Evidence: dated inventory at `doc/audits/w001_source_inventory.json`, initial
  decisions at `doc/audits/w001_acceptance.yaml`, checked policy cross-references.
- Destination: repository governance and later wave appendix; no new physics
  chapter or manuscript result. Tags: conventions, numerics, perturbation.
- Integration: documentation branch from existing local main; PR target main.
- Checks: parse JSON/YAML, validate local changed-document links, verify inventory
  blobs against recorded refs, check historical untracked hashes, git diff --check,
  project-local governance review; global sync preview only.
- Stop: reviewed policy and inventory commit; no whole-branch merge or deletion.

The local main baseline includes seven commits beyond the observed origin/main.
This task branch preserves that lineage. Do not open a main-targeted PR that
silently includes those inherited commits; T2/T3 must select or reconcile the
baseline before promotion. The policy commit can be isolated after dependency
review. Publishing the task branch is not acceptance or integration to main.

T1 completed: policy and preservation inventory passed focused validation and an
adversarial internal review after revisions. Evidence: `doc/audits/w001_verification.md`.
This does not mark any physics claim accepted or any branch integrated.

## Subsequent atomic tasks

| Task | Scope | Exit gate |
|---|---|---|
| T1 | Policy, source inventory, initial acceptance register | Internally reviewed records; no physics altered |
| T2 | Preserve untracked work; recover exact-perturbation history and source deltas | Explicit preservation checkpoints; original vs later claims distinguished |
| T3 | Create exploration and integrate source branches sequentially | Reconciled tree and records; disapproved material visibly separated |
| T4 | Recover perturbative contracts and review highest-risk foundational claims | Operator/projector/resolvent assumptions and actual implemented orders known |
| T5 | Define canonical labbook include map and incremental-note template | Logical main narrative and wave appendices cross-reference stable claims |
| T6 | Neutral manuscript filename, remove live SciPost labels, PRL format | All build/CI/Pages references updated; newly built draft clearly identified |
| T7 | First selected promotion to main | Bounded accepted claims, dependency checks, reviewed PR |

T2 must establish the safe inputs before T3. T4 does not promise to resolve the
prefactor issue. Later physics tasks require their own registered evidence.
Do not start all tasks in one turn. Allow the PI the requested reading interval.

## Classification and future scope

Review order: conventions/state space; mean-field/numerics; bulk/defects;
exact comparisons/perturbation; thermodynamics/correlations. These are search
facets, not unique bins. Model family, coupling/field domain, method, wave,
epistemic status, and claim dependencies are independent dimensions.

The main book can grow into parts for FM, AFM, extended models, and other lattices,
reusing the accepted foundations and adding explicit bridges. Its provisional
chapter outline is not frozen. Do not generalize FM [111] conclusions by taxonomy.

## References companion

Requested chat: `01a1b7f8-1799-7820-af9d-ed9dbbf3e0ee` in References. Pilot:
at most five necessary sources, immutable versions/hashes, source-faithful
Markdown, project collection, and paper cards. No canonical promotion based on
conversion success alone; no changes to kitaev-mf from that chat.

## Chronogram

Use reviewed sessions, not automatic deadlines: inventory/policy 1-2; branch
transition 1-2; foundational/contested-claim review 4-6; References pilot 1-2;
labbook integration 3-5; first promotions and manuscript selection 1-2. Re-estimate
from measured work. No new bulk physics benchmarks during preservation inventory.
