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

T1 baseline observation (operationally superseded by T3): local main included
seven commits beyond the observed origin/main.
This task branch preserves that lineage. Do not open a main-targeted PR that
silently includes those inherited commits; T2/T3 must select or reconcile the
baseline before promotion. The policy commit can be isolated after dependency
review. Publishing the task branch is not acceptance or integration to main.

T1 completed: policy and preservation inventory passed focused validation and an
adversarial internal review after revisions. Evidence: `doc/audits/w001_verification.md`.
This does not mark any physics claim accepted or any branch integrated.

## Task T2 intake: preservation, provenance, and version-control procedure

- Question: where is the earlier perturbation framework, and how will it actually
  constrain subsequent research and branch integration?
- PI correction: earlier work was ignored; alteration/corruption is not alleged
  or established. Inspect differences without assuming damaged derivations.
- Prerequisites: T1 inventory, original research ref/checkpoints, untouched Phase 3
  preservation inputs, and the latest user correction.
- Contract: preserve original tracked and untracked inputs without changing other
  worktrees; document chronology, content identity/differences, and missing use of
  earlier assumptions. No new scientific correctness or all-orders completion claim.
- Evidence: immutable archive tags, untracked snapshot commit, note blob comparison,
  `doc/audits/w001_perturbation_recovery.md`, and the version-control procedure.
- Destination: governance guide and chronological consolidation appendix. Tags:
  perturbation, conventions, physical-state-space; no publication result.
- Verification: compare snapshot tree to source commits and five recorded hashes;
  parse registers/state; inspect changed-file scope and git diff --check; internal review.
- Branch: temporary documentation branch stacked on T1 during the documented W001
  transition; no direct main push or main-targeted PR containing inherited commits.
- Stop: preservation/provenance and procedure commit; T3 owns branch integration.

T2 completed: seven original notes were identical across the compared tips;
six preservation tags were pushed and their remote targets inspected; Rule 10
and the version-control guide passed internal review. No physics was accepted.
Evidence: `doc/audits/w001_t2_verification.md`.

## Subsequent atomic tasks

| Task | Scope | Exit gate |
|---|---|---|
| T1 | Policy, source inventory, initial acceptance register | Internally reviewed records; no physics altered |
| T2 | Preserve untracked work; recover exact-perturbation history and previously ignored contracts | Explicit preservation checkpoints; original vs later claims distinguished |
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

## Continuous goal authorization

On 2026-11-07 Lu Freitas requested continuous execution through this chronogram,
subject to strategy/PI decision gates. Continue across atomic task turns after the
10-minute reading interval and steering check. Do not require a new permission
question for already authorized preservation, review, or reversible integration
preparation. Surface scientific decisions or material scope changes explicitly;
do not silently change criteria or publish unreviewed results. Keep the goal
active until its actual deliverables are complete or a user-directed pause applies.

## T3 integration intake and baseline

Intake, source accounting, conflict decisions, and verification are recorded in
`doc/audits/w001_integration.md`. The task integrates sources into exploration
only; main stays at its reviewed remote baseline. Original local-main history is
preserved under an archive branch/tag. Candidate workflows support exploration
and temporary branch PRs while requiring explicit manual publication from main.

T3 verification gate: original braiding test fails unchanged on its source tip.
PI chose to keep the claim unresolved and the integration PR in draft while
foundation review continues. T3 remains in progress for final merge; T4 may run
on the unified candidate first. No test threshold/xfail change is authorized by
this scheduling decision.

## T4 outcome

Completed the bounded foundation audit and reviewed conditional finite algebra.
Primary Kitaev source and finite diagnostics show why later completion claims
cannot certify full-Q resolvents, exact bulk prefactors, or truncation bounds.
The earlier notes remain unchanged and operative. Evidence and acceptance scope:
`doc/audits/w001_foundation_review.md`, `doc/audits/w001_acceptance.yaml`, and
`doc/audits/w001_t4_verification.md`. T3 merge stays draft-gated under PI steering;
T5 can integrate accepted statements and visibly qualified historical material.
