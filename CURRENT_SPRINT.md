# Current sprint — research and execution brief

**Owning milestone: M1.11** (completed Phase 1 closure in the current `PROJECT_STATE.yaml`).  
**Active execution sprint on `main`: none selected.** This is a truthful transition/handoff brief, not an authorization to begin Phase 2, reactivate M1.11 tasks, or change priority. The authoritative machine status is [PROJECT_STATE.yaml](PROJECT_STATE.yaml); accepted historical milestone summaries and the proposed next phases are in [the roadmap](doc/MILESTONE_ROADMAP.md).

## Baseline and current choice

The project state lists M1.11 as `completed` with no active tasks; the roadmap labels Phase 2 *active preparation*. A preparatory research question, open Issue, or draft PR must not be misreported as an approved sprint.

The natural next **candidate** is Phase 2 / M2.1 (bulk physical observables and predictive comparisons) under the same isotropic ferromagnetic Pauli-spin Kitaev model in a [111] field. **TODO(Lu):** confirm the selected milestone, its bounded first task and the exact experimental/physical acceptance criteria. The choice of a roadmap candidate here changes no status.

## Prerequisites and non-goals

- Local Majorana representation and projector sign inconsistencies: [K-019-A Issue #10](https://github.com/lufreitas0000/kitaev-mf/issues/10). The merged erratum on second-order Feshbach inference does *not* resolve every physical convention.
- Exact finite-torus gauge/parity, flux-sector thresholds and resolvent: [Issue #12](https://github.com/lufreitas0000/kitaev-mf/issues/12).
- Independent physical-spin finite-cluster references: [validation Issue #14](https://github.com/lufreitas0000/kitaev-mf/issues/14) and [kitaev-numerics Issue #1](https://github.com/lufreitas0000/kitaev-numerics/issues/1). Its competing Sprint 1A candidates require independent approval.
- Statistical experimental design: existing [Issue #7](https://github.com/lufreitas0000/kitaev-mf/issues/7) and [merged W1 design PR #16](https://github.com/lufreitas0000/kitaev-mf/pull/16) (documentation scope only; the superseded [draft PR #13](https://github.com/lufreitas0000/kitaev-mf/pull/13) was closed without merge); fit estimates and Wilson intervals are not analytical/physical guarantees.
- Geometry/KKT-reduction derivation: existing [Issue #8](https://github.com/lufreitas0000/kitaev-mf/issues/8).
- Vortex-seeded finite-field computations remain separately owner-gated: [Issue #15](https://github.com/lufreitas0000/kitaev-mf/issues/15).

**Non-goals:** no Phase 2 solver campaign, no new success tolerance, no registered-HC change, no numerical value assumed exact, no user-selected roadmap priority inferred from document creation.

## Proposed first bounded work, only after milestone selection

Prepare a reviewable theory specification of one spin-Hilbert observable and its Gaussian-MF counterpart on the *same* finite torus. Fix the operator ordering, Pauli/field units, BC/twist, physical projection, variational meaning and independent ED oracle before interpreting discrepancies. Link any future T1/T2 tasks in `PROJECT_STATE.yaml` instead of copying statuses here.

**Validation currently claimed:** only a planning-document source inspection, not a new physics or solver run. Normal software tests, independent W4 scientific review, explicit scope approval and Rule 06 atomic state/roadmap updates remain conditions of a subsequent implementation.

## Next decision / stop conditions

Ask Lu to select the next milestone and observable. Stop before activating physical calculations whenever the projector/sector definitions or statistical success criteria remain unresolved. The [idea backlog](BACKLOG.md) contains longer-horizon candidates; the [research journal](LOG_DISCOVERY_AND_PIVOTS.md) records obstructions and rationale without replacing registered evidence.
