# 2026-11-14 planning record migration crosswalk — kitaev-mf

**Scope:** documentation-only classification of the `kitaev-mf` portion (§2.1–2.5) of the owner-provided 2026-11-14 `maximalist_future_tasks_backlog.md`. The uploaded original stays in conversation provenance; the existing repository's roadmap and historical hypothesis records are untouched. Original user backlog entries are *candidate descriptions*, not independent evidence of correctness. See the [governance proposal](https://github.com/lufreitas0000/project-governance/pull/13).

| Original source item | Owning category / destination | Decision |
| --- | --- | --- |
| §2.1 existing Issues #6, #7 | GitHub Issues [#6](https://github.com/lufreitas0000/kitaev-mf/issues/6), [#7](https://github.com/lufreitas0000/kitaev-mf/issues/7) | Retain; no duplicate |
| §2.2 P1–P6 mass prefactor, mean-field analytic coefficient | BACKLOG IDEA-001, current claim/pilot evidence in `doc/rnd/`, Issue #6 | Open exploratory derivation; fit remains numerical |
| §2.2 continuum exact PT prefactor and static denominator | BACKLOG IDEA-002; [Issue #12](https://github.com/lufreitas0000/kitaev-mf/issues/12) | Prerequisite sector/resolvent contract |
| §2.2 L0 historical validity axes | Existing W2 pilot, K-019-A reviewed erratum and source hypotheses | Historical evidence review, not new sprint |
| §2.2 L1 statistics and paired solver benchmark | [Issue #7](https://github.com/lufreitas0000/kitaev-mf/issues/7), [draft PR #13](https://github.com/lufreitas0000/kitaev-mf/pull/13) | Independently gated implementation; no invented trial/sample thresholds |
| §2.2 L2 86/5 reduction, manifold, bound and T1–T7 | BACKLOG IDEA-003, IDEA-009; [Issue #8](https://github.com/lufreitas0000/kitaev-mf/issues/8) | No assumed D=5 closure or physics root count |
| §2.2 L3 S0–S5 exact/scMFT vortex constraints | BACKLOG IDEA-010; [Issues #12](https://github.com/lufreitas0000/kitaev-mf/issues/12), [#15](https://github.com/lufreitas0000/kitaev-mf/issues/15); roadmap M2.2 | Owner-controlled R3, no simulator activation |
| §2.2 L4 interval Newton and certified enumeration | BACKLOG IDEA-004 | Deferred pending compact domain and symmetry proof |
| §2.2 L5 physical predictive error budget | [Issue #14](https://github.com/lufreitas0000/kitaev-mf/issues/14), roadmap M2.1 | Blocked on independent ED |
| §2.3 exact PT Track A integration A-001..A-015 | Owning `doc/rnd/exact-perturbation-theory/BACKLOG.md` | Preserve original granular integration backlog, no new Issue per item |
| §2.3 exact PT Track B B-001..B-045 | Owning exact PT backlog and [Issue #12](https://github.com/lufreitas0000/kitaev-mf/issues/12) | Keep its internal mathematical task order and physical parity prerequisites |
| §2.4 Phase 2 M2.1 bulk atlas/flux/SDP and H1/H2/H3 | Existing `doc/MILESTONE_ROADMAP.md` M2.1; physical projection [Issue #14](https://github.com/lufreitas0000/kitaev-mf/issues/14) | Milestone-level, not idea-only; no unauthorized new sprint |
| §2.4 H4 spin-channel and BL-07 non-Gaussianity | BACKLOG IDEA-011, IDEA-005 | Exploratory physical-operator definitions |
| §2.4 M2.2 vortex, M2.3 thermodynamics, M2.4 static correlations | Existing `doc/MILESTONE_ROADMAP.md` | Future milestones; no task ID or numerical run inferred |
| §2.4 BL-03 spiral, BL-04 shear, BL-05 preconditioner | BACKLOG IDEA-008, IDEA-007, IDEA-006 | Deferred, source model-specific gates |
| §2.5 Phase 3 M3.0–M3.5 and H_DYN | Existing `doc/MILESTONE_ROADMAP.md` and `doc/rnd/phase_3_dynamical_spin_correlations/blueprint.md`; BACKLOG IDEA-012 for *alternative approximations only* | Future phase, not completed exact correlator reproduction |
| §2.5 Phase 4 FM/AFM, Phase 5 extended model, Phase 6 3D | Existing `doc/MILESTONE_ROADMAP.md` | Defer until selected; do not create fake active sprint |
| §1 cross-project methodology Waves 0–7 | `project-governance` coordination #10 and owning front Issues | Do not duplicate scientific governance in `kitaev-mf` |

## Review checklist

- [ ] Confirm that each original substantive `kitaev-mf` item has a traceable destination or justified deferment.
- [ ] Check any numerical or literature statements against primary sources before using them as acceptance targets.
- [ ] Confirm the roadmap and current YAML are still unchanged and any new `CURRENT_SPRINT.md` does not activate Phase 2.
- [ ] Independently review the new Rule 05/07 adoption and test `scripts/planning_records.py --root <this repo> check` from the governance repository.
- [ ] Preserve owner-controlled choice of next sprint, parity sector, R3 constraints, thresholds and physical publication approval.

**No code, registered hypothesis, source-roadmap milestone status, active task, benchmark, or acceptance criterion is changed in this migration.**
