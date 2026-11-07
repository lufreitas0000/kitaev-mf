# ADR-002: Consolidated Main, Exploration, and One Labbook

- Status: accepted by Lu Freitas in the 2026-11-07 consolidation discussion.
- Scope: project-local policy; does not modify other projects' governance.
- Decision: two long-lived branches named `main` and `exploration`; model-agnostic
  temporary task branches and explicit PR destinations; waves replace new phases.
- Research PR: exploration -> research/wNNN-tN-slug -> PR to exploration.
- Promotion PR: main -> promote/wNNN-tN-slug -> selected accepted changes -> PR to main.
- Maintenance/policy PR: main -> fix/ or doc/wNNN-tN-slug -> PR to main.
- Existing main is an inherited baseline pending claim review, not blanket approval.
- Old branch refs and untracked user work must be preserved before integration.
- One canonical labbook combines logical main chapters with historical technical
  appendices. Optional incremental task/wave notes are review artifacts.
- Multiple subject tags and independent model/method metadata support search;
  no exclusive subject-directory migration is authorized by this decision.
- Recent historical Phase 4 calculations and explanation are disapproved by the
  PI. Preserve them for audit and independently recover earlier perturbation work.
- References owns source archives/transcripts; project interpretation remains
  distinct. A separate References chat runs the bounded ingestion pilot.

Rules 01, 05, and 09 implement this decision. W001 sequences the transition;
branch integration, publication-source migration, scientific acceptance, and
release remain separate tasks. No server branch protections are claimed changed.
