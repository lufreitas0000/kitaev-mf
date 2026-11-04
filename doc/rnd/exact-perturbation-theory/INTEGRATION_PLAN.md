# Integration plan — perturbation research branch

## Objective

Bring the research branch into alignment with the current repository conventions on `main` without importing unrelated production history or weakening the independent exact-perturbation research track.

This is a selective shared-file integration, not a branch merge. The research branch currently diverges from `main`; shared files must be reviewed individually because some are normative, some encode active production architecture, and some contain build-system or workflow assumptions.

## Current audited baseline

- `main`: `bd4f177afa4c130fe71b3b69d71548de49991d8a`.
- `research/pure-kitaev-exact-perturbation-20261102`: 23 commits ahead and 19 behind `main`.
- Merge base: `897a28335386e0d20a9b955e2f8c80bb610aea4e`.
- The branch-local `AGENTS.md` predates the current modular governance and describes obsolete source paths (`src/domain`, `src/services`, `src/infra`).
- `main` now has `.agents/rules/00_workspace_hygiene.md` through `03_editorial_rigor.md`; none currently exists on the research branch.
- Root `README.md` is byte-identical on both branches, but its architecture paths are stale relative to the current `src/kitaev_mf/{domain,services,infra,app}` architecture. Copying it from `main` would not fix that.
- `PROJECT_STATE.yaml` on the research branch has stale M1.2 status and a research-track entry; `main` correctly records M1.2 as completed but has no research-track entry.
- `pyproject.toml` differs materially: main uses setuptools and main's Ruff policy; research uses flit, optional development dependencies, and a different Ruff ignore list. Do not overwrite it without checking the lock/install workflow.
- `.gitignore` differs: main includes ephemeral snapshot exclusions and additional reference-build ignores.
- Main has a Makefile and CI workflow absent or materially different on the research branch. Their commands assume requirements lockfiles and scripts; verify those inputs before adopting them.

## Integration principles

1. Keep the exact-perturbation work under `doc/rnd/exact-perturbation-theory/`; do not move or rewrite its derivations as part of governance synchronization.
2. Adopt main's modular rules as the repository-wide baseline, but retain a short research-specific section in root `AGENTS.md` covering analytical-first development, physical projection, branch isolation, and generated-artifact policy.
3. Treat `PROJECT_STATE.yaml` as shared project state: update the production milestone from current main, then preserve the research track as a separate named entry. Never replace one with the other.
4. Treat build and CI files as executable contracts. Reconcile dependencies, package discovery, lockfiles, scripts, and test paths before copying.
5. Use small, separately reviewable commits. No merge, rebase, force-push, or production-source integration is part of this plan.
6. Do not mark a task complete based only on a successful file write; run the documented checks on the branch.

## Shared-file integration matrix

| File or group | Main-branch state | Planned action | Risk / acceptance condition |
|---|---|---|---|
| `AGENTS.md` | Current modular governance and current architecture, but no perturbation-track detail | Replace branch version with a reconciled version based on main; add a concise research-track subsection | Must use actual `src/kitaev_mf/` paths and preserve research isolation / no-zero-shot-coding rule |
| `.agents/rules/*.md` | Four modular rules on main; absent on research branch | Bring all four over unchanged initially, then review rule 03 compatibility with existing R&D notes | Ensure claim-ID and bibliography requirements are either adopted with a migration plan or explicitly scoped; do not silently waive them |
| `README.md` | Identical on both branches and stale | Rewrite architecture and workflow sections against the actual current repository; add link to exact-perturbation R&D track | Do not copy stale paths; keep project-level README concise and make no claim that research code is production-ready |
| `PROJECT_STATE.yaml` | Main has M1.2 completed and verified K-001–K-005; research has stale M1.2 and research-track metadata | Three-way reconcile: current main production state + research-track block + actual branch/CI status | Preserve K-001–K-005 exactly unless independently re-audited; do not claim clean tree or passing CI without checking |
| `.gitignore` | Main has extra snapshot and reference-build exclusions | Merge missing relevant patterns into research version | Preserve research build path exclusion; do not ignore source notes or compact reviewed results |
| `pyproject.toml` | setuptools, Ruff ignores N802/N803/N806, strict mypy, pytest config | Reconcile with branch's Flit metadata and optional dev extras; decide one build backend and one dependency-install contract | Verify package layout, editable install, requirements files, and CI before changing backend |
| `Makefile`, `requirements*.in/txt` | Main has Makefile / lock workflow | Adopt only after all referenced files and scripts are present and commands pass | No dangling make targets or CI steps; keep WSL/Linux-compatible commands |
| `.github/workflows/ci.yml` | Main has lockfile + Ruff + mypy + pytest CI | Port workflow after dependency and path reconciliation | Must pass on research branch; ensure research branch is included in push triggers or PR validation |
| Other root configs / scripts | Potentially changed since merge base | Inventory from full tree diff before touching | Classify as shared infrastructure, production-only, or research-only; no bulk copy |

## Execution phases

### I0 — Freeze and inventory

- Record both branch HEADs and the merge base.
- Confirm no open local changes are at risk; preserve the research branch as the working target.
- Inventory all root-level and `.github/`, `.agents/`, `scripts/` paths on both refs.
- Classify each differing path: normative shared, build/CI shared, production-specific, research-specific, or generated.
- Produce a file-level decision table before editing.

**Exit:** complete inventory and no ambiguous file classified for blind copying.

### I1 — Governance foundation

- Port `.agents/rules/00_workspace_hygiene.md`–`03_editorial_rigor.md`.
- Reconcile root `AGENTS.md` with current main architecture and add research-specific constraints.
- Audit existing perturbation notes against the adopted editorial rule: claim IDs, self-contained derivations, notation, and verified references.
- Add a migration checklist for any pre-existing note that does not yet satisfy the rules.

**Exit:** one consistent contributor protocol; no contradictory directory or architecture instructions.

### I2 — Project documentation and state

- Rewrite root `README.md` to describe the actual package structure and supported workflows.
- Add a clearly scoped link to the perturbation R&D README.
- Reconcile `PROJECT_STATE.yaml` using main as the production-state source and research branch as the source for its research-track entry.
- Keep production milestones and research milestones separate.

**Exit:** README paths resolve; state accurately describes M1.2 and the perturbation track.

### I3 — Build, lint, and CI reconciliation

- Compare package metadata and dependency sources.
- Select one build backend only after checking packaging and editable-install behavior.
- Reconcile Ruff naming exceptions with mathematical notation policy.
- Ensure Makefile targets, requirements lockfiles, scripts, and CI commands refer to existing files.
- Port CI and test on the research branch.

**Exit:** clean editable install; Ruff, formatting, mypy, and pytest pass; CI workflow runs on the research branch.

### I4 — Final shared-file audit

- Re-run main-versus-research comparison.
- Inspect every remaining shared-file diff.
- Verify no production source was imported unintentionally.
- Confirm research notes, roadmap, backlog, and result policy remain intact.
- Commit the integration as small logical commits and report exact SHAs.

**Exit:** shared governance is current, research-specific state survives, and remaining divergence is intentional/documented.

## Suggested commit sequence

1. `docs(governance): port modular agent rules`
2. `docs(governance): reconcile research branch agent protocol`
3. `docs(project): reconcile README and project state`
4. `build: align package metadata and development dependencies`
5. `ci: align research branch validation workflow`
6. `chore: reconcile shared configuration and verify branch`

Do not combine all six into one large commit. If I3 reveals that main's build assumptions are not compatible, document the decision and keep the branch's working build setup rather than forcing parity.

## Explicit non-goals

- No wholesale merge or rebase from `main`.
- No import of unrelated mean-field production features into the perturbation solver.
- No change to the perturbation Hamiltonian, analytical derivations, or numerical conventions during governance synchronization.
- No claim that the research branch is production-ready.
