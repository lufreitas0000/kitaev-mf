# ADR-003: Linear `main` Git Topology, Kitaev 5-Repository Family Boundaries, and Build-Artifact Hygiene

- **Status:** Accepted (`RECON-20261116`, Stage 5)
- **Date:** 2026-11-18
- **Supersedes:** Unmerged `ADR-002` draft preserved in the pre-reconciliation index ([`doc/reconciliation/INVENTORY.json`](../reconciliation/INVENTORY.json))
- **Related:** [`ADR-001`](ADR-001-algebraic-gauge-and-solver-architecture.md), [`doc/reconciliation/PLAN.md`](../reconciliation/PLAN.md), [`doc/audits/w001_acceptance.yaml`](../audits/w001_acceptance.yaml)

## Context

Prior to November 2026 consolidation, `kitaev-mf` accumulated 17+ divergent refs, 5 concurrent Git worktrees, unintegrated Spin-ED and 2-RDM SDP prototypes, tracked binary PDFs, a mutable `latest` GitHub Release, and an untagged `1.0.0` changelog header alongside unresolved mathematical caveats (`K-019-A` second-order Feshbach–Schur erratum and open analytical evaluation `IDEA-001`). Maintaining parallel long-lived feature branches for a solo researcher created cognitive overhead, stale documentation pointers, and blurred boundaries across the five-repository Kitaev research family.

## Decisions

### 1. Minimal Linear Git Topology
- `kitaev-mf` maintains **one permanent branch (`main`)** and **one local worktree** (`/home/lucas/Projects/kitaev-mf`).
- Active changes use at most **one or two short-lived topic branches** (`<type>/<slug>` such as `feat/...`, `fix/...`, `doc/...`), integrated cleanly onto `main` via `--rebase` or `--squash --delete-branch`.
- Zero `archive/*` branches or tags are kept inside the live `kitaev-mf` repository. Historical pre-reconciliation refs were captured in an offline verified bundle outside the workspace (`/home/lucas/Projects/archive/kitaev-mf_backups/pre-reconciliation-20261118.bundle`) and indexed in [`doc/reconciliation/INVENTORY.json`](../reconciliation/INVENTORY.json).

### 2. Kitaev 5-Repository Family Boundaries
Each scientific layer owns a dedicated repository so `kitaev-mf` remains a pure, auditable 4-Majorana Gaussian self-consistent mean-field theory (`scMFT`) codebase:

| Repository | Scientific Scope & Ownership |
|---|---|
| `kitaev-mf` (this repo) | 4-Majorana Gaussian `scMFT` (`domain/` + `services/` + `infra/`), Grassmannian/DIIS/L-Broyden/Padé purification/KKT solvers, real-space vortex and defect clusters, RPA/linear response within `scMFT` |
| `kitaev-numerics` | Exact Spin-ED, Lanczos/KPM, OFTLM, 2-RDM SDP lower bounds, HTE, and `iPEPS`/tensor networks |
| `kitaev-walk-geometry` | Exact Majorana lattice walks, lattice Green's functions/resolvents, and full-$Q$ Feshbach–Schur perturbation theory |
| `ksl-sc-bilayer` | Kitaev spin liquid–superconductor heterostructures and proximity/Kondo bilayer couplings |
| `Julia/KitaevDerivation.jl` | Symbolic Clifford/CAR and Schrieffer–Wolff operator algebra CAS |

Non-`scMFT` prototypes (Spin-ED and 2-RDM SDP) are excluded from `src/kitaev_mf/` and handed off to `kitaev-numerics` (`/home/lucas/Projects/archive/kitaev-numerics-handoff/`). Exploratory pure-Kitaev perturbation material under `doc/rnd/exact-perturbation-theory/` stays strictly isolated from production `src/kitaev_mf/` modules (`AGENTS.md` Section 7).

### 3. Cross-Repository External Benchmark Data Contract
- External comparison datasets never enter `kitaev-mf` via copied solver code or ad-hoc arrays.
- Verified cross-repository benchmarks enter strictly under `data/external/<source_repo>/<dataset>.dat` paired with a sibling `<dataset>.meta.json` sidecar containing: `source_repo`, 40-character lowercase hex `git_sha`, `sha256` matching the `.dat` file, `unit_convention`, `geometry`, and `parameters`.
- Both CI (`.github/workflows/ci.yml`) and local tooling enforce this contract via `scripts/sync_external_data.py --verify`.
- Unverified historical comparison tables predating this contract remain quarantined in `data/external/historical_unverified/` and cannot be cited as live cross-repo validation.

### 4. Build-Artifact, SemVer, and Release Hygiene
- Compiled PDFs (`doc/labbook/labbook_kitaev_mf.pdf`, `doc/manuscript/manuscript_kitaev_mf.pdf`), build trees (`build/`), and `book-build.json` are excluded from Git tracking via `.gitignore` — Git stores deterministic `.tex` sources and `.dat` tables only.
- One consolidated workflow (`.github/workflows/docs.yml`) compiles both PDFs and the static Blueprint portal (`scripts/build_blueprint.py`), binds SHA-256 checksums into `book-build.json`, deploys Pages atomically on `main` pushes, and publishes GitHub Release PDFs **strictly on immutable `v*.*.*` SemVer tags** (no mutable `latest` tag or release).
- Package version remains `0.1.0` in `pyproject.toml` with pre-release changes tracked under `## [Unreleased]` in `CHANGELOG.md` until formal v1.0 milestone certification.

## Consequences
- Every pull request to `main` is linearly readable, verified by `CI`, and (when touching documentation) build-checked by `docs.yml`.
- Cross-repo physical comparisons (Gaussian `scMFT` vs. Spin-ED, 2-RDM SDP, and exact walk resolvents) gain cryptographic commit + SHA-256 provenance without coupling solvers across repositories.
- No binary PDF merge conflicts or mutable `latest` release drift can re-enter `kitaev-mf`.
