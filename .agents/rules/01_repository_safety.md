# Rule 01: Repository Safety, Version Control, and Four-Tier Data Integrity

## 1. Git branches and history

1. Inspect the current branch, status, and diff before making changes. Preserve unrelated working-tree content.
2. Never execute destructive Git operations (including `git reset --hard`, `git push --force`, or `git clean`) without explicit user authorization. The same applies to deleting branches, tags, or files.
3. Long-lived branches are `main` (consolidated) and `exploration` (research).
4. Temporary branches use `<kind>/wNNN-tN-<slug>`: `research/` for exploration,
   `feat/`, `fix/`, or `doc/` for implementation/maintenance, and `promote/` for
   accepted research. Names must not contain agent or model branding.
5. Research branches start from `exploration` and open PRs with base `exploration`.
   Consolidated maintenance and policy branches start from `main` and target `main`.
   Promotion branches start from `main`, contain only selected accepted changes and
   required dependencies from exploration, and target `main`.
6. For each PR state source base, target branch, writable scope, claim decisions,
   dependencies, and verification. A PR is not scientific acceptance by itself.
7. Never merge exploration wholesale into main. Never merge main into exploration
   without inspecting the effect on provisional work. Synchronize deliberately.
8. Use one active task branch at a time in this project. Temporary branches are
   short lived; old phase branches remain preservation inputs until accounted for.
   Branch/worktree deletion still requires explicit authorization.
9. Do not rewrite shared history or discard user files. Do not push routine work
   directly to main. Merge only after review and required checks; never infer CI.
10. Transition exception: W001 policy work starts from existing local main before
    exploration is created. Record that base in the dated audit. Creating the two
    branch names does not certify inherited main content.

For the step-by-step procedure and Git terms, see Rule 10 and
`doc/governance/VERSION_CONTROL_GUIDE.md`.

## 2. Four-Tier Version Control and Data Provenance Contract

Every repository separates source, compact benchmarks, heavy binary tensors, ephemeral build caches, and published deliverables across four tiers:

| Tier | Storage System | Permitted Contents | Rules & Provenance Contract |
| :--- | :--- | :--- | :--- |
| **Tier 1: Source & Compact Benchmarks** | **Git + GitHub** (`lufreitas0000`) | Code (`src/`, `test/`, `scripts/`), governance (`.agents/`, `AGENTS.md`, `PROJECT_STATE.yaml`, `LOG_DISCOVERY_AND_PIVOTS.md`), `doc/rnd/*.md`, `.tex` sources, `data/manifests/*.yaml`, and compact ASCII tables (`data/processed/benchmarks/*.dat` < 100 KB). | ASCII `.dat`/`.csv` files must start with `#` comment headers recording generator script, `git_sha`, `timestamp`, and model/mesh parameters so `pgfplots` figures are 100% traceable. |
| **Tier 2: Heavy Scientific Data** | **DVC + Google Drive Remote** (`data/raw/`) | Heavy binary arrays (`.h5`, `.npz`, `.pt`), Monte Carlo trajectories, wavefunction checkpoints, and raw reference PDFs. | Never commit raw binary arrays directly to Git. Track them via `.dvc` pointer files and SHA-256 manifests (`data/manifests/`) backed by Google Drive or designated remote storage. |
| **Tier 3: Ephemeral Build Cache** | **Git-Ignored Local (`build/`)** | Compiled `.aux`, `.log`, `.fls`, `cache_tikz/`, `.venv/`, `.pytest_cache/`, `.mypy_cache/`, `.ruff_cache/`, and local scratch PDFs. | Strictly git-ignored. Clean only via targeted `make clean-build` targets; never run unrestricted `git clean -fdx` on a research checkout. |
| **Tier 4: Certified Deliverables** | **GitHub Releases & GitHub Pages** | Compiled `labbook.pdf`, `manuscript.pdf`, `book-build.json` receipt, and the static **R&D Blueprint Portal** (`build/blueprint/`). | Published via CI workflows (`latex.yml`, `blueprint.yml`). Compiled PDFs must be accompanied by a `book-build.json` receipt binding `source_sha` to `pdf_sha256`. |

## 3. Atomic writes, security, and reproducibility

1. When writing important simulation output, avoid exposing partially written files as completed results. Write to a temporary file on the same filesystem, validate it, and atomically replace the destination with `os.replace()`.
2. Preserve sufficient provenance for committed numerical results: source revision, model and coupling parameters, lattice/geometry, boundary conditions, precision, solver settings, and random seeds when applicable.
3. Do not commit credentials, access tokens, private keys, personal data, `.env` files, or unreviewed machine-specific configuration.
4. Review staged changes with `git diff --cached` before committing.
5. Do not use generated output as the sole source of truth when the underlying inputs and reproducible procedure can be preserved.
