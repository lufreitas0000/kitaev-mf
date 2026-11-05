# Rule 01: Repository Safety, Version Control, and Four-Tier Data Integrity

## 1. Git branches and history

1. Inspect the current branch, status, and diff before making changes. Preserve unrelated working-tree content.
2. Never execute destructive Git operations (including `git reset --hard`, `git push --force`, or `git clean`) without explicit user authorization. The same applies to deleting branches, tags, or files.
3. Do not push routine development changes directly to `main`. Use a dedicated branch with a descriptive prefix:
   - `feat/<slug>` for functionality;
   - `fix/<slug>` for bug fixes;
   - `doc/<slug>` for documentation;
   - `research/<slug>` for isolated research;
   - `develop/<slug>` for a long-lived development integration branch;
   - `release/<version>` for release preparation.
4. Integrate through a pull request to `main`. Describe scope, compatibility impact, and verification performed.
5. A PR is ready to merge only after its diff has been reviewed, required checks have passed, and unresolved review comments and merge conflicts have been addressed. Do not claim checks passed without inspecting their actual results.
6. Do not rewrite shared branch history. Keep commits focused and avoid mixing scientific derivations, package-layout migrations, governance changes, and unrelated cleanup.
7. Ordinary feature and research branches should normally use `main` as their integration baseline. Long-lived branches should be synchronized with `main` at deliberate, documented checkpoints.
8. Do not merge, close, or delete branches/PRs unless authorized by the user or explicitly required by the agreed task.

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
