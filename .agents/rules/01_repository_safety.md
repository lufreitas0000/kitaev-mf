# Rule 01: Repository Safety and Data Integrity

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

## 2. Simulation data and generated artifacts

1. Treat raw simulation outputs and generated files (including HDF5, NPZ, tensor checkpoints, plots, and compiled PDFs) as generated artifacts, not source, by default. Commit them only when a reviewed, compact artifact is an explicit project deliverable.
2. Keep raw and processed simulation data under the repository's designated data locations or an explicitly configured external storage location. Do not assume that every machine has the same local data.
3. When writing important simulation output, avoid exposing partially written files as completed results. Prefer writing to a temporary file on the same filesystem, validating it, and atomically replacing the destination with `os.replace()`. Apply this where the file format and workflow support it.
4. Preserve sufficient provenance for committed numerical results: source revision, model and coupling parameters, lattice/geometry, boundary conditions, precision, solver settings, and random seeds when applicable.
5. Do not commit credentials, access tokens, private keys, personal data, or unreviewed machine-specific configuration.

## 3. Review and reproducibility

1. Review staged changes with `git diff --cached` before committing.
2. Do not use generated output as the sole source of truth when the underlying inputs and reproducible procedure can be preserved.
3. Keep large or transient outputs out of Git history; use documented external storage or release artifacts when persistence is required.
