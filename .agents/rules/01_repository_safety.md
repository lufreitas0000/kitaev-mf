# Rule 01: Repository Safety and Data Integrity

## 1. Git Branch & History Protection

1. Never execute destructive git commands (`git reset --hard`, `git push --force`, `git clean -fd`) without explicit user authorization.
2. Do not commit routine development changes directly to `main`. Use a dedicated branch:
   - `feat/<slug>` for new functionality;
   - `fix/<slug>` for bug fixes;
   - `doc/<slug>` for documentation and governance;
   - `research/<slug>` for isolated research work;
   - `release/<version>` for release preparation.
3. Open a pull request to `main` for integration. The PR must describe scope, compatibility impact, and verification performed.
4. A PR is merge-ready only after required CI checks pass, the diff has been reviewed, and no unresolved review comments or merge conflicts remain. Do not bypass failing checks merely to merge.
5. Prefer squash or merge commits according to repository policy; do not rewrite shared branch history. Integrate `main` into long-lived research branches at documented checkpoints.
6. Keep commits focused and atomic. Do not mix scientific derivations, package-layout migrations, governance changes, and unrelated cleanup in one commit.

## 2. Binary & Simulation Data Safeguards

1. Raw binary simulation outputs (`.h5`, `.npz`, `.pt`, generated `.pdf`, and generated `.png`) must not be committed directly to Git history unless a reviewed, compact artifact is explicitly required as a project deliverable.
2. All simulation outputs written to `data/raw/` or `data/processed/` must use atomic file replacement: write to a temporary file, flush and verify it, then use `os.replace()`.
3. Keep reproducibility metadata with any committed numerical result: code revision, model parameters, lattice/geometry, boundary conditions, precision, and relevant random seeds.
