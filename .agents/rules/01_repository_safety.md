# Rule 01: Repository Safety and Data Integrity

## 1. Git Branch & History Protection
1. Never execute destructive git commands (`git reset --hard`, `git push --force`, `git clean -fd`) without explicit user authorization.
2. Do not commit routine development changes directly to `main`. Use dedicated `feat/<slug>`, `fix/<slug>`, `doc/<slug>`, `research/<slug>`, or `release/<version>` branches.
3. Open a pull request to `main` describing scope, compatibility impact, and verification performed.
4. Merge only after required CI checks pass, the diff is reviewed, and no unresolved review comments or conflicts remain. Do not bypass failing checks.
5. Keep commits focused and atomic. Integrate `main` into long-lived research branches at documented checkpoints.

## 2. Binary & Simulation Data Safeguards
1. Raw binary simulation outputs (`.h5`, `.npz`, `.pt`, `.pdf`, `.png`) must never be committed directly to Git history.
2. All simulation outputs written to `data/raw/` or `data/processed/` must use atomic file replacement: write to `.tmp_<uuid>.h5`, flush, verify checksum, and execute `os.replace()`.
