# Rule 01: Repository Safety and Data Integrity

## 1. Git Branch & History Protection
1. Never execute destructive git commands (`git reset --hard`, `git push --force`, `git clean -fd`) without explicit user authorization.
2. Feature development must occur on dedicated branches (`feat/<slug>`, `fix/<slug>`, `doc/<slug>`) and pass `make lint` and `make test` prior to fast-forward merge into `main`.

## 2. Binary & Simulation Data Safeguards
1. Raw binary simulation outputs (`.h5`, `.npz`, `.pt`, `.pdf`, `.png`) must never be committed directly to Git history.
2. All simulation outputs written to `data/raw/` or `data/processed/` must use atomic file replacement: write to `.tmp_<uuid>.h5`, flush, verify checksum, and execute `os.replace()`.
