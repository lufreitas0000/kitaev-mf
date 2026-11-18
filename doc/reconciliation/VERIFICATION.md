# Reconciliation verification record

## 2026-11-16 First-task index
Base `770db3ead9d7834c75c9399efcb8e87a0a3451b6`, W001 input `0b75fdfabe3d09113c1615dac70545d73cf71e65`.
Catalog coverage and historical untracked inputs matched `INVENTORY.json`.

## 2026-11-18 Linear main absorption (`doc/absorb-historical-rnd`)
- Offline backup bundle verified at `/home/lucas/Projects/archive/kitaev-mf_backups/pre-reconciliation-20261118.bundle`.
- Historical W001 audit records, Phase 1–4 R&D notes (with explicit epistemic status headers), `data/external/historical_unverified/` pilot tables, and `CODE_SALVAGE_QUEUE.md` absorbed linearly onto `main`.
- `python3 -B doc/reconciliation/check_index.py`: verifies absorbed canonical paths and exact SHA-256 digests without requiring pruned side-branches or auxiliary worktrees.
