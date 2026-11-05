# Rule 07: Discovery and Obstruction Journaling

## 1. Requirement and Authority
Every technical repository in the workspace must maintain a `LOG_DISCOVERY_AND_PIVOTS.md` at its root. Any mathematical proof obstruction, algorithmic bottleneck, numerical instability, or design pivot must be recorded before pivoting definitions or discarding approaches.

## 2. Mandatory 5-Column Schema
All entries must conform to the strict five-column table format:

| Date | Module/Context | Obstruction | Pivot/Resolution | Status |
| :--- | :--- | :--- | :--- | :---: |
| YYYY-MM-DD | Exact path or symbol | Precise description of failure mode or obstruction | Rationale and technical direction taken | `Active` / `Resolved` / `Deferred` |

## 3. Operating Invariants
1. **Never Silently Pivot:** If a theorem, tactic, or solver fails, log the exact error or bottleneck first.
2. **Preserve Dead Ends:** Exploratory failures contain valuable physical and mathematical information preventing future agents from repeating identical traps.
3. **Commit Association:** Include the discovery log update in the same Git commit that implements the pivot or documents the blockage.
