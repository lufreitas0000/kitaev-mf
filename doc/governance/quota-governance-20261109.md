# Quota-governance integration candidate — 9 November 2026

Question: can bounded model routing and trustworthy quota observation be integrated without
altering scientific contracts, user customizations or active ledgers?

Prerequisites: release 2026-11.3 and the PI-authorized November 9 correction notes. This is
policy/tooling work; no physics result, hypothesis threshold, formal interface or research
acceptance is changed. Gemini-Task owns design; Codex implements and verifies.

Evidence destination: this note and the external verification handover at
`/home/lucas/scratch/quota-gov-implementation-20261109/HANDOVER.md`.

The snapshot implementation preserves established history on truncation or divergence.
Coverage is exact byte-prefix extension, not token magnitude, tuple ordering or hash rank.
All provider state remains read-only; derived telemetry smoke databases are disposable.

The integration is prepared in an isolated task worktree from the inspected local main.
Original dirty work is retained in its original checkout. Final merge and production
activation remain review gates. No research branch is promoted wholesale.

Verification: packaged 33/33; independent boundary/concurrency probes 9/9; installed
hub regressions 10/10. Syntax, bootstrap and scoped sync checks passed in the required
namespace sandbox. Real smoke used 61 private copied session logs and a disposable ledger.
The hub's full agentctl suite has 38 failures and 8 errors among 144 tests, exactly the
same failed-test set as the untouched local-main baseline: fake HTTP providers cannot
reach loopback inside the required network-disabled namespace. This is not a full-suite
pass. No sandbox weakening or live model call was used.
