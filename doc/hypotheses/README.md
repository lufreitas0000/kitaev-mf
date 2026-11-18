# Living Hypothesis Contracts (`doc/hypotheses/`)

This directory holds machine-verifiable, immutable falsifiable hypothesis contracts (`HC-XX.yaml`) governed by `project-governance` Rule 08 and protected against silent rewriting by `governance/profile.toml` (`protected_paths = ["doc/hypotheses/HC-*.yaml"]`).

## Distinction from Phase 1 Summary Registry

- [`doc/HYPOTHESIS_CONTRACTS.md`](../HYPOTHESIS_CONTRACTS.md) and Lab Book Chapter 11 ([`ch_m1_11_hypothesis_contracts.tex`](../labbook/chapters/ch_m1_11_hypothesis_contracts.tex)) summarize the historical Phase 1 contracts (`H1`–`H9`, `HC-P1-01`–`03`) authored prior to PG-10 `HC-*.yaml` adoption.
- All new quantitative scaling claims, thermodynamic extrapolations, and cross-repository falsification checks in Phase 2+ (including closure of `IDEA-001` analytical cubic prefactor evaluation and `K-019-A` full-$Q$ resolvent checks) MUST be pre-registered here as `HC-XX.yaml` files before promotional benchmarks are run.

## Rules for `HC-XX.yaml` Contracts

1. **Immutable Claim Core:** Once registered, `statement`, `assumptions`, `support_if`, and `refute_if` thresholds are human-owned and immutable.
2. **Append-Only Evidence History:** New numerical runs, bootstrap confidence intervals ($\mathrm{BC}_a$ $95\%$), or Bayesian model comparison metrics ($\Delta\mathrm{BIC}$) append to `history` / `evidence` without altering the original hypothesis statement.
3. **Refutations Are Discoveries:** Any refuted contract is marked `status: refuted` and logged in [`LOG_DISCOVERY_AND_PIVOTS.md`](../../LOG_DISCOVERY_AND_PIVOTS.md) with an explicit pivot or descendant hypothesis ID.
