---
id: RND-P1-M1.6-W4-PROTOCOL
title: "Wave 4 — Paired solver reliability and evaluation-cost design"
phase: 1
milestone: M1.6
status: draft
contracts:
  - K-012
  - K-013
---

# Wave 4 — Statistical benchmark specification (design only)

**Status:** Unregistered, nonconfirmatory protocol draft. **Owner:** \`kitaev-mf\`, related [Issue #7](https://github.com/lufreitas0000/kitaev-mf/issues/7) and [governance coordination #10](https://github.com/lufreitas0000/project-governance/issues/10). No old record, script, dataset or registered hypothesis is amended by this note. No tests or new benchmarks were executed.

## 1. What the historical script actually measures

Source: [\`statistical_solver_benchmark.py\`](https://github.com/lufreitas0000/kitaev-mf/blob/main/scripts/m1_6/statistical_solver_benchmark.py) and [\`solver.py\`](https://github.com/lufreitas0000/kitaev-mf/blob/main/src/kitaev_mf/services/solver.py), inspected 2026-11-14.

- The header and scorecard export claim **40 trials**, but \`n_trials = 30\` is executed in the loop.
- \`rng = default_rng(20261106)\` lies **outside** the nested method/delta loops. The different solver arms therefore receive successive, *different* random initial states: comparisons are **unpaired**.
- The script uses \`metrics.map_evaluations\` for convergent trials, but appends a **synthetic 60** for every nonconvergent trial, including a Picard run configured with cap 120. The solver reports its actual evaluation count separately; JFNK callback evaluations are included via \`extra_evals\` when reported by the accelerator.
- \`P_same\` is computed as \`same_count/max(1,conv_count)\`. Its natural estimand is **conditional among converged trials**, and is undefined if none converge; a reported zero then has a different meaning from zero matching states.
- The scorecard's \`r90\` is the **maximum sampled delta with an empirical rate at least 0.90**, not the largest statistically certified radius or evidence that every smaller perturbation is reliable. Wilson bounds are printed but not used to choose the reported radius.
- On the tested finite grid, a Wilson interval at each cell does not imply uniform coverage over all methods and delta values; multiplicity, dependencies and extrapolation remain separate design questions.

These are source-level observations, not conclusions about which solver is physically superior.

## 2. Proposed paired design and raw contract

Fix the physical Hamiltonian, gauge convention, reference branch and the exact definition of a successful solve before running. For each selected perturbation scale \(\delta\) and trial index \(i=1,\ldots,n\), generate **one** initial state \(X_0^{(i,\delta)}\) from a documented measure \(\mu_\delta\); run each solver \(m\) on an independently copied state with its own declared stopping and resource caps. Preserve the original random draw, seed, trial/arm IDs and source commit. Do not regenerate a different state inside each method loop.

Define \(C_{m,\delta,i}=\mathbf 1\{\text{converged under the fixed stopping contract}\}\). Define \(B_{m,\delta,i}=\mathbf 1\{\text{matching the predeclared branch/observable criterion}\}\) on **converged** outcomes, and \(J=C B\). Then the distinct quantities are

\[
\widehat P_{\mathrm{conv}}=\frac{1}{n}\sum_i C_i,\qquad
\widehat P_{\mathrm{same}\mid\mathrm{conv}}=
\frac{\sum_i J_i}{\sum_i C_i}\quad(\sum_i C_i>0),\qquad
\widehat P_{\mathrm{same\,and\,conv}}=\frac{1}{n}\sum_i J_i.
\]

Do not assign a numerical value to \(\widehat P_{\mathrm{same}\mid\mathrm{conv}}\) when the denominator vanishes. Paired arm comparison can use \(C_{m,i}-C_{m',i}\) and paired costs rather than treating solver outcomes on the shared draws as independent.

Record on **every** attempt: actual \(\Phi\)-call count from the solver, method-specific callback/Jacobian probes if available, iterations separately, elapsed wall time, outcome/failure reason, residual history (unpadded), branch descriptor and actual energy/observable values. The solver's reported \`map_evaluations\` should be audited against an independent call-counting wrapper before being called an exact oracle. A stopped run is a real recorded *failure* under its cap, not a missing observation whose cost should be replaced by that cap. Report unconditional and success-conditional cost distributions separately.

## 3. Reliability, intervals and radius

A two-sided 95% Wilson interval is one defensible *candidate* pointwise binomial interval only if individual trial outcomes correspond to a justified independent identically distributed sampling protocol. It is **not** a posterior, worst-case guarantee or proof of monotonicity in \(\delta\). For \(\widehat p=1\), its lower bound is \(n/(n+z^2)\); with \(z\approx 1.96\), \(n=30\) gives about \(0.8865\), below 0.90, while \(n=35\) exceeds 0.90. The **selection** of 95%, \(n\), threshold and multiplicity handling are owner-owned scientific decisions, not accepted here.

Before confirmatory benchmarking, specify one of:
(a) **empirical grid radius**, the largest tested \(\delta\) meeting an explicitly stated empirical rate criterion;
(b) **pointwise confidence-filtered grid**, a grid delta meeting an approved lower confidence bound criterion (without making claims for untested deltas);
(c) a **monotonicity-certified radius**, which additionally needs a mathematically or empirically justified monotonicity contract and appropriate simultaneous inference. These are not interchangeable.

## 4. Required owner and independent reviewer gates

**TODO(Lu):** define \(\mu_\delta\), its 11D block weights/metric and support; target physical/mean-field branch and \(B\) criterion; field and \(N_k\) scope; solver arms/caps; estimator and interval confidence/multiplicity; treatment of near-zero Dirac modes; success thresholds; hardware/resource budget; which outputs count as confirmatory. Do not silently retrofit these choices to historical runs.

**Next implementation (separate PR after design acceptance):** generate paired seed identifiers before method loops; record raw evaluation counts for all trials; distinguish nonconvergence/cost; fix the 30/40 header and inconsistent \`P_same\` undefined case; add deterministic fixtures that demonstrate pairing and JFNK callback-count handling. No numerical rerun, source weakening or HC claim transition is authorized by merging this specification.
