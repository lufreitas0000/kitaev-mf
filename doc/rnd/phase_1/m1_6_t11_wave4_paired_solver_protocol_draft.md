---
id: RND-P1-M1.6-T11-W4
title: "Wave 4 — Paired nonlinear-solver benchmark: statistical design draft"
phase: 1
milestone: M1.6
task: T11
status: draft
questions:
  - "What precisely is the sampling distribution and event measured by solver convergence and branch recovery?"
  - "How should paired seeds, failures, map-evaluation cost and uncertainty be recorded?"
  - "Which owner decisions are prerequisites to a confirmatory solver-ranking experiment?"
depends_on:
  - M1.6-T6
  - M1.6-T10
contracts:
  - K-011
  - K-012
  - K-013
artifacts:
  - scripts/m1_6/statistical_solver_benchmark.py
  - data/processed/benchmarks/m1_6/statistical_solver_scorecard.dat
  - data/processed/benchmarks/m1_6/pareto_speed_vs_radius.dat
---

# Wave 4 — Paired solver benchmark: W1 specification (draft only)

**Authority and limits:** This note is a *proposed* statistical design, not a pre-registered experiment, a reopened/completed M1.6 milestone, a modified historical hypothesis, or a numerical result. Wave 4 was initiated by Lu on 2026-11-14. Proposed T11 is a documentary identifier, **not an activated milestone/task in `PROJECT_STATE.yaml`**. Do not silently edit frozen history or use this draft to authorize confirmatory runs.

**Source baseline:** [kitaev-mf main `acb3fd1c288cb636ff18b996d7a690a5ca6261a0`](https://github.com/lufreitas0000/kitaev-mf/tree/acb3fd1c288cb636ff18b996d7a690a5ca6261a0). Read with the existing [M1.6 T6](https://github.com/lufreitas0000/kitaev-mf/blob/acb3fd1c288cb636ff18b996d7a690a5ca6261a0/doc/rnd/phase_1/m1_6_t6_statistical_scorecards_and_governance.md), [M1.6 T10](https://github.com/lufreitas0000/kitaev-mf/blob/acb3fd1c288cb636ff18b996d7a690a5ca6261a0/doc/rnd/phase_1/m1_6_t10_statistical_benchmark_evaluation.md), [historical benchmark script](https://github.com/lufreitas0000/kitaev-mf/blob/acb3fd1c288cb636ff18b996d7a690a5ca6261a0/scripts/m1_6/statistical_solver_benchmark.py), [scorecard](https://github.com/lufreitas0000/kitaev-mf/blob/acb3fd1c288cb636ff18b996d7a690a5ca6261a0/data/processed/benchmarks/m1_6/statistical_solver_scorecard.dat), [SolverConfig/Metrics](https://github.com/lufreitas0000/kitaev-mf/blob/acb3fd1c288cb636ff18b996d7a690a5ca6261a0/src/kitaev_mf/domain/solver.py), [accelerator dispatch](https://github.com/lufreitas0000/kitaev-mf/blob/acb3fd1c288cb636ff18b996d7a690a5ca6261a0/src/kitaev_mf/services/solver.py), and [seed generator](https://github.com/lufreitas0000/kitaev-mf/blob/acb3fd1c288cb636ff18b996d7a690a5ca6261a0/src/kitaev_mf/services/seeds.py). [Publication handoff Issue #7](https://github.com/lufreitas0000/kitaev-mf/issues/7) is related but does not authorize its own measurements.

## 1. What the legacy results do and do not establish

The historical script has `n_trials=30`, while its source docstring and exported `.dat` header state 40. The reported Wilson lower endpoint 0.886 when 30/30 converge is consistent with 30, not evidence of 40 runs. The actual scorecard contains five methods over eight amplitudes for one field `h=0.20`, with `16x16` Bloch mesh; do not infer `h=0.40` experiments from the script docstring.

- Its nested loops traverse *method first*, then perturbation magnitude, then trial. A single generator `default_rng(20261106)` is advanced continuously, yielding **different starting perturbations across solver methods**. Differences in reported success fractions are therefore **unpaired empirical comparisons**.
- The solver set is `diis`, `broyden`, `lbroyden`, `picard`, `jfnk`; source dispatch additionally provides `fp` (plain fixed-point), distinct from damped `picard`. An independent protocol should fix both labels and algorithmic hyperparameters.
- Converged runs contribute actual `metrics.map_evaluations`, whereas unsuccessful runs append the synthetic integer 60 to the median list, despite `picard` having a 120-iteration cap and JFNK using extra residual calls. This contaminates expected-cost and median-cost comparisons. The nonconvergent time and estimated FLOP totals are omitted entirely; displayed zero runtime/FLOPs at zero successes is a **missing/undefined conditional statistic**, not zero computational work.
- `P_same` is computed among successes only and uses tolerances on energy and onsite `V`, not a fully justified gauge/symmetry-invariant branch identification over `U,V,lambda`. In the existing table, `jfnk` at `delta=0.20` and `0.50` has `P_same=0.933` and `0.867` respectively, despite all runs converging. Thus even the historical table does not support universal target-branch recovery.
- The script computes `R90 = max{delta in scanned grid : observed P_conv(delta)>=0.90}`; it does **not** perform the logistic fit described in the older M1.6 T6 note, establish monotonicity, control multiplicity or give simultaneous confidence. This `R90` is an empirical scanned-grid diagnostic, not a confidence-certified radius or worst-case geometric basin.
- The two historical "curvature" arms rerun the same DIIS solver with no actual change in update rule; their reported C3 nematic diagnostic is **not** a Hessian-eigenvalue or independent saddle-classification oracle. Do not treat the label `ArmB_curvature_checked` as proof of curvature-check implementation.
- The symmetry-restricted `add_symmetric_perturbation` produces correlated perturbations of `U,V,lambda`; its claim of C3 and A/B preservation needs an actual group-action check, and representability of the entire perturbed state as an admissible Gaussian covariance is not established by this construction. Statistical conclusions are conditional on the *actual* sampling measure.

Retain the legacy tables and claims unchanged as dated exploratory artifacts. These defects are not evidence that any solver universally dominates another.

## 2. Mathematical estimands (proposed)

Fix model and numerical contract `theta=(K,h,cell,BZ_mesh,precision,branch,solver_configs,stopping_rule,resource_cap)`. Let `mu_{theta,delta}^{(E)}` be a **specified probability measure on one declared initial-state population** `E` at amplitude `delta`; a raw integer seed is a reproducibility identifier, not the measure itself. The algebraically valid, explicitly symmetry-restricted, and globally Gaussian-representable populations must be treated as **distinct ensembles**. Their probability measures and results cannot be pooled without a separately justified mixture estimand. For algorithm `a`, denote the completed/capped run `R_a(x0;theta)`.

Define convergence and the **three-way** independent branch decision:

$
C_a(x_0)=1\{\text{approved residuals meet tolerances within the approved budget}\},
\qquad
B_a(x_0)\in\{\mathrm{target},\mathrm{other},\mathrm{unknown}\}.
$

The oracle is applied **only to converged endpoints**. For a nonconverged run record `B=unknown` with a separate `converged=false` field; this is not an identified competing branch. When the oracle is absent, ambiguous, or lacks an independent reference, it must return `unknown`. The indicator of **confirmed target recovery** is `I_T=1{C=1 and B=target}`; it must not interpret unknown as other.

For each fixed stratum and population measure, let `n` be the number of trials, `n_C` the number converged, `n_T` confirmed-target converged, `n_O` confirmed-other converged and `n_U` converged but unclassified, so `n_C=n_T+n_O+n_U`. The identifiable empirical quantities are

$
\widehat P_{\mathrm{conv}}=\frac{n_C}{n},
\quad
\widehat P_{\mathrm{target}}\in\left[\frac{n_T}{n},\frac{n_T+n_U}{n}\right],
\quad
\widehat P_{\mathrm{same}\mid\mathrm{conv}}\in
\left[\frac{n_T}{n_C},\frac{n_T+n_U}{n_C}\right]\quad(n_C>0).
$

The intervals above are **partial-identification bounds**, not statistical confidence intervals. When `n_C=0`, the conditional probability is **undefined**, not zero. The exact product identity `P_target=P_conv P_same|conv` is meaningful for the **underlying well-defined events**; with unknown labels, report the bounds rather than substituting a point estimate. The observed unknown rate `n_U/n_C` is reported explicitly when defined.

The primary **resource** statistic is the actual `map_evaluations` and CPU/wall time recorded for **every run**, including failures, with a **common predeclared cap** or a justified algorithm-specific comparison budget. Report both total cost paid up to the cap and cost conditional on target success. A failure's observed cost is not a fictitious `K_max`. Iteration count is a separate metric and does not represent equal work for JFNK, where a single outer step may require multiple map evaluations. Time-to-success beyond an imposed cap is right-censored; cost *actually spent until stopping* is observed, not censored.

Comparisons between algorithms `a,b` use paired initial states `x0_j` and the joint outcomes `(C_{a,j},C_{b,j})`. Estimate difference `E_mu[C_a-C_b]` and paired cost differences on the same `j`; do not treat two arms sharing an initial condition as independent samples. Confidence intervals require an explicitly defined sampling design and suitable paired/cluster treatment.

## 3. Paired experimental unit and branch admissibility

**Unit:** one saved initial-state realization `(case_id, theta, delta, seed_id, initial_state_hash)`, run *unchanged* through each configured algorithm. Share the initial `U,V,lambda` **bytes/values**, not merely equal RNG seeds inside different methods. Distinct `case_id` values are independent only to the extent justified by the perturbation sampler; share their IDs across algorithm arms, and avoid artificial independence across repeated fields or amplitudes derived from the same direction.

**Population definitions (not approved sampling distributions):** Let $\mathcal X_{\mathrm{alg}}$ be the documented finite real coordinate arrays $(U,V,\lambda)$ satisfying the package's shape, finiteness and onsite skew relations. Let a *separately defined* group representation $\rho:G\to\mathrm{Aut}(\mathcal X_{\mathrm{alg}})$ determine $\mathcal X_{\mathrm{sym}}=\{x\in\mathcal X_{\mathrm{alg}}:\rho(g)x=x\ \forall g\in G\}$; this subset is not certified by naively applying a common $4\times4$ skew matrix across sites. Let $\mathcal X_{\mathrm{Gauss}}$ denote those correlator coordinates for which a **global CAR covariance** exists with consistent Fourier/bond/site marginals and the claimed Gaussian purity/mixed-state constraints. For a real antisymmetric covariance $\Gamma$, spectral admissibility requires $-\mathbb 1\preceq i\Gamma\preceq\mathbb 1$ in the compatible Majorana normalization (pure Gaussian states additionally obey the corresponding purity relation). Checking only selected two-site blocks is not a certificate of membership in $\mathcal X_{\mathrm{Gauss}}$; $\lambda$ are auxiliary multipliers rather than covariance entries.

A future sampling design must declare **one** population and an explicit distribution $\mu^{(E)}_{\theta,\delta}$ on it; a generator that only checks $\mathcal X_{\mathrm{alg}}$ must label its output **algebraic/unverified**, never symmetry-certified or Gaussian-certified. A separately checked invariant ensemble may be analyzed independently; broader algebraic robustness and physically admissible Gaussian robustness answer different questions. No physical population, global covariance construction or measure is chosen by this draft.

Before running a comparison, independently verify (or explicitly delimit) the intended symmetry orbit, physical/variational admissibility of the seed state, and the target branch. Target identification must specify allowable gauge/orientation transformations, boundary and momentum-space conventions, and an independent criterion beyond the solver's own `converged` flag. A stationary solution is not necessarily the global variational minimum or a physical projected spin state.

Treat `fp` (unmixed fixed point) and `picard` (damped Picard) as different algorithms. Compare `diis`, full `broyden`, `lbroyden`, and `jfnk` only with disclosed configurations and genuine implementations. Record any `skew_mode`, damping, restart, rank/history and JFNK tolerance, as well as shared stopping rules and solver-specific exceptions. Ensure equal semantics of convergence before comparing outcomes.

For fields near branch crossings and spinodals, stratify **by physical branch** and field segment; never average over different target branches or call a basin-crossing a numerical failure without a predeclared target criterion.

## 4. Radius and statistical uncertainty

For a chosen, **owner-approved** target probability `p_star`, distinguish:

1. *Observed grid radius* `R^{grid}_{p_star}`: largest tested amplitude with observed probability at least `p_star`; this need not be a contiguous interval and is descriptive only.
2. *Model-dependent interpolated radius*: extracted using an explicitly justified monotone/logistic or other fit, with diagnostics and uncertainty; a regression fit does not certify the functional form.
3. *Confidence-backed grid radius*: requires a predeclared simultaneous lower confidence procedure across tested amplitudes (including multiplicity and any shared directions). A conservative definition requires **all** grid points at or below the reported radius to meet the simultaneous lower-bound criterion. This is statistical coverage *under mu*, never a deterministic geometric worst-case radius.

For each Bernoulli probability, Wilson intervals can be useful *pointwise* when underlying trials are iid Bernoulli at a fixed design point. They do not automatically address paired differences, adaptive amplitude scanning, simultaneous grid claims, shared random directions, model misspecification or a posterior probability of a claim. Specify the confidence interpretation and multiplicity adjustment **before** confirmatory tests.

Do not manufacture a sample count, confidence level, `p_star`, field grid, sigma level, numerical threshold, or failure decision to suit historical success. Put all quantitative choices in an explicit owner-reviewed protocol version. Pilot runs for implementation validation and confirmatory runs must be separated, with independent fresh cases after the design is frozen.

## 5. Proposed traceable output contract (not implemented)

At minimum, each raw trial record should contain `protocol_id`, exact `git_sha` and dirty-state indicator, script/config hash, `case_id`, independent `initial_state_hash`, generator identifier/version and seed/direction, `K`, full field vector and units, lattice and boundary conditions, mesh/precision, algorithm and all hyperparameters, resource cap, stop reason, `converged`, individual residuals, `iterations`, **actual** `map_evaluations`, actual measured wall/CPU time (where available), final-state hash, branch decision (`target`/`other`/`unknown`), branch-oracle provenance, raw unpadded residual history, numerical warnings and solver failures. Make `unknown` and `not measured` distinct from zero.

Separate a durable trial/metadata manifest from computed summary tables. Generate numerical tables only from the raw records under validated, atomic output transactions, and retain enough evidence to reproduce a scorecard without rerunning the physics. Avoid storing enormous raw arrays in Git; follow the repository's DVC/manifests policy. Preserve the legacy `.dat` files rather than in-place rewriting. Pin CPU/GPU hardware, thread count, BLAS backend and environment when reporting time; count map-evaluations consistently and independently test JFNK's inner residual callback accounting.

**Independent tests for a later W2 implementation:** shared `initial_state_hash` across methods; demonstrably different hashes across genuinely independent cases; exact success and failure counts; undefined conditional quantities at denominator zero; real nonconvergence cost; JFNK inner-map count; all required residuals; branch-oracle invariance under allowed symmetries; seeded repeatability; malformed/partial-record rejection; raw-to-summary reconstruction; simultaneous uncertainty coverage via a synthetic Bernoulli simulation with known probability. A numerical pilot must not overwrite historical benchmarks or be presented as confirmatory.

## 6. Gates, owner decisions and work decomposition

**TODO(Lu) — Physical/mathematical:** choose model parameter/field points including whether to sample near the first-order transition and the spinodal; select QSL and polarized reference branches; define separate measures on `X_alg`, `X_sym` or `X_Gauss` and their actual certification routines; define the group action, allowable gauge/branch equivalences and independent three-way target/other/unknown oracle. **Wave 4 approval alone authorizes software infrastructure, not these scientific selections.**

**TODO(Lu) — Computational:** approve which solver variants count as separate algorithms, comparable resource cap (iteration vs actual `Phi` work), residuals and numerical tolerances, mesh and precision, hyperparameter selection/tuning policy, warm-start handling, runtime platform and memory/time budget.

**TODO(Lu) — Statistical:** specify primary estimand(s), `p_star`, confidence/error interpretation, trial count, field/perturbation sampling and dependence, multiplicity plan, unknown-branch policy, target success condition and precommitted ranking/acceptance criteria. The historical 0.90 and 95% choices are *prior claims*, not automatically approvals for this new confirmatory design.

**Proposed execution order:**
1. **W1 (this PR):** inspect source/legacy data; define events, measures and confounders; expose owner decisions. Deliver a draft protocol only.
2. **W4 independent review:** try to falsify the design, especially branch oracle, sampling independence, fairness of JFNK cost, radius interpretation and plausible negative controls. Submit exact objections and required fixes.
3. **Owner freeze:** select/approve actual numerical parameters and statistical decision thresholds. Only then assign/activate appropriate milestone and update `PROJECT_STATE.yaml` with `doc/MILESTONE_ROADMAP.md` atomically if its scope/tasks change.
4. **W2 infrastructure (separately authorized by Lu after W4 review):** implement paired-case identity and hashing, a three-way branch-classification interface, raw trial recording with actual costs, and descriptive/identified-bound analysis with synthetic unit tests. **Do not implement a physical perturbation sampler, branch oracle, numerical runner or confirmatory inference until the corresponding scientific contract is approved.**
5. **Bounded exploratory pilot:** test reproducibility, logs, branch classes and resource accounting; keep its seed cohort separate from confirmatory samples. Repair mechanical bugs without retroactively loosening scientific criteria.
6. **Confirmatory benchmark and W4:** run new approved paired cohort; report paired uncertainty and raw reproducibility; independently review and decide which claims are supported, unresolved or contradicted. Scientific approval and any publication handoff remain separate.

**STOP:** one documentary research note, one draft PR. No modifications to source, tests, state, legacy hypotheses, numerical data, thresholds, canonical governance or site files. Do not merge or initiate a benchmark run without separate authorization.
