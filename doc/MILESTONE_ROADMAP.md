# Milestone Roadmap — kitaev-mf

Branch: develop/m1-foundations-review
Phase: Phase 1 — Analytical Rigor & Domain Foundations
Active batch: M1.4.3b — Reproducible K-006 quadrature
Progress policy: every commit atomically updates this file and PROJECT_STATE.yaml; see .agents/rules/04_milestone_progress.md.
Last known branch head before this update: 50d8ce87591fa76289a841890314f8d626355829

## 0. Resume protocol
A fresh agent reads AGENTS.md, all .agents/rules/, PROJECT_STATE.yaml, and this roadmap; inspects git status --short --branch and remote head; fetches and fast-forwards only develop/m1-foundations-review using git fetch origin then git merge --ff-only origin/develop/m1-foundations-review; reconciles validation from user reports; and includes both progress files in every commit. Never infer tests passed.

## 1. Global constraints
- Work directly on develop/m1-foundations-review. No PR, main merge, branch deletion, or branch switch unless requested.
- GitHub Actions excludes develop/**; local WSL validation is required.
- Strict Red-Green-Refactor. Define mathematical contracts and failing tests before implementation.
- Keep services pure; isolate CLI/I/O.
- Do not implement production observables/energetics before contracts and independent tests are accepted.
- Distinguish exact algebra, derivation, literature comparison, numerical convergence, automated tests, and user-run validation.
- Keep files focused, preferably below 150–200 lines. No scratch artifacts.

## 2. Completed foundations
| Milestone | Scope | Status | Evidence |
|---|---|---|---|
| M1.1 | Local Majorana algebra, representation, projection, quartic Hamiltonian | completed; user-reported WSL pass | Notes 00a–00e; algebra/projection/exchange tests |
| M1.2 | Generalized decoupling and Bloch matrices | completed | Note 02; test_hamiltonian |
| M1.3 | Project state and invariant provenance | completed | Note 00f; K-006 formula |
| M1.4.1 | Gaussian observable/energetics contracts | completed | Note 04 |
| M1.4.2 | Independent finite-Fock observable tests | completed; user-reported WSL pass | Note 04a; test_observables |
| M1.4.3 | Gaussian exchange and double-counting tests | awaiting local WSL validation | Note 04b; test_energetics |
| M1.4.3a | Pure-Kitaev energy normalization audit | analytical audit; note correction recorded | Note 04c |

M1.4.3 remains incomplete until user reports local validation of test_energetics.py and broader checks.

## 3. Active batch — M1.4.3b: Reproducible K-006 quadrature
Objective: evaluate Ubar = (1/3) integral over [0,1)^2 of |1 + exp(-2 pi i k1) + exp(-2 pi i k2)| with deterministic midpoint quadrature, bounded memory, and convergence diagnostics.

The absolute-value integrand is continuous with conical/non-differentiable behavior at isolated Dirac points. This avoids the 0/0 expression in the alternative normalized integrand. Use midpoint meshes N divisible by three: half-shifted coordinates do not coincide with fractional Dirac coordinates 1/3 or 2/3. This avoids sampling the Dirac points, but does not remove the cusp's effect on convergence.

### Ordered sub-batches
| ID | Small task | Deliverables | Status |
|---|---|---|---|
| M1.4.3b.0 | Correct spin-coupling map and audited midpoint values | Note 04c | in_progress |
| M1.4.3b.1 | Specify API; failing tests for mesh restrictions, chunk invariance, direct small-mesh agreement, convergence order | test_k006_quadrature.py | planned |
| M1.4.3b.2 | Implement chunked NumPy midpoint integral and convergence estimator | services/quadrature.py | planned |
| M1.4.3b.3 | Add reproducible CLI, mesh sequence, observed order, empirical Richardson estimate, limitations | scripts/k006_quadrature.py; Note 04d | planned |
| M1.4.3b.4 | Focused test, full suite, Ruff, mypy, diff hygiene; user reports WSL results | local validation record | awaiting_local_validation |

Acceptance:
1. Integrate continuous |f(k)|, not the ratio with vanishing denominator.
2. Require positive mesh sizes divisible by three; convergence meshes double exactly.
3. Chunk size changes memory only; result invariant within floating-point summation tolerance.
4. Report observed order from at least three meshes. Richardson estimate is conditional on asymptotic power-law behavior, not a rigorous bound.
5. CLI prints meshes, values, successive differences, observed order, and empirical fine-mesh error.
6. Cross-check U_b^00 = -Ubar and gauge link +1.

Immediate next action: M1.4.3b.1 only—add failing tests, inspect red result, and update both progress files in that same commit. No implementation in the red commit.

## 4. Remaining M1.4 execution order

### M1.4.3c — Close Gaussian energetics validation
- Reconcile T04b-03 precision with accepted K-006 output; keep regression tolerance distinct from quadrature precision.
- Run pytest -q test/unit/test_energetics.py.
- Run full unit suite, ruff check ., ruff format --check ., mypy src/, git diff --check.
- Record exact outcomes and user confirmation in both progress files.
- Complete M1.4.3 only after user-reported WSL success.

### M1.4.4 — Production observables and energetics (after M1.4.3)
1. M1.4.4a: contracts/tests for one-/two-point observables, plaquette Wilson loop/flux, exchange energy, fields, constraints, scalar double-counting.
2. M1.4.4b: implement pure observable service against domain types/backend ports.
3. M1.4.4c: implement plaquette geometry and ordered-edge Wilson loop; test orientation, gauge invariance, pure-Kitaev zero flux.
4. M1.4.4d: define energetics contracts/tests, including per-bond/cell/site normalization.
5. M1.4.4e: implement pure Gaussian energetics service without duplicated Wick algebra.
6. M1.4.4f: integrate quadratic spectral energy with accepted band/BZ counting.
7. M1.4.4g: integration tests across observables, correlators, solver, energy functional.
Each sub-batch receives a focused R&D note/section and atomic progress update.

### M1.4.5 — Consolidated verification
- Full unit/integration suites; Ruff lint/format; strict mypy.
- Revalidate K-001–K-006 and normalization.
- Audit interfaces, errors, CPU/CUDA parity, import cycles.
- User WSL validation required before completion.

## 5. M1.5 — Research report
1. M1.5a: freeze conventions and claim inventory.
2. M1.5b: derive local algebra, projection, exchange decoupling, Bloch spectrum, correlators, observables, energetics.
3. M1.5c: numerical provenance, convergence, tolerances, reproducibility commands.
4. M1.5d: verify primary literature and bibliography metadata.
5. M1.5e: independent review, report build, final consistency audit.

## 6. Decisions and blockers
- Spin mapping: matching Kp sum sigma_i sigma_j to -Ks sum S_i S_j with S=sigma/2 requires Ks=-4Kp. Do not invert.
- K-006 prior 0.5248645 is not a high-precision target. Independent midpoint recomputation is near 0.5248657458506; reproducible implementation and empirical error estimate pending.
- Observed N^-3 behavior is numerical evidence, not a rigorous quadrature error bound.
- M1.4.2 tests do not establish plaquette Wilson-loop/flux behavior; add in M1.4.4.
- M1.4.3 tests await local WSL validation.

## 7. Change log
- Roadmap created; M1.4.3b active; mandatory per-commit progress synchronization introduced.
