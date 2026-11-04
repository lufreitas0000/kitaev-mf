# Milestone Roadmap — kitaev-mf

Branch: develop/m1-foundations-review
Phase: Phase 1 — Analytical Rigor & Domain Foundations
Active batch: M1.4.3b.5 — Explicit SO(4) matrices and type-safe constants
Progress policy: every commit atomically updates this file and PROJECT_STATE.yaml; see .agents/rules/04_milestone_progress.md.
Base commit before this update: a813f6ca886dc77f4647032630f40bea4efd143a

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
| M1.4.3b.0 | Correct spin-coupling map and audited midpoint values | Note 04c | completed |
| M1.4.3b.1 | Specify API; tests for mesh restrictions, chunk invariance, direct small-mesh agreement, convergence order | test/unit/test_k006_quadrature.py | completed; 9 user-reported WSL tests passed |
| M1.4.3b.2 | Implement chunked NumPy midpoint integral and convergence estimator | src/kitaev_mf/services/quadrature.py | completed; focused and full unit tests passed in user WSL |
| M1.4.3b.3 | Add reproducible CLI, mesh sequence, observed order, empirical Richardson estimate, limitations | scripts/k006_quadrature.py; Note 04d | completed |
| M1.4.3b.4 | Focused test, full suite, Ruff, mypy, diff hygiene; user reports WSL results | local validation record | awaiting_local_validation |

Acceptance:
1. Integrate continuous |f(k)|, not the ratio with vanishing denominator.
2. Require positive mesh sizes divisible by three; convergence meshes double exactly.
3. Chunk size changes memory only; result invariant within floating-point summation tolerance.
4. Report observed order from at least three meshes. Richardson estimate is conditional on asymptotic power-law behavior, not a rigorous bound.
5. CLI prints meshes, values, successive differences, observed order, and empirical fine-mesh error.
6. Cross-check U_b^00 = -Ubar and gauge link +1.

M1.4.3b.1 tests, b.2 implementation, and b.3 CLI/provenance are committed and have now been executed successfully in user WSL. Reported results: 9 focused quadrature tests passed, 3 energetics tests passed, 70 tests passed in the full unit suite, and the quadrature CLI completed with values converging to approximately 0.5248657458506 and observed order near three. Mypy on src/ passed. Ruff identified one B905 issue and formatting differences; git diff --check identified one extra blank line at EOF. This commit fixes the B905 issue and EOF whitespace, formats the active-batch source/CLI/tests, and introduces a unified validation target. The formatter check intentionally covers all production source and scripts plus the two tests introduced/modified in this batch; unrelated legacy test formatting is not churned. The full-repository formatter remains available through make format.

After the syntax correction, the next `make verify` run reported 12 failed and 58 passed. All 12 failures shared a `NameError` caused by missing dictionary-comprehension iteration clauses in `_local_operators` in `test_local_majorana_representations.py` and `test_physical_projection.py`. The affected spin and gauge dictionaries are corrected in the current commit. The user subsequently reported that Ruff lint passed and the full-tree formatter found 33 files already formatted, but Mypy reported nine dtype-annotation errors in duplicated `_bilinear` helpers across four test modules. Production constants are now explicit `float64` arrays, and independent Kronecker definitions plus `N=(M+G)/2` are checked in `test_local_majorana_algebra.py`. This change is committed but has not yet been run locally.

### M1.4.3b.5 — Explicit SO(4) matrices and dtype contract
- Replace production `np.kron`-constructed `G_MATRICES` and `M_MATRICES`, and arithmetic-generated `N_MATRICES`, with explicit `np.float64` 4x4 arrays.
- Annotate each public dictionary as `dict[int, NDArray[np.float64]]`; the intended dtype contract is real float64.
- Preserve independent verification: reconstruct G and M from Kronecker-product definitions in the algebra test and verify `N^gamma=(M^gamma+G^gamma)/2` for every component. Keep the explicit N Kronecker identities.
- This improves transparency and type stability, not runtime performance: these small module-level matrices are not a numerical hot path.
- No physics convention, matrix ordering, sign, or tolerance should change.
- Local acceptance: `mypy src test scripts`; focused local-algebra, representation, projection, and exchange tests; then `make verify`.
- Do not mark complete until the user reports exact WSL output. No post-change tests have been executed by the assistant.

## 4. Remaining M1.4 execution order

### Immediate handover — M1.4.3b.5 local validation
1. Read `AGENTS.md`, `.agents/rules/04_milestone_progress.md`, `PROJECT_STATE.yaml`, and this roadmap.
2. Inspect `git status --short --branch`; fetch and fast-forward only `develop/m1-foundations-review` using the documented protocol.
3. Run `mypy src test scripts`. If dtype errors persist, inspect NumPy typing behavior and constant annotations before changing helper signatures or adding casts.
4. Run `pytest -q test/unit/test_local_majorana_algebra.py test/unit/test_local_majorana_representations.py test/unit/test_physical_projection.py test/unit/test_exchange_hamiltonian.py`.
5. Run `make verify`; capture each stage and final Git status. Do not infer a clean tree from test success.
6. Report exact output. Only then update both progress files atomically with the user-reported validation and close the batch if all acceptance criteria pass.


### M1.4.3c — Close Gaussian energetics validation
- Reconcile T04b-03 precision with accepted K-006 output; keep regression tolerance distinct from quadrature precision.
- Run pytest -q test/unit/test_energetics.py.
- Run the unified make verify target: coverage-enabled pytest, K-006 quadrature CLI, full-tree Ruff lint and formatting, Mypy across src/test/scripts, committed-branch and local Git whitespace checks, and a final working-tree report.
- Confirm the Git status output and report any remaining local modifications.
- Record exact outcomes and user confirmation in both progress files.
- Complete M1.4.3 only after unified verification is clean.

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
- K-006 prior 0.5248645 is not a high-precision target. The reproducible midpoint implementation and CLI have now been executed in user WSL; the numerical result is near 0.5248657458506.
- Observed N^-3 behavior is numerical evidence, not a rigorous quadrature error bound.
- M1.4.2 tests do not establish plaquette Wilson-loop/flux behavior; add in M1.4.4.
- M1.4.3 unit tests passed in user WSL; full static/Git verification remains pending.

## 7. Change log
- Corrected the machine-readable M1.4.3b definition of done and recorded completed normalization subtask M1.4.3b.0.
- Corrected the CLI invocation formatting in Note 04d; no scientific content changed.
- M1.4.3b.0 normalization mapping and midpoint table corrected; b.1 tests, b.2 implementation, and b.3 CLI/note committed. Awaiting local validation.
- M1.4.3b.1 contract tests committed; M1.4.3b.2 implementation and M1.4.3b.3 provenance/CLI added. Local WSL validation pending.
- M1.4.3b.1 test contract is being added; no local test execution is claimed.
- Roadmap created; M1.4.3b active; mandatory per-commit progress synchronization introduced.

- Added `make verify` to run tests, static analysis, whitespace checks, and report working-tree state in sequence.
- Replaced Kronecker-generated production SO(4) matrices with explicit float64 arrays; retained independent Kronecker and N=(M+G)/2 checks. Pending user WSL Mypy/test validation.

- Completed the active-batch formatting adjustment in `test_energetics.py`; no numerical or physical expressions were changed.

- Extended `make verify` to execute the reproducible K-006 quadrature CLI between pytest and static checks; this is a numerical diagnostic, not a unit test.

- Expanded the formatter check back to the full `src/`, `test/`, and `scripts/` tree after applying Ruff's reported formatting-only changes to the remaining legacy test modules; no test logic or physics expressions were changed.

- Fixed the missing comma in the `np.allclose` argument list of `test_local_majorana_algebra.py` that caused pytest collection to stop; this is a syntax-only correction. Awaiting a complete rerun of `make verify`.


- Fixed missing dictionary-comprehension iteration clauses in both finite-Fock-space local-operator helpers after the user-reported WSL run returned 12 failed and 58 passed. All 12 failures shared the same `NameError`; the affected algebraic assertions were not reached. Awaiting focused and unified WSL revalidation.
