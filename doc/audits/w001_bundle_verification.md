# W001-T3 bundle integration verification — 2026-11-08 UTC

- Read archive without executing bytecode; preserved text members byte for byte
  and recorded archive/member hashes in the import manifest.
- Inspected GitHub run 42423051113/job 128555418191 logs at da17ffe: missing
  matplotlib during two test imports; lint/type/lock stages had passed.
- Python 3.12 `.venv/bin/pip-compile --resolver=backtracking --quiet -o
  requirements.txt requirements.in` and equivalent dev-lock command completed.
  Existing locked versions retained; only missing dependency closure added.
- `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/pytest -q --tb=short
  test/unit/test_phase2_chiral_vortex_transport.py`: 5 passed in 3.76 seconds.
- Equivalent focused invocation of `test/unit/test_w001_gaussian_quadrature.py`:
  2 passed in 0.20 seconds, including independent 4x4 Fock exponentials.
- `.venv/bin/ruff check src test scripts
  doc/rnd/exact-perturbation-theory/code/w001_gaussian_quadrature.py`: PASS.
- `.venv/bin/ruff format --check` over the same scope: 173 files formatted.
- `.venv/bin/mypy` over the same scope: PASS, 173 source files.
- Isolated candidate CLI at L=4/orders(24,10) completed; reported total
  13.5611188552 remains an unaccepted unprojected numerical candidate.
- Full required suite with coverage in ignored `build/w001/bundle-check`,
  PYTHONPATH bound to that copy and single-threaded BLAS: **354 passed in
  126.77 seconds**, 96% coverage. This run contains the five chiral tests; the
  two subsequently added Fock/guard tests passed separately as recorded above.
  The final remote suite should contain 356 required tests. No current-head
  remote result is asserted by this pre-push receipt.
- An initial suite invocation was interrupted after 83 passes because its cwd
  was the primary checkout rather than the disposable copy. It had not reached
  the figure-writing tests; tracked benchmark sources were unchanged. The
  complete run above used the correct isolated cwd.
- Explicit opt-in research check: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
  .venv/bin/pytest -q --tb=short
  test/research_checks/check_legacy_mean_field_braiding.py`: **1 failed** in
  1.40 seconds, enclosing phase 0.5904, distance to pi 2.55114 > 0.50. This is
  the preserved unsupported guarantee, deliberately retired from required CI
  by explicit PI decision. No xfail or tolerance relaxation.
- Terra 6.1 adversarial referee: PASS for bounded revised integration; reviewer
  did not run tests. This is not physical coefficient or many-body braid approval.
- `git diff --check`, manifest hashes, YAML parsing, preservation of the three
  original assertion ASTs, and staged-scope inspection checked before commit.

Remote verification after push will be appended with the exact source SHA and
run outcomes. No main merge, Pages publication, or release occurred.
