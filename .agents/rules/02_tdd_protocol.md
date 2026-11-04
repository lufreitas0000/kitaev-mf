# Rule 02: Test-Driven Development (TDD) Protocol

## 1. Strict Red-Green-Refactor Cycle
1. **Red:** Before implementing any domain entity in `src/kitaev_mf/domain/` or physics algorithm in `src/kitaev_mf/services/`, write a unit test in `test/unit/` encoding the exact analytical invariant or limiting case. Verify that the test fails.
2. **Green:** Implement the minimal pure mathematical code required to satisfy the test contract to machine precision (`atol <= 1e-14` for double-precision algebraic and spectral identities).
3. **Refactor:** Enforce strict static typing (`mypy --strict`) and formatting (`ruff check` and `ruff format --check`) without altering observable behavior.

## 2. Mandatory Physics Verification Benchmarks
- Every matrix representation must pass explicit Hermiticity/skew-symmetry, trace orthogonality, and Lie bracket commutation checks.
- Every momentum-space Hamiltonian builder must reproduce the exact analytical Kitaev dispersion $E(\mathbf{k}) = \pm \frac{\vert{}K\vert{}}{2}\vert{}f(\mathbf{k})\vert{}$ for $c^0$ and flat bands $\pm \frac{\vert{}K\vert{}}{2} U_\gamma^{00}$ for $c^\gamma$ at $J = \Gamma = \Gamma' = \mathbf{h} = 0$, and satisfy supercell Brillouin-zone folding parity to `< 1e-14`.
