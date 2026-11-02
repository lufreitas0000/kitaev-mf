# Numerical implementation boundary

No solver module is implemented yet. The first implementation should be a
small NumPy/SciPy fixed-gauge solver, independent of the existing mean-field
backend.

## Proposed modules

- `geometry.py`: periodic honeycomb cell indexing and bond translations.
- `gauge.py`: link configurations, plaquette fluxes, Wilson-loop/boundary
  twist bookkeeping, and gauge-string constructors.
- `majorana.py`: construction of the bipartite matrix \(T\) and
  antisymmetric matrix \(A\).
- `spectrum.py`: singular values, Gaussian vacuum energy, and Fock-space
  spectral bounds.
- `parity.py`: physical parity constraint, implemented only after its
  finite-size convention is derived.
- `spin_ed.py`: small-cluster spin-space reference calculation.

The first solver contract is
\[
u^\alpha_{\mathbf r}\in\{\pm1\},\quad K>0
\quad\Longrightarrow\quad
s_n(T)\geq0,\quad
E_{\mathrm{vac}}=-\frac12\sum_n s_n(T).
\]

Geometry and gauge inputs should be immutable data. Numerical routines
should be pure functions; file output and plotting belong in separate
scripts.

## Tests before implementation

1. \(A\) is real antisymmetric.
2. Positive eigenvalues of \(iA\) equal singular values of \(T\).
3. Flux-free finite-torus spectrum agrees with the Bloch dispersion.
4. One link flip changes precisely its two adjacent plaquette fluxes.
5. Gauge-equivalent link configurations have identical singular values.
6. Physical parity and sector energies agree with spin-space ED on small
   clusters.

Keep tests in the existing `test/unit/` and `test/integration/` layout.
Do not create a second test framework under this research directory.
