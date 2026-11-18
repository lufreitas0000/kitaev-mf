# PI-supplied W001 bundle: preserved source, not an accepted derivation

The manifest binds the original ZIP and its members to SHA-256 hashes. Text
members are preserved byte for byte with `.txt` suffixes so they cannot silently
replace live modules. Bundled Python bytecode is recorded but not imported.
The ZIP remains in the user's Downloads folder; no heavy archive was committed.

The braiding proposal is adapted as a **separate fixed-kappa chiral-matter**
service, with validation and isolation/overlap/closure guards. It does not replace
the four-flavor finite-field mean-field API or certify its old braiding claim.

The numerical perturbation prototype is adapted at
`../../code/w001_gaussian_quadrature.py`. Run from repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  doc/rnd/exact-perturbation-theory/code/w001_gaussian_quadrature.py \
  --length 4 --second-order 24 --third-order 10
```

It prints explicitly unaccepted **unprojected** candidate results and does not
import `kitaev-walk-geometry` or write to a machine-specific output location.
The original all-orders note is preserved unchanged; its anti-Hermitian local
constraint, unproven physical projection, exactness, stability, and all-orders
claims are not adopted. See `doc/audits/w001_bundle_review.md` from repository root.
