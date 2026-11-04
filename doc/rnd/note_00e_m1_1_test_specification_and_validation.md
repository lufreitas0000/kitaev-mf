# Note 00e: M1.1 Test Specification and Validation Evidence

## 1. Scope and validation policy

This note consolidates the finite-dimensional test contracts for M1.1. It maps the proposed A01--A20 inventory to the implemented tests in Notes 00a--00d and records what the tests establish.

The tests use explicit finite-dimensional Jordan--Wigner matrix representations of the Majorana algebra. These are independent computational realizations used to check the identities; they are not substitutes for the analytical derivations in the preceding notes. Floating-point residuals establish numerical consistency at the stated tolerance, not formal proofs.

Unless a test states otherwise, matrix comparisons use absolute tolerance $10^{-12}$ and zero relative tolerance. The tested identities concern finite local Fock spaces and finite graphs; no thermodynamic-limit statement is inferred.

## 2. Test inventory and implementation map

| ID | Contract | Implemented test(s) | Coverage |
|---|---|---|---|
| A01 | Majorana CAR and normalization | T00a-01 | Direct anticommutators in an independent Jordan--Wigner representation |
| A02 | Explicit gauge/spin bilinear contractions | T00a-04 | Component-wise contractions against the stated $G^\alpha$ and $M^\alpha$ matrices |
| A03 | $[\mathcal O[A],\mathcal O[B]]=i\mathcal O[[A,B]]$ | T00a-03 | Hermiticity and commutator identity over a basis of real skew matrices |
| A04 | Gauge Lie algebra and sign | T00a-05; T00b-05 | Matrix commutators and the conventional gauge-generator sign |
| A05 | Spin Lie algebra and normalization | T00a-05; T00b-05 | Matrix commutators and conventional spin-generator normalization |
| A06 | Gauge--spin commutation | T00a-05; T00b-05 | Matrix and represented-generator mixed commutators |
| A07 | Trace orthogonality and normalization | T00a-06 | Gauge, spin, and mixed trace Gram matrices |
| A08 | Parity Hermiticity and involution | T00b-01 | $\mathcal D^\dagger=\mathcal D$ and $\mathcal D^2=I$ |
| A09 | Casimir--parity identities | T00b-03 | Both spin and gauge Casimirs as parity polynomials |
| A10 | Local representation decomposition | T00b-04 plus Note 00b, Sections 3--4 | Sector dimensions, Casimir eigenvalues, and analytical identification of the two irreducible doublets |
| A11 | Physical projector | T00c-01; T00c-02 | Local orthogonal projector and finite two-site product projector; distinguishes local constraints from total parity |
| A12 | Gauge-singlet characterization | T00b-04 plus Note 00b, Section 4 | Gauge generators annihilate the physical sector; the complementary sector has nonzero gauge Casimir, establishing the converse |
| A13 | Unprojected Pauli-bilinear Clifford relations | T00c-03 | Squares and anticommutators on the full four-dimensional local Fock space |
| A14 | Projected Pauli algebra | T00c-04; T00c-05; T00c-06 | Full-space product correction, local projected multiplication, and two-site tensor-product commutation |
| A15 | Physical spin normalization | T00c-07 | $S^\alpha_{\rm phys}=\sigma^\alpha/2$, the projected $\mathfrak{su}(2)$ commutator, and $S^2=3I/4$ |
| A16 | Exchange-tensor symmetry and bond convention | T00d-01; T00d-02 | Symmetric/antisymmetric tensor components, DM orientation, and reverse-edge transpose convention |
| A17 | Quartic Hamiltonian prefactors | T00d-03; T00d-04 | Exchange and field substitutions with their exact numerical prefactors |
| A18 | Constraint multiplier normalization and projection | T00d-05; T00d-06 | Multiplier normalization and vanishing after physical projection |
| A19 | Ordered-edge counting | T00d-02; T00d-06 | Ordered half-sum equals one directed-edge contribution |
| A20 | Small-system projected Hamiltonian equivalence | T00d-06 | Direct two-site matrix comparison for a generic real tensor, fields, and nonzero multipliers |

The inventory contains 20 contracts. They are not 20 separate pytest functions: several contracts are tested jointly by the same independent matrix calculations.

## 3. Analytical statements versus numerical evidence

The following are analytical finite-dimensional statements established by the derivations in Notes 00a--00d:

1. The Majorana CAR and the bilinear commutator identity follow algebraically from the defining CAR.
2. The gauge and spin Lie algebras, their mutual commutation, and the parity--Casimir identities follow from the explicit matrix definitions and the bilinear map.
3. The local Fock representation decomposes as one spin doublet with gauge singlet quantum numbers plus one gauge doublet with spin singlet quantum numbers. The dimensions and Casimir values identify the irreducible sectors.
4. The physical projector is the product of local projectors, not a projector onto a single global fermion-parity sector.
5. The Pauli-bilinear product identity on the extended space contains a gauge-generator correction. The correction vanishes on the physical sector.
6. The Majorana exchange and field prefactors follow by direct substitution of $i c_i^T N^\alpha c_i=2\sigma_i^\alpha$.
7. The multiplier term is proportional to the gauge generators and vanishes when sandwiched by the physical projector.
8. The ordered-edge identity follows from $J_{ji}=J_{ij}^{T}$ and commutation of even operators on distinct sites.

The tests provide independent finite-matrix checks of these statements, including signs, normalizations, ranks, and small-system operator equivalence. They do not establish the identities for arbitrary system size by numerical induction; that generality comes from the analytical arguments.

## 4. Test implementation and reproducibility

The test modules are:

- `test/unit/test_local_majorana_algebra.py`: T00a-01 through T00a-07.
- `test/unit/test_local_majorana_representations.py`: T00b-01 through T00b-05.
- `test/unit/test_physical_projection.py`: T00c-01 through T00c-07.
- `test/unit/test_exchange_hamiltonian.py`: T00d-01 through T00d-06.

The Jordan--Wigner construction is implemented in the test modules rather than delegated to the production Majorana builders. Production matrices are imported only where the test explicitly compares the implementation's fixed $G$, $M$, or $N$ conventions with the independent operator construction.

The principal numerical comparison tolerance is $10^{-12}$ with `rtol=0.0`. Tests involving projector ranks use traces of finite-dimensional projectors; the analytical rank statements are separately derived in Notes 00b and 00c.

## 5. Local validation record

**Validation status at note creation: pending user-side WSL execution against the commit containing this note and T00c-07.**

The user previously reported successful local runs for M1.1a--M1.1c. The M1.1d run exposed a missing local `identity` binding in T00d-03; that source defect was corrected on the development branch. The corrected M1.1d suite and the newly added T00c-07 have not yet been confirmed by a post-update local run. No claim of a fully passing M1.1 suite is made here.

Recommended validation from the repository root:

```sh
pytest -q \
  test/unit/test_local_majorana_algebra.py \
  test/unit/test_local_majorana_representations.py \
  test/unit/test_physical_projection.py \
  test/unit/test_exchange_hamiltonian.py

ruff check \
  test/unit/test_local_majorana_algebra.py \
  test/unit/test_local_majorana_representations.py \
  test/unit/test_physical_projection.py \
  test/unit/test_exchange_hamiltonian.py

mypy \
  test/unit/test_local_majorana_algebra.py \
  test/unit/test_local_majorana_representations.py \
  test/unit/test_physical_projection.py \
  test/unit/test_exchange_hamiltonian.py
```

A successful run validates the current test implementation under the local environment. It does not replace the analytical arguments above.

## 6. Explicit scope boundary

M1.1 establishes the local Majorana algebra, local representation and physical projection, and the exact finite-system spin-to-Majorana Hamiltonian correspondence. It does not establish a mean-field decoupling, a self-consistency equation, a saddle-point approximation, an observable formula, or a thermodynamic phase diagram. Those belong to subsequent milestones and require their own contracts and validation evidence.
