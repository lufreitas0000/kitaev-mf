# W001-T4: Foundation review and claim dispositions

Reviewed source baseline: `91d59bde6bf3b96e8ac221134d1a64dd659379ae`.
Scope: original conventions/projectors, later Feshbach/third-order claims, and
finite diagnostic witnesses. This is not a complete audit of all 350 tests or
a proof of an exact bulk prefactor. Historical Phase 4 remains disapproved;
braiding remains unresolved under the PI's explicit steering.

## 1. Original requirements reinstated

The seven original exact-perturbation notes are unchanged, not missing or believed
corrupted. Their requirements must govern subsequent implementation:

- Use one declared Pauli Hamiltonian, coupling/field map, site order, directed
  bonds, and Fourier convention. A positive rescaling of a one-particle matrix
  does not establish a physical many-body gap correspondence.
- Work in the physical spin Hilbert space or explicitly derive the gauge/parity
  projection. A Gaussian vacuum need not survive projection.
- Distinguish a joint zero-plaquette-flux sector from a ground eigenspace, a
  restricted external energy window, and a one-particle occupied-band projector.
- Control the actual eliminated-block resolvent at the external energy. An
  unperturbed threshold does not establish uniform convergence at nonzero field.

The original STATUS.md still records open prerequisites while later BACKLOG
checkboxes declare completion. Neither is treated as current proof. Their actual
source/test coverage is assessed below; no history is overwritten.

### Explicit operator-level unit substitution

For the same Hamiltonian, substituting S=sigma/2 into
`-J_spin sum S_i^alpha S_j^alpha - sum h_spin^alpha S_i^alpha`
gives `K_Pauli=J_spin/4` and `h_Pauli=h_spin/2`. Thus
`J_spin=4 K_Pauli`, `h_spin=2 h_Pauli`. The production parameter convention
`+K_prod sum sigma_i sigma_j` corresponds to `K_prod=-K_Pauli` for FM coupling.
This exact coefficient substitution says nothing about spectral projectors,
physical parity, or which gap is being compared.

## 2. Four different projectors

Write these with different symbols in the canonical exposition:

1. `P_ff`: joint +1 eigenspace of all physical commuting plaquette operators,
   including global Wilson sectors unless a restriction is explicitly specified.
2. `P_0`: spectral ground projector of the chosen unperturbed finite Hamiltonian.
3. `P_Lambda`: external energy window inside Ran(P_ff), not a globally gapped
   complement in a gapless matter system.
4. `Pi_occ(k)`: occupied one-particle band projector, not a many-body flux projector.

`compute_spin_resolvent_third_order` recomputes its intermediate complement from
all excited eigenvectors of H_K, independently of its passed P. It can serve as
a fixed-energy ground-subspace calculation under matching hypotheses; it does
not generally implement the P_ff/Q_flux decomposition in the original note.
Its third-order formula also needs the appropriate perturbation convention
(e.g. PVP=0) before it is called the full energy-independent degenerate effective
Hamiltonian; energy dependence/folded terms cannot silently be ignored.

## 3. Conditional finite-dimensional result that can be retained

Let H0 and V be self-adjoint on a finite-dimensional Hilbert space, P an
orthogonal projector with [P,H0]=0, and Q=I-P. All compressed spectra below mean
operators restricted to Ran(P) or Ran(Q); extraneous zeros from the ambient space
are not sector eigenvalues. Set H(lambda)=H0+lambda V.

When E-QH(lambda)Q is invertible on Ran(Q), the Q component of an eigenvector
satisfies

`Q psi = (E-QH(lambda)Q)^(-1) QH(lambda)P psi`.

Substitution into its P equation gives the exact energy-dependent Schur/Feshbach
map

`F(E,lambda)=PH(lambda)P + PH(lambda)Q (E-QH(lambda)Q)^(-1) QH(lambda)P`.

This proves the eigenvalue reduction for nonzero P psi. Pure-Q eigenvectors are
not reconstructed by it; the invertibility hypothesis and chosen spectral range
must be retained. Write R0(E)=(E-QH0Q)^(-1) on Ran(Q). Factorization gives

`(E-QH0Q-lambda QVQ)^(-1) = R0(E) (I-lambda QVQ R0(E))^(-1)`.

If `abs(lambda) ||QVQ R0(E)|| < 1`, the finite-operator Neumann series converges
in norm. With [P,H0]=0 this yields

`F(E,lambda)=PH0P+lambda PVP`
` + sum_(n>=0) lambda^(n+2) PVQ R0(E) [QVQ R0(E)]^n QVP`.

This is a conditional fixed-energy operator series, not an evaluation of its
coefficients for the physical honeycomb problem, and not a thermodynamic theorem.
For real E below the true lower edge b of QH0Q,

`||R0(E)||=1/(b-E)`.

This follows by diagonalizing the finite self-adjoint compressed operator. If
E<=E0+Lambda and b=E0+Delta_flux, Lambda<Delta_flux gives the upper bound
`||R0(E)||<=1/(Delta_flux-Lambda)`. For the perturbed block one must check its
own invertibility or a sufficient Neumann condition. No volume-uniform estimate
follows from the unperturbed bound alone.

These statements retain their explicit hypotheses. They do not certify the
current services, which do not enforce the complete contract.

## 4. Implementation and evidence gaps

| Claim or component | What the source/test actually establishes | Current disposition |
|---|---|---|
| Single-W spin thresholds | Compares two signs of one supplied involution, then relabels the lower sector as reference | Needs revision for a fixed joint flux-free contract |
| Removal of projected zeros | Removes every small eigenvalue, including genuine sector zeros | Energy-reference covariance fails in a 2x2 witness |
| Torus eliminated-block minimum | Samples one adjacent pair pattern across four twists | A candidate sector minimum, not the infimum over all nonzero-flux configurations |
| `1/delta_flux` torus resolvent certificate | Reciprocal of the sampled candidate gap | Not an established upper bound for the full-Q resolvent |
| Cubic third-order test | An expression with three factors of a linear V is homogeneous of degree three | Finite homogeneity regression, not an independent coefficient evaluation |
| Scalar-gap coefficient | Implements the assumed h_x h_y h_z / Delta_eff^2 expression | Approximation; exact coefficient and K-019 reconciliation remain unsupported |
| Fifth-order and 2/L bounds | Returns assumed ratios and cutoffs; no higher-order coefficient or uniform remainder estimate is calculated | Heuristic indicators, not certified bounds |
| One-spin Zeeman Y component | Returns +sigma_y for a positive y field, rather than -sigma_y in the declared standard basis | Convention mismatch; preserve and repair only in a separate reviewed task |
| Exact-matter geometry/parity tests | Finite geometry identities and selected finite oracle comparisons | Useful bounded evidence; no blanket thermodynamic or zero-mode-state classification certificate |

For the torus bound, let C be the sampled subset of Ran(Q_flux). Then
`min spectrum(H0|C) >= inf spectrum(H0|Ran(Q_flux))` when the sectors are comparable.
Thus the sampled gap is an upper estimate of the full gap, not a certified lower
estimate. Its reciprocal can underestimate the true resolvent norm. A search over
four twists does not remedy missing flux patterns or separations.

The existing physicality predicate also returns True for every occupation vector
when a zero singular value is detected. A degenerate zero-mode energy can admit
physical states without making every arbitrarily labeled vector physical. This
requires a separate state-versus-energy contract, not a blanket acceptance claim.

## 5. Primary-source reconciliation

The local References pilot supplies original source packages, not a trusted
finished transcript. I verified the archive and original TeX hashes directly;
the separate literal text companion was used only for navigation. Version and
hashes are pinned in `doc/reference/literature_manifest.yaml`.

**Kitaev, cond-mat/0506438v3, hexagon.tex, Section 6:**

- Under “Derivation of an effective Hamiltonian,” original lines 1637-1644 retain
  a nonzero time-reversal-even second-order term. The unconditional statement
  in the later reconciliation note that first-order vanishing makes second
  order vanish is not supported and is contradicted by this source.
- Original lines 1645-1665 describe the gauge-sector Green function and then
  introduce a qualitative scalar denominator. Label `Hmageff` uses scaling, not
  a calculated exact scalar coefficient. This fails to establish the project
  coefficient; it does not prove that no separate exact derivation could exist.
- Both bilinear and four-fermion third-order configurations appear. The latter
  cannot be removed from the full effective operator without an explicit scope.
- Labels `spectrum1` and `maggap` give signed eigenvalues +/-sqrt(|f|^2+Delta^2)
  and mass Delta=6 sqrt(3) kappa. The positive excitation at a cone is |Delta|;
  the separation of the two signed bands there is 2|Delta|. Many-body physical
  spin gaps and spectroscopy thresholds are additional quantities needing parity,
  sector, and matrix-element analysis.

The production gap benchmark explicitly computes e5(K)-e4(K), a signed-band
separation. Multiplying by a matrix-normalization factor preserves that meaning;
it does not alone turn it into the lowest physical spin excitation or an exact
coefficient. The named “global gap” is a minimum along a sampled high-symmetry
path, not a certified minimum over the full Brillouin zone.

**Pedrocchi–Chesi–Loss, 1105.4573v2, label `eq:P0`:** the physicality condition
depends on lattice labeling, topology, gauge configuration, and eigenmode parity.
The rectangular code/Clifford tests must be assessed within their fixed ordering;
this source does not license transferring a parity sign between conventions.
Selected source-span inspection is not full-paper or transcript certification.

## 6. Finite diagnostic witnesses

Reproducer: `scripts/audits/w001_foundation_checks.py`. Compact observations:
`doc/audits/w001_foundation_observations.json`; raw generated output stays under
`build/rnd/exact-perturbation-theory/`. Production source hashes are included.

- H=diag(-3,-1), W=diag(1,-1): both reference and service return a gap of 2.
  Adding 3I leaves the exact sector gap 2 but the service raises ValueError after
  removing the genuine zero reference eigenvalue.
- H=diag(-4,-3,-2,-1) with commuting W1=(+,+,-,-), W2=(+,-,+,-): the joint
  flux-free projector has complement gap 1; the single-W call returns 2. This
  illustrates scope, not a claimed bulk-honeycomb counterexample.
- Reversing the sign-sector labeling makes the service swap P/Q to retain a
  positive gap. This is not the fixed-projector contract.
- The one-spin y-field call exactly equals +sigma_y, while the documented
  standard Zeeman operator is -sigma_y. No floating-point tolerance is involved.
- Independent tensor-product construction on the isolated six-spin ring gives
  a four-dimensional ground eigenspace, gap about 0.5359, PVP norm 2.53e-16,
  and second-order norm 0.038124 at field components (0.02,0.02,0.02). Its
  non-scalar residual is 0.0235253. The ground grouping scale 1e-10 is more than
  nine orders below the next energy and projector residual is 6.10e-16.
  This is a finite open-ring witness, not a periodic bulk result.

Independently of that numerical example, for a true ground projector P0 and
R0=Q(E0-QH0Q)^(-1)Q with a positive gap on Q, second order satisfies

`<psi|P0 V R0 V P0|psi> = -||(QH0Q-E0)^(-1/2) QV psi||^2`.

It is negative semidefinite and need not vanish when P0VP0=0. Taking the exact
2x2 example H0=diag(0,Delta), P0=diag(1,0), V=h sigma_x gives
`H_eff^(2)=-(h^2/Delta)P0`, Delta>0. Thus the claimed inference is false even
before invoking any bulk approximation. Time reversal does not forbid every
even-order effective term; it constrains their transformation properties.

## 7. What may enter the book, and what must wait

The logical chapters may use the explicit convention map, distinct-projector
definitions, conditional finite-dimensional Schur/Feshbach theorem, and the
second-order sign identity with their proofs and assumptions. The chronological
appendix should retain the diagnostics, primary-source locators, unsupported
completion claims, and repair backlog. No exact bulk prefactor, linked-cluster
convergence, physical ED gap identification, or non-Abelian braid certificate is
admitted by this review.

Required later implementation work: fixed-sector compressed eigenspaces including
true zeros; explicit joint projectors where that is the requested scope; declared
ground-versus-flux resolvents and field conventions; independently evaluated
coefficient/observable bridges; actual truncation estimates if claimed. Do not
silently implement those repairs during a documentation audit.

## 8. Independent review and next handoff

Adversarial review passed for the conditional finite algebra and bounded audit:
Schur factorization, Neumann ordering, fixed-energy expansion, spectral norm bound,
and the second-order sign identity were checked against their stated assumptions.
The reviewer read the primary TeX spans and diagnostic construction but did not
independently execute the numerical script. This is not integration merge approval.

The next labbook task must retain the separation between source-supported facts,
mathematical conditional results, observed finite witnesses, and open physical
questions. Other foundational book claims, such as Gaussian upper bounds in an
enlarged versus projected physical state space and Bloch particle-hole coupling,
are not accepted merely by this perturbation audit; review them before reuse.
No new bulk computations, source-code fixes, or manuscript result are required
by this audit's completion.
