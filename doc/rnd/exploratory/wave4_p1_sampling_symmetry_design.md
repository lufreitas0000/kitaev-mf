---
id: RND-W4-P1-SAMPLING-SYMMETRY
title: "Wave 4 P1 — Sampling populations, covariance admissibility, and symmetry actions"
phase: 1
status: draft
questions:
  - "Which population of initial mean-field inputs is the solver-basin probability conditioned on?"
  - "Which covariance, symmetry, and parity constraints can be established independently?"
  - "What would falsify the legacy symmetric-sampler claim before any numerical pilot?"
depends_on:
  - "kitaev-mf PR #16 (merged statistical specification)"
  - "kitaev-mf PR #17 (merged paired-case infrastructure)"
  - "kitaev-mf Issue #7"
  - "kitaev-mf Issue #8"
contracts: []
artifacts:
  - src/kitaev_mf/domain/state.py
  - src/kitaev_mf/domain/lattice.py
  - src/kitaev_mf/services/seeds.py
  - src/kitaev_mf/domain/wave4_contract.py
  - doc/rnd/phase_1/m1_6_t11_wave4_paired_solver_protocol_draft.md
---

# Wave 4 P1 — Initial-state populations and symmetry: theory/design only

**Status and authority.** Initiated by Lu as a bounded W4-P1 design task. This is a **noncanonical, unregistered research proposal** on a separate draft PR; not a completed symmetry derivation, certified physical sampler, numerical result, new hypothesis registration, or approval of a physical pilot. Do not edit the historical benchmark, production solver, canonical rules, \`PROJECT_STATE.yaml\`, or the merged W4 W1/W2 implementation under this PR.

**Inspected source snapshot:** [\`kitaev-mf/main\` at \`770db3ead9d7834c75c9399efcb8e87a0a3451b6\`](https://github.com/lufreitas0000/kitaev-mf/tree/770db3ead9d7834c75c9399efcb8e87a0a3451b6). Coordinating issues: [#7 (paired solver experiment)](https://github.com/lufreitas0000/kitaev-mf/issues/7), [#8 (symmetry and parameter counts)](https://github.com/lufreitas0000/kitaev-mf/issues/8). The accepted [W4 specification](https://github.com/lufreitas0000/kitaev-mf/blob/770db3ead9d7834c75c9399efcb8e87a0a3451b6/doc/rnd/phase_1/m1_6_t11_wave4_paired_solver_protocol_draft.md) defines only generic populations and statistical estimands. The source-derived claims below are explicitly distinguished from proposed mathematical constructions.

## 1. Model and coordinate bookkeeping

The intended family is the ferromagnetic Kitaev honeycomb model with field \`h || [111]\`, described in the repository by \`K=-1\` and the chosen spin/Zeeman conventions. The **actual field magnitudes, unit conventions, finite torus, mesh, boundary conditions and branch references have not been approved for W4**. They are parameters of a future \`theta\`, not implicit defaults.

Let \`n\` be the number of sites (sublattice labels) in a chosen cell and \`m\` its **directed** bond entries, each with its own cell offset \`R_b\` and type \`gamma_b\`. The package stores

$$
x=(U,V,\lambda),\quad
U\in\mathbb R^{m\times4\times4},\quad
V\in\mathfrak{so}(4)^n,\quad
\lambda\in\mathbb R^{n\times3}.
$$

The directed \`U_b\` is a **general real 4-by-4 cross-site block**; the reverse bond obeys \`U_{ji}(-R_b)=-U_{ij}(R_b)^T\`, but this does *not* imply \`U_b^T=-U_b\`. The \`lambda\` variables are gauge Lagrange multipliers, **not CAR covariance entries**.

For the actual 2-site primitive, 4-site rectangular and 8-site \`make_honeycomb_supercell(2,2)\` geometry, \`m=3n/2\`. The following counts are exact **storage/algebraic** dimensions, with no physical-representability or symmetry claim:

| Cell | \`n\` | \`m\` | Stored real slots \`16m+19n\` | Independent algebraic real coordinates \`16m+9n\` | Symmetry-invariant dimension | Gaussian-representable dimension |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| Primitive 2-site | 2 | 3 | 86 | 66 | Not derived | Not determined |
| Rectangular 4-site | 4 | 6 | 172 | 132 | Not derived | Not determined |
| 2x2 8-site supercell | 8 | 12 | 344 | 264 | Not derived | Not determined |

Here \`V_s\` contributes \`4*3/2=6\` independent entries per site, not 16. \`state.py\` checks dimensions and onsite skew symmetry numerically (absolute deviation \`<=1e-12\`), but **does not itself validate all dtypes, finiteness, global CAR positivity or a spatial symmetry representation**. The present [Wave 4 case constructor](https://github.com/lufreitas0000/kitaev-mf/blob/770db3ead9d7834c75c9399efcb8e87a0a3451b6/src/kitaev_mf/domain/wave4_trials.py) adds finite float64 snapshot checks and labels cases \`algebraic_unverified\`; that label is intentional.

Do not identify these counts with a pure-Gaussian manifold dimension, physical spin Hilbert-space dimension, or the coordinate dimension of a symmetry-restricted self-consistent map. A finite one-cell torus can also *alias distinct translated bonds* unless its periodic geometry and momentum representation are specified carefully.

## 2. Three different populations; no mixing of estimands

Let \`S_alg(theta)\` be the finite-dimensional set of real \`(U,V,lambda)\` with the shape, skew and finiteness constraints, plus any explicitly stated magnitude/normalization bound. The code's constructor-level tolerance and a strict mathematical definition are distinct contracts.

**E_alg — algebraic robustness:** A candidate explicit measure is the push-forward of a standard Gaussian \`Z~N(0,I_d)\` by

$$
x_\delta=x_\star+\delta L_\theta Z,\qquad L_\theta:\mathbb R^d\longrightarrow
\mathbb R^{16m}\oplus\mathfrak{so}(4)^n\oplus\mathbb R^{3n}.
$$

This defines a probability measure \`mu_alg(theta,delta)\` once \`x_star\`, an ordered basis, \`L_theta\` (including separate scales for \`U,V,lambda\`), and the RNG are fixed. **It does not have bounded support**; \`delta\` is a scale, not a deterministic ball radius. A bounded radial/directional distribution is a *different* candidate with a different estimand. These points require predeclared treatment when interpreting any \`R90\`.

**E_sym — explicitly invariant coordinate robustness:** Choose a finite group \`G\` and a **verified** action \`rho(g)\` on the directed-bond/site/multiplier coordinates compatible with the chosen finite torus, bond orientation, gauge convention and field. Define

$$
S_{\rm sym}=\{x\in S_{\rm alg}: \rho(g)x=x\quad\forall g\in G\}.
$$

The sampling measure must have its support in \`S_sym\` and state whether it is a conditional law, an average/projected push-forward, or a distribution built directly in a basis of the invariant subspace. These generally differ. A claimed \`C3\` action cannot be justified by merely tying sites' onsite \`V\` entries or setting bond diagonals equal.

**E_Gauss — globally CAR-representable robustness:** Fix **one finite torus** with \`N\` actual sites and 4 Majoranas per site (4N modes), a consistent set of directed bond images and boundary conditions. For the convention

$$
\omega(c_a c_b)=\delta_{ab}+i\Omega_{ab},
\qquad \Omega^T=-\Omega\in\mathbb R^{4N\times4N},
$$

mixed quasifree states are characterized by the exact spectral condition

$$
-\mathbb 1\preceq i\Omega\preceq\mathbb 1
\quad\Longleftrightarrow\quad \|\Omega\|_2\le1.
$$

For *pure* quasifree states additionally \`Omega^2=-I\`. With the repository's \`U,V\` definitions and \`C_ij\` convention, the relevant covariance blocks are \`Omega_ii=-V_i\` and \`Omega_(i,0),(j,R_b)=-U_b\` for the correctly embedded directed bond. Define the observation map \`Pi_cell(Omega)=(U(Omega),V(Omega))\` only after bond-cell translations and normalization have been fixed. Then \`S_Gauss\` is the image under \`Pi_cell\` of admissible **global** covariances, optionally crossed with separately chosen \`lambda\` coordinates. It need not coincide with \`S_alg\` or \`S_sym\`; \`lambda\` cannot be certified by a fermionic covariance bound. Being Gaussian in the **enlarged Majorana Fock space** also does not alone certify projection into the physical-spin parity/gauge sector.

**Exact counterexample to local-only certificates.** Take three distinct sites, one Majorana flavor on each, and embed the real skew block

$$
\Omega_3=\begin{pmatrix}
0&t&t\\ -t&0&0\\ -t&0&0
\end{pmatrix},\qquad 1/\sqrt 2<t\le1,
$$

into a larger Majorana covariance with other entries zero. Each two-site principal block corresponding to pairs (0,1) or (0,2) has spectral norm \`|t|<=1\`, yet the full covariance has spectral norm \`sqrt(2)|t|>1\`. Thus positivity of each tested bond's 8-by-8 two-site correlation matrix cannot certify a single globally representable state. This is a **finite-dimensional algebraic fact**, independent of mean-field numerics.

For a **candidate finite-torus pure-Gaussian push-forward**: start from an *independently certified* reference \`Omega_*\` with \`Omega_*^2=-I\`, sample an antisymmetric \`A\` from a declared normalized law, and set

$$
Q_\delta=\exp(\delta A),\qquad
\Omega_\delta=Q_\delta\Omega_*Q_\delta^T.
$$

Then \`Q_delta\` is special orthogonal and \`Omega_delta\` is pure/quasifree by construction; project through \`Pi_cell\` and sample \`lambda\` separately. This construction is **proposed, not implemented, not uniform Haar**, and may break translation/sublattice/\`C3\` symmetries. \`Omega_*\` must truly exist on the chosen torus; the nominal \`make_qsl_seed\` or \`make_polarized_seed\` arrays alone do not constitute that global construction.

If a verified orthogonal group action \`R_g\` and \`G\`-invariant \`Omega_*\` exist, one may draw \`A\` from a declared distribution and set \`A_G=|G|^{-1}sum_g R_g A R_g^T\`. Then \`exp(delta A_G) Omega_* exp(-delta A_G)\` is invariant **provided the representation and baseline invariance really hold**. This supplies a *candidate* intersection population \`E_Gauss ∩ E_sym\`, not an established current-code feature.

**Rule:** Report \`P_conv\`, \`P_target\` and costs separately for each fully specified population and reference branch. Never pool \`E_alg\`, \`E_sym\`, and \`E_Gauss\` samples as if their measures or physical meanings matched.

## 3. Source audit: current "symmetric" perturbation is not a certificate

The existing [\`services/seeds.py\`](https://github.com/lufreitas0000/kitaev-mf/blob/770db3ead9d7834c75c9399efcb8e87a0a3451b6/src/kitaev_mf/services/seeds.py) contains three independent issues:

1. In \`add_symmetric_perturbation\`, \`gamma=(b % 3)+1\` ignores \`cell.bonds[b].bond_type\` and the optional \`cell\` argument is not used. The 2-site primitive bond list is ordered \`[z,x,y]=[3,1,2]\`, so the resulting \`U_b[gamma,gamma]\` changes **different Majorana flavors than those selected by the actual bond types**. The 4-site rectangular bond order is \`[z,z,x,x,y,y]\`, also incompatible with that index rule. This is a source-code counterexample to the advertised bond-type construction, independent of any physical rotation test.
2. The shared \`dv_skew\` is an **arbitrary** skew 4-by-4 onsite matrix. Even under the simplest illustrative cyclic permutation \`P:1->2->3->1\` of the three gauge flavors (matter flavor 0 fixed), the skew basis element \`E_01-E_10\` maps to \`E_02-E_20\` and is not invariant. This simple action is a **test candidate**, not yet a proven physical \`C3\` implementation with bond/gauge transformations.
3. Applying the same perturbation to both sublattices or setting the added \`lambda\` components equal does not establish their transformation under **actual** site exchange, oriented-bond reversal or a local gauge action. The old \`check_seed_representability\` only checks individual 8-by-8 bond blocks. \`add_controlled_perturbation\` is algebraic noise, despite the phrase "Gaussian perturbation" in its docstring; it does not enforce a global quasifree covariance.

These defects are **negative design findings**, not evidence that the legacy numerical solver returned a physically incorrect phase diagram. No historical data or claim labels are altered in this PR.

## 4. Symmetry implementation obligations, coordinated with Issue #8

For each actual finite cell/torus choose the geometric map \`g\` on site labels, Bravais offsets and **directed** bond entries, including when the transformed bond is represented with reversed orientation. Define its action on Majorana generators \`c_{i,mu}\` via orthogonal matrices and appropriate site relabeling, respecting \`{c,c}=2I\`. Derive the induced action on \`U\`, \`V\` and \`lambda\` rather than guessing equality of blocks. Verify these axioms **before claiming \`C3\` or sublattice exchange**:

- Representation/group law: \`R_e=I\` and \`R_g R_h=R_(gh)\` (or a precisely defined projective/gauge-equivalent law), together with orthogonality and CAR preservation.
- Source geometry: transformed bonds, types and offsets match the finite lattice, including orientation-reversal signs and momentum/gauge conventions; choose a torus admitting the claimed rotation.
- Physical action: the transformed Hamiltonian at the declared \`theta\` is equivalent to the original, including Zeeman field and gauge conventions. Separate \`h=0\` from \`h || [111]\`.
- Coordinate consistency: projection of a transformed global covariance onto \`(U,V)\` commutes with the derived coordinate action. Multipliers transform as the correct constraint covectors, not automatically like a bare spin vector.
- Subspace closure: compute \`rho(g) Phi(x)\` versus \`Phi(rho(g)x)\` using a suitably instrumented future map; agreement at finitely many samples is **numerical evidence**, not a group-theoretic proof. A symmetry-constrained solver is also a *different algorithm* from an unrestricted solver.
- **No unsupported dimension:** determine the rank of the group projector/invariant subspace in real coordinates before claiming \`D=5\`, \`D<=19\`, or any speedup. Store counts separately from the raw \`16m+9n\` algebraic dimension above.

Issue #8 remains the owner of the **full parameter-count and 2/4/8-cell representation-theoretic derivation**; W4-P1 only supplies the sampling-relevant proof obligations and specific source counterexamples, without duplicating or closing #8.

## 5. Candidate operational sampling designs (owner choice pending)

| Question/ensemble | Constructible object | Mathematical certification required | What it would measure | Status |
| --- | --- | --- | --- | --- |
| A. Algebraic-coordinate perturbations | \`x_*+delta L Z\` on \`U,V,lambda\` | Dimensions/onsite skew/finite values; exact law, metric and support | Robustness of a numerical solver to arbitrary admissible *coordinates* | Candidate, physical validity unverified |
| B. Explicitly \`C3\`-restricted coordinates | Random coordinates in computed fixed subspace of \`rho(G)\` | Representation, Hamiltonian equivariance, group fixed-subspace basis | Restricted-symmetry basin recovery | Blocked by group-action derivation |
| C. Finite-torus pure Gaussian covariances | \`Omega_delta=Q_delta Omega_* Q_delta^T\`, then \`Pi_cell\` | Independently valid \`Omega_*\`, torus embedding, CAR and purity, declared law | Robustness within enlarged-Fock quasifree states | Theoretical candidate, no sampler |
| D. Symmetric Gaussian covariances | Group-averaged Lie-algebra perturbations of invariant \`Omega_*\` | All of B/C, including nonempty intersection and field-dependent group | Robustness *conditional on both constraints* | Blocked by representation and baseline |
| E. Physical-spin-projected family | Explicit normalized physical-sector construction, if possible | Actual gauge/parity projector, nonzero projected norm and observable consistency | Robustness for a **physical spin** family | Separate research question, not interchangeable with C |

This table is a **menu of distinct scientific experiments**, not five approved solver arms. \`P_target\` also needs an independent branch oracle; none of these population choices supplies it automatically. A proposed measure on a null or unknown-reference population is not a meaningful experiment.

## 6. Falsification checks and evidence gates (no runs in this PR)

**S1 — Symmetry counterexample/recovery.** On the actual 2-site and 4-site cells, inspect \`bond_type\` vs index-selected \`gamma\`; construct a non-invariant skew onsite direction for a declared \`R_{C3}\`. Require a corrected group action and a test showing invariant seed and perturbation before any \`C3\` label.

**S2 — Local/global covariance counterexample.** Use the analytic 3-mode star example above. A future certification routine must reject the globally inadmissible covariance even though each selected pair is admissible.

**S3 — Global Gaussian sampler.** For a finite, specified torus, verify skew/realness, \`||Omega||_2<=1\`, purity when claimed, and trace/antisymmetry of each projected onsite/bond block. Verify deterministic reconstruction of identical samples from complete generator settings *or* record exact initial arrays through W4's merged archive. Check nonempty support rather than condition on rare invalid draws silently.

**S4 — Invariance and distribution.** Prove exact group-action closure for the defined ensemble (or explicitly label tests as finite-sample checks); measure violations under full coordinates as negative controls. Document any projection-induced change of the distribution and any dependence between amplitudes built from one random direction.

**S5 — Physical meanings.** Distinguish (i) stationary solution; (ii) stable local minimum; (iii) globally lowest energy among found branches; (iv) projected physical spin state. The three-way target/other/unknown classifier remains independent and its reference must be fixed before recording pilot outcomes.

No tolerance, sample count, field window, Gaussian parity sector, statistical confidence level, \`R90\` probability target, or acceptance threshold is adopted here.

## 7. Decisions and next independently gated tasks

**TODO(Lu) — Scientific target:** Select whether W4 initially tests numerical **algebraic-coordinate robustness** (A) or globally Gaussian-admissible robustness (C). B/D are subsequent distinct symmetry-conditional questions. Decide whether pure/mixed quasifree states and a fixed global parity component matter; specify which reference branch and whether \`lambda\` is held fixed or perturbed separately.

**TODO(Lu) — Domain and measure:** Choose finite cell/torus with a non-aliased interpretation of \`R_b\`; choose the base state \`Omega_*\` if Gaussian sampling is desired; define the push-forward measure, norm/metric, bounds and amplitude scaling; fix the [111] field convention. Issue #8 must derive the genuine group action before any invariance certificate.

**TODO(Lu) — Next implementation scope:** A later W4-P1b PR may implement **only the chosen certified sampler and its mathematical negative controls**, not physical pilot runs. W4-P2 must independently implement branch-equivalence oracles; W4-P3 must instrument solver/callback costs and expected-cohort completeness; statistical and numerical pilot choices require a separately frozen protocol.

**Independent review request:** Challenge the finite-torus construction, \`C3\` actions, site/bond tensor signs, strict population distinctions and whether a proposed measure can be sampled without hidden rejection/conditioning. Check the explicit bond-order source counterexample. Do not promote this design note to a validated sampler or merge without reviewed scope.

**STOP:** this PR is **one exploratory document only**, with no code, tests, numerical data, legacy record, scientific verdict, \`PROJECT_STATE.yaml\` or roadmap mutations. W4 remains scientifically open.
