---
id: RND-P1-M1.10-T1-W2
title: "Wave 2 — K-010/K-019 cubic-gap scientific-claim evidence pilot"
phase: 1
milestone: M1.10
task: T1
status: draft
questions:
  - "What does the finite-mesh Dirac-gap fit establish, and what is still asymptotic or predictive?"
  - "Which perturbative terms and operator identities are established under which projectors and conventions?"
  - "What independent evidence would distinguish the exact-spin and mean-field cubic-gap coefficients?"
depends_on:
  - M1.6-T5
  - M1.10-T1
contracts:
  - K-010
  - K-019
artifacts:
  - scripts/m1_6/benchmark_majorana_gap.py
  - data/processed/benchmarks/m1_6/majorana_gap_scaling.dat
  - doc/rnd/phase_1/m1_10_t1_exact_perturbation_and_normalization_reconciliation.md
---

# Wave 2: cubic-gap evidence pilot — K-010 and K-019

**Status: UNREGISTERED METHODOLOGICAL DRAFT; INDEPENDENT REVIEW REQUIRED.** This is neither an amended hypothesis contract nor a scientific verdict. P1–P6 below are *candidate propositions/questions*, not six accepted claims. The current K-010/K-019 statuses, historical findings and acceptance thresholds are not modified.

**Source snapshot:** `lufreitas0000/kitaev-mf` `main` at full Git commit `4c2d629ec7ee924c731c566fc00d9dc0411217c6`, inspected 2026-11-14. All [S01]–[S26] links below are pinned to this revision and include line locators. The [Wave 1 noncanonical methodology](https://github.com/lufreitas0000/project-governance/blob/9ef392d6295f560a00de8ce41c00258d9c905d82/docs/SCIENTIFIC_CLAIM_METHODOLOGY_WAVE1_DRAFT.md) supplies questions, **not** a new registry or validator. No benchmark was rerun, no numerical data were generated and no analytic perturbation coefficient was derived for this draft.

## 1. Scientific question and the two historical records

Can a field-induced Majorana Dirac mass in the **physical** isotropic spin-1/2 Kitaev model be related, with a controlled normalization and approximation, to the **self-consistent four-flavor Gaussian mean-field** Dirac gap?

- **K-010** records the `numerical_power_law_fit_with_finite_mesh_regression` status and a reported regression with prefactor `2.5921 ± 0.0034`, exponent `2.9936 ± 0.0006`, `R²=0.99999`, interval `0.005≤h≤0.150` and `24×24` momentum mesh [S01]. The legacy hypothesis table labels its combined weak-field claim `PASSED` [S03]. Those are historical record labels, not an independent proof of an asymptotic exponent.
- **K-019** records `verified_operator_derivation`, an asserted projected three-spin/matter-hopping identity, a fourfold spectral conversion, and a static-denominator comparison [S02]. The legacy table labels `HC-P1-01` `PASSED` [S03]. The independent exact-perturbation track explicitly reports **no implemented perturbative coefficients** and lists physical parity and resolvent domains as pending [S08,S09]. Thus the existing K-019 status cannot certify every analytic or physical inference in its accompanying note.

Neither a passing arithmetic unit test nor the historical status converts a fitted mean-field prefactor into a physical-spin-model perturbation coefficient. The unresolved issues below are **new audit observations**, not retroactive reclassifications.

## 2. Hamiltonians, units, and spectral quantities

The two tracks use *different signs for the coupling symbol*, but describe matching ferromagnetic Pauli-spin Hamiltonians under the following map, derived in the exact-track normalization note [S10,S24]:

$$
\begin{aligned}
H_K^{\mathrm{P}} &= K_{\mathrm P}\sum_{\langle ij\rangle_\gamma}
\sigma_i^\gamma\sigma_j^\gamma,\quad K_{\mathrm P}<0,\\
H_K^{\mathrm R} &=-K_{\mathrm R}\sum_{\langle ij\rangle_\gamma}
\sigma_i^\gamma\sigma_j^\gamma,\quad K_{\mathrm R}>0,\qquad
K_{\mathrm P}=-K_{\mathrm R}.
\end{aligned}
$$

Each directed `A→B` physical bond is counted once. For conventional `S^\gamma=\sigma^\gamma/2`, the spin-exchange map is `K_{\mathrm R}=J/4` and the field map between `-\mathbf h_S\cdot\mathbf S` and `-\mathbf h_\sigma\cdot\boldsymbol\sigma` is `\mathbf h_\sigma=\mathbf h_S/2`. This Pauli-to-spin change is **distinct** from the Bloch-matrix scaling [S07,S10]. The production benchmark fixes `K_{\mathrm P}=-1` and `\mathbf h_\sigma=(h/\sqrt3)(1,1,1)`, where `h=|\mathbf h_\sigma|` is the scanned scalar [S04]. In general, report a dimensional gap as `|K|\,c(h/|K|)^p`; write `c h^3/|K|^2` only **after** establishing or assuming `p=3`.

The production matrix is `H_{\mathrm B}(\mathbf k)=i A(\mathbf k)/4`, an 8×8 Hermitian matrix in the two-site, four-flavor cell [S12]. Let its ascending eigenvalues be `w_1(\mathbf k),\ldots,w_8(\mathbf k)`. The *actual regression observable* is the **particle-hole band separation at the chosen Dirac corner**:

$$
\Delta_{\mathrm B,K}(h)
:=w_5(\mathbf K;h)-w_4(\mathbf K;h),\quad
\mathbf K=(1/3,2/3).
$$

This is exactly the script's `gap_k` [S04]. It must not silently be identified with the *minimum positive quasiparticle energy*, a many-body spin excitation threshold, or a minimum over the full Brillouin zone. In a symmetric ± pair, `w_5-w_4=2w_5`, so the definition of a "gap" already contains a potential factor of two relative to a minimum positive frequency.

## 3. K-010: what the saved finite calculation establishes

**Executed protocol reconstructed from source (not rerun):**

| Item | Verified recorded implementation | Audit implication |
| --- | --- | --- |
| Input | 17 field values including zero, `0.005,\ldots,0.350`, `K_{\mathrm P}=-1`, field components `h/\sqrt3` [S04] | There is no finite-field data below `h=0.005` |
| Branch | Forward continuation from the QSL seed; state is carried forward **only on convergence** [S04] | Branch identification depends on continuation protocol |
| Solver | CPU backend; DIIS; `24×24` BZ integration; at most 40 iterations; `tol_correlator=tol_gauge=10^{-8}`; mixing `0.6`; skew retraction [S04] | Convergence is for this discretized solver and its implemented residual, not a continuum certificate |
| Dirac observable | Eight eigenvalues at exactly `\mathbf K`, `w_5-w_4` [S04] | Direct algebraic gap at a fixed momentum |
| So-called global gap | Minimum of `w_5-w_4` over **75 points** on `Γ→K→M→Γ`, not a 2D BZ minimization [S04] | Equality to `gap_k` in saved rows does **not** prove a global BZ minimum |
| Fit data | The eleven positive field values in `[0.005,0.150]`; zero field excluded [S04,S05] | Empirical finite-window scaling only |
| Fit model | Unweighted `np.polyfit(np.log(h),np.log(gap_k),deg=1)` with `p` and `\log c` [S04] | The script does **not** generate `R²`, `prefactor_stderr` or `exponent_stderr` |
| Output | Committed `majorana_gap_scaling.dat` with per-point gap, convergence flag, evaluation count and wall time [S05] | A source-linked dataset exists, but its header gives an abbreviated generating commit `910bdc6`, not a fully pinned environment/dependency/seed record |

All 17 saved rows show `conv=True`; the tabulation displays equality of the fixed-K and sampled-path gaps **at its printed precision** [S05]. The fitted values and their claimed uncertainties are recorded in `PROJECT_STATE.yaml` [S01] and the earlier note [S06], but the script has **no uncertainty-estimation code** or model/window-sensitivity comparison. Nor does it validate the implicit finite-mesh projection near Dirac zero modes, whose implemented spectral projector uses a fixed numerical `zero_tol=10^{-12}` and half-weight convention [S14]. A fit can be informative while its error bars remain untraceable.

**Critical mismatch with existing stronger prose:** the earlier R&D note calls a finite fit confirmation of "**exact cubic scaling**" [S06], and the labbook claims the "true global minimum" and a theorem-level cubic selection rule [S21]. The benchmark cannot verify either as stated. Moreover, **time-reversal oddness and cyclic `C_3` invariance alone do not select `h_xh_yh_z` uniquely**: `h_x+h_y+h_z` is also odd and `C_3`-invariant. Excluding a linear mass requires additional *operator representation, symmetry and/or flux selection rules*; an asserted sublattice argument must be checked in the actual mean-field representation. Even if a signed mass starts at cubic order, an unsigned gap behaves as `|m(h)|` and does not itself transform as a TR-odd scalar. The finite regression does not establish a convergent expansion or an `O(h^5)` remainder [S21].

## 4. K-019: perturbative projection and the second-order error

In the exact-spin convention `H_0=H_K^{\mathrm R}`, `V=-\sum_{i,\gamma}h_\gamma\sigma_i^\gamma`. A physical spin insertion flips fluxes adjacent to its bond. **Conditionally** on a properly defined flux-free subspace `\mathcal H_{\mathrm{ff}}` and sector projector `P_{\mathrm{ff}}`, the selection rule is

$$
P_{\mathrm{ff}}\,V\,P_{\mathrm{ff}}=0.
$$

This is an **intra-sector** statement, not a statement about all virtual processes. The historical note writes `P` both as "zero-flux physical ground-state subspace" and as the zero-flux projector, without resolving whether it contains *all* flux-free matter excitations or only the degenerate `E_0` eigenspace [S07]. These two projectors have different resolvent and perturbative contracts, especially for gapless matter [S08,S09].

For an isolated finite-volume eigenspace with `H_0P=E_0P`, `Q=1-P`, and a well-defined inverse `(E_0-QH_0Q)^{-1}` on `Q\mathcal H`, a standard energy-dependent Feshbach expansion has

$$
H_{\mathrm{eff}}^{(1)}=PVP=0,\qquad
H_{\mathrm{eff}}^{(2)}=
PVQ(E_0-QH_0Q)^{-1}QVP.
$$

The **source's inference** `PVP=0\Rightarrow H_{\mathrm{eff}}^{(1)}=H_{\mathrm{eff}}^{(2)}=0` [S07] is not mathematically justified. Under the additional assumption `QH_0Q-E_0>0` and nonzero `QV|\psi\rangle`, for normalized `|\psi\rangle\in\operatorname{Ran}P` one obtains

$$
\langle\psi|H_{\mathrm{eff}}^{(2)}|\psi\rangle
=-\bigl\langle QV\psi,(QH_0Q-E_0)^{-1}QV\psi\bigr\rangle<0.
$$

Thus a **nonzero scalar energy shift** (or other flux-preserving even-order structure) may occur even if no TR-odd chiral mass appears at second order. This strict inequality is **conditional** on the stated spectral and nonzero-coupling assumptions, not a numerical evaluation of the Kitaev model. For `P=P_{\mathrm{ff}}` containing an interval of matter energies, one cannot generally substitute a common `E_0` into the entire retained block; the energy-dependent Feshbach map and its actual resolvent domain must be specified.

The substantive question is therefore **not** whether the *whole* second order vanishes. It is which second-order operators are allowed, which reduce to a constant or coupling renormalization, and whether the particular **chiral, TR-breaking mass operator** first arises at third order. That requires an explicit selection-rule derivation and careful treatment of intermediate matter/flux states, denominators, degeneracy and gauge/physical projection. The exact track explicitly lists these as unfinished tasks [S08,S09,S23].

### 4.1 Three-spin to matter-hopping identity: conditional, not yet globally certified

K-019 reports a local triad correspondence, schematically,

$$
P_{\mathrm{ff}}\,
\sigma_j^\alpha\sigma_k^\gamma\sigma_l^\beta
\,P_{\mathrm{ff}}
\stackrel{\text{fixed gauge}}{\longleftrightarrow}
\epsilon^{\alpha\beta\gamma} i c_j^0c_l^0,
$$

where two intervening directed bond-link operators are set to `u=+1` [S07]. This is meaningful only after the **local constraint, site and flavor ordering, intermediate link signs, orientation and gauge choice** are fixed. The expression on the right is *not by itself* a gauge-invariant physical spin operator; prior to gauge fixing it must carry the corresponding Wilson/link string. A projected operator identity does **not** derive its perturbative coefficient: different ordered virtual paths acquire different resolvent weights.

Two concrete convention discrepancies must be resolved, not silently reconciled:

1. K-019 defines `\mathcal D_k=-c_k^0c_k^xc_k^yc_k^z=b_k^xb_k^yb_k^zc_k^0=+1` [S07]; the M1.1 physical-projection note defines `\mathcal D_k=+c_k^0c_k^1c_k^2c_k^3=+1` and uses `P_{G,k}=(1+\mathcal D_k)/2` [S18,S19]. The explicit independent local projection tests also select the **positive eigenvalue of the latter product** [S25]. These are opposite parity signs under an identical displayed flavor ordering, unless a further change of basis/operator convention is made explicit. This affects the triad identity's sign and physical-sector interpretation.
2. K-019 states `N^\gamma_{0\gamma}=-1` in its Fock/KKT discussion [S07], while implemented `N^\gamma_{0\gamma}=+1` for every flavor [S20]; the independent normalization note also uses `+1` [S10]. This is a direct matrix-entry mismatch, not a verified dressing mechanism.

These discrepancies should be resolved by a separately reviewed operator-level calculation; **no sign or source is changed in this wave**.

### 4.2 Third-order structure versus coefficient

A projected path of three spin insertions can generate a TR-breaking next-nearest-neighbor matter hopping in an appropriate gauge-fixed representation. The note writes the chiral term with a coefficient `\kappa` and then replaces its full intermediate Feshbach resolvents by a **single static denominator** `\Delta E` [S07]. Such a replacement is a **physical approximation**, not a consequence of the exact triad identity. It also ignores the dispersion and matrix elements of the intermediate two-flux-plus-matter states.

In that approximation, with `\mathbf h=(h/\sqrt3)(1,1,1)`, the *note's stated convention* yields

$$
\kappa_{\mathrm{static}}=
\frac{h_xh_yh_z}{(\Delta E)^2}
=\frac{h^3}{3\sqrt3(\Delta E)^2},
\qquad
\Delta_{\mathrm{static}}=12\sqrt3\kappa_{\mathrm{static}}
=\frac{4h^3}{(\Delta E)^2}.
$$

Taking `\Delta E=0.264|K|` gives the note's approximate reference number `57.3921\,h^3/|K|^2` [S07,S16]. **This is arithmetic conditional on the chosen approximate denominator and mass/gap convention.** It is not a computed continuum coefficient, a rigorous bound, or an uncertainty estimate. The research roadmap separately lists approximate *nearest-neighbor* and *asymptotically separated* vison-pair reference energies, and requires physical-parity-corrected thresholds on matched finite systems [S09]. In particular, "isolated two-vison threshold", adjacent-pair energy and far-pair threshold must not be used interchangeably.

The exact track explicitly says "No perturbative coefficients have been implemented" and that physical matter parity, vortex-pair energies and resolvent margins remain pending [S08]. The historical note's wording that the continuum-integrated coefficient "is evaluated explicitly" [S07] therefore cannot be read as a demonstrated completed calculation at the inspected revision.

## 5. What the factor of four establishes — and does not

For the **same** antisymmetric matrix `A(\mathbf k)`, with the production convention `H_{\mathrm B}=iA/4`,

$$
\lambda_j(iA(\mathbf k))=4\,w_j(\mathbf k),\qquad
[\lambda_5-\lambda_4](iA)=4\Delta_{\mathrm B}(\mathbf k).
$$

This is a **linear-algebra identity** (positive rescaling of a Hermitian matrix), supported directly by the builder [S12] and independently derived using the real-space quadratic operator `H=(i/4)c^TAc` and the two-Majorana example [S10,S11]. The unperturbed matter dispersion test instead compares the builder eigenvalues with `\pm(|K|/2)|f(\mathbf k)|`, whereas the exact canonical positive matter frequency is `2|K||f(\mathbf k)|` in the matched convention [S11,S17].

For **the same empirical fit**, this gives `4\times2.5921=10.3684` for the corresponding **two-sided matrix eigenvalue separation**. Retaining its observed exponent, the exact rescaling statement is

$$
4\Delta_{\mathrm B,K}^{\mathrm{fit}}(h)
=4\,c_{\mathrm{fit}}|K|\,(h/|K|)^{p_{\mathrm{fit}}}.
$$

The expression `10.3684\,h^3/|K|^2` additionally **replaces** the fitted `p=2.9936` by the hypothesized asymptotic exponent 3. No analytic coefficient is obtained by that replacement.

The matrix rescaling neither (i) proves a valid nonzero-field **mean-field KKT/Fock response expansion** nor (ii) identifies a many-body *physical-spin* excitation threshold, including projection and parity. The exact-track normalization note explicitly states that its positive scaling leaves the ideal occupied projector invariant but may change numerical zero-mode classification [S11,S14]. The label "physical excitation gap" must define whether it means minimum positive canonical mode, ± band separation, fixed-K Dirac splitting, global momentum minimum, or a physical spin spectrum. A factor of four alone does not resolve these distinctions.

The bare squared energy-scale comparison `(1.0497315/0.2640)^2=15.8106` in K-019 [S07,S16] is likewise an **arithmetic ratio of model-dependent denominators**, not a determined KKT vertex correction. The note *proposes* that induced `\lambda^\gamma=O(h)` and off-diagonal Fock terms generate a chiral mass [S07,S22]; a closed-form response coefficient, analytic regularity at the gapless Dirac point, physical projection and independent comparison remain unproved. Fitting a "dressing factor" to reproduce a target ratio would be circular.

## 6. P1–P6: six **candidate, unregistered** propositions

The following partition is provisional and intentionally not the historical registry. Each candidate states a distinct type of evidence obligation; none carries an owner-approved `support_if` or `refute_if`.

| Candidate and logical form | Domain / assumptions | Existing evidence and epistemic classification | Independent challenge, potential falsifier, limitation; owner decision |
| --- | --- | --- | --- |
| **P1 — Hamiltonian convention map (identity):** `K_{\mathrm P}=-K_{\mathrm R}`, `K_{\mathrm R}=J/4`, with consistent field rescaling | Pauli `\sigma` versus spin `S=\sigma/2`, same oriented bonds, same physical spin Hilbert space | **Exact algebra**, with derivation in [S10] and independent implementation convention [S24] | Independently expand the bond operator and Zeeman term; a mismatched coefficient on one bond refutes the proposed map. **TODO(Lu):** endorse the standard comparison convention |
| **P2 — Flux-sector selection (conditional identity):** `P_{\mathrm{ff}}VP_{\mathrm{ff}}=0`; no inference that `H_{\mathrm{eff}}^{(2)}=0` | Physical Hilbert space; well-defined plaquette and Wilson sectors; separate `P_{\mathrm{ff}}` from `P_{E_0}` | **Conditional selection-rule argument**, reported in [S07], with second-order whole-operator vanishing **contradicted as an inference** | Independent finite-cluster spin operator block calculation and resolvent; a nonzero retained single-spin matrix element refutes the selection claim under the defined sector. **TODO(Lu):** specify `P`, parity, boundary/Wilson sectors |
| **P3 — Projected triad operator (conditional identity):** specified three-spin flux-return path maps to an oriented gauge-string matter hopping, with coefficient sign derived from ordering | Three-site path; nontrivial physical projector; fixed `u_{jk},u_{kl}` and sign convention | **Reported formal mapping**, not independently certified; incompatible local parity conventions appear in [S07,S18,S25] | Symbolically or exactly diagonalize a finite local/cluster CAR representation and test both sides with gauge strings; disagreement refutes the *precise* claimed identity. **TODO(Lu):** choose and verify projector/sign/orientation convention |
| **P4 — Chiral perturbation (asymptotic/conditional):** leading allowed TR-breaking effective hopping can have `\kappa\propto h_xh_yh_z`, but continuum coefficient remains to be derived | Isolated finite-volume spectral subspace or controlled Feshbach energy window; field small relative to justified sector threshold; physical parity and all intermediate states | **Approximate/partly analytical structure**, **unresolved exact coefficient and remainder** [S07,S08,S09]; static denominator is a model assumption | Enumerate all low-order projected paths and their resolvents; a nonzero lower-order *mass* or missing third-order path falsifies the proposed selection rule. **TODO(Lu):** regulator, physical sector, resolvent error/limiting regime |
| **P5 — Spectral normalization (matrix identity, interpretation conditional):** `iA=4H_{\mathrm B}` and hence fourfold scaling of *the same* spectral difference | Matched four-Majorana quadratic operator, Fourier basis and matrix `A`, unchanged self-consistent state | **Exact matrix rescaling** [S11,S12]; broader physical identification **open**; arithmetic test [S16] is not a spectrum test | Independently form real-space `A`, compare canonical frequencies and Bloch matrices at generic `k` and spin-parity reference; failure of eigenvalue rescaling refutes algebra, mismatch of *observables* refutes an interpretation. **TODO(Lu):** which gap quantity the comparison should target |
| **P6 — MF cubic gap and physical relevance (empirical / predictive):** finite-window `\Delta_{\mathrm B,K}\simeq c h^p`; separate conjecture `\lim_{h\to0}\Delta_{\mathrm B,K}/(h^3/|K|^2)=c_0>0` | `K_{\mathrm P}=-1`, QSL continuation, sampled `24×24` mesh, `h\in[0.005,0.150]` for the regression; *additional* regularity/physical-sector assumptions for asymptotic or spin prediction | **Saved numerics + fitted parameters** [S01,S04,S05]; **not** proved asymptotics, global BZ gap or predictive validity | Independent full-BZ minimization, fit-window/mesh/solver sensitivity, independent small-system ED/controlled exact-spin perturbation and variational checks; a lower off-path gap or changed fit limits a stronger claim. **TODO(Lu):** observables, metrics, fit uncertainty, decisive comparison criterion |

**Nonvacuity and relevance:** P1 requires an actual common Hamiltonian map; P2–P4 require an admissible physical finite system and a nonempty retained subspace with a justified resolvent; P5 requires the *same* assembled antisymmetric matrix; P6 requires a converged QSL branch and a specified nonzero field window. A plausible symbolic formula in an empty sector or a regression over a selected path is not self-validating. The target scientific relevance is specifically the low-field **physical Kitaev spin Hamiltonian**, not merely spectral correctness of an enlarged Gaussian mean-field model.

## 7. Three independent validity axes for this pilot

| Axis | K-010 / MF evidence | K-019 / perturbative evidence | What is still required |
| --- | --- | --- | --- |
| **Computational validity** | Recorded 17-point dataset, converged flags and a fixed-point routine; spectral builder and direct regression observable are identifiable [S04,S05,S12,S13] | One unit test checks only hard-coded denominator arithmetic and `4c_{\mathrm{fit}}` [S16] | Reproduce fit from documented data with honest uncertainty; full-BZ and finite-mesh/precision study; independent operator-level test, exact-resolvent calculation |
| **Admissibility / stationarity** | The solver measures `\Delta U,\Delta V` and averaged gauge residual and a rebuilt-matrix commutator [S13]; Gaussian energy functional exists [S15,S26] | Physical parity and flux projectors remain unresolved for the exact-spin comparison [S08,S09] | Verify covariance representability, physical-projection assumptions, full required stationarity/KKT conditions, and distinguish stationary from lowest-energy admissible state |
| **Predictive validity** | No direct implication for the physical spin gap from a self-consistent Gaussian band fit | Static denominator is approximate, exact continuum coefficient and parity constraints absent [S07,S08] | Independent physical spin ED and matched exact limits/resolvent; error budget separating numerical from model-form discrepancy |

In particular, the implemented `compute_commutator_stationarity` computes a **BZ-averaged Frobenius norm** of a projector formed from the input Hamiltonian against a Hamiltonian rebuilt from output correlators [S13]. The code's projector routine actually returns `2\Pi^-` (with half occupations at zero modes), not a raw rank-`N` projector [S14]. A small commutator can be a useful necessary diagnostic, but it neither proves `\Phi(X)=X` in every coordinate nor a globally minimizing, physically projected state. The energy implementation includes Lagrange-multiplier terms [S15]; one must not confuse a stationary point of that extended functional with a variational bound in the original physical spin Hilbert space.

## 8. Candidate independent tasks, ranked by evidential leverage

These are **future proposals only**. No numerical/solver rerun, code/test changes, hypothesis registration or scientific promotion is authorized by this note.

| Priority | Separately verifiable task | Input and independent oracle | Scope / feasibility and expected new evidence |
| --- | --- | --- | --- |
| **1 — blocking** | Fix a single physical projector, parity convention, directed-bond signs and spin/Pauli mapping; derive the triad identity independently | Local CAR matrices, spin-space small cluster and gauge-string tracking [S18,S19,S20,S23] | Finite-volume symbolic/ED check is feasible; resolves opposite-parity and `N` sign assumptions without guessing a coefficient |
| **2 — blocking** | Derive `PVP`, `H_{\mathrm{eff}}^{(2)}` and all allowed third-order *operators* with stated `P/Q` and resolvent domain | Finite-cluster physical-spin ED and independent resolvent block calculation [S08,S09,S23] | Separates a nonzero scalar second-order shift from the absence of a mass; prerequisite to coefficient claims |
| **3 — high** | Reconcile gap definitions and the factor-four mapping with the canonical real-space quadratic generator, including ± splitting versus positive frequency | Independent finite-cell antisymmetric `A`/singular-value comparison [S10,S11,S17] | Narrow linear-algebra test; distinguishes actual frequency normalization from projected many-body gap |
| **4 — high** | Analyze zero-field gaplessness and the response/analyticity needed for a leading cubic MF mass; derive KKT/Fock terms with independent signs | Explicit response equations, representation-level symmetry and projected comparator [S07,S13,S15,S20,S21] | Nontrivial analytic problem; cannot infer `c_{\mathrm{MF}}` from numerical agreement |
| **5 — high** | Compute the *physical* finite-system flux thresholds and a controlled two-flux-plus-matter continuum resolvent | Finite physical sector and parity-aware exact-spin ED [S08,S09,S23] | Demanding but decisive for exact-spin coefficient; only later discuss thermodynamic control |
| **6 — medium** | Audit saved fit error bars, branch selection, full BZ minimum, fit-window and mesh scaling independently | Archived `.dat` plus separately specified numerical protocol [S01,S04,S05] | Feasible independent numerical check after authorization; cannot turn an empirical fit into proof |
| **7 — eventual** | Compare matched physical-spin and projected/unprojected MF observables; constrain model-form error | Spin-space ED/analytic limits with explicit sectors and normalization | Required for predictive claims; no blanket equivalence from matrix scaling |

**Practical sequencing:** resolve *P2–P3 projection/selection* and *P5 normalization* as small independent packages before attempting P4 exact continuum work or promoting P6 to a physical prediction. This ordering follows logical dependencies, not the order of the repository's milestone labels.

## 9. Decisions reserved to the scientific owner

All items below remain open; the label explicitly denotes **no authorization or scientific acceptance**.

1. **TODO(Lu):** The exact physical Hilbert space, flux/Wilson sector and definition of `P` (flux-free full sector versus degenerate ground-energy subspace), boundary conditions and matter parity for perturbation theory.
2. **TODO(Lu):** Choose and reconcile the physical projector and local Majorana orientation (opposite `\mathcal D` signs), the `N^\gamma_{0\gamma}` sign conflict, and gauge-string/link orientation before treating P3 as an operator theorem.
3. **TODO(Lu):** Which observable is the intended "gap": fixed-K Bloch ± splitting, canonical positive quasiparticle threshold, global BZ minimum or projected physical-spin excitation gap; how its finite-size limit is defined.
4. **TODO(Lu):** Which exact third-order *operator*, energy-dependent Feshbach map, resolvent domain and allowed intermediate sectors constitute the physical target; how scalar second-order shifts are separated from the chiral mass and what remainder is required.
5. **TODO(Lu):** Whether and under which independent evidence conditions the original symmetry-based `O(h^3)` and `O(h^5)` assertions may be promoted; the analytic regularity and gauge/matter spectral assumptions in the gapless regime.
6. **TODO(Lu):** Statistical interpretation/source of the stated fit stderr, finite-window uncertainty, solver and mesh tolerances, full BZ minimization criterion, independent branch selection and later `support_if/refute_if` rules.
7. **TODO(Lu):** Which `\Delta_{\mathrm{flux}}` (adjacent pair, far pair, physical parity-corrected, finite or thermodynamic) supplies an admissible static approximation, if any; no inferred bare-to-dressed factor.
8. **TODO(Lu):** Scope and cost of an independent finite-cluster spin ED comparison; thresholds for calling a mean-field gap *predictive* rather than only computationally consistent.
9. **TODO(Lu):** Which, if any, P1–P6 formulations should eventually become **new registered** hypotheses after independent review; no v2 adoption or retroactive record updates here.

## 10. Evidence-level conclusion and W4 challenge

**Supported at source level:** the existing scripts define a particular finite-mesh, fixed-K Bloch eigenvalue separation and produce a committed set of converged output rows; a two-sided Bloch spectral separation rescales by four when comparing `iA/4` to `iA` for the same matrix. An exact Pauli/spin coupling conversion is documented. These are much narrower than an exact physical-spin-model cubic-gap coefficient.

**Not established by the reviewed evidence:** a full-BZ gap minimum, rigorous fit standard errors, a controlled `h\to0` asymptotic exponent or prefactor, vanishing of the entire second-order projected operator, a uniquely reconciled physical projector/sign, a continuum-resolvent coefficient, a closed-form MF KKT/Fock coefficient, or quantitative predictive agreement with the physically projected spin Hamiltonian.

**Independent W4 reviewer questions:** (i) Can every operator and coefficient be reconstructed from the *same* convention and projection? (ii) Does the source really prove an asserted vanishing term rather than only an absent *chiral* contribution? (iii) Does the regression's quantity match the reported physical observable and global/minimum qualifier? (iv) Which numeric uncertainties were computed rather than copied? (v) Could an independent test have exposed a false result? (vi) Has any conditional structure been presented as a concluded scientific result? (vii) Would this evidence analysis change a real acceptance decision compared to the existing governance, or is it redundant prose?

This note remains an unregistered, independently reviewable research-methodology pilot. **Stop after Wave 2 draft PR and coordinator handoff.**

## Source index — fixed Git revision

[S01]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/PROJECT_STATE.yaml#L290-L306
[S02]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/PROJECT_STATE.yaml#L416-L426
[S03]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/HYPOTHESIS_CONTRACTS.md#L19-L28
[S04]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/scripts/m1_6/benchmark_majorana_gap.py#L27-L164
[S05]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/data/processed/benchmarks/m1_6/majorana_gap_scaling.dat#L1-L24
[S06]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/rnd/phase_1/m1_6_t5_critical_zoom_and_supercells.md#L58-L88
[S07]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/rnd/phase_1/m1_10_t1_exact_perturbation_and_normalization_reconciliation.md#L25-L116
[S08]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/rnd/exact-perturbation-theory/STATUS.md#L3-L113
[S09]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/rnd/exact-perturbation-theory/ROADMAP.md#L17-L79
[S10]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/rnd/exact-perturbation-theory/notes/02_normalization_and_bloch_conventions.md#L1-L89
[S11]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/rnd/exact-perturbation-theory/notes/02_normalization_and_bloch_conventions.md#L247-L329
[S12]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/src/kitaev_mf/services/hamiltonian.py#L118-L192
[S13]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/src/kitaev_mf/services/solver.py#L62-L160
[S14]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/src/kitaev_mf/services/correlators.py#L11-L47
[S15]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/src/kitaev_mf/services/energetics.py#L13-L81
[S16]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/test/unit/test_post_spinodal_kkt_audit.py#L25-L38
[S17]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/test/unit/test_hamiltonian.py#L78-L123
[S18]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/rnd/phase_1/m1_1_t3_physical_projection_pauli_algebra.md#L24-L47
[S19]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/rnd/phase_1/m1_1_t3_physical_projection_pauli_algebra.md#L135-L171
[S20]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/src/kitaev_mf/domain/constants.py#L81-L109
[S21]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/labbook/chapters/ch_m1_6_solvers_and_basins.tex#L245-L283
[S22]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/rnd/phase_1/t01_epistemic_demarcation_and_gap_reconciliation.md#L24-L70
[S23]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/rnd/exact-perturbation-theory/BACKLOG.md#L28-L63
[S24]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/src/kitaev_mf/domain/parameters.py#L19-L74
[S25]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/test/unit/test_physical_projection.py#L54-L100
[S26]: https://github.com/lufreitas0000/kitaev-mf/blob/4c2d629ec7ee924c731c566fc00d9dc0411217c6/doc/rnd/phase_1/m1_2_t1_majorana_mean_field_decoupling.md#L19-L85
