# KKT Majorana scMFT — physical validation checklist

**Status:** proposed evidence audit, not scientific acceptance.  
**Baseline:** [2971ad21ce72](https://github.com/lufreitas0000/kitaev-mf/tree/2971ad21ce7207bcd1be5a4e44456be9a9d28886).  
**AI contribution:** ChatGPT (OpenAI), 2026-11-25; static repository review and planning only; no tests executed or independent approval implied.

## Objective and decision rule

Determine whether existing KKT-constrained four-Majorana Gaussian scMFT describes the physical Kitaev spin liquid in a field, quantitatively and with the correct fractionalization. Retain KKT until a controlled benchmark reveals a physical deficiency. A nonzero gauge Casimir of the auxiliary Gaussian parent is not, by itself, a failure: even the exact zero-field fixed-gauge parent need not be physical before projection. Exact local projection removes gauge-variant bilinear expectations but need not remove gauge-invariant four-point correlations. Do not feed the projected covariance into the unchanged Wick update.

A proposed change must identify the failed benchmark, root cause (normalization, finite-size/solver error, extension ambiguity, variational manifold or observable evaluation), and demonstrate improvement on matched controls without degrading accepted benchmarks. No new penalty or projection optimization is authorized by this checklist.

## Model and metrology

Start with the pure isotropic FM honeycomb model in Pauli units, K=-1, J=Gamma=Gamma'=0, T=0; field magnitude h along [111] with components h/sqrt(3), Zeeman term minus the field dotted into Pauli operators. Verify exact repository sign and normalization bridge in #25 before comparing historical results. Start with h=0, 0.001, 0.002, 0.005, 0.01; extend to the existing HC-01 grid where necessary. The susceptibility at zero field is generally **nonzero**, although spontaneous magnetization vanishes.

Primary translationally invariant control: L=20 (800 spins); L=21 (882 spins) tests the Dirac-commensurate mesh separately. Converge over L=24,36,48,72 or larger and noncommensurate controls as justified. Evaluate the Dirac point directly rather than inferring a tiny O(h^3) mass from a coarse mesh minimum. Distinguish a mass, a positive eigenenergy and a band separation; report canonical units. Match geometries, field definitions, twists, parity and flux sectors for direct finite-torus comparisons. Use the existing smaller ED tori for matched oracle tests, not as substitutes for bulk results. At h=0.01 the nominal h^3/|K|^2 scale is 1e-6, not the actual gap coefficient.

Each box is unchecked pending the new audit. For each ID record PASS / FAIL / PARTIAL / NOT TESTED / OUT OF SCOPE, exact scope, existing or new artifact path, pinned commit, command, inputs, tolerances and result. Separate exact algebra/proof, auxiliary Gaussian results, projected physical finite-size evidence, perturbative approximations and extrapolations. A passed auxiliary test is not automatically a passed physical test. Preserve old evidence and failures. Pre-register new Phase 2+ quantitative contracts per [registry](../../HYPOTHESIS_CONTRACTS.md). Follow [ADR-003](../../adr/ADR-003-linear-main-and-family-boundaries.md): exact spin ED/SDP and exact perturbative/dynamical oracles belong to sibling repos; ingest their verified data, not their implementations.

## A — Algebra, constraints and representation

- [ ] A01 CAR and physical spin algebra in locked conventions.
- [ ] A02 Local gauge-generator algebra, signs and normalization.
- [ ] A03 Casimir identity and equivalence of zero Casimir with physical parity; compare independent CAR oracle.
- [ ] A04 Equivalence of enlarged-space spin extensions on the physical subspace, including exchanges.
- [ ] A05 Zeeman-extension invariance: for fixed exchange extension, unrestricted KKT multipliers absorb any added linear gauge-generator term. Check restricted ansatz/multipliers, matched branches and observable definitions; do not assume arbitrary exchange-extension invariance.
- [ ] A06 Gauge covariance of ansatz and self-consistency; distinguish local gauge transformations from physical lattice symmetries.
- [ ] A07 Constrained energy stationarity, correct tangent gradient, dual residual, distinction between stationarity, sign-map fixed point and stable minimum.
- [ ] A08 Wick contraction, energy normalization and double counting against independent Fock oracle.

## B — Exact zero-field limit

- [ ] B01 Ground-state energy per site with consistent normalization and sectors.
- [ ] B02 Nearest-neighbor spin correlations.
- [ ] B03 Vanishing spin correlations beyond nearest neighbors in the exact model.
- [ ] B04 Zero spontaneous magnetization (not zero susceptibility).
- [ ] B05 Matter-Majorana dispersion over the Brillouin zone.
- [ ] B06 Dirac-point positions and velocity.
- [ ] B07 Flux-free ground-state sector, full plaquette operator versus factorized proxy.
- [ ] B08 Physical two-vison excitation energy.
- [ ] B09 Vison interaction versus separation.
- [ ] B10 Wilson-loop sectors and boundary dependence.
- [ ] B11 Physical matter-fermion parity, zero-projection sectors and correct finite-size ground state.
- [ ] B12 No spurious symmetry breaking; distinguish ansatz restrictions from independent tests.

## C — Weak-field response

- [ ] C01 Physical linear susceptibility tensor and absolute coefficient from independent zero-field spectral response.
- [ ] C02 Second-order energy correction, consistent with susceptibility.
- [ ] C03 Hellmann–Feynman/thermodynamic consistency m=-de/dh with numerical uncertainty.
- [ ] C04 Cubic leading chiral mass scaling, distinguish auxiliary and physical effective mass.
- [ ] C05 Absolute coefficient against full virtual-flux resolvent and physical parity; do not use constant-denominator approximation as exact.
- [ ] C06 Angular dependence proportional to h_x h_y h_z at leading order, including component-zero controls.
- [ ] C07 Signed mass/chirality reversal under field reversal.
- [ ] C08 Leading field renormalization of Dirac positions/velocities where symmetry permits.
- [ ] C09 Leading change in physical plaquette flux expectation.
- [ ] C10 Field-induced changes in spin correlations.
- [ ] C11 Absence of lower-order chiral mass or other symmetry-forbidden contributions.

## D — Topology and fractionalization

- [ ] D01 Auxiliary Majorana Chern number in a gapped band regime.
- [ ] D02 Berry-curvature integration and mesh convergence; no Chern assignment at an unresolved gap closing.
- [ ] D03 Chiral Majorana edge spectrum with geometry and boundary controls.
- [ ] D04 Vison-bound Majorana modes in the weak-field non-Abelian phase.
- [ ] D05 Localization length and finite-size effects.
- [ ] D06 Splitting/hybridization of separated vison-bound modes.
- [ ] D07 Physical torus ground-state sectors and parity; thermodynamic threefold degeneracy in the gapped Ising phase.
- [ ] D08 Physical topological entanglement entropy, gamma=ln 2 in the gapped Ising phase; no inference from Gaussian parent alone.
- [ ] D09 Chiral central charge |c_-|=1/2 and low-T thermal Hall coefficient in the gapped phase; distinguish auxiliary prediction from physical response.
- [ ] D10 Stability and phase boundaries under increasing field; benchmark against independent physical evidence.

## E — Dynamical response (later; observable evaluation may need extension)

- [ ] E01 Exact zero-field dynamical spin structure factor.
- [ ] E02 Physical spin threshold versus gapless zero-field matter spectrum; resolve sector selection.
- [ ] E03 Spectral-weight redistribution under field.
- [ ] E04 Equal-time and frequency-integrated sum rules with correct normalization.
- [ ] E05 Static susceptibility from dynamical response.
- [ ] E06 Distinguish matter single-particle gaps from physical excitation gaps.
- [ ] E07 Field-induced vison dynamics.
- [ ] E08 Long-distance field-induced spin correlations.

## F — Numerical reliability and stability

- [ ] F01 Self-consistency residual, independent tolerances.
- [ ] F02 Local gauge first-moment residual.
- [ ] F03 Eigensolver versus matrix-sign covariance.
- [ ] F04 Covariance antisymmetry, CAR representability and purity.
- [ ] F05 Seed dependence within matched branch, including symmetry-unrestricted seeds where relevant.
- [ ] F06 Constrained stability/Hessian versus merely stationary solutions.
- [ ] F07 Competing metastable branches and spinodal status.
- [ ] F08 Momentum mesh convergence of physical and auxiliary quantities.
- [ ] F09 Dirac-commensurate versus noncommensurate meshes.
- [ ] F10 Direct Dirac-point evaluation.
- [ ] F11 Precision/residual sensitivity at cubic mass scale.
- [ ] F12 Real-space/Bloch agreement at same finite size and twists.
- [ ] F13 Independent implementation and repeatability, failures retained.
- [ ] F14 Computational scaling with system size and solver resources.

## G — Independent many-body evidence

- [ ] G01 Matched spin ED ground-state energies; no general upper-bound claim for unprojected Gaussian energies.
- [ ] G02 Matched magnetization/susceptibility.
- [ ] G03 Matched physical spin correlations.
- [ ] G04 Matched plaquette flux distribution.
- [ ] G05 Matched excitation spectra and sector quantum numbers.
- [ ] G06 Compare global projected Gaussian states with the same optimized parent and controls.
- [ ] G07 Physical overlap, energy variance and absolute projection weight; account for statistical errors and zero-projection sectors.
- [ ] G08 Compare independent DMRG/tensor-network evidence with geometry/extrapolation caveats.
- [ ] G09 High-field transition against independent physical data; distinguish KKT branch crossings from physical transitions.
- [ ] G10 High-field asymptotic energy and response.

## Existing evidence crosswalk (initial reconnaissance, not verdicts)

| Tests | Existing sources to inspect and re-run | Qualification |
|---|---|---|
| A01–A04, A08, F04 | [spin/gauge note](../phase_1/m1_1_t2_spin_gauge_representation.md), [observable tests](../phase_1/m1_4_t2_observable_test_specification.md), test/unit/test_local_majorana_algebra.py, test/unit/test_local_majorana_representations.py, test/unit/test_physical_projection.py, test/unit/test_quasifree_states.py | Existing Gaussian Fock-oracle tests explicitly do not establish projected physical equivalence; #25 and #27 govern convention/foundational re-audit. |
| A07, C02–C03, F01–F07 | [KKT stationarity note](../phase_1/m1_4_5_t2_kkt_stationarity_and_response_asymptotics.md), test/integration/test_variational_response.py, test/unit/test_response_and_thermodynamic_invariants.py, test/unit/test_post_spinodal_kkt_audit.py | Recheck dimensions, tangent projection, branch stability, sign/normalization and derivative conventions; historical solver evidence is not independent physical validation. |
| B01–B06, B12 | test/unit/test_solvable_limits.py, test/unit/test_k006_quadrature.py, test/unit/test_symmetry_selection_rules.py, data/processed/benchmarks/dispersion_path.dat, doc/HYPOTHESIS_CONTRACTS.md | Candidate translationally invariant evidence; audit test scope, exact reference and finite-size parity. |
| B07–B11, D04–D06 | [HC-03 vison scaling](../phase_2/wave1_hc03_multi_l_vison_scaling.md), test/unit/test_wave1_visons.py, test/unit/test_real_space_vortices.py | #12, #15, #27; distinguish physical flux sectors, gauge-fixed proxies, constrained relaxation and physical excitations. |
| C01–C07, F08–F11 | [HC-01 cubic gap](../phase_2/wave1_hc01_cubic_gap_kkt_reconciliation.md), scripts/wave1/run_hc01_cubic_gap.py, data/processed/benchmarks/susceptibility_curve.dat, data/processed/benchmarks/hellmann_feynman_scaling.dat | HC-01 reports canonical auxiliary gap/h^3=10.641666 and exponent 2.999694, not an independently established physical prefactor; h=0.002 at N_k=48 flagged unconverged. #30 addresses independent coefficients. |
| F12–F14 | test/unit/test_real_space_and_purification.py, test/unit/test_hardware_parity.py, test/integration/test_m1_system_convergence.py, test/unit/test_wave4_trials.py | Existing scaling and paired trials; audit matched conditions and raw evidence. |
| G01–G07, B11 | [HC-02 projected/ED/SDP sandwich](../phase_2/wave1_hc02_ed_sdp_gutzwiller_sandwich.md), scripts/wave1/run_hc02_variational_sandwich.py, test/unit/test_wave1_cubic_and_gutzwiller.py | Reported N=8,12,18 projected results exist; #14 and #27 address matched finite-torus and projection qualifications. VMC estimates are statistical, not deterministic bounds. |
| D01–D03, D07–D10, E01–E08, G08–G10 | [milestone roadmap](../../MILESTONE_ROADMAP.md), [hypothesis contracts](../../HYPOTHESIS_CONTRACTS.md), doc/rnd/phase_3_dynamical_spin_correlations/theory/09_exact_vs_mean_field_benchmark.md | Inventory existing work before new tasks; no pass implied. |

## Execution order and issue ownership

1. [#29](https://github.com/lufreitas0000/kitaev-mf/issues/29): inventory and re-audit translationally invariant evidence (A/B/F and existing C/G evidence); report gaps before implementation.
2. [#30](https://github.com/lufreitas0000/kitaev-mf/issues/30): independent weak-field susceptibility and cubic mass coefficient after #29; coordinate with #12 for exact perturbation reference.
3. Existing [#25](https://github.com/lufreitas0000/kitaev-mf/issues/25) and [#27](https://github.com/lufreitas0000/kitaev-mf/issues/27): conventions, projection and foundational claims; no duplicated Casimir project.
4. Existing [#14](https://github.com/lufreitas0000/kitaev-mf/issues/14) and [#15](https://github.com/lufreitas0000/kitaev-mf/issues/15): matched ED and flux/defect sectors. D/E/G follow only after evidence inventory and appropriate reference availability.

The three initial scientific questions are Zeeman representation invariance, absolute susceptibility, and absolute cubic mass coefficient. Defer new theory until an actual deficiency survives controls. This is a planning record, not a change to PROJECT_STATE, milestone commitments or locked hypothesis thresholds.
