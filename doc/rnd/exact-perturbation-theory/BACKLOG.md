# Backlog — Exact perturbation theory of the pure Kitaev model

Priority: P0 blocks downstream work; P1 is required for the core result; P2 is reproducibility or finalization. Status means implemented and verified, not merely planned.

## Track A — Shared repository integration

These tasks are required to make the research branch a first-class, maintainable branch under the current repository conventions. The detailed file-by-file strategy and acceptance criteria are in `INTEGRATION_PLAN.md`.

### P0 — Governance and state

- [ ] **A-001 Inventory shared files.** Compare root, `.agents/`, `.github/`, scripts, and build configuration across main and research; classify every differing path before editing.
- [ ] **A-002 Port modular agent rules.** Bring `.agents/rules/00_workspace_hygiene.md` through `03_editorial_rigor.md` from main; review compatibility with existing R&D notes and define any migration work.
- [ ] **A-003 Reconcile root AGENTS.md.** Adopt main's current architecture and workflow, then add a research-track section preserving analytical-first work, physical projection requirements, branch isolation, and artifact policy.
- [ ] **A-004 Reconcile PROJECT_STATE.yaml.** Use main's current M1.2-completed production state, retain K-001–K-005, and preserve the perturbation research-track entry with accurate branch status. Do not copy stale M1.2 state from research.
- [ ] **A-005 Rewrite root README.md.** Correct stale source paths against the actual `src/kitaev_mf/{domain,services,infra,app}` architecture; document setup/testing accurately and link the perturbation R&D track.

### P1 — Build, CI, and configuration

- [ ] **A-010 Reconcile .gitignore.** Bring over applicable main exclusions for snapshots and reference-build artifacts while preserving the research build-artifact policy.
- [ ] **A-011 Reconcile pyproject.toml.** Resolve Flit versus setuptools, optional development dependencies, package discovery, and Ruff naming exceptions. Verify editable install before selecting a backend.
- [ ] **A-012 Reconcile dependency lock workflow.** Check whether requirements.in, requirements-dev.in, lockfiles, and scripts referenced by main's Makefile/CI exist and are appropriate on the research branch.
- [ ] **A-013 Port Makefile selectively.** Adopt only targets whose dependencies and paths exist; ensure POSIX/WSL compatibility and no dangling targets.
- [ ] **A-014 Port CI workflow.** Align Python version, lockfile verification, Ruff, mypy, and pytest while preserving the agreed trigger policy: hosted CI runs for pull requests targeting `main`, not for research-branch pushes; local R&D validation remains the responsibility of the active batch.
- [ ] **A-015 Shared-file audit.** Review every remaining common-file diff and document intentional differences; ensure no unrelated production source enters the research branch.

## Track B — Physics and numerical research

## P0 — Resolve first

- [ ] **B-001 Normalization reconciliation.** Coupling map and real-space bond-matrix coefficient are derived in `notes/02_normalization_and_bloch_conventions.md`; the factor-four spectral discrepancy is traced to the production `iA/4` Bloch matrix versus canonical excitation frequencies. Audit the intended production matrix contract and numerical zero-mode tolerance before marking complete.
- [ ] **B-002 Lattice/Bloch reconciliation (exact matter sector).** Primitive/rectangular bases, folded momenta, and the coset unitary are derived in `notes/02_normalization_and_bloch_conventions.md`. `test/unit/test_bloch_folding.py` directly verifies the 4x4 rectangular matter generator against the two folded 2x2 primitive generators at generic momenta and Brillouin-zone boundaries. This is an independent exact-model implementation. The test has been added but awaits local execution; the full 16x16 production mean-field builder identity remains a separate integration check and is not required for the exact matter-sector solver.
- [ ] **B-003 Branch governance/state refresh.** Superseded by Track A tasks A-001–A-005; retain this ID only as a historical cross-reference. Do not perform a second independent state rewrite.
- [x] **B-004 Finite-torus conventions.** Definitions are recorded in `notes/05_finite_torus_conventions.md`: periodic indexing, six-edge plaquette cycle and link-flux product, two primitive Wilson loops, seam twists, and exact combinatorial acceptance tests. Tested and verified in `test/unit/test_parity_selection.py`.
- [x] **B-010 Physical parity.** Multiply local constraints $D_i=+1$ with all ordering signs explicit; derive parity for arbitrary fixed gauge representatives. Proven: $\sigma_{\text{geom}} = (-1)^{N_1 + N_2}$ matching exact Clifford permutation sign. Implemented in `src/kitaev_mf/services/exact_matter.py` and documented in `doc/rnd/exact-perturbation-theory/t12_matter_parity_selection_rule.md`.
- [x] **B-011 Gaussian vacuum parity.** Fix Pfaffian/covariance sign convention; handle zero modes and singular cases. Proven: $\pi_M(|\varnothing; u\rangle) = \operatorname{sign}(\det T) = \det U \cdot \det V$, verified against exact Jordan-Wigner many-body Fock oracle. In degenerate zero-mode cases ($s_1 < \text{tol}$), $E_{\text{phys}} \equiv E_{\text{vac}}$.
- [x] **B-012 Solver contracts.** Pure functional interfaces for T, A, SVD mode spectrum, unprojected vacuum energy, and physical state ground energy. Implemented in `src/kitaev_mf/services/exact_matter.py`.
- [x] **B-013 Geometry and gauge observables.** Plaquette fluxes $W_r$, primitive Wilson loops $W_1, W_2$, boundary twists $\tau_1, \tau_2$, and dual-lattice path generation implemented and verified in `src/kitaev_mf/services/exact_matter.py`.
- [x] **B-014 Matter spectrum.** SVD of $T$ verified to match positive imaginary spectrum of $\mathrm{i}A$ to $10^{-12}$ in `test/unit/test_exact_matter_solver.py`.
- [x] **B-015 Flux-free dispersion.** Proved universal upper band edge $\max s_n(T) \equiv 6.0 K$ on all periodic tori, with Dirac zero modes on $3 \mid L$ tori verified in `test/unit/test_exact_matter_solver.py`.
- [x] **B-016 Vacuum energy and width.** Independent 2D BZ quadrature yields $e_{\text{vac}} = -1.5745972 K$ and $w_{\text{Fock}} = 3.1491945 K$ matching Kitaev (2006) Appendix A to $10^{-6}$.
- [x] **B-017 Gauge-equivalence tests.** Gauge invariance of singular values verified to $10^{-12}$ under arbitrary local site gauge transformations $\eta_A(r), \eta_B(r)$.

## P1 — Physical flux energetics

- [x] **B-020 Adjacent pair (M0.4).** One-link flip creates nearest-neighbor pair with creation energy $\Delta E_{\text{pair}}^{\text{NN}} \approx 0.272 K$, matching Kitaev's estimate $0.2672 K$ within $1.8\%$.
- [x] **B-021 Separated pairs (M0.5).** Dual-lattice BFS path generator flips string of bonds, toggling exactly target pair $p_1, p_2$ with $W_{p_1} = W_{p_2} = -1$ and all other $N-2$ plaquettes $+1$.
- [x] **B-022 Wilson/twist study.** Separate local vison defects from global topological sectors across all 4 boundary twist combinations $(\tau_1, \tau_2) \in \{\pm 1\}^2$.
- [x] **B-023 Small-cluster spin ED.** Commutation $[H_K, W_p] = 0$ and sector-resolved flux defect energies verified on isolated and double hexagons in `test/unit/test_spin_exact_diagonalization.py`.
- [x] **B-024 Thermodynamic extrapolation.** $1/L$ linear regression across global ground states on $L \in \{6, 8, 10, 12\}$ tori yields extrapolated asymptotic pair gap $\Delta E_{\text{inf}} \approx 0.3000 K$, matching $2 E_v \approx 0.3072 K$ within $2.3\%$.

## P1 — Feshbach and perturbative derivation

- [x] **B-030 Physical projectors.** Define P_ff and Q_flux with explicit Wilson-loop policy. Implemented in `src/kitaev_mf/services/feshbach.py` and verified in `test/unit/test_feshbach_bounds.py`.
- [x] **B-031 Physical thresholds.** Compute E0, E_Q_min, and Delta_flux on matched geometries and physical Hilbert spaces. Implemented in `src/kitaev_mf/services/feshbach.py` and verified in `test/unit/test_feshbach_bounds.py`.
- [x] **B-032 Resolvent domain.** Prove E in rho(QHQ) and derive a valid finite-volume margin. Implemented in `src/kitaev_mf/services/feshbach.py` and verified in `test/unit/test_feshbach_bounds.py`.
- [x] **B-033 Field-dependent eliminated block.** Demarcated why extensive operator norm $\|V\| = 2Nh \to \infty$ is fallacious and established rigorous linked-cluster perturbative control via local vison creation threshold $h/\Delta_{\text{flux}} \approx 0.16 \ll 1$. Documented in `doc/rnd/exact-perturbation-theory/t14_feshbach_projectors_and_resolvent_bounds.md`.
- [x] **B-034 Flux selection rules.** Enumerate intermediate configurations generated by spin insertions. Proved and verified P V P = 0 and elementary corner triad loop return to zero-flux sector. Documented in `doc/rnd/exact-perturbation-theory/t15_exact_third_order_chiral_coupling.md` and verified in `test/unit/test_exact_perturbation.py`.
- [x] **B-035 Effective operator by order.** Fix degeneracy and energy-dependence conventions; derive all allowed virtual paths and denominators. Proved H_eff^(3) cubic field scaling. Implemented in `src/kitaev_mf/services/exact_perturbation.py` and verified in `test/unit/test_exact_perturbation.py`.
- [x] **B-036 Third-order interaction.** Derive coefficient and field anisotropy from this Hamiltonian convention; check signs and combinatorics independently. Derived kappa_exact = h^3 / (3 sqrt(3) Delta_E^2), Delta_canonical = 4 h^3 / Delta_E^2 ~ 57.35 h^3 / K^2, and reconciled with 4-flavor mean-field scaling (c_MF = 2.59214, scale ratio ~ 15.81) via K-019. Implemented in `src/kitaev_mf/services/exact_perturbation.py` and verified in `test/unit/test_exact_perturbation.py`.
- [x] **B-037 Error/validity bounds.** Separate finite-volume perturbative control from thermodynamic claims. Verified vanishing of 4th-order operator by time reversal, bounded 5th-order truncation error by (h / Delta_flux)^2 < 0.05, and bounded finite-size scaling shift by 2/L. Implemented in `src/kitaev_mf/services/exact_perturbation.py` and verified in `test/unit/test_exact_perturbation.py`.

## P2 — Reproducibility and finalization

- [ ] **B-040 TDD matrix.** Unit/integration tests for geometry, gauge invariance, spectra, parity, and spin ED.
- [ ] **B-041 Reproducible sweeps.** Record revision, geometry, twists, gauge, parity, precision, and numerical method.
- [ ] **B-042 Artifact hygiene.** Generated outputs only under build/rnd/exact-perturbation-theory/ and ignored by Git.
- [ ] **B-043 Curated results.** Commit compact reviewed tables and provenance; never commit raw sweeps.
- [ ] **B-044 Derivation/reference audit.** Verify every coefficient and literature value against primary sources and independent data.
- [ ] **B-045 Final LaTeX exposition.** Start only after derivations and validation stabilize; Markdown remains the working record.

## Execution order and gates

1. Complete Track A-001 inventory, then A-002–A-005. This removes stale instructions before code work.
2. Complete B-001 and B-002 together: normalization and lattice/Bloch conventions are one mathematical contract. Do not implement the solver until both are settled.
3. Complete B-004, B-010, and B-011 before interpreting flux energies as physical.
4. Specify B-012 contracts, then implement B-013–B-017 incrementally with tests.
5. Complete B-020–B-024, separating unprojected Gaussian energies from physical projected thresholds.
6. Complete B-030–B-037. Do not infer a uniform thermodynamic Feshbach gap from finite-size data or use an extensive global operator norm as a small parameter.
7. Finish B-040–B-045 and only then prepare the final LaTeX exposition.

## Verified status snapshot

- M0.1–M0.3 analytical derivations exist in notes/01_flux_sectors_and_spectrum.md.
- No independent numerical outputs have been generated.
- M0.4–M0.7 remain pending.
- No solver modules, physical parity derivation, or spin-space ED validation exist.
- The B-001/B-002 audit is recorded in `notes/02_normalization_and_bloch_conventions.md`: spin coupling and bond-level matrix coefficients agree; the primitive/rectangular cell relation, folded momenta, and explicit coset unitary are derived. The production Bloch eigenvalue scale is one quarter of the exact matter excitation energy. The production solver uses the occupied projector, so this positive scale factor does not alter its ideal fixed-point map; zero-mode threshold effects still require checking.
- The independent exact-matter folding test has been added to check the full 4x4 matrix identity (not merely sorted eigenvalues) at generic momenta and Brillouin-zone boundaries; local execution is pending. The 16x16 production mean-field builder check is intentionally outside the exact matter-sector scope.
- Track A is planned; shared-file synchronization has not yet been applied.
