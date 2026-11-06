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
- [ ] **B-004 Finite-torus conventions.** Definitions are recorded in `notes/05_finite_torus_conventions.md`: periodic indexing, six-edge plaquette cycle and link-flux product, two primitive Wilson loops, seam twists, and exact combinatorial acceptance tests. Implement and execute the geometry tests before closing.

## P1 — Projection and baseline solver

- [ ] **B-010 Physical parity.** Multiply local constraints $D_i=+1$ with all ordering signs explicit; derive parity for arbitrary fixed gauge representatives.
- [ ] **B-011 Gaussian vacuum parity.** Fix Pfaffian/covariance sign convention; handle zero modes and singular cases.
- [ ] **B-012 Solver contracts.** Specify immutable geometry/link inputs and pure outputs for T, A, singular values, vacuum energy, and Fock bounds.
- [ ] **B-013 Geometry and gauge observables.** Implement bond indexing, plaquette fluxes, Wilson loops, twists, and gauge transformations.
- [ ] **B-014 Matter spectrum.** Construct T; test A antisymmetry and equality of singular values of T with positive eigenvalues of iA.
- [ ] **B-015 Flux-free dispersion.** Compare finite-torus spectrum with Bloch dispersion, Dirac-point sampling, and modulo-three effects.
- [ ] **B-016 Vacuum energy and width.** Independent quadrature and finite-size convergence of E_vac/N, W_Fock/N, and parity endpoint corrections.
- [ ] **B-017 Gauge-equivalence tests.** Gauge-related link configurations have identical spectra and energies.

## P1 — Physical flux energetics

- [ ] **B-020 Adjacent pair (M0.4).** One-link flip; verify exactly two adjacent plaquette fluxes change; compute Gaussian energy, then physical parity-corrected energy.
- [ ] **B-021 Separated pairs (M0.5).** Construct strings with controlled endpoint separation and compare paths with identical endpoints.
- [ ] **B-022 Wilson/twist study.** Separate local flux defects from global topological sectors.
- [ ] **B-023 Small-cluster spin ED.** Construct spin Hamiltonian and plaquette operators; resolve flux sectors inside degenerate eigenspaces.
- [ ] **B-024 Thermodynamic extrapolation.** Study size, shape, modulo-three class, and twists; report uncertainty and distinguish Gaussian from physical thresholds.

## P1 — Feshbach and perturbative derivation

- [ ] **B-030 Physical projectors.** Define P_ff and Q_flux with explicit Wilson-loop policy.
- [ ] **B-031 Physical thresholds.** Compute E0, E_Q_min, and Delta_flux on matched geometries and physical Hilbert spaces.
- [ ] **B-032 Resolvent domain.** Prove E in rho(QHQ) and derive a valid finite-volume margin.
- [ ] **B-033 Field-dependent eliminated block.** Analyze QVQ; do not use an unjustified extensive operator-norm smallness condition.
- [ ] **B-034 Flux selection rules.** Enumerate intermediate configurations generated by spin insertions.
- [ ] **B-035 Effective operator by order.** Fix degeneracy and energy-dependence conventions; derive all allowed virtual paths and denominators.
- [ ] **B-036 Third-order interaction.** Derive coefficient and field anisotropy from this Hamiltonian convention; check signs and combinatorics independently.
- [ ] **B-037 Error/validity bounds.** Separate finite-volume perturbative control from thermodynamic claims.

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
