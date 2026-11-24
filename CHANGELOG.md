# Changelog

All notable user-visible changes to `kitaev-mf` are documented here.

The package version (`0.1.0`) is maintained in `pyproject.toml`. See
[Rule 04: Versioning, Changelog, and Releases](.agents/rules/04_versioning_and_releases.md)
for the SemVer policy and tag-based release procedure.

## [Unreleased]

*(Pre-`1.0.0` development tracking under package version `0.1.0` per Rule 04; formal releases will be tagged `v*.*.*` after scientific milestone closure.)*

### [2026-11-24] [M2.2-WAVE2A / UOW-01] [Actor: blue-lead] — Cross-QSL `(L^a, R^a)` / `W_p` / `(u_ij, v_ij)` convention bridge port & exact Fock parity tests (`kitaev-mf#25`)
#### Ideas & Contracts Evolved
- Codified the opt-in cross-QSL `conventions-v2.0` boundary contract (`CONVENTIONS.md`, `KNOWLEDGE_MAP.md`) preserving `kitaev-mf`'s locked internal `(-c_j, b_j^x, b_j^y, b_j^z)` basis while proving exact lossless equivalence to `ksl-sc-bilayer` (`(+c_j, b_j^x, b_j^y, b_j^z)`), `kitaev-numerics`, and `kitaev-walk-geometry` via the parity involution `P_0 = diag(-1, 1, 1, 1)` (`P_0 L^a P_0 = R^a`, `+0.5 M^a = R^a`, `-0.5 G^a = L^a`, `lambda_v2 = -2 lambda_kmf`, and sign flip on off-diagonal matter–gauge blocks `U_b^{0,a}`).
- Verified on the exact 16-state two-site Fock space (`< 1e-13`) and the exact 4096-state 6-site sparse Fock hexagon space (`< 1e-13`) that `P_+ W_p^{spin} P_+ = +P_+ prod_{b in partial p} (i b_{A_b}^{alpha_b} b_{B_b}^{alpha_b}) P_+` and that directed pure-Kitaev bond scalars `(u_ij^alpha, v_ij^alpha)` are `P_0`-invariant.

#### Code & Artifacts Changed
- Added `src/kitaev_mf/domain/qsl_bridge.py` (`SO4_L`, `SO4_R`, `P0_PARITY_BRIDGE`, `DefaultQSLConventionBridge`, state/Lagrange/coupling/chiral-`kappa` adapters, and `(u_ij, v_ij)`/`W_p` extractors), added `QSLConventionBridgePort` to `src/kitaev_mf/domain/ports.py`, and re-exported public bridge symbols from `src/kitaev_mf/domain/__init__.py`.
- Added independent exact 16-state & 4096-state Fock unit tests in `test/unit/test_qsl_convention_bridge.py`.
- Added root `CONVENTIONS.md`, living `KNOWLEDGE_MAP.md`, and passive Pre-/Post-UoW observatory entry in `governance/info_telemetry.jsonl`.

### [2026-11-19] [M2.1-WAVE1] [Actor: antigravity] — Wave 1 `HC-01..03` sweeps, Gutzwiller sandwich, and stratified two-vison bridge
#### Ideas & Contracts Evolved
- Quantified `HC-01` Feshbach–Schur cubic Dirac gap (`c_MF_cont = 10.641666`, `p = 2.999694`), isolating KKT $\lambda^{*,\gamma} G^\gamma$ inside $Q$-space (`<= 4.04e-15` without $\lambda^*$) alongside 95% Hartree $PQ$-screening.
- Verified `HC-02` `24/24` $E_{\mathrm{SDP}}^{\mathrm{cert}} \le E_0^{\mathrm{ED}} \le E_{\mathrm{proj}}$ sandwich bounds across $N \in \{8, 12, 18\}$ with Zero-Variance $N=12$ $h=0$ exact projection (`-0.8166902250`, $e_{\mathrm{loc}}\text{-std} = 3.08 \times 10^{-16}$) and 50-chain $\mathrm{BC}_a$ 95% CIs.
- Reconciled `HC-03` normalization (`kappa_static = kappa_KWG / 2` matching `kitaev-walk-geometry` to `< 4.4e-9`) across 94 chiral runs ($L \in \{8..20\}$), discovering period-3 intervalley $\cos(\mathbf{K}\cdot\mathbf{R})$ modulation (`Delta_BIC in [2.39, 4.68]` for raw monotone fits vs $\xi_M \approx \{3.03, 1.36, 0.76\}$ on fixed-phase $R=3 \to 6$).

#### Code & Artifacts Changed
- Added `src/kitaev_mf/services/{cubic_gap_perturbation,gutzwiller_projection,vortex_pair_scaling}.py`, `scripts/wave1/run_hc0{1,2,3}_*.py`, unit suites `test/unit/test_wave1_*.py`, external sync under `data/external/{kitaev_numerics,kitaev_walk_geometry}/`, Front-1 notes `doc/rnd/phase_2/wave1_hc0*.md`, Wave 1 `.dat`/`.meta.json` tables, and `learning_log` updates in `doc/hypotheses/HC-01..03.yaml`.

### [2026-11-19] [UOW-20261119-01] [Actor: antigravity] — Add HORIZON.md, rename thesis entrypoint, unify uv venv, and stage Wave 1 HC-01..03.yaml
#### Ideas & Contracts Evolved
- Staged `doc/hypotheses/HC-01.yaml`–`HC-03.yaml` (`status: draft`, `TODO(owner)` thresholds preserved) linking `IDEA-001` / `K-019-A`, `kitaev-numerics` ED/Clarabel SDP sandwich comparison, and `kitaev-walk-geometry:HC-11` two-vison FSS scaling ahead of Wave 1 sweeps.
- Authored `HORIZON.md` and updated sibling table from archived `KitaevDerivation.jl` to pure-Python `kitaev-derivation`.

#### Code & Artifacts Changed
- Renamed `doc/reference/thesis/main.tex` to `doc/reference/thesis/thesis_kitaev_mf_ref.tex`, switched `Makefile` `venv`/`lock`/`sync` targets to `uv`, and updated `PROJECT_STATE.yaml` (`doc/labbook/labbook_kitaev_mf.tex`) and stale `/home/lucas/Projects/` documentation references.

### [2026-11-18] [RECON-20261116-T5] [Actor: antigravity-impl] — Enforce `pipefail` in `Makefile` `verify` gate
#### Ideas & Contracts Evolved
- None (build-harness reliability fix; closes the `pytest | grep ... || true` exit-code masking caveat noted at `RECON-20261116` closure).

#### Code & Artifacts Changed
- Updated `Makefile` (`SHELL := /bin/bash -o pipefail`, dropped `|| true` from step `[4/4]`) so any `pytest` failure fails `make verify` cleanly.
- Added UoW `Goal`/`Scope`/`Definition of Done (DoD)`/`Next Steps` markers to `README.md` and `Owning milestone: RECON-20261116` to `CURRENT_SPRINT.md` so `scripts/uow.py check` passes cleanly.

### Added (`RECON-20261116` Consolidation & Pure-`scMFT` Core Kernels)
- **Linear Git Topology & 5-Repo Family Boundaries (`ADR-003`):** Formalized one permanent `main` branch and one worktree workflow, documented repository boundaries across `kitaev-mf` (pure `scMFT`), `kitaev-numerics`, `kitaev-walk-geometry`, `ksl-sc-bilayer`, and `Julia/KitaevDerivation.jl`, and added PG-10 `native-v2` profile (`governance/profile.toml`) plus `doc/hypotheses/README.md`.
- **Cross-Repository External Data Contract (`scripts/sync_external_data.py`):** Implemented CLI `--ingest` / `--verify` checks requiring sibling `.meta.json` sidecars (`source_repo`, 40-hex `git_sha`, matching `sha256`, `unit_convention`, `geometry`, `parameters`) for active datasets under `data/external/<source_repo>/`, while quarantining legacy comparison tables in `data/external/historical_unverified/`.
- **Track 1 Pure-`scMFT` Real-Space & Rational Kernels (`PR #23`, `W4` reviewed):** Ported `InhomogeneousMeanFieldState` (`domain/state.py`), `SolverConfig.precondition_dual` and SPD-floored Lehmann dual Hessian preconditioning (`services/solver.py`), Kenney–Laub Padé rational matrix sign purification (`purify_pade_sign`, `services/purification.py`), real skew-symmetric canonical Schur decomposition (`canonical_skew_schur`, `services/canonical_skew.py`), and Parlett–Reid 12-Majorana local plaquette Pfaffian flux metrology (`compute_skew_pfaffian`, `RealSpaceTorusCluster.compute_plaquette_pfaffian_flux`, and `compute_all_plaquette_pfaffian_fluxes`, `services/real_space.py`) with corresponding unit tests.
- **Historical R&D, W001 Audit & Wave 4 Absorption (`PR #20`, `PR #22`):** Integrated Wave 4 statistical sampling/symmetry contracts, W001 foundational audit reports (`doc/audits/`), Phase 1–4 R&D notes (`doc/rnd/`), `doc/rnd/exact-perturbation-theory/`, and deduplicated `LOG_DISCOVERY_AND_PIVOTS.md`.
- **CI, Docs/Pages & Website Feed Planning:** Consolidated `.github/workflows/latex.yml` and `blueprint.yml` into `.github/workflows/docs.yml` (atomic GitHub Pages PDFs + Blueprint portal + `book-build.json`, tag-only `v*.*.*` release PDFs, no mutable `latest` release), added Lab Book appendix [`app_architecture_and_foundations.tex`](doc/labbook/appendices/app_architecture_and_foundations.tex), and drafted [`doc/editorial/PERSONAL_SITE_FEED_PLAN.md`](doc/editorial/PERSONAL_SITE_FEED_PLAN.md).

### Changed & Fixed
- **SemVer & Artifact Hygiene:** Replaced the premature `[1.0.0]` heading with `[Unreleased]`, renamed `doc/manuscript/scipost_kitaev_mf.tex` to `doc/manuscript/manuscript_kitaev_mf.tex`, removed stale `doc/labbook/preamble.tex`, and excluded compiled PDFs and `book-build.json` from Git tracking.
- **Epistemic Qualifications (`H6`, `K-019-A`, Parity Sign):** Clarified empirical finite-mesh regressions ($\Delta_K^{\mathrm{MF}} \approx 2.5921 h^3$, $\|\boldsymbol{\lambda}^*\|_2 \approx \sqrt{5} h$) vs. open analytical evaluation (`IDEA-001`), documented the `K-019-A` second-order Feshbach–Schur erratum (`m1_10_t1_k019a_second_order_erratum.md`), and harmonized the local parity sign convention $\mathcal{D}_i = +c_i^0 c_i^1 c_i^2 c_i^3 = -(-1)^{n_1+n_2}$ across Lab Book and manuscript sections.

### Prior Phase 1 Foundation Work (2026-11-06)
- **Algebraic Foundations & CAR Algebra (`M1.1`):** Verified local 4-Majorana CAR Clifford algebra, Casimir invariants, Dirac fermion representation, and single-generator parity (`K-001`–`K-003`).
- **Bloch Momentum & Supercell Folding (`M1.2`):** Explicit 2-site and 4-site reciprocal Bloch Hamiltonian matrices $\mathcal{H}(\mathbf{k})$, exact pure-Kitaev matter and flat gauge dispersions (`K-004`), and unitary folding equivalence matrix $W(\mathbf{q})$ (`K-005`).
- **Spectral Projector & DIIS Solver (`M1.3`–`M1.4.5`):** Zero-temperature spectral projector, Borwein–Pearson high-precision quadrature `K-006`, scalar double-counting subtraction `K-007`, Joint Anderson Acceleration, Parseval/semidefinite representability gates, and Klein duality boundary splitting `K-009`.
- **Statistical Solver Scorecards, Supercells & Real-Space Purification (`M1.6`–`M1.11`):** Grassmannian fixed-point solvers (Picard, DIIS, Full Broyden, damped L-Broyden, JFNK), cubic Dirac gap fits (`K-010`), spinodal $h_{c2}^+ = 0.400$ and KKT post-spinodal audit (`K-011`, `K-023`), supercell scaling (`K-015`, `K-016`), real skew-symmetric arithmetic (`K-017`), eigenvector-free Newton–Schulz matrix sign purification (`K-018`, `K-021`), vacancy/disorder conditioning (`K-022`), and two-vison real-space vortex branch cuts (`K-020`).
- **Monographs & Epistemic Portal:** Authored the Living Lab Book (`doc/labbook/labbook_kitaev_mf.tex`), preprint manuscript (`doc/manuscript/manuscript_kitaev_mf.tex`), Hypothesis Contracts registry (`doc/HYPOTHESIS_CONTRACTS.md`), and static Blueprint portal generator (`scripts/build_blueprint.py`).
