# Changelog

All notable user-visible changes to `kitaev-mf` are documented here.

The package version (`0.1.0`) is maintained in `pyproject.toml`. See
[Rule 04: Versioning, Changelog, and Releases](.agents/rules/04_versioning_and_releases.md)
for the SemVer policy and tag-based release procedure.

## [Unreleased]

*(Pre-`1.0.0` development tracking under package version `0.1.0` per Rule 04; formal releases will be tagged `v*.*.*` after scientific milestone closure.)*

### Added (`RECON-20261116` Consolidation & Pure-`scMFT` Core Kernels)
- **Linear Git Topology & 5-Repo Family Boundaries (`ADR-003`):** Formalized one permanent `main` branch and one worktree workflow, documented repository boundaries across `kitaev-mf` (pure `scMFT`), `kitaev-numerics`, `kitaev-walk-geometry`, `ksl-sc-bilayer`, and `Julia/KitaevDerivation.jl`, and added PG-10 `native-v2` profile (`governance/profile.toml`) plus `doc/hypotheses/README.md`.
- **Cross-Repository External Data Contract (`scripts/sync_external_data.py`):** Implemented CLI `--ingest` / `--verify` checks requiring sibling `.meta.json` sidecars (`source_repo`, 40-hex `git_sha`, matching `sha256`, `unit_convention`, `geometry`, `parameters`) for active datasets under `data/external/<source_repo>/`, while quarantining legacy comparison tables in `data/external/historical_unverified/`.
- **Track 1 Pure-`scMFT` Real-Space & Rational Kernels (`PR #23`, `W4` reviewed):** Ported `InhomogeneousMeanFieldState` (`domain/state.py`), `SolverConfig.precondition_dual` and SPD-floored Lehmann dual Hessian preconditioning (`services/solver.py`), Kenney–Laub Padé rational matrix sign purification (`purify_pade_sign`, `services/purification.py`), real skew-symmetric canonical Schur decomposition (`canonical_skew_schur`, `services/canonical_skew.py`), and Parlett–Reid 12-Majorana local plaquette Pfaffian flux metrology (`compute_skew_pfaffian`, `compute_inhomogeneous_plaquette_fluxes_pfaffian`, `services/real_space.py`) with corresponding unit tests.
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
