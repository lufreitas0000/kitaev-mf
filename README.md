# kitaev-mf: Majorana Mean-Field Solver for Extended Kitaev Systems

[![CI](https://github.com/lufreitas0000/kitaev-mf/actions/workflows/ci.yml/badge.svg)](https://github.com/lufreitas0000/kitaev-mf/actions/workflows/ci.yml)
[![Build & Publish Docs](https://github.com/lufreitas0000/kitaev-mf/actions/workflows/docs.yml/badge.svg)](https://github.com/lufreitas0000/kitaev-mf/actions/workflows/docs.yml)
[![Hypothesis Contracts](https://img.shields.io/badge/Hypothesis_Contracts-Phase_1_%2B_PG--10-1A365D)](doc/HYPOTHESIS_CONTRACTS.md)
[![Living Lab Book PDF](https://img.shields.io/badge/Living_Lab_Book-Pages_PDF-8B2635)](https://lufreitas0000.github.io/kitaev-mf/kitaev-mf-labbook.pdf)
[![Manuscript PDF](https://img.shields.io/badge/Manuscript-Pages_PDF-2E6F40)](https://lufreitas0000.github.io/kitaev-mf/kitaev-mf-manuscript.pdf)

> ### 📄 Quick Access to Compiled PDFs & Key Deliverables
> - **Living Lab Book Monograph (11 chapters + architecture appendix):** [**Read Online (GitHub Pages)**](https://lufreitas0000.github.io/kitaev-mf/kitaev-mf-labbook.pdf) &nbsp;|&nbsp; [LaTeX Source](doc/labbook/labbook_kitaev_mf.tex) &nbsp;|&nbsp; Local build: `make labbook`
> - **Technical Manuscript (`manuscript_kitaev_mf`):** [**Read Online (GitHub Pages)**](https://lufreitas0000.github.io/kitaev-mf/kitaev-mf-manuscript.pdf) &nbsp;|&nbsp; [LaTeX Source](doc/manuscript/manuscript_kitaev_mf.tex) &nbsp;|&nbsp; Local build: `make manuscript`
> - **Centralized Hypothesis Contracts Registry & PG-10 Contracts:** [**View Contracts Registry**](doc/HYPOTHESIS_CONTRACTS.md) &nbsp;|&nbsp; [PG-10 YAML Guide](doc/hypotheses/README.md) &nbsp;|&nbsp; [View on GitHub Pages](https://lufreitas0000.github.io/kitaev-mf/HYPOTHESIS_CONTRACTS.md)
> - **Interactive R&D Epistemic Portal:** [**Browse Live Portal**](https://lufreitas0000.github.io/kitaev-mf/) (Mermaid epistemic DAG, R&D notes, verified invariants, `.dat` benchmark tables, and SHA-256 `book-build.json` receipt)

`kitaev-mf` is a computational and analytical research library dedicated strictly to **4-Majorana Gaussian self-consistent mean-field theory (`scMFT`)** for Kitaev quantum magnets on honeycomb lattices, formulated on the orthogonal Grassmannian manifold bundle $\mathrm{Gr}_{\mathbb{C}}(2N_{\mathrm{sub}}, 4N_{\mathrm{sub}})$ and managed under a linear single-`main` Git workflow ([`ADR-003`](doc/adr/ADR-003-linear-main-and-family-boundaries.md)).

**Goal:** Deliver a self-contained, algebraically verified, and numerically audited 4-Majorana Gaussian `scMFT` solver, Living Lab Book, and publication manuscript for extended Kitaev spin liquids, cleanly decoupled from exact many-body and walk-resolvent sibling repositories (`ADR-003`).
**Scope:** Owns pure-`scMFT` core theory/code (`src/kitaev_mf/`), benchmark tables (`data/processed/benchmarks/`), `doc/labbook/`, `doc/manuscript/`, `doc/rnd/`, and SHA-256 external data ingestion (`scripts/sync_external_data.py`); excludes Spin-ED, 2-RDM SDP (`kitaev-numerics`), and full-$Q$ resolvents (`kitaev-walk-geometry`).
**Definition of Done (DoD):** Direct `pytest` (`make verify` with `pipefail`), `ruff`, `mypy`, `check_index.py`, `sync_external_data.py --verify`, `planning_records.py check`, `uow.py check`, `make labbook`, `make manuscript`, and CI/Pages builds pass with zero unmasked failures and independent W4 review receipts.
**Next Steps:** See [`CURRENT_SPRINT.md`](CURRENT_SPRINT.md) for human-gated `IDEA-001` / `K-019-A` closure or `doc/hypotheses/HC-01.yaml` Phase 2 registration.

---

## Compiled Monographs & Key Deliverables

| Deliverable | GitHub Pages (Direct Browser View) | Repository Source | Description |
| :--- | :--- | :--- | :--- |
| **Living Lab Book (PDF)** | [**Open PDF (Pages)**](https://lufreitas0000.github.io/kitaev-mf/kitaev-mf-labbook.pdf) | [`doc/labbook/labbook_kitaev_mf.tex`](doc/labbook/labbook_kitaev_mf.tex) | Comprehensive R&D monograph (`memoir`, 11 chapters + architecture appendix, TikZ schematics, and PGFPlots benchmarks). |
| **Manuscript (PDF)** | [**Open PDF (Pages)**](https://lufreitas0000.github.io/kitaev-mf/kitaev-mf-manuscript.pdf) | [`doc/manuscript/manuscript_kitaev_mf.tex`](doc/manuscript/manuscript_kitaev_mf.tex) | Technical working paper synthesizing Phase 1 `scMFT` derivations, solver geometry, and empirical benchmarks. |
| **R&D Epistemic Portal** | [**Interactive Portal**](https://lufreitas0000.github.io/kitaev-mf/) | [`scripts/build_blueprint.py`](scripts/build_blueprint.py) | Live GitHub Pages portal with PDF links, Mermaid epistemic DAG, R&D notes, `.dat` tables, and `book-build.json`. |
| **Hypothesis Contracts Registry** | [**View on Pages**](https://lufreitas0000.github.io/kitaev-mf/HYPOTHESIS_CONTRACTS.md) | [`doc/HYPOTHESIS_CONTRACTS.md`](doc/HYPOTHESIS_CONTRACTS.md) | Historical Phase 1 contract audit table (`PH-*`, `AH-*`, `HC-P1-*`, `C-M1.9-*`) plus [`doc/hypotheses/`](doc/hypotheses/README.md) rules. |
| **Linear Topology & Boundaries ADR** | — | [`doc/adr/ADR-003-linear-main-and-family-boundaries.md`](doc/adr/ADR-003-linear-main-and-family-boundaries.md) | Single-`main` workflow, Kitaev 5-repository division, `data/external/` contract, and CI/Pages release policy. |
| **Master Solver Resource Table** | — | [`master_solver_resource_table.dat`](data/processed/benchmarks/m1_9/master_solver_resource_table.dat) | Solver selection rules, complexity, RAM footprint, iteration counts, and wall times across Phase 1 benchmarks. |
| **Verified Physics Invariants** | — | [`PROJECT_STATE.yaml`](PROJECT_STATE.yaml) | Single source of truth indexing verified invariants (`K-001`–`K-023`), milestone history (`M1.1`–`M1.11`), and `RECON-20261116` tasks. |
| **Milestone Roadmap & Sprint Brief** | — | [`doc/MILESTONE_ROADMAP.md`](doc/MILESTONE_ROADMAP.md) / [`CURRENT_SPRINT.md`](CURRENT_SPRINT.md) | Strategic roadmap and owner-gated sprint status. |
| **Discovery & Pivot Journal** | — | [`LOG_DISCOVERY_AND_PIVOTS.md`](LOG_DISCOVERY_AND_PIVOTS.md) | Chronological log of mathematical obstructions, refuted conjectures, errata, and algorithmic pivots. |

> **Artifact Policy (`ADR-003`):** Compiled `.pdf` binaries are never committed to Git. `.github/workflows/docs.yml` builds and publishes both PDFs atomically to GitHub Pages on pushes to `main` and attaches release PDFs strictly to immutable `v*.*.*` SemVer tags.

---

## Kitaev 5-Repository Family & External Data Gate (`ADR-003`)

To prevent scope creep and keep `kitaev-mf` auditable as a pure 4-Majorana Gaussian `scMFT` engine, exact non-mean-field and complementary layers reside in dedicated sibling repositories:

| Repository | Owned Responsibility | Interface into `kitaev-mf` |
|---|---|---|
| **`kitaev-mf`** | 4-Majorana Gaussian `scMFT` (`domain/` + `services/` + `infra/`), Grassmannian/DIIS/L-Broyden/Padé purification/KKT solvers, real-space defects/vortices, `scMFT` response | Native source & tests |
| **`kitaev-numerics`** | Exact Spin-ED, Lanczos/KPM, OFTLM, 2-RDM SDP lower bounds, HTE, `iPEPS`/tensor networks | `data/external/kitaev_numerics/*.dat` + `.meta.json` |
| **`kitaev-walk-geometry`** | Exact Majorana lattice walks, lattice Green's functions/resolvents, full-$Q$ Feshbach–Schur perturbation theory | `data/external/kitaev_walk_geometry/*.dat` + `.meta.json` |
| **`ksl-sc-bilayer`** | Kitaev spin liquid–superconductor heterostructures and proximity/Kondo bilayer couplings | Independent repository |
| **`kitaev-derivation`** | Pure-Python `SymPy` Clifford/CAR & Schrieffer–Wolff projector CAS + 66-state sparse Fock oracle (`lufreitas0000/kitaev-derivation`, replacing archived `KitaevDerivation.jl`) | Referenced symbolic derivations & Wave 3 ring exchange |

Cross-repository comparison datasets are verified in CI via `python3 -B scripts/sync_external_data.py --verify`; unverified historical tables remain quarantined under [`data/external/historical_unverified/`](data/external/historical_unverified/README.md).

---

## Strategic Multi-Phase Architecture

The research program is organized into dedicated model-family phases:

### 1. Phase 1: Pure Ferromagnetic Kitaev Honeycomb ($K = -1.0$) with $[111]$ Field — Numerical Quantification & Algorithmic Foundations (**Core `scMFT` Foundation; `IDEA-001` & `K-019-A` Open**)
- **Exclusive Focus:** Exhaustive numerical quantification of all fixed-point solver algorithms on the orthogonal Grassmannian $\mathrm{Gr}(4, 8)$ and real-space torus clusters: Damped Picard, Joint DIIS (Anderson Acceleration), Full Good Broyden, Limited-Memory Broyden (`L-Broyden`, $\mathcal{O}(mD)$), Jacobian-Free Newton-Krylov (`JFNK`), and eigenvector-free Matrix Sign Purification ($\Omega = \mathrm{i}\operatorname{sgn}(\mathrm{i}\mathbf{A})$).
- **Mastery of Numerical Nuances:**
  - *Solver Pareto Frontier (`K-012`, `K-013`):* Full Broyden and Joint DIIS expand the 90% Wilson survival basin radius by $>35\times$ ($R_{90} = 0.350\text{--}0.500$ vs. $0.010$ for Picard) and accelerate convergence by $4.26\times$. Damping the L-Broyden base step ($\alpha = 0.5$, `K-021`) stabilizes the $(D-m)$-dimensional orthogonal complement on high-dimensional real-space manifolds ($D \sim 10^4$).
  - *First-Order Spinodal & KKT Admissibility (`K-011`, `K-023`):* Thermodynamic crossing at $h_c^* = 0.3325 \pm 0.0005$ with coexistence domain-wall surface tension $\sigma_{\mathrm{DW}} = +0.042|K| > 0$, forward spinodal at $h_{c2}^+ = 0.400 \pm 0.005$, and complete elimination of post-spinodal unconverged saddle excursions via KKT dual stationarity enforcement ($\epsilon_{\mathrm{gauge}} \le 5.6 \times 10^{-15}$).
  - *Cubic Dirac Gap Channel Decomposition & Open Analytical Work (`K-010`, `K-019`, `K-019-A`):* Finite-mesh 4-flavor `scMFT` yields empirical regressions $\Delta_K^{\mathrm{MF}}(h) \approx (2.5921 \pm 0.0034) h^3$ ($\Delta_{\mathrm{canonical}}^{\mathrm{MF}} = 4\Delta_K^{\mathrm{MF}} \approx 10.3684 h^3$) and $\|\boldsymbol{\lambda}^*\|_2 \approx \sqrt{5} h$; closed-form BZ evaluation remains open under `IDEA-001`. Quantitative comparison against pure-Kitaev perturbation theory accounts for canonical normalization ($4\times$), the two-gauge vs. two-vison denominator scale, self-consistent KKT counter-field vertex dressing, and the [`K-019-A` second-order Feshbach–Schur erratum](doc/rnd/phase_1/m1_10_t1_k019a_second_order_erratum.md) (full-$Q$ resolvent open under Issues `#10`, `#12`).
  - *Real-Space Matrix Sign Purification (`K-017`, `K-018`, `K-022`):* Track A real skew-symmetric arithmetic halves memory ($2.0\times$ with $\frac{26}{3}M^3$ real FLOPs), while Track B cubic Newton-Schulz purification bypasses $3 \mid L$ Dirac resonance condition-number squaring ($\kappa(\mathbf{S}) \approx 3.6 \times 10^{31}$) and LAPACK `zheevd` secular deflation stalling, achieving a **$2.31\times$ statistical mean speedup** ($109.7\text{--}114.6\text{ GFLOP/s}$ sustained `dgemm` throughput at $4.01 \times 10^{-15}$ accuracy, with up to $63.53\times$ against single-run deflation stalls) and 100% convergence under bond disorder and vacancy zero modes ($\kappa(\mathbf{A}) = 3.82 \times 10^4$).

#### Phase 1 Certified Master Solver Selection & Resource Table (`HC-P1-03`)

| Physical Regime | System Size | State Dim $D$ | Certified Optimal Solver | Step Complexity | Peak RAM | Typical Iters | Wall Time (i5-13500) |
| :--- | :--- | :---: | :--- | :---: | :---: | :---: | :---: |
| **Bulk 2-Site ($T=0$, Weak/Int. $h$)** | $N_k = 24^2$ ($N=1152$) | $86$ | **Joint DIIS ($m=6$)** | $\mathcal{O}(N_k n_c^3)$ | $2.1\text{ MB}$ | $7\text{--}10$ | $0.18\text{ s}$ |
| **Bulk 2-Site (Near Spinodal $h_{c2}^+$)** | $N_k = 24^2$ ($N=1152$) | $86$ | **Full Broyden ($\alpha=0.5$)** | $\mathcal{O}(D^2 + N_k n_c^3)$ | $2.4\text{ MB}$ | $10\text{--}14$ | $0.27\text{ s}$ |
| **Bloch Supercells ($n_c \in [4, 32]$)** | $N_k = 12^2$, $n_c \le 32$ | $\le 1376$ | **Joint DIIS ($m=8$)** | $\mathcal{O}(N_k n_c^3)$ | $18.5\text{ MB}$ | $18\text{--}23$ | $1.42\text{ s}$ |
| **Real-Space Clean Torus ($T=0$)** | $L=12$ ($M=1152$) | $11{,}520$ | **Purif NS + DIIS ($m=8$)** | $\mathcal{O}(k_{\mathrm{NS}} M^3)$ | $32.0\text{ MB}$ | $7\text{--}9$ | $5.60\text{ s}$ |
| **Real-Space Vortex Pair / Defect** | $L \in [8, 12]$ ($M \le 1152$) | $\le 11{,}520$ | **Purif NS + L-Broyden ($m=10, \alpha=0.5$)** | $\mathcal{O}(k_{\mathrm{NS}} M^3 + mD)$ | $12.4\text{ MB}$ | $10\text{--}13$ | $7.80\text{ s}$ |
| **Large Real-Space Defect Torus** | $L=24$ ($M=4608$) | $46{,}080$ | **Outer-Warm Purif + L-Broyden ($\alpha=0.5$)** | $\mathcal{O}(k_{\mathrm{NS}} M^3)$ | $170\text{ MB}$ | $10\text{--}15$ | $48\text{ s}$ ($4\text{ s}$ GPU) |
| **Finite-Temperature Bulk ($T > 0$)** | $N_k = 36^2$, $n_c=2$ | $86$ | **Bloch FD `eigh` + Joint DIIS** | $\mathcal{O}(N_k n_c^3)$ | $4.5\text{ MB}$ | $6\text{--}12$ | $0.35\text{ s}$ |

---

### 2. Phase 2: Low-Energy Physics of the Pure Ferromagnetic Kitaev $[111]$ Model (**ACTIVE PREPARATION**)
Phase 2 is the direct physical continuation of Phase 1 in the **identical model regime** ($K = -1.0$, $\mathbf{h} \parallel [111]$), conducting calculations primarily on real-space clusters with Bloch momentum space used for reality checks:
1. **$T = 0$ Bulk Ground State (`M2.1`):** High-resolution phase diagram ($h \in [0, 0.60]$), magnetic susceptibility tensor $\chi^{\alpha\beta}(h)$, topological Chern number $\nu \in \{-1, 0\}$, factored 6-link vs. exact 12-Majorana Pfaffian flux $\mathcal{W}_p = \operatorname{Pf}(\Omega|_{\partial p})$, and Gutzwiller / 2-RDM semidefinite lower bounds.
2. **$T = 0$ State with Isolated Pairs of $\mathbb{Z}_2$ Vortices (`M2.2`):** Self-consistent two-vison interaction potential $V(R, h)$ on large real-space tori ($L \in \{12, 16, 24\}$), non-Abelian Majorana bound-mode hybridization splitting $\delta\epsilon_0(R) \propto e^{-R/\xi(h)}$, and local spin/gauge screening clouds ($\xi_{\mathrm{gauge}} \approx 1.1 a$).
3. **Finite-Temperature ($T > 0$) States (`M2.3`):** Self-consistent Fermi–Dirac thermal covariance $\Omega_T = \mathrm{i}\tanh(\mathrm{i}\beta\mathbf{A}/2)$, two-peak specific heat $C_v(T)$, entropy release $S(T)$ from $\ln 2$, and thermal melting of the first-order spinodal hysteresis loop.
4. **Static Spin Correlations (`M2.4`):** Equal-time structure factor $S^{\alpha\beta}(\mathbf{q})$ and real-space correlation functions $\langle \sigma_i^\alpha \sigma_j^\beta \rangle$ across the ground state, two-vortex state, and thermal ensemble.

---

### 3. Phase 3: Dynamical Spin Correlations & Majorana Zero Mode Spectroscopy (`scMFT` Response & Cross-Repo Baselines)
Phase 3 evaluates dynamical spin correlations and STM tunneling signatures of Majorana zero modes within 4-Majorana `scMFT`, benchmarked against imported external pure-Kitaev datasets (`ADR-003`):
1. **External Exact Quench Baselines (`kitaev-walk-geometry` / `kitaev-numerics` $\to$ `data/external/`):** Exact gauge-quench dynamics (Knolle et al. 2014, Baskaran et al. 2007) and arXiv:2209.12930 MZM reference spectra are produced in sibling repositories and imported under SHA-256 provenance manifests.
2. **Majorana Mean-Field Generalization (`M3.4`):** Dynamical spin structure factor $S^{\alpha\beta}(\mathbf{q}, \omega)$ within the 4-Majorana `scMFT` framework of `kitaev-mf`, benchmarked against the imported pure-Kitaev reference tables (`M3.5`).
3. **Scanning Tunneling Spectroscopy of MZMs:** Local tunneling density of states (LDOS) and in-gap Majorana bound-state signatures in isolated `scMFT` vortex pair sectors.

---

### 4. Subsequent Phases (Dedicated Model Families)
Following Phase 3, subsequent model families are developed sequentially:
- **Phase 4:** Antiferromagnetic Kitaev Honeycomb Model ($K = +1.0$) in $[111]$ Field (intermediate gapless U(1) spin liquid candidate, emergent spinon Fermi surfaces).
- **Phase 5:** Extended Non-Integrable Honeycomb Exchange ($J\text{--}K\text{--}\Gamma\text{--}\Gamma'\text{--}D$) for $\alpha\text{-RuCl}_3$ and $\text{Na}_2\text{IrO}_3$.
- **Phase 6:** Three-Dimensional Tricoordinated Kitaev Lattices (hyperhoneycomb and stripyhoneycomb).

---

## Repository Architecture

The production Python package follows strict hexagonal (ports-and-adapters) separation:

- `src/kitaev_mf/domain/`: mathematical entities, Clifford CAR operators, lattice definitions, model parameters, inhomogeneous real-space state (`InhomogeneousMeanFieldState`), solver configuration, and port protocols.
- `src/kitaev_mf/services/`: pure-`scMFT` algorithms (Bloch Hamiltonians, spectral projectors, Grassmannian fixed-point solvers with SPD-floored Lehmann dual Hessian preconditioning, real-space skew-symmetric operators, Newton–Schulz and Kenney–Laub Padé `[3/2]` matrix sign purification `purify_pade_sign`, `canonical_skew_schur` real canonical Schur blocks, 12-Majorana Parlett–Reid Pfaffian flux metrology `compute_skew_pfaffian`, observables, and energetics).
- `src/kitaev_mf/infra/`: concrete linear-algebra, hardware FLOP profiling (`ResourceProfiler`), and persistence adapters.
- `test/unit/` and `test/integration/`: automated verification suite enforcing algebraic, kinematic, external-data (`test_sync_external_data.py`), and numerical invariants.
- `doc/labbook/`: `memoir` LaTeX source, chapters, architecture appendix, and PGFPlots figures for the Living Lab Book (`make labbook`).
- `doc/manuscript/`: LaTeX source for `manuscript_kitaev_mf.pdf` (`make manuscript`).
- `doc/rnd/`: Front-1 R&D technical notes (`phase_1/`), exact perturbation theory exploration notes (`exact-perturbation-theory/`), and peer-review records (`doc/reviews/`).

## Development & Verification Setup

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel pip-tools
pip-sync requirements.txt requirements-dev.txt
python -m pip install -e . --no-deps
```

Run the verification and document build suite locally:

```bash
make verify                                          # ruff check, ruff format --check, mypy, and pytest
python3 -B scripts/sync_external_data.py --verify    # enforce ADR-003 data/external/ provenance sidecars
python3 -B doc/reconciliation/check_index.py         # verify absorbed reconciliation catalog and hashes
make labbook                                         # compiles gitignored doc/labbook/labbook_kitaev_mf.pdf
make manuscript                                      # compiles gitignored doc/manuscript/manuscript_kitaev_mf.pdf
```

## Governance and Releases

The pre-release package version (`0.1.0`) is declared in `pyproject.toml`, with ongoing changes documented under `## [Unreleased]` in [`CHANGELOG.md`](CHANGELOG.md). See [`ADR-003`](doc/adr/ADR-003-linear-main-and-family-boundaries.md), [`governance/profile.toml`](governance/profile.toml), [`AGENTS.md`](AGENTS.md), and [`.agents/rules/04_versioning_and_releases.md`](.agents/rules/04_versioning_and_releases.md) for SemVer tags (`v*.*.*`), Hypothesis Contracts, and Pages publishing.
