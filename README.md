# kitaev-mf: Majorana Mean-Field Solver for Extended Kitaev Systems

[![CI](https://github.com/lufreitas0000/kitaev-mf/actions/workflows/ci.yml/badge.svg)](https://github.com/lufreitas0000/kitaev-mf/actions/workflows/ci.yml)
[![Build LaTeX Living Lab Book & Manuscript](https://github.com/lufreitas0000/kitaev-mf/actions/workflows/latex.yml/badge.svg)](https://github.com/lufreitas0000/kitaev-mf/actions/workflows/latex.yml)

`kitaev-mf` is a high-performance computational and analytical research library for Majorana parton mean-field theories of Kitaev quantum magnets on honeycomb lattices, formulated as a nonlinear fixed-point problem on the orthogonal Grassmannian manifold bundle $\mathrm{Gr}_{\mathbb{C}}(2N_{\mathrm{sub}}, 4N_{\mathrm{sub}})$.

---

## Compiled Monographs & Key Deliverables

| Deliverable | GitHub Pages (Direct Browser View) | Latest CI Release Download | Repository Path | Description |
| :--- | :--- | :--- | :--- | :--- |
| **Living Lab Book (PDF)** | [**Open PDF (Pages)**](https://lufreitas0000.github.io/kitaev-mf/kitaev-mf-labbook.pdf) | [`kitaev-mf-labbook.pdf`](https://github.com/lufreitas0000/kitaev-mf/releases/download/latest/kitaev-mf-labbook.pdf) | [`doc/labbook/labbook_kitaev_mf.pdf`](doc/labbook/labbook_kitaev_mf.pdf) | **74-page** comprehensive R&D monograph (11 chapters, full proofs, TikZ schematics, and PGFPlots benchmarks). |
| **SciPost Monograph (PDF)** | [**Open PDF (Pages)**](https://lufreitas0000.github.io/kitaev-mf/kitaev-mf-manuscript.pdf) | [`kitaev-mf-manuscript.pdf`](https://github.com/lufreitas0000/kitaev-mf/releases/download/latest/kitaev-mf-manuscript.pdf) | [`doc/manuscript/scipost_kitaev_mf.pdf`](doc/manuscript/scipost_kitaev_mf.pdf) | **9-page** executive publication paper synthesizing the core breakthroughs of Phase 1. |
| **R&D Epistemic Portal** | [**Interactive Portal**](https://lufreitas0000.github.io/kitaev-mf/) | — | [`scripts/build_blueprint.py`](scripts/build_blueprint.py) | Live GitHub Pages portal with PDF links, Mermaid epistemic DAG, 70 R&D notes, and 26 `.dat` benchmarks. |
| **Hypothesis Contracts Registry** | [**View on Pages**](https://lufreitas0000.github.io/kitaev-mf/HYPOTHESIS_CONTRACTS.md) | — | [`doc/HYPOTHESIS_CONTRACTS.md`](doc/HYPOTHESIS_CONTRACTS.md) | Centralized registry of all **22 Hypothesis Contracts** (`PH-01`–`06`, `AH-01`–`05`, `HC-P1-01`–`03`, `C-M1.9-T1`–`T6`, `H1`–`H11` / `HC-P2-01`–`05`) with pass/refuted verdicts and post-mortems. |
| **Master Solver Resource Table** | — | — | [`master_solver_resource_table.dat`](data/processed/benchmarks/m1_9/master_solver_resource_table.dat) | Certified solver selection rules, complexity, RAM footprint, iteration counts, and wall times across all physical regimes. |
| **Verified Physics Invariants** | — | — | [`PROJECT_STATE.yaml`](PROJECT_STATE.yaml) | Single source of truth indexing **23 verified invariants** (`K-001` through `K-023`) and milestone provenance (`M1.1`–`M1.11`). |
| **Milestone Roadmap** | — | — | [`doc/MILESTONE_ROADMAP.md`](doc/MILESTONE_ROADMAP.md) | Multi-phase strategic roadmap (Phase 1 Completed $\to$ Phase 2 FM Kitaev `[111]` Physics $\to$ Phase 3+ Model Families). |
| **Discovery & Pivot Journal** | — | — | [`LOG_DISCOVERY_AND_PIVOTS.md`](LOG_DISCOVERY_AND_PIVOTS.md) | Chronological log of all mathematical obstructions, refuted conjectures, and algorithmic pivots. |

---

## Strategic Multi-Phase Architecture

The research program is organized into dedicated model-family phases:

### 1. Phase 1: Pure Ferromagnetic Kitaev Honeycomb ($K = -1.0$) with $[111]$ Field — Numerical Quantification & Algorithmic Foundations (**COMPLETED & CERTIFIED**)
- **Exclusive Focus:** Exhaustive numerical quantification of all fixed-point solver algorithms on the orthogonal Grassmannian $\mathrm{Gr}(4, 8)$ and real-space torus clusters: Damped Picard, Joint DIIS (Anderson Acceleration), Full Good Broyden, Limited-Memory Broyden (`L-Broyden`, $\mathcal{O}(mD)$), Jacobian-Free Newton-Krylov (`JFNK`), and eigenvector-free Matrix Sign Purification ($\Omega = \mathrm{i}\operatorname{sgn}(\mathrm{i}\mathbf{A})$).
- **Mastery of Numerical Nuances:**
  - *Solver Pareto Frontier (`K-012`, `K-013`):* Full Broyden and Joint DIIS expand the 90% Wilson survival basin radius by $>35\times$ ($R_{90} = 0.350\text{--}0.500$ vs. $0.010$ for Picard) and accelerate convergence by $4.26\times$. Damping the L-Broyden base step ($\alpha = 0.5$, `K-021`) stabilizes the $(D-m)$-dimensional orthogonal complement on high-dimensional real-space manifolds ($D \sim 10^4$).
  - *First-Order Spinodal & KKT Admissibility (`K-011`, `K-023`):* Thermodynamic crossing at $h_c^* = 0.3325 \pm 0.0005$ with coexistence domain-wall surface tension $\sigma_{\mathrm{DW}} = +0.042|K| > 0$, forward spinodal at $h_{c2}^+ = 0.400 \pm 0.005$, and complete elimination of post-spinodal unconverged saddle excursions via KKT dual stationarity enforcement ($\epsilon_{\mathrm{gauge}} \le 5.6 \times 10^{-15}$).
  - *Cubic Dirac Gap Reconciliation (`K-010`, `K-019`):* Exact reconciliation between the 4-flavor mean-field cubic gap $\Delta_K^{\mathrm{MF}}(h) = 2.5921 h^3$ ($\Delta_{\mathrm{canonical}}^{\mathrm{MF}} = 4\Delta_K^{\mathrm{MF}} = 10.3684 h^3$) and Kitaev's 3rd-order perturbative vison gap ($57.50 h^3$) via the squared energy-denominator ratio $(2\Delta_{\mathrm{gauge}}/\Delta_{\mathrm{flux}})^2 = 15.81\times$ and linear KKT counter-field screening ($\|\boldsymbol{\lambda}^*\|_2 \approx \sqrt{5}h$, $2.85\times$).
  - *Real-Space Matrix Sign Purification (`K-017`, `K-018`, `K-022`):* Track A real skew-symmetric arithmetic halves memory ($2.0\times$), while Track B cubic Newton-Schulz purification bypasses $3 \mid L$ Dirac resonance condition-number squaring ($\kappa(\mathbf{S}) \approx 3.6 \times 10^{31}$) and LAPACK `zheevd` secular deflation stalling, achieving a **$63.53\times$ speedup** ($104.5\text{--}114.6\text{ GFLOP/s}$ sustained `dgemm` throughput at $4.2 \times 10^{-15}$ accuracy) and 100% convergence under bond disorder and vacancy zero modes ($\kappa(\mathbf{A}) = 3.82 \times 10^4$).

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
Phase 2 is the direct physical continuation of Phase 1 in the **identical model regime** ($K = -1.0$, $\mathbf{h} \parallel [111]$). Using the optimal solvers certified in Phase 1, Phase 2 maps the complete phase diagram and computes all physical observables across four regimes (`HC-P2-01`–`05` in [`doc/HYPOTHESIS_CONTRACTS.md`](doc/HYPOTHESIS_CONTRACTS.md)):
1. **$T = 0$ Bulk Ground State (`M2.1`):** High-resolution phase diagram ($h \in [0, 0.60]$), magnetic susceptibility tensor $\chi^{\alpha\beta}(h)$, topological Chern number $\nu \in \{-1, 0\}$, factored 6-link vs. exact 12-Majorana Pfaffian flux $\mathcal{W}_p = \operatorname{Pf}(\Omega|_{\partial p})$, and Gutzwiller / 2-RDM semidefinite lower bounds.
2. **$T = 0$ State with Isolated Pairs of $\mathbb{Z}_2$ Vortices (`M2.2`):** Self-consistent two-vison interaction potential $V(R, h)$ on large real-space tori ($L \in \{12, 16, 24\}$), non-Abelian Majorana bound-mode hybridization splitting $\delta\epsilon_0(R) \propto e^{-R/\xi(h)}$, and local spin/gauge screening clouds ($\xi_{\mathrm{gauge}} \approx 1.1 a$).
3. **Finite-Temperature ($T > 0$) States (`M2.3`):** Self-consistent Fermi–Dirac thermal covariance $\Omega_T = \mathrm{i}\tanh(\mathrm{i}\beta\mathbf{A}/2)$, two-peak specific heat $C_v(T)$, entropy release $S(T)$ from $\ln 2$, and thermal melting of the first-order spinodal hysteresis loop.
4. **Static & Dynamical Spin Correlations (`M2.4`):** Equal-time structure factor $S^{\alpha\beta}(\mathbf{q})$ and dynamical spin structure factor $S^{\alpha\beta}(\mathbf{q}, \omega)$ via four-Majorana Lehmann spectral convolution across the ground state, two-vortex state, and thermal ensemble.

---

### 3. Subsequent Phases (Dedicated Model Families)
Following Phase 2, each distinct model family is developed in its own dedicated phase:
- **Phase 3:** Antiferromagnetic Kitaev Honeycomb Model ($K = +1.0$) in $[111]$ Field (intermediate gapless U(1) spin liquid candidate, emergent spinon Fermi surfaces).
- **Phase 4:** Extended Non-Integrable Honeycomb Exchange ($J\text{--}K\text{--}\Gamma\text{--}\Gamma'\text{--}D$) for $\alpha\text{-RuCl}_3$ and $\text{Na}_2\text{IrO}_3$.
- **Phase 5:** Three-Dimensional Tricoordinated Kitaev Lattices (hyperhoneycomb and stripyhoneycomb).

---

## Repository Architecture

The production Python package follows strict hexagonal (ports-and-adapters) separation:

- `src/kitaev_mf/domain/`: mathematical entities, Clifford CAR operators, lattice definitions, model parameters, and port protocols.
- `src/kitaev_mf/services/`: physics and numerical algorithms (Bloch Hamiltonians, spectral projectors, Grassmannian fixed-point solvers, real-space skew-symmetric operators, Newton–Schulz purification, observables, and energetics).
- `src/kitaev_mf/infra/`: concrete linear-algebra, hardware FLOP profiling (`ResourceProfiler`), and persistence adapters.
- `test/unit/` and `test/integration/`: **153 automated verification tests** (95% coverage) enforcing algebraic, kinematic, and numerical invariants.
- `doc/labbook/`: `memoir` LaTeX source, chapters, and PGFPlots figures for the 74-page Living Lab Book (`make labbook`).
- `doc/manuscript/`: SciPost LaTeX source for the 9-page Executive Monograph (`make manuscript`).
- `doc/rnd/`: Front-1 R&D technical notes (`phase_1/`), exact perturbation theory derivations (`exact-perturbation-theory/`), and peer-review records (`peer_review/`).

## Development & Verification Setup

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel pip-tools
pip-sync requirements.txt requirements-dev.txt
python -m pip install -e . --no-deps
```

Run the full verification and document build suite locally:

```bash
make verify      # ruff check, ruff format --check, mypy, and pytest (153 tests)
make labbook     # compiles doc/labbook/labbook_kitaev_mf.pdf (74 pages)
make manuscript  # compiles doc/manuscript/scipost_kitaev_mf.pdf (9 pages)
```

## Governance and Releases

The package version (`1.0.0`) is declared in `pyproject.toml`. See [`AGENTS.md`](AGENTS.md), [`CHANGELOG.md`](CHANGELOG.md), and [`.agents/rules/04_versioning_and_releases.md`](.agents/rules/04_versioning_and_releases.md) for the SemVer policy, Hypothesis Contract requirements, and release procedure.
