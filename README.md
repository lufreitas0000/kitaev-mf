# kitaev-mf: Majorana Mean-Field Solver for Extended Kitaev Systems

`kitaev-mf` is a high-performance computational research library for analytical and numerical studies of Majorana parton mean-field theories of Kitaev quantum magnets on honeycomb lattices.

---

## Strategic Multi-Phase Architecture

The scientific roadmap is structured across dedicated model families:

- **Phase 1: Pure Ferromagnetic Kitaev Honeycomb ($K = -1.0$) with $[111]$ Field — Numerical Quantification & Algorithmic Foundations (COMPLETED)**
  - Exhaustive characterization of nonlinear fixed-point solvers on the orthogonal Grassmannian $\Gr(4, 8)$ and real-space torus clusters: Damped Picard, Joint DIIS (Anderson Acceleration), Full Broyden, Limited-Memory Broyden ($\mathcal{O}(mD)$), JFNK, and eigenvector-free Matrix Sign Purification ($\Omega = \ii \operatorname{sgn}(\ii A)$).
  - Theoretical and empirical mastery of all numerical nuances: basin survival radii ($R_{90}$ with Wilson score CIs), critical slowing down at the spinodal $h_{c2}^+ = 0.400$, precision loss at Dirac mesh resonances ($3 \mid L$), and hardware resource throughput ($104.5\text{ GFLOP/s}$ via BLAS-3 `dgemm`).
  - Formal documentation in a 72-page Living Lab Book (`doc/labbook/labbook_kitaev_mf.pdf`), 11 Hypothesis Contracts (`ch_m1_11`), and 20 verified invariants (`K-001` through `K-020`).

- **Phase 2: Physics of the Ferromagnetic Kitaev $[111]$ Regime (ACTIVE / IN PREPARATION)**
  - Direct physical continuation using the optimal solvers established in Phase 1 (Joint DIIS and warm-started Newton-Schulz purification).
  - Mapping the comprehensive low-energy phase diagram and calculating all physical observables:
    - **$T = 0$ Ground State:** Order parameters, Wilson loop continuum, magnetization curves, and topological mass gap $\Delta_K(h) \propto h^3$.
    - **$T = 0$ Topological Defect State:** Real-space clusters with isolated pairs of $\mathbb{Z}_2$ vortices ($W_p = -1$), two-vison interaction potential $V(R)$, in-gap Majorana zero modes, and local electric polarization textures.
    - **Finite-Temperature States:** Thermal flux proliferation, vison crossover, and specific heat.
    - **Dynamical Correlations:** Dynamical structure factor $S^{\alpha\beta}(\mathbf{q}, \omega)$ and response functions.

- **Subsequent Phases:** Sequential investigation of distinct model families (e.g. Antiferromagnetic Kitaev $K > 0$, full $J$--$K$--$\Gamma$--$\Gamma'$--$D$ exchange, and 3D tricoordinated lattices).

---

## Architecture

The production Python package follows ports-and-adapters separation:

- `src/kitaev_mf/domain/`: mathematical entities, lattice definitions, model parameters, and port protocols.
- `src/kitaev_mf/services/`: physics algorithms operating on domain types and interfaces.
- `src/kitaev_mf/infra/`: concrete linear-algebra, hardware, and persistence adapters.
- `src/kitaev_mf/app/`: reserved for command-line and workflow orchestration if and when this layer is established.
- `test/unit/` and `test/integration/`: focused verification suites.
- `doc/adr/`: architecture decision records.
- `doc/rnd/`: exploratory derivations and numerical research notes.

The exact pure-Kitaev perturbation project is a separate research track. Its conventions and results are not automatically production interfaces or validated mean-field results.

## Development setup

Use the checked-in requirements specifications and lockfiles:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip setuptools wheel pip-tools
pip-sync requirements.txt requirements-dev.txt
python -m pip install -e . --no-deps
```

The canonical CI checks are defined in [`.github/workflows/ci.yml`](.github/workflows/ci.yml). Run the corresponding checks locally before opening a pull request:

```bash
ruff check src test
ruff format --check src test
mypy src test
pytest --cov=kitaev_mf --cov-report=term-missing test/
```

## Governance and releases

The package version is declared once in `pyproject.toml`. See
[`CHANGELOG.md`](CHANGELOG.md) and
[`.agents/rules/04_versioning_and_releases.md`](.agents/rules/04_versioning_and_releases.md)
for the SemVer policy, changelog format, and release procedure. All changes are integrated through reviewed pull requests; releases are made from `main`.
