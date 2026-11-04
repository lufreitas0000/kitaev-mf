# kitaev-mf: Majorana Mean-Field Solver for Extended Kitaev Systems

`kitaev-mf` is a research codebase for numerical and analytical work on Majorana parton mean-field descriptions of extended Kitaev models on honeycomb lattices.

## Architecture

The production Python package follows ports-and-adapters separation:

- `src/kitaev_mf/domain/`: mathematical entities, lattice definitions, model parameters, and port protocols.
- `src/kitaev_mf/services/`: physics algorithms operating on domain types and interfaces.
- `src/kitaev_mf/infra/`: concrete linear-algebra, hardware, and persistence adapters.
- `src/kitaev_mf/app/`: command-line and workflow orchestration.
- `test/unit/` and `test/integration/`: focused verification suites.
- `doc/adr/`: architecture decision records.
- `doc/rnd/`: exploratory derivations and numerical research notes.

The exact pure-Kitaev perturbation project is a separate research track under
[`doc/rnd/exact-perturbation-theory/`](doc/rnd/exact-perturbation-theory/README.md).
Its conventions and results are not automatically production interfaces or validated mean-field results.

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
