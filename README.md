# kitaev-mf

A Python research codebase for Majorana-parton mean-field calculations in
extended Kitaev models on the honeycomb lattice. The repository also hosts a
separate analytical research track on perturbation theory around the pure
Kitaev spin liquid; that track is documented independently and is not yet a
production solver.

## Repository layout

The current layout on this branch is:

- `src/domain/`: domain-level mathematical objects and model parameters.
- `src/infra/`: numerical infrastructure, including the linear-algebra backend.
- `test/`: unit tests.
- `doc/`: project documentation, derivations, and research notes.
- `doc/rnd/exact-perturbation-theory/`: isolated research on the pure Kitaev
  model, physical flux sectors, and low-energy perturbation theory.

This branch retains the earlier `src/domain` / `src/infra` layout. The
production architecture on `main` has since moved to
`src/kitaev_mf/{domain,services,infra,app}`; synchronizing that architecture
with this research branch is a separate integration task. Do not assume the
new package paths exist on this branch yet.

## Development setup

Python 3.10 or newer is required. Install the project and its development
dependencies in a virtual environment:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ".[dev]"
```

The development extra provides pytest, pytest-cov, Ruff, and mypy. The current
CI checks are:

```bash
ruff check .
mypy src test
pytest --cov=src --cov-report=term-missing
```

Run the checks from the repository root. CI configuration is in
`.github/workflows/ci.yml`.

## Exact pure-Kitaev perturbation research

The independent research track is maintained under
[`doc/rnd/exact-perturbation-theory/`](doc/rnd/exact-perturbation-theory/README.md).

Its scope includes the unperturbed matter spectrum, flux-sector energetics,
physical gauge projection and matter parity, and the resolvent domain for
low-energy effective descriptions. Analytical derivations are primary;
numerical calculations are intended as independent finite-size validation.
No production Python modules for this research track have been introduced.

Start with its [research status](doc/rnd/exact-perturbation-theory/STATUS.md),
[roadmap](doc/rnd/exact-perturbation-theory/ROADMAP.md), and
[backlog](doc/rnd/exact-perturbation-theory/BACKLOG.md).

## Scientific conventions

The pure-Kitaev research notes state their Hamiltonian normalization,
Majorana convention, bond orientation, and physical constraint explicitly.
Those conventions have not yet been mapped at operator level to the
production mean-field Bloch convention. Do not compare coupling prefactors or
spectral scales across the two tracks until that mapping is established.
