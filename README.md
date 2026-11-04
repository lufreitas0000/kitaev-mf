# kitaev-mf

Research code for analytical and numerical work on Majorana-parton mean-field descriptions of the extended Kitaev model on the honeycomb lattice.

The project separates exact finite-dimensional algebraic identities from the mean-field approximation used for lattice calculations. Local Majorana CAR, the SO(4) spin/gauge representation, physical projection identities, and finite-site exchange-operator relations are tested in finite Fock space. Lattice spectra, correlators, and self-consistent mean-field states are numerical outputs of the chosen approximation; they are not exact solutions of the interacting extended model.

## Current scope

- Local four-Majorana representation, SO(4) spin and gauge generators, and physical-sector projection.
- Exchange tensors and exact finite-site Majorana bilinear/quartic operator identities.
- Momentum-space Bloch Hamiltonians, spectral projectors, correlators, and saddle-point infrastructure.
- Reproducible K-006 pure-Kitaev correlator quadrature and convergence diagnostics.
- Gaussian observable and energetics contracts, with production integration developed in milestone order.

The pure-Kitaev reference calculation uses the zero-flux background and conventions documented under `doc/rnd/`. Numerical convergence estimates are empirical unless a rigorous error analysis is explicitly supplied.

## Repository layout

- `src/kitaev_mf/domain/`: mathematical constants, model data, and domain-level types.
- `src/kitaev_mf/services/`: numerical and analytical algorithms.
- `src/kitaev_mf/`: package implementation and adapters.
- `test/unit/`: unit and finite-Fock-space algebra tests.
- `scripts/`: reproducible command-line calculations.
- `doc/rnd/`: derivations, contracts, numerical provenance, and research notes.
- `PROJECT_STATE.yaml` and `doc/MILESTONE_ROADMAP.md`: machine-readable and human-readable development handover state.

## Development setup

Install the package with development dependencies in an isolated Python environment:

```sh
python -m pip install -e ".[dev]"
```

Inspect available targets:

```sh
make help
```

Run the configured local verification sequence:

```sh
make verify
```

This target runs the configured test suite, the K-006 numerical diagnostic, Ruff lint/format checks, Mypy, Git whitespace checks, and a final working-tree report. A successful run is local evidence for that environment; it does not replace independent analytical validation.

## Research conventions

Use the derivation notes and invariant provenance as the authority for signs, normalization, Fourier conventions, and physical-sector definitions. Distinguish per-bond, per-unit-cell, and per-site energies, and distinguish exact algebraic statements from numerical convergence evidence and mean-field approximations.
