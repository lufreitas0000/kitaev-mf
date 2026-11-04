# kitaev-mf

Majorana parton mean-field solver and observable framework for extended Kitaev models ($J\text{-}K\text{-}\Gamma\text{-}\Gamma'\text{-}D\text{-}\mathbf{h}$) on honeycomb lattices.

## Overview
This repository provides a modular, verified computational framework rederiving and generalizing the theoretical models of:
> **L. R. D. Freitas, T. Bauer, R. Egger, and R. G. Pereira**, *Electric polarization near vortices in the extended Kitaev model*, npj Quantum Materials **9**, 33 (2024). [DOI: 10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5).

## Current Status
- **Active Sprint:** `M1.4.4` (Production Observables and Energetics Services)
- **Branch:** `develop/m1-foundations-review`
- **Quality Gate:** Strict TDD, Ruff linting/formatting, and Mypy strict (`--strict`).

## Development & Verification
```bash
make verify
```
Runs the full test suite, high-resolution K-006 quadrature, Ruff format/lint checks, Mypy type-checking, and Git whitespace validation.
