# kitaev-mf: Majorana Mean-Field Solver for Extended Kitaev Systems

`kitaev-mf` is an exact, high-performance numerical and analytical laboratory designed to rederive and generalize the results of **Freitas, Bauer, Egger, and Pereira, *npj Quantum Mater.* 9, 33 (2024)** [arXiv:2311.15750].

The engine computes self-consistent Majorana parton mean-field states for the extended $J\text{-}K\text{-}\Gamma\text{-}\mathbf{h}$ Kitaev model on honeycomb lattices with:
1. Rigorous $SO(4)$ gauge constraint projection ($c_i^T G^\gamma c_i = 0$).
2. Competing magnetic order multi-sublattice unit cells (Zigzag, Stripy, Néel, Polarized, KSL).
3. Inhomogeneous real-space geometries with mobile/metastable $\mathbb{Z}_2$ vison excitations.
4. Local electrostatic tip modulation ($V_0$) and induced Dzyaloshinskii-Moriya ($D(V_0)$) interactions.

## Architecture
The repository adheres to **Hexagonal Architecture (Ports and Adapters)**, strict **SOLID** design, and **TDD**:
- `src/domain/`: Pure mathematical types (lattice geometry, exchange couplings, mean-field states).
- `src/services/`: Pure physical algorithms (Hamiltonian builders, spectral projectors, Anderson Acceleration).
- `src/infra/`: Hardware adapters (PyTorch GPU/CPU dispatch, atomic HDF5 storage, plotting).
- `doc/notes/`: Formal, compilable LaTeX monographs matching publication standards.
- `doc/rnd/`: Markdown scratchpads for exploratory derivations and algorithm convergence logs.

## Setup & Testing
```bash
# Install dependencies in editable mode
pip install -e ".[dev]"

# Run test suite
pytest
```
