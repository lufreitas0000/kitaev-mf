# Rule 00: Workspace Hygiene and Directory Taxonomy

## 1. Zero Untracked Artifacts
1. The working tree must remain free of untracked temporary scripts, scratch plots, or debug dumps at the end of every task.
2. Ephemeral inspection files (such as `REPO_SNAPSHOT.md` generated via `make snapshot`) must be listed in `.gitignore` and never staged or committed.
3. All code generation or bulk file updates must be delivered via idempotent POSIX bash `cat << 'EOF'` heredocs.

## 2. Canonical Directory Taxonomy
- `.agents/rules/`: Normative operational, TDD, safety, and scientific editorial rules.
- `doc/adr/`: Immutable Architecture Decision Records (`ADR-NNN-<slug>.md`).
- `doc/rnd/`: Self-contained theoretical and numerical R&D notes (`note_NN_<slug>.md`).
- `doc/reference/`: Read-only reference literature and thesis LaTeX sources.
- `src/kitaev_mf/domain/`: Pure mathematical entities, lattice graphs, algebraic constants, and Port protocols. Zero dependencies on `services/` or `infra/`.
- `src/kitaev_mf/services/`: Pure physics algorithms (Hamiltonian builders, correlators, saddle-point solvers, observables). Depends exclusively on `domain/`.
- `src/kitaev_mf/infra/`: Driven hardware and I/O adapters (PyTorch CPU/CUDA backend, atomic HDF5 storage).
- `src/kitaev_mf/app/`: Driving CLI and simulation workflow orchestrators.
- `test/unit/` & `test/integration/`: Mirroring `domain/`, `services/`, and `infra/`.
