# Rule 00: Workspace Hygiene, Directory Taxonomy, and File Granularity

## 1. Zero Untracked Artifacts
1. The working tree must remain free of untracked temporary scripts, scratch plots, or debug dumps at the end of every task.
2. Ephemeral inspection files (such as `REPO_SNAPSHOT.md` generated via `make snapshot`) must be listed in `.gitignore` and never staged or committed.
3. All code generation or bulk file updates must be delivered via idempotent POSIX bash `cat << 'EOF'` heredocs.

## 2. Granularity and Small-File Discipline (Anti-Monolith Principle)
1. **Single Responsibility:** Every file—source code, test suite, script, or documentation note—must have a single, clearly delineated purpose.
2. **File Size Target:** Files should remain concise and modular (ideally under 150–200 lines). Monolithic all-in-one scripts or mega-notes are prohibited.
3. **Granular Verification:** Do not create monolithic multi-topic audit scripts in `scripts/`. Distinct tests belong in individual, focused test modules inside `test/unit/`.

## 3. Canonical Directory Taxonomy
- `.agents/rules/`: Normative operational, TDD, safety, and scientific editorial rules.
- `doc/adr/`: Immutable Architecture Decision Records (`ADR-NNN-<slug>.md`).
- `doc/rnd/`: Self-contained theoretical and numerical R&D notes (`note_NN_<slug>.md`).
- `doc/reference/`: Read-only reference literature and thesis LaTeX sources.
- `src/kitaev_mf/domain/`: Pure mathematical entities, lattice graphs, algebraic constants, and Port protocols. Zero dependencies on `services/` or `infra/`.
- `src/kitaev_mf/services/`: Pure physics algorithms. Depends exclusively on `domain/`.
- `src/kitaev_mf/infra/`: Driven hardware and I/O adapters.
- `src/kitaev_mf/app/`: Driving CLI and simulation workflow orchestrators.
- `test/unit/` & `test/integration/`: Mirroring `domain/`, `services/`, and `infra/`.
