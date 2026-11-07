# Rule 00: Workspace Hygiene and Repository Layout

> W001 transition: Rules 05 and 09 govern new wave/task names, note paths,
> research acceptance, and labbook ownership. Phase/milestone examples below
> describe legacy records; preserve their identifiers rather than bulk-renaming.

## 1. Working-tree hygiene

1. Inspect `git status --short --branch` and the existing diff before editing. Preserve unrelated user changes.
2. At task completion, remove only temporary files created by the task when they are no longer needed. Do not require a globally clean working tree: intentional untracked research files, user data, and pre-existing changes must be preserved.
3. Generated artifacts belong in ignored build/cache locations where practical. Do not add broad ignore rules merely to hide unexplained files; inspect and classify each artifact.
4. Never stage or commit generated snapshots, debug dumps, temporary scripts, caches, or raw simulation output unless a reviewed deliverable explicitly requires them.
5. Do not run cleanup commands that may remove user data or pre-existing files without explicit authorization. See Rule 01.

## 2. File and module granularity

1. Give each source, test, script, and documentation file a coherent primary responsibility.
2. Arbitrary or unindexed note filenames (such as `note_04g.md`) are strictly prohibited. New notes must carry a deterministic wave/task index (`doc/rnd/waves/w001_t1_<slug>.md`). Existing phase/milestone filenames (e.g. `doc/rnd/phase_1/m1_4_5_t2_kkt_stationarity_and_response_asymptotics.md`) remain stable historical paths.
3. Keep tests focused on a behavior or contract. Use separate test modules when their subject matter or fixtures differ materially.
4. Avoid monolithic audit scripts. Reuse the repository's existing scripts and test structure before adding new tooling.

## 3. Repository layout

The repository layout is organized by the three-tier documentation and hexagonal package architecture:

- `.agents/rules/`: Normative operational, engineering, testing, safety, and scientific-documentation policies.
- `data/processed/benchmarks/`: Plain text `.dat` tables with SHA-256 provenance headers, serving as the single source of truth for benchmark figures and tests.
- `build/`: Ephemeral intermediate logs, run trajectories (`runs.jsonl`), quick-look diagnostic plots, and build products (strictly git-ignored).
- `data/raw/` and `data/manifests/`: Optional storage for external heavy array checkpoints and scan manifests (not used for routine milestones).
- `doc/adr/`: Architecture decision records, named `ADR-NNN-<slug>.md`.
- `doc/governance/`: Governance architecture and repository-level policy records.
- `doc/rnd/phase_<N>/`: Tier-1 Context of Discovery (machine-first R&D logs indexed as `m<phase>_<milestone>_t<task>_<slug>.md` with mandatory YAML frontmatter).
- `doc/rnd/exploratory/`: General surveys, literature comparisons, and exploratory calculations.
- `doc/labbook/`: Tier-2 Context of Justification (living human-first `memoir` LaTeX document with `#FAF6EE` background, `#2B2B2B` charcoal text, compiling into `doc/labbook/main.pdf`).
- `doc/report/`: Tier-3 Phase Monograph (formal publication deliverable compiled into `doc/report/main.pdf` at the close of a phase).
- `doc/reference/`: Reference material, including imported thesis sources.
- `src/kitaev_mf/domain/`: Mathematical entities, lattice structures, parameters, state, and port protocols. Keep independent of concrete services and infrastructure.
- `src/kitaev_mf/services/`: Physics algorithms and computational services operating on domain abstractions.
- `src/kitaev_mf/infra/`: Concrete numerical backends, hardware adapters, and tabular data I/O.
- `src/kitaev_mf/app/`: Command-line interface and trajectory orchestration.
- `test/unit/`: Focused unit tests verifying mathematical and domain contracts.
- `test/integration/`: System convergence, response scaling, and duality verification suites.
- `scripts/`: Repository utilities such as snapshot generation, benchmark export, and state management.

Do not create directories solely to make the tree conform to this description. Follow existing package and test conventions, and update this rule when the architecture materially changes.
