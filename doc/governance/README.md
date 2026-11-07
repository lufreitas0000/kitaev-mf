# Repository Governance Architecture

**Status:** Operational baseline
**Scope:** Repository organization, document authority, and engineering workflow

This document defines how repository guidance, executable configuration, scientific records, and project planning relate to one another.

## W001 accepted transition

ADR-002 and Rules 05/09 govern new waves, multi-tag classification, branch targets,
claim acceptance, and publication ownership. Phase-indexed areas below remain
historical locations. New R&D notes use `doc/rnd/waves/`; their creation is task
scoped. The canonical existing labbook root is
`doc/labbook/labbook_kitaev_mf.tex`; `main.pdf` below is a stale legacy example.
Current publication sources live in `doc/manuscript/`; `doc/report/` below is
historical/planned, not a second required book. Manuscript filename/PRL migration
is W001-T6 and has not yet run. See `doc/plans/w001_consolidation.md`.

These rules are project-local. The workspace synchronizer has no project selector;
only a dry-run preview is appropriate here, with no propagation to other projects.

## 1. Document authority and precedence

| Artifact | Authoritative for | Not authoritative for |
|---|---|---|
| GitHub repository settings | Server permissions and branch protections | Local development workflows |
| CI workflows and Makefile | Executable checks and validation sequences | Requirements not encoded in scripts |
| `pyproject.toml` | Package version and tool configuration | Scientific derivations |
| `.agents/rules/` | Agent and contributor operating policy | Proof of general mathematical claims |
| `PROJECT_STATE.yaml` | Concise durable project state and active sprint tasks | Live CI status or transient git state |
| `doc/MILESTONE_ROADMAP.md` | Strategic milestone roadmap and backlog | Ephemeral daily task bookkeeping |
| `doc/rnd/phase_<N>/` | Tier-1 machine-first R&D notes (Context of Discovery) | Polished external publications |
| `doc/labbook/` | Tier-2 human-first Living Lab Book (Context of Justification) | Raw unverified machine logs |
| `doc/report/` | Tier-3 publication monograph (Context of Exposition) | In-progress daily exploration |
| `data/processed/benchmarks/` | Authoritative numerical benchmark tables (`.dat`) | Raw heavy binary arrays |

## 2. Canonical repository areas

| Path | Responsibility |
|---|---|
| `data/processed/benchmarks/` | Single source of truth for ASCII `.dat` benchmark tables with provenance headers |
| `data/raw/` | Heavy binary array checkpoints and raw wavefunction tensors (DVC-tracked) |
| `data/manifests/` | Scan trajectory configurations and run catalog manifests |
| `doc/labbook/` | Living human-readable `memoir` LaTeX labbook (`main.pdf`) |
| `doc/rnd/phase_<N>/` | Phase-indexed machine-first agent R&D notes with YAML metadata |
| `doc/rnd/exploratory/` | General exploratory surveys and literature calculations |
| `doc/report/` | Formal phase monograph deliverable |
| `src/kitaev_mf/domain/` | Mathematical entities, parameters, lattice graphs, and port protocols |
| `src/kitaev_mf/services/` | Pure physics algorithms and computational services |
| `src/kitaev_mf/infra/` | Hardware adapters (Torch CPU/CUDA) and tabular `.dat` file I/O |
| `src/kitaev_mf/app/` | CLI trajectory driver and catalog query interface |
| `test/unit/` | Mathematical and algebraic unit verification suites |
| `test/integration/` | Variational response, convergence, and duality integration suites |
| `scripts/` | Deterministic automation tools (benchmarks, state management, snapshots) |
