# Rule 03: Scientific Documentation, Mathematical Rigor, and Three-Tier Architecture

Apply these standards to research notes, living lab books, architecture records, code documentation, and publication deliverables.

## 1. Three-tier documentation architecture

The repository maintains an explicit separation between the Context of Discovery and the Context of Justification:

### Tier 1: Context of Discovery (Machine-First Agent Logs)
- **Location**: `doc/rnd/phase_<N>/m<phase>_<milestone>_t<task>_<slug>.md`
- **Audience**: Automated agents and LLM thought partners.
- **Contract**: Every note must open with strict YAML frontmatter:
  ```yaml
  ---
  id: RND-P<phase>-M<milestone>-T<task>
  title: "<Descriptive Title>"
  phase: <int>
  milestone: "<milestone_id>"
  task: "<task_id>"
  status: "<planned|in_progress|completed>"
  contracts: ["<TST-...>"]
  claims: ["<CLM-...>"]
  artifacts: ["<data/processed/benchmarks/...>"]
  ---
  ```
- **Style**: Mathematically formal, explicit, and operational. Free of rhetorical introductions, personal development framing, or conversational filler. Must state: analytical contracts, exact parameter domains, observed contraction rates $\mu$, failed seeds, and numerical bounds.

### Tier 2: Context of Justification (Living Human-First Lab Book)
- **Location**: `doc/labbook/` compiling to `doc/labbook/main.pdf`
- **Audience**: The human investigator and project owner.
- **Contract**:
  - Class: `memoir` with soft parchment beige background (`#FAF6EE`) and dark charcoal font (`#2B2B2B`).
  - Chapter files are indexed by milestone: `chapters/ch_m<phase>_<milestone>_<slug>.tex`.
  - Updated incrementally at the completion of every milestone sprint.
  - Focuses on narrative pedagogical exposition, physical motivation, phase stability, and failure analysis.
  - All numerical plots and tables must load directly from `data/processed/benchmarks/*.dat` via `pgfplotstable` or deterministic figure generation—hardcoded numerical arrays in TeX are strictly forbidden.

### Tier 3: Phase Monograph (Context of Exposition)
- **Location**: `doc/report/` compiling to `doc/report/main.pdf`
- **Audience**: External peer review, academic archive, or thesis publication.
- **Contract**: Produced at the conclusion of a macro Phase (e.g. M1.9). Distills the verified results into a concise, definition-theorem-proof document linking claims directly to verified invariant tags (`CLM-M1-*`).

## 2. Mathematical claims and provenance

1. State definitions, domains, assumptions, conventions, and boundary conditions before relying on them.
2. Give proofs for nontrivial claims when the note presents them as established results. Identify omitted arguments and cite a verifiable source when relying on an external result.
3. Distinguish exact algebraic statements, analytical results, finite-volume statements, numerical evidence, approximations, extrapolations, and conjectures.
4. Use persistent claim identifiers (`CLM-M<phase>-...`) when a result is tracked across research notes, tests, project-state invariants, or publication sources.
5. A test passing establishes only its tested computational contract and parameter domain; it is not a proof of a theorem.

## 3. Conventions and notation

State conventions explicitly whenever they affect signs, prefactors, symmetry representations, or physical interpretation:

- lattice graph, sublattice and bond orientation;
- Fourier-transform and Brillouin-zone conventions;
- Majorana or parton representation and normalization;
- gauge constraints and physical-state projection;
- index ranges, summation conventions, and matrix ordering;
- boundary conditions and finite-size geometry.

For the production four-flavor representation, retain the established project notation:
- $\Lambda = (\mathcal{V},\mathcal{E})$ denotes the lattice graph;
- Latin indices $i, j$ label vertices; Greek indices $\alpha, \beta, \gamma \in \{1, 2, 3\}$ label Cartesian spin/bond directions;
- $c_i = (c_i^0, c_i^1, c_i^2, c_i^3)^T$ denotes the four-flavor Majorana multiplet;
- $\sigma_i^\gamma = \frac{i}{2}c_i^T\mathrm{N}^\gamma c_i = i c_i^0c_i^\gamma$ with $\mathrm{N}^\gamma = \frac{1}{2}(\mathrm{M}^\gamma + \mathrm{G}^\gamma)$;
- local gauge constraint: $T_i^\gamma = \frac{i}{8}c_i^T\mathrm{G}^\gamma c_i = 0$.

## 4. Single responsibility for numerical benchmark data

1. Numerical benchmark tables belong exclusively in `data/processed/benchmarks/<slug>.dat`.
2. Every `.dat` file must begin with comment headers stating the benchmark title, ISO timestamp, Git commit SHA, description, and column schema.
3. Plotting scripts, LaTeX tables (`pgfplotstable`), and regression tests must consume these `.dat` files directly. Do not duplicate floating-point data tables across multiple repositories or source files.
