# Rule 03: Scientific Documentation and Mathematical Rigor

Apply these standards to research notes, architecture records, code documentation, and publication-facing material. The intended level is rigorous and readable scientific exposition; journal-specific formatting is governed by the target venue when a manuscript is prepared.

## 1. Research notes and publication-facing documents

1. `doc/rnd/` is the exploratory record for derivations, conventions, algorithmic edge cases, numerical experiments, failed approaches, and research decisions. Notes may develop non-linearly and should be self-contained within their stated scope.
2. Publication-facing documents must provide a coherent narrative for an external reader. They must not depend on unpublished scratch notes for definitions, proofs, or essential reasoning.
3. A publication document may cite or acknowledge a research note only when that is appropriate for its intended audience and the note is an accessible, stable reference. Do not impose a blanket rule that all internal notes must be cited or that all such references are forbidden.
4. Do not assume a canonical monograph directory exists. Use the location established by the relevant project deliverable and update the repository documentation when that structure is adopted.
5. Keep numbering of definitions, lemmas, equations, and theorems local to the document or explicitly defined document series.

## 2. Mathematical claims and provenance

1. State definitions, domains, assumptions, conventions, and boundary conditions before relying on them.
2. Give proofs for nontrivial claims when the note presents them as established results. Identify omitted arguments and cite a verifiable source when relying on an external result.
3. Distinguish exact algebraic statements, analytical results, finite-volume statements, numerical evidence, approximations, extrapolations, and conjectures.
4. Use persistent claim identifiers when a result is tracked across research notes, tests, project-state invariants, or publication sources. Do not require an identifier for every routine definition or statement.
5. When a tracked claim is transferred to a publication source, record its identifier and relevant source/test provenance in comments or project metadata where this improves traceability. A test is evidence for its tested computational contract; it is not a proof of a theorem.
6. Record limiting-case checks and independent validations when they materially support the result. Explain what those checks do and do not establish.

## 3. Conventions and notation

State conventions explicitly whenever they affect signs, prefactors, symmetry representations, or physical interpretation. In particular, document as relevant:

- lattice graph, sublattice and bond orientation;
- Fourier-transform and Brillouin-zone conventions;
- Majorana or parton representation and normalization;
- gauge constraints and physical-state projection;
- index ranges, summation conventions, and matrix ordering;
- boundary conditions and finite-size geometry.

For the production four-flavor representation, retain the established project notation unless an explicit, reviewed convention change is made. In the notation used by the current R&D notes:

- `\Lambda = (\mathcal{V},\mathcal{E})` denotes the lattice graph;
- Latin indices label vertices; Greek indices may label Cartesian spin/bond directions, with their ranges stated locally;
- `c_i=(c_i^0,c_i^1,c_i^2,c_i^3)^T` denotes the four-flavor Majorana vector;
- the representation `\sigma_i^\gamma = \frac{i}{2}c_i^T\mathrm{N}^\gamma c_i = i c_i^0c_i^\gamma`, with `\mathrm{N}^\gamma=\frac12(\mathrm{M}^\gamma+\mathrm{G}^\gamma)`, and the associated constraint `T_i^\gamma=\frac{i}{8}c_i^T\mathrm{G}^\gamma c_i=0`, are convention-specific. State the convention and verify it against the defining algebra before reuse.
- reserve `\eta` for the established project uses (Bogoliubov quasiparticles or the Hund-to-Hubbard ratio) unless a local notation declaration explicitly overrides it.

Use upright Roman type for named operators and labels where conventional (for example, `\operatorname{tr}`, `\operatorname{Pf}`, `\mathrm{d}k`, `\mathrm{MF}`, `\mathrm{BdG}`) and blackboard bold for number fields and identity operators (for example, `\mathbb{R}`, `\mathbb{C}`, `\mathbb{Z}_2`, `\mathbb{1}`). Punctuate displayed equations as part of the surrounding sentence.

## 4. Code and mathematical notation

1. Use descriptive English identifiers by default. Conventional uppercase names for matrices, Hamiltonians, exchange tensors, and operator representations are acceptable when they improve correspondence with the mathematical notation.
2. Ruff naming exceptions are defined in `pyproject.toml`; do not assume additional exceptions. Do not suppress lint rules locally without a concrete reason.
3. Avoid semicolon-compressed statements and minified code. Preserve readability and explicit control flow.

## 5. References and bibliographic integrity

1. Cite primary literature, peer-reviewed work, authoritative books, and original software documentation as appropriate to the claim.
2. Verify bibliographic metadata against publisher records, DOI registries, arXiv records, or the cited source itself. Never invent or infer a DOI, arXiv identifier, title, author list, page range, or publication year.
3. Use a consistent citation style appropriate to the document. Author-year citations are suitable for Markdown research notes; follow the target journal or publisher style for manuscripts.
4. Include a references section in research notes and architecture records when they rely on external literature or documentation. Do not add an empty bibliography to documents with no external references.
5. The thesis bibliography at `doc/reference/thesis/refs.bib` is a useful local reference source, not a substitute for verifying the metadata of a new citation.
