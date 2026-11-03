# Rule 03: Scientific Editorial, Provenance, and Mathematical Rigor Standards

All documentation in `doc/rnd/`, `doc/adr/`, and final LaTeX deliverables must adhere to the editorial and mathematical standards of *Physical Review Letters* (APS) and *Nature Physics*.

---

## 1. Dual-Track Epistemology: R&D Lab Notebook vs. Canonical Monograph

A strict separation of concerns is maintained between exploratory derivations and final publications:

1. **Track 1: Exploratory R&D Notes (`doc/rnd/`):**
   - Act as an internal, modular laboratory notebook capturing derivations, algorithmic edge cases, and numerical experiments.
   - Progression is non-linear and modular. Notes do **not** map 1-to-1 to chapters of the final report.
   - Internal numbering (Definitions, Lemmas, Theorems) is strictly local to each note. Cross-note structural coupling (e.g., numbering Note 02 items to match future Chapter 2) is forbidden.

2. **Track 2: Canonical LaTeX Monograph Deliverable (`doc/report/` or `doc/monograph/`):**
   - The final milestone deliverable is a self-contained, publication-grade scientific document.
   - It follows an optimal, linear pedagogical narrative designed for an external reader.
   - It contains zero references to intermediate scratchpad notes (`doc/rnd/`).
   - All environments (`\begin{definition}`, `\begin{theorem}`) follow standard LaTeX counter sequences governed exclusively by the monograph's own outline.

---

## 2. Claim Provenance Protocol (Persistent Claim IDs)

To maintain end-to-end traceability between laboratory exploration, numerical test suites, and canonical publications without structural coupling:

1. **Claim ID Tagging in R&D Notes:**
   Every formal mathematical statement (Definition, Lemma, Theorem, Proposition, Corollary) in `doc/rnd/` must carry a persistent, immutable claim identifier in an HTML comment tag immediately following the statement header:
   ```markdown
   **Theorem (Casimir–Parity Identity)** <!-- claim: CLM-M1-ALG-002 -->
   Let $\mathcal{Z}_i := c_i^0 c_i^1 c_i^2 c_i^3$...
   ```
   - Schema: `CLM-<MILESTONE>-<DOMAIN>-<ID>` (e.g., `CLM-M1-ALG-001`, `CLM-M1-BLO-003`, `CLM-M1-COR-002`).
   - Claim IDs are immutable: they never change even if the note is revised, reordered, or split.

2. **Provenance Metadata in LaTeX Sources:**
   When claims are synthesized into the final LaTeX deliverable, their Claim ID, R&D origin, and unit test verification suite are recorded in comment metadata directly preceding the LaTeX environment:
   ```latex
   % provenance: CLM-M1-ALG-002
   % source: doc/rnd/note_01_algebraic_foundations.md
   % verified_by: test/unit/test_algebra.py::test_casimir_parity_and_pauli
   \begin{theorem}[Casimir--Parity Identity]
   \label{thm:casimir_parity}
   Let $\mathcal{Z}_i = c_i^0 c_i^1 c_i^2 c_i^3$...
   \end{theorem}
   ```

3. **Project State Tracking:**
   Key benchmark invariants in `PROJECT_STATE.yaml` (e.g., `K-001`, `K-002`) link directly to the corresponding `CLM` identifiers and pytest targets.

---

## 3. Deductive Mathematical Structure in R&D Notes
Every theoretical R&D note (`doc/rnd/note_NN_<slug>.md`) must remain self-contained within its specific scope and follow the deductive progression:
1. **Context & Conventions:** Explicitly state index sets, graph conventions, and Fourier definitions.
2. **Lemmas & Intermediate Steps:** Document intermediate algebraic contractions and Lie brackets.
3. **Theorems & Proofs:** Provide step-by-step analytical derivations.
4. **Limiting Invariants:** Explicitly verify how general expressions reduce to solvable limits (such as the pure Kitaev limit $J = \Gamma = \Gamma' = D = \mathbf{h} = 0$).

---

## 4. Canonical Index and Notation Conventions
- **Lattice Graph:** $\Lambda = (\mathcal{V}, \mathcal{E})$ with vertex set $\mathcal{V}$ and directed edge set $\mathcal{E}$.
- **Vertex Indices:** Latin letters $i, j, k, l \in \mathcal{V}$. Nearest-neighbor set of vertex $i$ is $\mathcal{V}_i$.
- **Directed Bonds:** $\langle ij \rangle_\gamma \in \mathcal{E}$ denotes a nearest-neighbor bond of type $\gamma \in \{1, 2, 3\} \equiv \{x, y, z\}$ oriented from sublattice $A$ (site $i$) to sublattice $B$ (site $j$).
- **Flavor & Spin Indices:** Greek indices $\alpha, \beta, \gamma \in \{1, 2, 3\}$ label Cartesian spin/bond directions; $\mu, \nu, \rho \in \{0, 1, 2, 3\}$ label all four Majorana flavors $c_i = (c_i^0, c_i^1, c_i^2, c_i^3)^T$.
- **Fixed Parton Representation:** Fix the Kitaev representation $\sigma_i^\gamma = \frac{i}{2} c_i^T \mathrm{N}^\gamma c_i = i c_i^0 c_i^\gamma$ where $\mathrm{N}^\gamma = \frac{1}{2}(\mathrm{M}^\gamma + \mathrm{G}^\gamma)$, paired with the local gauge constraint $T_i^\gamma = \frac{i}{8} c_i^T \mathrm{G}^\gamma c_i = 0$. Reserve $\eta$ exclusively for Bogoliubov quasiparticles ($\eta_{n,\mathbf{k}}$) or the Hund-to-Hubbard ratio ($J_H/U$).
- **Typography:** Punctuate all display equations as part of the surrounding grammatical sentence. Use upright Roman font for operators/labels ($\operatorname{tr}$, $\operatorname{Pf}$, $\mathrm{d}k$, $\mathrm{MF}$, $\mathrm{BdG}$) and Blackboard bold for number fields and identity operators ($\mathbb{R}, \mathbb{C}, \mathbb{Z}_2, \mathbb{1}$).

---

## 5. Python Code Naming for Mathematical Physics
- In compliance with mathematical physics conventions, linear algebra matrices, operator representations, and exchange tensors MAY be named using uppercase single- or multi-letter symbols ($H, U, V, A, B, J, M, G, N$) in both argument lists and local scopes.
- Ruff rules `N802` (function names), `N803` (argument names), and `N806` (variable names) are explicitly ignored in `pyproject.toml` to prevent conflicts with standard Hamiltonian and covariance matrix notation.
- One-line compound statements with semicolons (`;`) and non-standard minified Python remain strictly prohibited (`E701`, `E702`).

---

## 6. Citation and Bibliographic Verification Protocol
1. **In-Text Citations:** Always use the strict `Author (Year)` or `(Author et al., Year)` format in prose (e.g., `Kitaev (2006)`, `(Freitas et al., 2024)`, `(Rau et al., 2014)`). Never use bare numbered brackets `[1]` in Markdown notes.
2. **Terminal Bibliography:** Every Markdown document in `doc/rnd/` and `doc/adr/` must end with a `## References` section containing the complete, verified bibliographic details cross-checked against publisher DOIs and `doc/reference/thesis/refs.bib`.
