# Rule 03: Scientific Editorial and Mathematical Rigor Standards

All documentation in `doc/rnd/`, `doc/adr/`, and module docstrings must adhere to the editorial and mathematical standards of *Physical Review Letters* (APS) and *Nature Physics*.

## 1. Deductive Mathematical Structure
Every theoretical R&D note (`doc/rnd/note_NN_<slug>.md`) must be **strictly self-contained** and follow the deductive progression:
1. **Context & Physical Motivation:** Define the physical regime, symmetries, and model assumptions.
2. **Definitions & Lattice/Operator Conventions:** Explicitly state the graph $\Lambda = (\mathcal{V}, \mathcal{E})$, Hilbert space, index sets, and Fourier conventions before writing derived equations.
3. **Lemmas & Propositions:** State intermediate algebraic identities (e.g., Clifford contractions, Casimir invariants, block symmetries).
4. **Theorems & Proofs:** Provide step-by-step analytical derivations without skipping non-trivial index contractions or sign permutations.
5. **Corollaries & Limiting Cases:** Verify how general expressions reduce to known exact limits (e.g., pure Kitaev limit $J = \Gamma = \Gamma' = \mathbf{h} = 0$).

## 2. Canonical Index and Notation Conventions
- **Lattice Graph:** $\Lambda = (\mathcal{V}, \mathcal{E})$ with vertex set $\mathcal{V}$ and directed edge set $\mathcal{E}$.
- **Vertex Indices:** Latin letters $i, j, k, l \in \mathcal{V}$ (never $u, v$). Nearest-neighbor set of vertex $i$ is $\mathcal{V}_i$.
- **Directed Bonds:** $\langle ij \rangle_\gamma \in \mathcal{E}$ denotes a nearest-neighbor bond of direction $\gamma \in \{1, 2, 3\} \equiv \{x, y, z\}$ oriented from sublattice $A$ (site $i$) to sublattice $B$ (site $j$).
- **Flavor & Spin Indices:** Greek indices $\alpha, \beta, \gamma \in \{1, 2, 3\}$ label Cartesian spin/bond directions; $\mu, \nu, \rho \in \{0, 1, 2, 3\}$ label all four Majorana flavors $c_i = (c_i^0, c_i^1, c_i^2, c_i^3)^T$.
- **Fixed Parton Representation:** We fix the Kitaev representation $\sigma_i^\gamma = \frac{i}{2} c_i^T \mathrm{N}^\gamma c_i = i c_i^0 c_i^\gamma$ where $\mathrm{N}^\gamma = \frac{1}{2}(\mathrm{M}^\gamma + \mathrm{G}^\gamma)$, paired with the local gauge constraint $T_i^\gamma = \frac{i}{8} c_i^T \mathrm{G}^\gamma c_i = 0$. Do not introduce artificial unphysical interpolation parameters ($\eta$) into the spin representation; reserve $\eta = J_H / U$ exclusively for the Hund-to-Hubbard ratio in the microscopic Schrieffer-Wolff expansion.
- **Typography:** Punctuate all display equations as part of the surrounding grammatical sentence. Use upright Roman font for operators/labels ($\mathrm{tr}$, $\mathrm{Pf}$, $\mathrm{d}k$, $\mathrm{MF}$, $\mathrm{BdG}$) and Blackboard bold for number fields and identity operators ($\mathbb{R}, \mathbb{C}, \mathbb{Z}_2, \mathbb{1}$).

## 3. Python Code Naming for Mathematical Physics
- In compliance with mathematical physics conventions, linear algebra matrices, operator representations, and exchange tensors MAY be named using uppercase single- or multi-letter symbols ($H, U, V, A, B, J, M, G, N$) in both argument lists and local scopes.
- Ruff rules `N802` (function names), `N803` (argument names), and `N806` (variable names) are explicitly ignored to prevent conflicts with standard Hamiltonian and covariance matrix notation.
- One-line compound statements with semicolons (`;`) and non-standard minified Python remain strictly prohibited (`E701`, `E702`).

## 4. Citation and Bibliographic Verification Protocol
1. **In-Text Citations:** Always use the strict `Author (Year)` or `(Author et al., Year)` format in prose (e.g., `Kitaev (2006)`, `(Freitas et al., 2024)`, `(Rau et al., 2014)`). Never use bare numbered brackets `[1]` in Markdown notes.
2. **Terminal Bibliography:** Every Markdown document in `doc/rnd/` and `doc/adr/` must end with a `## References` section containing the complete, verified bibliographic details:
   - Full author list, *Article or Book Title*, Journal Name **Volume**, Article ID / Page Range (Year), DOI URL, and arXiv identifier.
   - Every reference must be cross-checked against `doc/reference/thesis/refs.bib` or verified publisher metadata. Hallucinated titles, authors, or DOIs are strictly prohibited.
