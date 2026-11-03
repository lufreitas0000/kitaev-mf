# Operating Protocol for AI Agents and Collaborators

Welcome to `kitaev-mf`. This repository is a research laboratory for Majorana parton mean-field theory (MMFT) on the extended Kitaev model ($J\text{-}K\text{-}\Gamma\text{-}\Gamma'\text{-}D\text{-}\mathbf{h}$), rederiving and generalizing **Freitas, Bauer, Egger, and Pereira, *npj Quantum Mater.* 9, 33 (2024)** [arXiv:2311.15750].

All human and AI contributors must strictly adhere to the modular rule sets in `.agents/rules/`:
- `.agents/rules/00_workspace_hygiene.md`: Directory taxonomy, zero untracked artifacts, and POSIX heredoc standards.
- `.agents/rules/01_repository_safety.md`: Branch discipline and atomic HDF5 data integrity.
- `.agents/rules/02_tdd_protocol.md`: Strict Red-Green-Refactor TDD workflow and `< 1e-14` precision gates.
- `.agents/rules/03_editorial_rigor.md`: APS *PRL*/*PRB* & *Nature* mathematical style, vertex/edge notation, and verified `Author (Year)` bibliography rules.

---

## 1. Core Engineering & Architectural Principles
1. **Zero-Shot Coding is Prohibited:** Formulate the mathematical contracts, agree on the Hexagonal ports/adapters design, and write failing tests before implementing services.
2. **Hexagonal Architecture (Ports & Adapters):**
   - `src/kitaev_mf/domain/`: Pure mathematical entities (`UnitCell`, `Bond`, `ModelParameters`, `MeanFieldState`), algebraic constants, and Port protocols.
   - `src/kitaev_mf/services/`: Pure physics functions operating on `domain/` entities and Port interfaces.
   - `src/kitaev_mf/infra/`: Concrete hardware and storage adapters (`TorchBackend`, HDF5 writers).
3. **Hardware Backend Discipline:** Linear algebra operations must pass through the `LinalgBackend` protocol supporting both CPU and CUDA (`complex128` / `float64`) without tensor allocation leakage.

---

## 2. Agent Handoff & Context Snapshot Protocol
At the beginning of any session:
1. Run `make snapshot` to generate an up-to-date, uncomitted `REPO_SNAPSHOT.md` containing the repository tree, `PROJECT_STATE.yaml`, ADRs, R&D notes, and source code for AI context ingestion.
2. Review `PROJECT_STATE.yaml` and `doc/adr/` to identify active milestones and verified invariants.

At the end of any session:
1. Run `make lint` and `make test` to ensure 100% test and type-checking compliance.
2. Update `PROJECT_STATE.yaml` and remove any temporary snapshot files (`make clean`).

---

## 3. Canonical Mathematical Notation Summary
* **Lattice Graph:** $\Lambda = (\mathcal{V}, \mathcal{E})$ with vertices $i, j, k, l \in \mathcal{V}$ and directed bonds $\langle ij \rangle_\gamma \in \mathcal{E}$ ($\gamma \in \{1, 2, 3\} \equiv \{x, y, z\}$). Primitive lattice constant $a = \vert{}\mathbf{a}_1\vert{} = \vert{}\mathbf{a}_2\vert{} = 1$ ($d_{\mathrm{NN}} = 1/\sqrt{3}$).
* **Majorana CAR Algebra:** $\{c_i^\mu, c_j^\nu\} = 2\delta_{ij}\delta^{\mu\nu}\mathbb{1}$ with $\mu, \nu \in \{0, 1, 2, 3\}$.
* **Fixed Kitaev Spin Representation:** $\sigma_i^\gamma = \frac{i}{2} c_i^T \mathrm{N}^\gamma c_i = i c_i^0 c_i^\gamma$, where $\mathrm{N}^\gamma = \frac{1}{2}(\mathrm{M}^\gamma + \mathrm{G}^\gamma)$.
* **Gauge Constraint Generators:** $T_i^\gamma = \frac{i}{4} c_i^T \mathrm{G}^\gamma c_i = 0$ on $\mathcal{H}_{\mathrm{phys}}$ ($D_i = c_i^0 c_i^1 c_i^2 c_i^3 = +1$).
* **Ultralocal Mean-Field State:** Bond correlators $\mathrm{U}_{ij}^{\mu\nu} = \langle i c_i^\mu c_j^\nu \rangle$ (`shape: (N_bonds, 4, 4)`), on-site correlators $\mathrm{V}_i^{\mu\nu} = \langle i c_i^\mu c_i^\nu \rangle$ (`shape: (N_sites, 4, 4)`), and Lagrange multipliers $\lambda_i^\gamma$ (`shape: (N_sites, 3)`).
