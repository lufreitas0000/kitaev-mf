# Operating Protocol for AI Agents and Collaborators

Welcome to `kitaev-mf`. This repository is a research laboratory for Majorana parton mean-field theory (MMFT) on the extended Kitaev model ($J\text{-}K\text{-}\Gamma\text{-}\Gamma'\text{-}D\text{-}\mathbf{h}$), rederiving and generalizing **Freitas, Bauer, Egger, and Pereira, *npj Quantum Mater.* 9, 33 (2024)** [arXiv:2311.15750].

All human and AI contributors must strictly adhere to the modular rule sets in `.agents/rules/`:
- `.agents/rules/00_workspace_hygiene.md`: Directory taxonomy, zero untracked artifacts, and POSIX heredoc standards.
- `.agents/rules/01_repository_safety.md`: Branch discipline and atomic HDF5 data integrity.
- `.agents/rules/02_tdd_protocol.md`: Strict Red-Green-Refactor TDD workflow and `< 1e-14` precision gates.
- `.agents/rules/03_editorial_rigor.md`: APS *PRL*/*PRB* & *Nature* mathematical style, vertex/edge notation, and verified `Author (Year)` bibliography rules.
- `.agents/rules/04_milestone_progress.md`: mandatory atomic updates to PROJECT_STATE.yaml and doc/MILESTONE_ROADMAP.md in every commit.

---

## 1. Core Engineering & Architectural Principles
1. **Zero-Shot Coding is Prohibited:** Formulate the mathematical contracts, agree on the Hexagonal ports/adapters design, and write failing tests before implementing services.
2. **Hexagonal Architecture (Ports & Adapters):**
   - `src/kitaev_mf/domain/`: Pure mathematical entities (`UnitCell`, `Bond`, `ModelParameters`, `MeanFieldState`), algebraic constants, and Port protocols (`LinalgBackendPort`, `HamiltonianBuilderPort`). Zero dependencies on `services/` or `infra/`.
   - `src/kitaev_mf/services/`: Pure physics algorithms operating on `domain/` entities and Port interfaces.
   - `src/kitaev_mf/infra/`: Concrete hardware and storage adapters (`TorchBackend`, HDF5 writers).
3. **Hardware Backend Discipline:** Linear algebra operations pass through the `LinalgBackendPort` protocol supporting both CPU and CUDA (`complex128` / `float64`) without tensor allocation leakage.

---

## 2. Canonical Mathematical & Lattice Conventions
* **Lattice Graph:** $\Lambda = (\mathcal{V}, \mathcal{E})$ with vertices $i, j \in \mathcal{V}$ and directed bonds $\langle ij \rangle_\gamma \in \mathcal{E}$ oriented from sublattice $A$ to sublattice $B$. Primitive lattice constant $a = 1$ ($d_{\mathrm{NN}} = 1/\sqrt{3}$).
* **Majorana CAR Algebra:** $\{c_i^\mu, c_j^\nu\} = 2\delta_{ij}\delta^{\mu\nu}\mathbb{1}$ with $\mu, \nu \in \{0, 1, 2, 3\}$.
* **Fixed Kitaev Spin Representation:** $\sigma_i^\gamma = \frac{i}{2} c_i^T \mathrm{N}^\gamma c_i = i c_i^0 c_i^\gamma$, where $\mathrm{N}^\gamma = \frac{1}{2}(\mathrm{M}^\gamma + \mathrm{G}^\gamma)$.
* **Gauge Link Convention:** On directed bonds $i \in A \to j \in B$, the canonical zero-flux pure Kitaev state sets $\mathrm{U}_{\langle ij \rangle_\gamma}^{\gamma\gamma} = \langle i c_i^\gamma c_j^\gamma \rangle = +1.0$.
* **Momentum Phase Factors:** For fractional coordinates $(k_1, k_2) \in [0, 1)^2$ and bond cell offset $(n_1, n_2) \in \mathbb{Z}^2$, the Bloch phase is explicitly $\exp(i 2\pi (k_1 n_1 + k_2 n_2))$.
* **Bloch-Majorana Prefactor:** $\mathcal{H}(\mathbf{k}) = \frac{i}{4}\mathcal{M}(\mathbf{k})$, reproducing exact matter dispersion $E(\mathbf{k}) = \pm \frac{|K|}{2}|1 + e^{-i 2\pi k_1} + e^{-i 2\pi k_2}|$ ($E(\mathbf{0}) = 1.5|K|$) and flat gauge bands $\pm \frac{|K|}{2}|u_{\mathrm{matter}}|$.
