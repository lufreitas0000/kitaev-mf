# Operating Protocol for AI Agents and Collaborators

Welcome to `kitaev-mf`. This repository is a research laboratory for Majorana parton mean-field theory on the extended Kitaev model ($J\text{-}K\text{-}\Gamma\text{-}\mathbf{h}$). All human and AI contributors must strictly adhere to this operating agreement.

---

## 1. Core Engineering Principles
1. **Zero-Shot Coding is Prohibited:** Before generating new modules, architectures, or non-trivial implementations, formulate the design, specify the mathematical contracts, and agree on the approach.
2. **Test-Driven Development (TDD):** Every domain rule or mathematical function must have a corresponding test in `test/unit/` or `test/integration/`. Never mark a milestone task complete without a passing test.
3. **Immutability & Pure Functions:** Numerical domain models in `src/domain/` must be frozen dataclasses or pure type representations. Physics operations in `src/services/` must be pure functions with zero hidden I/O side effects.
4. **Hardware Backend Discipline:** Linear algebra operations must pass through the `LinalgBackend` protocol (`src/infra/linalg_backend.py`) supporting both CPU and GPU execution without tensor allocation leakage.
5. **Atomic Data Storage:** Any write operation saving simulation results (`.h5`) must write to a temporary file (`.tmp_<hash>.h5`) and perform an atomic rename (`os.replace`) only upon successful completion.

---

## 2. Agent Handoff Protocol
At the beginning of any session:
1. Read `PROJECT_STATE.yaml` to identify the active sprint, current milestone, and verified equations.
2. Check `doc/adr/` for recent architectural decision records to avoid repeating discarded approaches.

At the end of any session:
1. Update `PROJECT_STATE.yaml` with the current status, new benchmarks, or blockers.
2. If an approach or algorithm fails, document the hypothesis, parameters, and failure mode in a new markdown file inside `doc/rnd/`.
3. If an architectural decision or mathematical convention is finalized, create a new numbered record in `doc/adr/ADR-XXX-<name>.md`.

---

## 3. Mathematical Notation Standards
* Spin Pauli matrices: $\sigma_i^\alpha$ with $\alpha \in \{1, 2, 3\} \equiv \{x, y, z\}$.
* Majorana generators: $\{c_i^\mu\}$ with $\mu \in \{0, 1, 2, 3\}$ satisfying $\{c_i^\mu, c_j^\nu\} = 2\delta_{ij}\delta^{\mu\nu}\mathbb{1}$.
* Gauge generators: $G^\gamma$ satisfying $T_i^\gamma = \frac{i}{4} c_i^T G^\gamma c_i$.
* Mean-field parameters: $U_{ij}^{\mu\nu} = \langle i c_i^\mu c_j^\nu \rangle$, $V_i^{\mu\nu} = \langle i c_i^\mu c_i^\nu \rangle$.
* Citations: Always use author-year format (e.g., Freitas et al. (2024), Kitaev (2006), Rau et al. (2014)).
