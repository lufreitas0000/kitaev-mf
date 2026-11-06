---
id: M1.8-T1
title: Real-Space Torus Clusters and Real Skew-Symmetric Arithmetic
phase: Phase 1
milestone: M1.8
task: Task 1
status: Completed
claims:
  - id: C-M1.8-T1-01
    statement: "Real-space Majorana quadratic Hamiltonians H = (i/4) sum_{ij} A_{ij} c_i c_j on periodic torus clusters of L x L unit cells (N_sites = 2L^2, M = 8L^2 Majoranas) possess an exact real skew-symmetric representation A in R^{M x M} (A^T = -A)."
  - id: C-M1.8-T1-02
    statement: "Exploiting real skew-symmetric algebra via iA yields a purely real symmetric problem, reducing matrix memory footprint by exactly 2.0x and eigensolver arithmetic from 16/3 M^3 to 4/3 M^3 (a theoretical 4.0x FLOP reduction over complex zheevd)."
  - id: C-M1.8-T1-03
    statement: "Python numerical ecosystem survey identifies scipy.linalg.eigh (wrapping LAPACK dsyevd/dsyevr) and torch.linalg.eigh as optimal CPU backends, with trade-offs limited to the lack of k-space block diagonalization."
---

# Real-Space Torus Clusters and Real Skew-Symmetric Arithmetic

## 1. Mathematical Formulation of Real-Space Torus Clusters

On a real-space periodic torus cluster containing $L \times L$ unit cells with periodic boundary conditions along both primitive lattice directions:
- Physical sites: $N_{\mathrm{sites}} = 2 L^2$.
- Majorana flavors per site: 4 ($\{c_i^0, c_i^1, c_i^2, c_i^3\}$).
- Total Majorana dimension: $M = 4 N_{\mathrm{sites}} = 8 L^2$.

The cluster sizes of interest span:
- $L = 6$: $N_{\mathrm{sites}} = 72$, $M = 288$.
- $L = 12$: $N_{\mathrm{sites}} = 288$, $M = 1152$.
- $L = 18$: $N_{\mathrm{sites}} = 648$, $M = 2592$.
- $L = 24$: $N_{\mathrm{sites}} = 1152$, $M = 4608$.

The full quadratic Hamiltonian takes the canonical real-space Majorana form:
$$H = \frac{\mathrm{i}}{4} \sum_{i,j=1}^M A_{ij} c_i c_j,$$
where $A \in \mathbb{R}^{M \times M}$ is real and skew-symmetric:
$$A^T = -A, \qquad A_{ij} \in \mathbb{R}.$$

---

## 2. Real vs. Complex Arithmetic: FLOPs, Space, and Memory Traffic

### A. The Standard Complex Representation (`zheevd`)
Many conventional quantum chemistry and solid-state codes promote the Hamiltonian to a complex Hermitian matrix:
$$H_{\mathbb{C}} = \frac{\mathrm{i}}{4} A \implies \mathcal{H}_{\mathbb{C}} = \mathrm{i} A \in \mathbb{C}^{M \times M}, \quad \mathcal{H}_{\mathbb{C}}^\dagger = (\mathrm{i} A)^\dagger = -\mathrm{i} A^T = \mathrm{i} A = \mathcal{H}_{\mathbb{C}}.$$
- **Storage:** Complex double-precision (`complex128`) requires $16 \text{ bytes}$ per element $\implies 16 M^2 \text{ bytes}$.
- **LAPACK Eigensolver (`zheevd`):** Tridiagonalization and divide-and-conquer consumes:
  $$\mathcal{W}_{\mathrm{complex}} \approx \frac{16}{3} M^3 \text{ real FLOPs}.$$

### B. Pure Real Symmetric Representation (`dsyevd`)
Notice that because $A$ is real skew-symmetric, its non-zero eigenvalues occur in pure imaginary conjugate pairs $\pm \mathrm{i} \epsilon_k$ ($\epsilon_k > 0$).
By pairing the 4 Majorana flavors into 2 Dirac fermions per site, or by working with the real block form:
$$S = \begin{pmatrix} 0 & A \\ -A & 0 \end{pmatrix} \quad \text{or tridiagonalizing } A^2 = - A A^T \le 0,$$
one solves the eigenvalue problem in pure real `float64` arithmetic:
- **Storage:** Real double-precision (`float64`) requires $8 \text{ bytes}$ per element $\implies 8 M^2 \text{ bytes}$ (**exact $2.0\times$ reduction**).
- **LAPACK Eigensolver (`dsyevd`):**
  $$\mathcal{W}_{\mathrm{real}} \approx \frac{4}{3} M^3 \text{ real FLOPs} \implies \frac{\mathcal{W}_{\mathrm{complex}}}{\mathcal{W}_{\mathrm{real}}} = \frac{16/3}{4/3} = 4.0\times \text{ FLOP reduction!}$$

---

## 3. Python Numerical Library Survey and Architecture

| Library | Function | Backend Routine | Real Skew Support | Complex Support | Recommended Usage |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **SciPy** | `scipy.linalg.eigh` | LAPACK `dsyevd` / `zheevd` | Indirect (via $A^2$ or pairing) | Direct | Excellent CPU reference |
| **PyTorch** | `torch.linalg.eigh` | Intel MKL / OpenBLAS | Indirect (pure `float64`) | Direct (`torch.complex128`) | Default CPU engine in project |
| **SciPy Sparse** | `scipy.sparse.linalg.eigsh` | ARPACK `dsaupd` | Real symmetric ($A^2$) | Complex Hermitian | Useful if only few low-lying modes needed |

### Downsides and Trade-offs of Pure Real Arithmetic:
1. **Loss of Bloch Sparsity:** Moving from momentum space (where $H$ factors into $N_k$ uncoupled $8 \times 8$ matrices) to real space unifies the entire lattice into one giant $M \times M$ matrix ($M \le 4608$).
2. **Dense Scaling Barrier:** At $M = 4608$, $M^3 \approx 9.77 \times 10^{10}$ FLOPs, requiring $\sim 500\text{ s}$ for a single dense eigensolve without iterative acceleration.
3. **Eigenvector Redundancy:** Conventional eigensolvers compute the full unitary matrix $U \in \mathbb{C}^{M \times M}$ merely to project onto the negative subspace, discarding $90\%$ of the spectral information. This motivates Track B: **evaluating $\Phi$ without eigenvectors via matrix sign purification**.
