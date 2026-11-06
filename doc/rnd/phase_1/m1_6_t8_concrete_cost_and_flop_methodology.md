---
id: RND-P1-M1.6-T8
title: "Concrete Computational Complexity, FLOP Formulas, and Memory Budgets"
phase: 1
milestone: M1.6
task: T8
status: completed
claims:
  - CLM-M1-FLOP-001
  - CLM-M1-FLOP-002
  - CLM-M1-FLOP-003
---

# Concrete Computational Complexity, FLOP Formulas, and Memory Budgets

## 1. Concrete FLOP Formulas for the Spectral Map $\Phi$

Rather than relying on vague asymptotic $\mathcal{O}(\cdot)$ notations, we derive the exact non-asymptotic floating-point operation (FLOP) count for evaluating $\Phi(\Omega, \boldsymbol{\lambda})$.

### 1.1 Bloch Momentum Space Evaluation
For a primitive cell with $n_c$ physical sites and $M_c = 4n_c$ Majorana flavors integrated over a 2D Brillouin zone mesh of $N_k = L_k^2$ points:

1. **Hamiltonian Assembly:**
   Building the $M_c \times M_c$ complex skew-Hermitian Bloch matrix $\mathcal{H}(\mathbf{k}) = \mathrm{i} A(\mathbf{k})$ for $N_k$ points requires evaluating matrix elements for $N_b = \frac{3}{2}n_c$ bonds:
   $$
   \mathrm{FLOPs}_{\mathrm{build}} = N_k \left( 36 N_b + 12 N_s \right) = N_k (66 n_c) \text{ FLOPs}.
   $$

2. **Complex Hermitian Eigendecomposition (`zheevd`):**
   Standard divide-and-conquer eigensolvers require:
   $$
   \mathrm{FLOPs}_{\mathrm{eig}} \approx \frac{16}{3} M_c^3 \cdot N_k = \frac{16}{3} (4 n_c)^3 N_k = \frac{1024}{3} n_c^3 N_k \text{ FLOPs}.
   $$
   For the 2-site cell ($n_c = 2, M_c = 8$):
   $$
   \mathrm{FLOPs}_{\mathrm{eig}}(2) = \frac{16}{3} (8^3) N_k \approx 2730.7 \cdot N_k \text{ FLOPs}.
   $$

3. **Fermi Spectral Projection / Covariance Construction:**
   Constructing the negative-energy projector $P(\mathbf{k}) = \sum_{n=1}^{M_c/2} |u_n(\mathbf{k})\rangle \langle u_n(\mathbf{k})|$ via complex GEMM:
   $$
   \mathrm{FLOPs}_{\mathrm{proj}} \approx 8 \cdot M_c \cdot M_c \cdot (M_c/2) \cdot N_k = 4 M_c^3 N_k = 256 n_c^3 N_k \text{ FLOPs}.
   $$

4. **Brillouin Zone Numerical Quadrature:**
   Summing projected covariances and Fourier factors over $N_k$ points for each bond:
   $$
   \mathrm{FLOPs}_{\mathrm{quad}} \approx 2 N_k N_b M_c^2 \text{ FLOPs}.
   $$

**Total FLOPs per $\Phi$ Evaluation (Bloch):**
$$
\mathrm{FLOPs}(\Phi_{\mathrm{Bloch}}) \approx \left( \frac{1792}{3} n_c^3 + 48 n_c^3 + 66 n_c \right) N_k \approx \frac{1936}{3} n_c^3 N_k \text{ FLOPs}.
$$
For $n_c = 2$ and $N_k = 24 \times 24 = 576$:
$$
\mathrm{FLOPs}(\Phi_{\mathrm{Bloch}}) \approx 5162.7 \times 576 \approx 2.97 \times 10^6 \text{ FLOPs}.
$$

---

## 2. Solver Overhead per Step

Each solver introduces an algebraic overhead on top of the map evaluation $\Phi$:

| Solver Paradigm | Additional Operations per Step | Overhead FLOPs | Storage Requirement (Bytes) |
| :--- | :--- | :--- | :--- |
| **Plain Fixed-Point** | Vector copy | $0$ | $8 D$ |
| **Damped Picard** | Vector addition $x + \alpha r$ | $2 D$ | $16 D$ |
| **Joint DIIS ($m$)** | Inner products $B_{ij} = \langle r_i, r_j \rangle$ + solve $(m+1)$ KKT system + linear combination | $2 m D + \frac{2}{3}m^3 + 2 m D$ | $8(2m+1)D + 8(m+1)^2$ |
| **Good Broyden (Full)** | Dense matrix-vector product $H_k F_k$ + outer product update of $H_{k+1}$ | $4 D^2 + 4 D$ | $8 D^2 + 32 D$ |
| **L-Broyden ($m$)** | Two-loop vector recursions ($m$ pairs $(s_j, y_j)$) | $4 m D$ | $8(2m+3)D$ |
| **JFNK** | $k_{\mathrm{GMRES}}$ inner evaluations of $\Phi$ + GMRES Krylov Arnoldi orthogonalization | $k_{\mathrm{GMRES}} \cdot \mathrm{FLOPs}(\Phi) + 2 k_{\mathrm{GMRES}}^2 D$ | $8(k_{\mathrm{GMRES}}+4) D$ |

### Dominance of $\Phi$ Evaluation
Because state dimension for $n_c=2$ is $D = 30$, the linear algebra overhead of DIIS ($m=5$) is:
$$
\mathrm{FLOPs}_{\mathrm{DIIS}} \approx 4(5)(30) + \frac{2}{3}(125) \approx 683 \text{ FLOPs},
$$
which is only $0.023\%$ of the $2.97 \times 10^6$ FLOPs required by $\Phi$. Thus, $\Phi$ accounts for $> 99.9\%$ of total FLOPs in the 2-site model.

### The Full Broyden Memory Wall
Notice the contrast between Full Broyden and L-Broyden:
- For $n_c = 2$ ($D = 30$): Full Broyden storage is $8(900) = 7.2\text{ KB}$ (trivial).
- For an $L=24$ real-space cluster ($N_s = 1152$, $D = 17{,}280$):
  - **Full Broyden:** Storing dense $H$ requires $8 \times (17{,}280)^2 \approx 2.39\text{ GB}$ of RAM, and each step requires $4 \times (17{,}280)^2 \approx 1.2 \times 10^9$ FLOPs just for the quasi-Newton update.
  - **L-Broyden ($m=10$):** Requires only $8(23)(17{,}280) \approx 3.18\text{ MB}$ of RAM and $4(10)(17{,}280) \approx 6.9 \times 10^5$ FLOPs per step.
This proves that Full Broyden is fundamentally unviable for large supercells, whereas L-Broyden and DIIS scale effortlessly.
