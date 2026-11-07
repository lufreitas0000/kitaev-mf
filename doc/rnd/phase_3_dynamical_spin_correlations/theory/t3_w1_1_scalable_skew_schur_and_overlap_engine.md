---
title: "Task T3.W1.1: Scalable Canonical Skew-Schur Solver and SVD-Accelerated Singular Thouless Overlap Engine (L <= 37)"
status: completed
author: "Lu Freitas"
date: "2026-11-07"
phase: "Phase 3 (Wave 1)"
milestone: "M3.1 / M3.2"
task: "T3.W1.1"
contracts: ["H_SCALE_01", "H_SCALE_02", "H_SCALE_03", "H_SCALE_04"]
primary_model: "Pure FM Kitaev (K > 0) with chiral 3-spin term kappa (arXiv:2209.12930)"
geometry: "Real-space honeycomb lattice clusters (torus and cylinder geometries, L up to 37, N <= 2738 sites)"
guiding_questions:
  - "How can the real Schur decomposition of skew-symmetric Hamiltonians A in R^{N_sites x N_sites} be stably standardized into canonical paired 2x2 blocks for large lattices up to L=37 (N_sites = 2738, dim A = 2738)?"
  - "Why is the naive O(N) determinant loop in singular Bogoliubov resolution a bottleneck at N=1369, and how does SVD singular-space projection identify the optimal particle-hole interchange mode in O(1) post-SVD operations?"
  - "How can the memory footprint and FLOP throughput of the Thouless pairing matrix Z_b = -(X^dagger)^{-1} Y^T be optimized for clusters hosting thousands of Majoranas?"
  - "Does the accelerated overlap engine preserve the canonical Bogoliubov unitarity relations ||XX^dagger + YY^dagger - I|| < 1e-12 and Majorana matrix element completeness sum rule sum |M|^2 = 1.0000 on large clusters?"
---

# Task T3.W1.1: Scalable Canonical Skew-Schur Solver and SVD-Accelerated Singular Thouless Overlap Engine ($L \le 37$)

## 1. Scope & Scientific Objectives
In arXiv:2209.12930 (*Scanning tunneling spectroscopy of Majorana zero modes in a Kitaev spin liquid*, Phys. Rev. B **107**, 054432), thermodynamic limit spectroscopy requires honeycomb clusters of linear dimension $L = 37$ ($N_{\mathrm{sites}} = 2 \times 37^2 = 2738$ sites, $M = 2738$ Majoranas, $N = 1369$ complex Dirac modes) with vortex separation $R \ge 10a$ to suppress finite-size Majorana zero mode (MZM) hybridization splitting below the spectroscopic resolution ($\delta\epsilon_0 < 10^{-6} |K|$).

While the prototype Phase 3 engine confirmed the mathematical formalism on small lattices ($L \in \{4, 6, 8\}$, $N_{\mathrm{sites}} \le 128$), scaling to $L = 37$ encounters two algorithmic bottlenecks:
1. **Dense Eigensystem Standardization & Sorting:**
   Standard `scipy.linalg.schur(A, output="real")` on $2738 \times 2738$ skew-symmetric matrices requires efficient, cache-friendly block-diagonal sorting of $N = 1369$ pairs $[[0, \epsilon_k], [-\epsilon_k, 0]]$ without memory replication or loss of orthogonality in $Q \in \mathrm{O}(2738)$.
2. **Singular Thouless Resolution Bottleneck:**
   When the Bogoliubov transformation between unperturbed sector $G$ and quenched sector $G' = G \cdot \sigma_j$ has $\det W = -1$ (opposite fermion parity), the vacuum overlap vanishes ($\det X = 0$). The prototype resolved this by performing an exhaustive loop over all $N$ modes:
   $$\max_{m \in \{0, \dots, N-1\}} |\det X_m|, \qquad X_m = \text{replace row } m \text{ of } X \text{ with row } m \text{ of } Y^*.$$
   For $N = 1369$, computing 1369 distinct $1369 \times 1369$ determinants requires $1369 \times \mathcal{O}(N^3) \approx 1369 \times 3.5\text{ GFLOP} \approx 4.8\text{ TFLOP}$, taking several minutes per spin-flip site.

The scientific objectives of **Task T3.W1.1** are:
1. **High-Performance Canonical Skew Solver (`H_SCALE_01`):**
   Implement an optimized real skew-Schur solver with standardized column reordering that constructs canonical orthogonal transformations $Q \in \mathrm{O}(N_{\mathrm{sites}})$ satisfying $Q^T A Q = \bigoplus_{k=1}^N \begin{pmatrix} 0 & \epsilon_k \\ -\epsilon_k & 0 \end{pmatrix}$ with $\epsilon_1 \le \dots \le \epsilon_N$, ensuring $Q^T Q = I$ to machine precision ($< 10^{-14}$) and $\det Q = \operatorname{sgn}(\operatorname{Pf}(A)) \in \{-1, +1\}$ for $N_{\mathrm{sites}} \le 2738$.
2. **SVD-Projected Particle-Hole Mode Selector (`H_SCALE_02`):**
   Formulate and implement an $\mathcal{O}(1)$ singular mode identification method based on the singular value decomposition of $X$:
   $$X = U_X \Sigma_X V_X^\dagger, \qquad \Sigma_X = \operatorname{diag}(s_1, \dots, s_N), \quad 0 \le s_1 \le \dots \le s_N \le 1.$$
   Because the zero singular value mode $v_1$ spans the kernel $\ker(X)$, the singular row replacement that maximizes $|\det X_m|$ is uniquely determined by the left-singular projection $m^* = \arg\max_m |(U_X)_{m, 1}|$, eliminating the $\mathcal{O}(N^4)$ determinant loop.
3. **Memory & Scaling Benchmarking (`H_SCALE_03`):**
   Profile runtime and memory scaling across cluster sizes $L \in \{6, 12, 18, 24, 37\}$, demonstrating stable execution within host RAM bounds without paging and achieving speedup $\ge 10\times$ for $N \ge 64$ under single-threaded BLAS execution.
4. **Unitarity & Sum-Rule Invariants on Large Clusters (`H_SCALE_04`):**
   Verify that canonical Bogoliubov identities $\|XX^\dagger + YY^\dagger - I\|_\infty < 10^{-12}$ and the local transition completeness sum rule $\sum_\lambda |\langle \lambda | c_j | \phi_0 \rangle|^2 = 1.000 \pm 0.005$ remain strictly preserved for $L \ge 12$ on genuine singular quench configurations.

---

## 2. Mathematical Derivations

### 2.1 SVD Analysis of Singular Bogoliubov Transformations
Let $X, Y \in \mathbb{C}^{N \times N}$ satisfy the canonical Bogoliubov identities:
$$X X^\dagger + Y Y^\dagger = I_N, \qquad X Y^T + Y X^T = 0 \iff X^* Y^\dagger + Y^* X^\dagger = 0.$$
From $X X^\dagger + Y Y^\dagger = I$, the singular values of $X$ satisfy $s_k \in [0, 1]$, and the singular values of $Y$ are $t_k = \sqrt{1 - s_k^2}$.

When $\det X = 0$, exactly one singular value vanishes (for a single bond-flip quench): $s_1 = 0$, while $s_k > 0$ for $k \ge 2$.
Let $X = U \Sigma V^\dagger$ be the SVD of $X$. Then:
$$\ker(X^\dagger) = \operatorname{span}\{u_1\}, \qquad \ker(X) = \operatorname{span}\{v_1\},$$
where $u_1$ is the first column of $U$ (left singular vector with singular value 0) and $v_1$ is the first column of $V$.

To determine the relation between $u_1$ and $v_1$ under $Y$, apply the second canonical identity $X^* Y^\dagger + Y^* X^\dagger = 0$ to $u_1$:
$$(X^* Y^\dagger + Y^* X^\dagger) u_1 = X^* (Y^\dagger u_1) + Y^* (X^\dagger u_1) = X^* (Y^\dagger u_1) + 0 = 0.$$
Therefore, $Y^\dagger u_1 \in \ker(X^*)$.
Since $\ker(X^*) = \operatorname{span}\{v_1^*\}$, there exists a scalar $\alpha \in \mathbb{C}$ such that:
$$Y^\dagger u_1 = \alpha v_1^*.$$
Evaluating the Euclidean norm squared:
$$\|Y^\dagger u_1\|^2 = u_1^\dagger (Y Y^\dagger) u_1 = u_1^\dagger (I - X X^\dagger) u_1 = u_1^\dagger u_1 - 0 = 1.$$
Since $\|v_1^*\| = 1$, we must have $|\alpha| = 1$, i.e., $\alpha = \mathrm{e}^{\mathrm{i}\theta}$.
Taking the complex conjugate of $Y^\dagger u_1 = \alpha v_1^*$ yields $Y^T u_1^* = \alpha^* v_1$.
By the dual Bogoliubov identity $Y^\dagger X + X^T Y^* = 0$, applying it to $v_1 \in \ker(X)$ gives $X^T (Y^* v_1) = 0$, so $Y^* v_1 \in \ker(X^T) = \operatorname{span}\{u_1^*\}$.
Therefore, $Y^* v_1 = \alpha u_1^*$, and taking components along coordinate $m$:
$$(Y^* v_1)_m = \alpha (u_1^*)_m \implies |(Y^* v_1)_m| = |\alpha| |(u_1^*)_m| \equiv |(u_1)_m| \quad \text{identically for all } m \in \{0, \dots, N-1\}.$$

Now consider replacing row $m$ of $X$ with row $m$ of $Y^*$:
$$X_m = X + e_m (y_m^* - x_m),$$
where $e_m$ is the $m$-th standard basis vector, $x_m = e_m^T X$ is the $m$-th row of $X$, and $y_m^* = e_m^T Y^*$ is the $m$-th row of $Y^*$.
This is an exact **rank-1 modification**:
$$X_m = X + e_m w_m^T, \qquad w_m = y_m^* - x_m.$$

By the cofactor expansion along row $m$:
$$\det X_m = \sum_{k=1}^N (Y^*)_{mk} C_{mk},$$
where $C_{mk}$ is the cofactor of $X_{mk}$. Since $s_1 = 0$, the adjugate matrix $\operatorname{adj}(X)$ has rank 1 and is proportional to $v_1 u_1^\dagger$:
$$\operatorname{adj}(X) = \left( \prod_{k=2}^N s_k \right) v_1 u_1^\dagger.$$
Therefore:
$$|\det X_m| = \left| \sum_{k=1}^N (Y^*)_{mk} [\operatorname{adj}(X)]_{km} \right| = \left( \prod_{k=2}^N s_k \right) |(u_1)_m| \cdot |(Y^* v_1)_m| = \left( \prod_{k=2}^N s_k \right) |(u_1)_m|^2.$$

Because $\prod_{k=2}^N s_k > 0$ is a strictly positive constant independent of $m$, the mode $m$ that maximizes $|\det X_m|$ is **identically and uniquely**:
$$m^* = \arg\max_{m \in \{0, \dots, N-1\}} |(u_1)_m|.$$

This rigorously proves that **a single SVD of $X$ determines both the optimal mode index and the entire vector of determinants $|\det X_m|$ analytically in $\mathcal{O}(N^3)$ operations**, completely eliminating the $\mathcal{O}(N^4)$ determinant loop!

### 2.2 Skew-Schur Transformation Determinant and the Pfaffian Identity
For any real skew-symmetric matrix $A \in \mathbb{R}^{2N \times 2N}$ and orthogonal matrix $Q \in \mathrm{O}(2N)$, the Jacobi–Cayley Pfaffian congruence identity states:
$$\operatorname{Pf}(Q^T A Q) = \det(Q) \operatorname{Pf}(A).$$
When canonical real Schur decomposition standardizes $A$ into block-diagonal form $\bigoplus_{k=1}^N \begin{pmatrix} 0 & \epsilon_k \\ -\epsilon_k & 0 \end{pmatrix}$ with all $\epsilon_k \ge 0$, its Pfaffian is:
$$\operatorname{Pf}\left(\bigoplus_{k=1}^N \begin{pmatrix} 0 & \epsilon_k \\ -\epsilon_k & 0 \end{pmatrix}\right) = \prod_{k=1}^N \epsilon_k \ge 0.$$
Therefore:
$$\det(Q) = \operatorname{sgn}(\operatorname{Pf}(A)) \in \{-1, +1\}.$$
Whenever $\operatorname{Pf}(A) < 0$ (such as in odd-parity topological vortex sectors on periodic tori), $\det(Q) = -1$ strictly. Thus $Q \in \mathrm{O}(2N)$, and the determinant $\det(Q) \in \{-1, +1\}$ is an intrinsic topological diagnostic rather than an arbitrary gauge choice.

---

## 3. Falsification Contracts (`H_SCALE_01` – `H_SCALE_04`)

| Contract ID | Invariant / Property | Numerical Criterion | Falsification Condition |
| :--- | :--- | :--- | :--- |
| **`H_SCALE_01`** | Skew-Schur Canonical Orthogonality | $\|Q^T Q - I\|_\infty < 10^{-14}$, $\|Q^T A Q - \bigoplus \epsilon_k J\|_\infty < 10^{-13}$, $\det Q \in \{-1, 1\}$ | Orthogonality residual $\ge 10^{-12}$ |
| **`H_SCALE_02`** | SVD Particle-Hole Mode Exactness | $m^*_{\mathrm{SVD}} == m^*_{\mathrm{det\_loop}}$ and $|\det X_{m^*}^{\mathrm{SVD}}| \ge 0.999 \max_m |\det X_m|$ | Discrepancy in optimal mode index or $|\det X'| < 10^{-6}$ |
| **`H_SCALE_03`** | Computational Speedup | $\tau_{\mathrm{SVD}} / \tau_{\mathrm{det\_loop}} \le \frac{1}{10}$ for $N \ge 64$ under single-threaded BLAS | Speedup $< 10\times$ for $N \ge 64$ |
| **`H_SCALE_04`** | Large-Cluster Unitarity & Sum Rule | $\|XX^\dagger + YY^\dagger - I\|_\infty < 10^{-12}$, $\sum_\lambda |\langle \lambda | c_j | \phi_0 \rangle|^2 = 1.000 \pm 0.005$ on $L \ge 12$ under singular quench | Unitarity error $\ge 10^{-11}$ or sum rule defect $> 0.01$ |
