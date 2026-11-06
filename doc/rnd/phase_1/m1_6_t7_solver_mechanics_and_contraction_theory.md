---
id: RND-P1-M1.6-T7
title: "Nonlinear Solver Mechanics, Contraction Theory, and Norm Invariance on the Grassmannian"
phase: 1
milestone: M1.6
task: T7
status: completed
claims:
  - CLM-M1-SOLV-001
  - CLM-M1-SOLV-002
  - CLM-M1-SOLV-003
---

# Nonlinear Solver Mechanics, Contraction Theory, and Norm Invariance on the Grassmannian

## 1. Mathematical Structure of the Self-Consistent Map

The mean-field problem is formulated as finding the fixed point of a nonlinear map $\Phi: \mathcal{M} \to \mathcal{M}$ on the manifold:
$$
\mathcal{M} = \bigoplus_{b=1}^{N_b} \mathfrak{so}(4) \oplus \bigoplus_{s=1}^{N_s} \mathbb{R}^3,
$$
where each state vector $x = (\Omega, \boldsymbol{\lambda})$ consists of skew-symmetric bond covariance matrices $V_b \in \mathfrak{so}(4)$ and local Lagrange multipliers $\boldsymbol{\lambda}_s \in \mathbb{R}^3$. The state vector has dimension:
$$
D = 6 N_b + 3 N_s.
$$
For a 2-site primitive cell with 3 bonds, $D = 6(3) + 3(2) = 24$. For general $n_c$-site cells ($N_b = \frac{3}{2}n_c$, $N_s = n_c$), $D = 15 n_c$.

The fixed-point condition $x^* = \Phi(x^*)$ is equivalent to vanishing nonlinear residual:
$$
F(x) \equiv x - \Phi(x) = 0.
$$

---

## 2. Solver Paradigms: Algorithmic Mechanics

We systematically evaluate six distinct solver paradigms:

### 2.1 Plain Fixed-Point Iteration (FP)
$$
x_{k+1} = \Phi(x_k).
$$
- **Convergence Order:** Strictly linear with asymptotic factor $\rho = \rho(J_\Phi(x^*))$, where $\rho(J_\Phi)$ is the spectral radius of the Jacobian of $\Phi$.
- **Failure Mode:** Diverges or oscillates whenever $\rho(J_\Phi) \ge 1$.

### 2.2 Damped Picard Iteration (Relaxation)
$$
x_{k+1} = (1 - \alpha) x_k + \alpha \Phi(x_k) = x_k - \alpha F(x_k), \qquad \alpha \in (0, 1].
$$
- **Convergence Order:** Strictly linear with effective Jacobian $J_{\alpha} = I - \alpha(I - J_\Phi)$.
- **Benefit:** Stabilizes mildly oscillatory eigenvalues ($\mathrm{Re}(\mu) < 0$) and contracts the spectral radius into the unit circle when $\alpha < 2 / (1 + \rho(J_\Phi))$.

### 2.3 Direct Inversion in the Iterative Subspace (DIIS / Anderson Acceleration)
Given history $\{x_i, r_i\}_{i=k-m+1}^k$ where $r_i = \Phi(x_i) - x_i$, DIIS finds the linear combination $x_{\mathrm{opt}} = \sum_{i} c_i x_i$ that minimizes the squared residual norm:
$$
\min_{\sum c_i = 1} \left\| \sum_{i=k-m+1}^k c_i r_i \right\|_2^2,
$$
leading to the small $(m+1) \times (m+1)$ linear KKT system:
$$
\begin{pmatrix}
B & -\mathbf{1} \\
-\mathbf{1}^T & 0
\end{pmatrix}
\begin{pmatrix}
\mathbf{c} \\
\mu
\end{pmatrix}
=
\begin{pmatrix}
\mathbf{0} \\
-1
\end{pmatrix}, \qquad B_{ij} = \langle r_i, r_j \rangle.
$$
The next iterate is $x_{k+1} = \sum_i c_i \Phi(x_i)$.
- **Convergence Order:** Asymptotically superlinear. Eliminates dominant oscillatory error components.
- **Complexity:** Subspace solve is $\mathcal{O}(m^2 D + m^3)$. For $m \le 8$, overhead is entirely negligible compared to $\Phi$.

### 2.4 Good Broyden (Quasi-Newton Rank-1 Secant)
Maintains an approximation to the inverse Jacobian $H_k \approx -[\nabla F(x_k)]^{-1}$:
$$
s_k = x_{k+1} - x_k, \qquad y_k = F(x_{k+1}) - F(x_k),
$$
$$
H_{k+1} = H_k + \frac{s_k - H_k y_k}{s_k^T H_k y_k} s_k^T H_k.
$$
- **Convergence Order:** Superlinear ($q_k \to 0$).
- **Complexity:** Full Broyden explicitly stores the $D \times D$ dense matrix $H_k$, requiring $\mathcal{O}(D^2)$ storage and $\mathcal{O}(D^2)$ matrix-vector products. For large supercells ($D > 10^4$), this encounters a severe memory and cache wall.

### 2.5 Limited-Memory Broyden (`L-Broyden`)
Instead of storing $H_k \in \mathbb{R}^{D \times D}$, `L-Broyden` stores only the last $m$ secant pairs $\{s_j, y_j\}_{j=k-m}^k$. The action $H_k F(x_k)$ is evaluated via recursive vector dot products:
- **Storage:** Exactly $2 m D$ floating-point numbers ($\mathcal{O}(mD)$ space).
- **Arithmetic:** $4 m D$ FLOPs per step ($\mathcal{O}(mD)$ work).
- **Significance:** Retains the superlinear acceleration of Broyden while scaling to millions of degrees of freedom.

### 2.6 Jacobian-Free Newton–Krylov (JFNK)
Solves the true Newton equation $\nabla F(x_k) \Delta x_k = -F(x_k)$ using GMRES, evaluating the directional derivative matrix-free via forward finite differences:
$$
\nabla F(x_k) v \approx \frac{F(x_k + \epsilon v) - F(x_k)}{\epsilon}, \qquad \epsilon = \frac{\sqrt{\epsilon_{\mathrm{mach}}}}{\|v\|_2}.
$$
- **Convergence Order:** Quadratic outer convergence ($e_{k+1} \le C e_k^2$).
- **Cost Trade-Off:** Each GMRES inner step requires 1 evaluation of $\Phi$. An outer step taking 15 GMRES iterations consumes 16 calls to $\Phi$, drastically increasing total work unless preconditioned.

---

## 3. Contraction Classification and Norm Equivalence Invariance

Let $e_k = x_k - x^*$. Convergence is classified by the limiting ratio:
$$
\lim_{k \to \infty} \frac{\|e_{k+1}\|}{\|e_k\|^p} = C.
$$
- $p = 1, C \in (0, 1)$: **Linear convergence** (Fixed-point, Damped Picard).
- $p \in (1, 2)$ or $q_k = \|e_{k+1}\|/\|e_k\| \to 0$: **Superlinear convergence** (DIIS, Broyden, L-Broyden).
- $p = 2$: **Quadratic convergence** (Exact Newton, JFNK).

### Norm Equivalence Theorem
On finite-dimensional vector space $\mathbb{R}^D$, any two norms $\|\cdot\|_a$ and $\|\cdot\|_b$ are equivalent:
$$
c_1 \|v\|_b \le \|v\|_a \le c_2 \|v\|_b \quad \forall v \in \mathbb{R}^D.
$$
Taking the $k$-th root of the error ratio as $k \to \infty$:
$$
\rho_{\mathrm{eff}} \equiv \lim_{k\to\infty} \|e_k\|^{1/k}.
$$
Because $(c_1)^{1/k} \to 1$ and $(c_2)^{1/k} \to 1$ as $k \to \infty$, the asymptotic convergence rate $\rho_{\mathrm{eff}}$ is **strictly invariant under the choice of norm** ($\ell_1, \ell_2, \ell_\infty$, Frobenius). Only the pre-asymptotic constant is norm-dependent.
