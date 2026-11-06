---
id: RND-P1-M1.6-T9
title: "Critical Slowing Down vs. Metastable Nematic Saddles: Theory and Curvature Diagnostics"
phase: 1
milestone: M1.6
task: T9
status: completed
claims:
  - CLM-M1-SADD-001
  - CLM-M1-SADD-002
---

# Critical Slowing Down vs. Metastable Nematic Saddles: Theory and Curvature Diagnostics

## 1. Physical and Mathematical Distinctions

In nonlinear fixed-point iterations, solver stalling (residual failing to decrease below tolerance) can arise from two mathematically and physically distinct mechanisms:

### 1.1 Physical Spinodal Breakdown (Critical Slowing Down)
Near a second-order critical point or at the edge of metastability (spinodal threshold $h_{c2}^+ = 0.400$), the physical ground state energy surface $E_{\mathrm{MF}}(x)$ flattens along the soft mode direction:
$$
\lambda_{\min}(\nabla^2 E_{\mathrm{MF}}(x^*)) \to 0^+.
$$
By the stationarity condition, the fixed-point map $\Phi(x)$ is related to the gradient of the variational energy. Linearizing near $x^*$:
$$
J_\Phi(x^*) \approx I - M^{-1} \nabla^2 E_{\mathrm{MF}}(x^*).
$$
As the smallest eigenvalue of the Hessian approaches zero, the spectral radius of the fixed-point Jacobian approaches unity from below:
$$
\rho(J_\Phi) = 1 - \frac{\lambda_{\min}}{\|M\|} \to 1^-.
$$
The contraction ratio $q_k = \|r_{k+1}\|/\|r_k\| \to 1$, producing **critical slowing down**. Crucially:
- The energy curvature along all physical directions remains **strictly positive semi-definite** ($v^T \nabla^2 E v \ge 0$).
- Magnetization and bond correlators remain strictly $C_3$-symmetric.

### 1.2 Trapping Near an Unstable Nematic Saddle
At intermediate fields ($h \in [0.15, 0.30]$), large random initial perturbations can drag the state vector into the basin of attraction of an **unstable nematic saddle point** $x_{\mathrm{saddle}}$.
At this saddle point:
- The residual appears small ($\|F(x_{\mathrm{saddle}})\| \ll 1$), giving a false appearance of convergence.
- The state breaks discrete $C_3$ rotation symmetry ($m_x \neq m_y \neq m_z$, $U_x \neq U_y$).
- The variational Hessian possesses **negative eigenvalues**:
$$
\exists v \in T_{x}\mathcal{M}: \quad v^T \nabla^2 E_{\mathrm{MF}} v < 0.
$$
An unconstrained quasi-Newton solver (such as plain Broyden or undamped DIIS) that does not test second-order curvature can remain trapped on the flat plateau of this saddle for dozens of iterations or report false convergence.

---

## 2. A/B Curvature Diagnostic Protocol

To distinguish true spinodal slowing down from saddle trapping, we implement a directional finite-difference Hessian curvature test:

### Directional Curvature Estimator
Given recent step direction $v = \Delta x_k / \|\Delta x_k\|_2$:
$$
\kappa_v \equiv v^T \nabla^2 E_{\mathrm{MF}}(x_k) v \approx \frac{E(x_k + \epsilon v) - 2 E(x_k) + E(x_k - \epsilon v)}{\epsilon^2}, \qquad \epsilon = 10^{-4}.
$$

### Classification Rule:
1. **True Ground State / Spinodal Slowing Down:** $\kappa_v > 0$ and $N_{\mathrm{nem}} < 10^{-5}$. Continue solver with Anderson acceleration or Picard damping.
2. **Nematic Saddle Trapping:** $\kappa_v < -10^{-4}$ or ($N_{\mathrm{nem}} > 10^{-2}$ with $\|r_k\|_\infty > 10^{-6}$). Apply an escape perturbation along the negative curvature eigenvector $v_{\mathrm{neg}}$ or re-project onto the $C_3$-symmetric manifold.
