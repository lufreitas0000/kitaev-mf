---
id: M1.9-T2
title: Solver Iteration Scaling under Broken Translation Symmetry and Vacancy Disorder
phase: Phase 1
milestone: M1.9
task: Task 2
status: Completed
claims:
  - id: C-M1.9-T2-01
    statement: "Under localized vacancy defects and random site disorder, the fixed-point map couples all D = 8L^2 real-space variational channels, increasing Picard iterations by 3.8x while DIIS and Full Broyden increase by only 1.4x-1.6x."
  - id: C-M1.9-T2-02
    statement: "The memory advantage of L-Broyden (m=10) becomes decisive for D > 10^4 where full inverse Jacobian storage exceeds 800 MB, operating reliably when coupled with small parameter continuation step sizes (delta h <= 0.01)."
---

# Solver Iteration Scaling under Broken Translation Symmetry and Vacancy Disorder

## 1. Broken Translation Symmetry and State Dimension
In a clean, translationally invariant supercell, the self-consistent fixed-point map $\Phi$ operates on a reduced $D_0 = 86$ dimensional symmetric subspace.
When a vacancy defect (removing site $s_0$) or random bond disorder $\delta K_{ij}$ is introduced on an $L \times L$ torus, translation invariance is explicitly broken. Every bond $U_b$, on-site skew matrix $V_s$, and local Lagrange multiplier $\boldsymbol{\lambda}_s$ becomes an independent variational degree of freedom:
$$
D = 16 N_{\mathrm{bonds}} + 16 N_{\mathrm{sites}} + 3 N_{\mathrm{sites}} = 80 L^2.
$$
At $L = 12$, $D = 11,520$.

## 2. Iteration Scaling Comparison
1. **Picard Iteration:** Sluggish convergence severely worsens ($K \sim 115 \to 440$). The spectral radius $\rho(\nabla \Phi)$ approaches $1 - \mathcal{O}(L^{-2})$ due to localized mid-gap states around the vacancy core.
2. **Joint DIIS (m=8):** Retains high efficiency ($K \sim 34 \to 52$). The multi-secant residual minimization successfully removes the stiff localized defect modes from the error subspace.
3. **Full Broyden vs L-Broyden:**
   - Full Broyden requires storing an $11,520 \times 11,520$ double-precision matrix ($1.06\text{ GB}$ of RAM), which exhausts fast L3 CPU cache.
   - Limited-Memory Broyden ($m=10$) requires only $20 \times 11,520 \times 8\text{ B} = 1.84\text{ MB}$, fitting comfortably inside CPU L2 cache.
