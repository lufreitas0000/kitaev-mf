---
id: M1.9-T2
title: "Solver Iteration Scaling under Broken Translation Symmetry and Vacancy Disorder"
status: Completed
date: "2026-11-06"
questions:
  - "How does breaking translation symmetry via dilute site vacancies and bond disorder alter the real-space variational dimension and solver conditioning?"
  - "Can Joint DIIS and damped L-Broyden maintain rapid convergence when localized vacancy zero modes spike the condition number to kappa(A) ~ 3.8e4?"
depends_on:
  - M1.8-T1
  - M1.9-T1
contracts:
  - K-021
  - K-022
artifacts:
  - src/kitaev_mf/services/real_space.py
  - scripts/m1_9/benchmark_wave3_advanced_solvers.py
  - data/processed/benchmarks/m1_9/wave3_disorder_scaling.dat
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
