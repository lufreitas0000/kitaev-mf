---
id: M1.9-T3
title: "L-Broyden Damping Reconciliation and Multi-Secant Operator Theory"
status: Completed
date: "2026-11-06"
questions:
  - "What is the mathematical mechanism underlying the narrow basin of attraction (R_90 = 0.010) in undamped unit-step Limited-Memory Broyden?"
  - "Why does base damping alpha = 0.5 expand the L-Broyden survival basin by more than 20x in real-space fixed-point iterations?"
depends_on:
  - M1.6-T7
  - M1.6-T10
contracts:
  - K-013
  - K-021
artifacts:
  - src/kitaev_mf/services/fixed_point.py
  - scripts/m1_9/benchmark_wave3_advanced_solvers.py
  - data/processed/benchmarks/m1_9/wave3_purification_and_solvers.dat
---

# L-Broyden Damping Reconciliation and Multi-Secant Operator Theory

## 1. Diagnosing the L-Broyden Basin Narrowing
In milestone M1.6, full `BroydenSolver` achieved $R_{90} \ge 0.50$ while `LimitedMemoryBroydenSolver` ($m=10$) collapsed to $R_{90} = 0.01$.
Mathematical auditing reveals that:
$$
H_m v = -\alpha v + \sum_{j=0}^{m-1} u_j (w_j^T v).
$$
For any vector $v \in \mathrm{span}\{s_0, \dots, s_{m-1}\}^\perp$, the update acts as $H_m v = -\alpha v$.
In the initial `LimitedMemoryBroydenSolver` implementation, $\alpha$ was set to $1.0$, meaning that in all $D - m$ unconstrained directions, the solver took raw, un-damped Picard steps! Since Picard with $\alpha = 1.0$ is violently unstable for the Majorana mean-field map ($\rho > 1$), errors along the orthogonal complement blew up.

## 2. Damping Fix and Verification
Setting $\alpha = 0.5$ (matching the base damping of the successful Full Broyden and Damped Picard methods) stabilizes the orthogonal complement, expanding the L-Broyden basin to $R_{90} = 0.22$ and $R_{50} = 0.40$.
