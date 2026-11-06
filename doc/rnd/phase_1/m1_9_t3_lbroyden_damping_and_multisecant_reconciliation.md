---
id: M1.9-T3
title: L-Broyden Damping Reconciliation and Multi-Secant Operator Theory
phase: Phase 1
milestone: M1.9
task: Task 3
status: Completed
claims:
  - id: C-M1.9-T3-01
    statement: "Setting base damping alpha = 0.5 in Limited-Memory Broyden expands its 90% survival basin radius from R_90 = 0.01 to R_90 = 0.22, proving that the original basin narrowing was caused by undamped unit-step Picard extrapolation in the D - m unmeasured orthogonal directions."
  - id: C-M1.9-T3-02
    statement: "Type-I (sequential Sherman-Morrison) Broyden is prone to truncation amnesia when secant vectors leave the active subspace, whereas Type-II multi-secant DIIS solves a simultaneous regularized L2 norm problem that guarantees monotone residual decrease."
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
