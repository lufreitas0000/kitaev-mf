---
id: M1.9-T1
title: "Self-Consistent Field Iteration with Matrix Sign Purification and Warm-Start Gradient Corrections"
status: Completed
date: "2026-11-06"
questions:
  - "Why does direct inner warm-starting Y_0 = Omega_{k-1} freeze Newton-Schulz matrix sign purification?"
  - "How do outer-accelerated fixed-point updates (DIIS / damped L-Broyden) circumvent the idempotency trap in real-space self-consistent loops?"
depends_on:
  - M1.8-T2
  - M1.8-T3
contracts:
  - K-018
  - K-021
  - HC-P1-03
artifacts:
  - src/kitaev_mf/services/purification.py
  - src/kitaev_mf/services/real_space.py
  - test/unit/test_real_space_and_purification.py
---

# Self-Consistent Field Iteration with Matrix Sign Purification and Warm-Start Gradient Corrections

## 1. The Fixed-Point Pitfall of Direct Idempotency
In milestone M1.8, we established that the single-particle ground-state covariance matrix $\Omega \in \mathbb{R}^{M \times M}$ satisfying:
$$
\Omega = \mathrm{i}\operatorname{sgn}(\mathrm{i}\mathbf{A}), \qquad \Omega^T = -\Omega, \qquad \Omega^2 = -\mathbf{I}
$$
can be evaluated without calculating eigenvectors via the cubic Newton-Schulz iteration:
$$
\mathbf{Y}_{n+1} = \frac{1}{2}\mathbf{Y}_n(3\mathbf{I} + \mathbf{Y}_n^2).
$$
When initialized from the cold-start scaled Hamiltonian $\mathbf{Y}_0 = -\mathbf{A} / \|\mathbf{A}\|_2$, the iteration converges in $n \approx 11\text{--}13$ steps.

In a self-consistent field (SCF) loop where the Hamiltonian $\mathbf{A}^{(k)} = \mathbf{A}[\Omega^{(k-1)}]$ changes by a small increment $\delta\mathbf{A} = \mathbf{A}^{(k)} - \mathbf{A}^{(k-1)}$, one desires to *warm-start* the purification from $\Omega^{(k-1)}$.
However, naive substitution $\mathbf{Y}_0 = \Omega^{(k-1)}$ fails completely:
Because $(\Omega^{(k-1)})^2 = -\mathbf{I}$, we have:
$$
\mathbf{Y}_1 = \frac{1}{2}\Omega^{(k-1)}(3\mathbf{I} - \mathbf{I}) = \Omega^{(k-1)}.
$$
The iteration does not sample the new Hamiltonian $\mathbf{A}^{(k)}$ at all!

## 2. Gradient-Corrected Purification Warm Starts
To couple the new Hamiltonian $\mathbf{A}^{(k)}$ into the projector without cold-starting, we employ the first-order density-matrix perturbation (commutator flow):
$$
\mathbf{Y}_0 = \Omega^{(k-1)} - \tau \left( \mathbf{A}^{(k)} + \Omega^{(k-1)} \mathbf{A}^{(k)} \Omega^{(k-1)} \right),
$$
where $\tau = \frac{1}{2\|\mathbf{A}^{(k)}\|_2}$.
Because $\Omega^{(k-1)} \mathbf{A} \Omega^{(k-1)}$ projects out the occupied-occupied and virtual-virtual blocks, $\mathbf{Y}_0$ possesses small off-diagonal corrections of order $\|\delta\mathbf{A}\|$.
Iterating Newton-Schulz from this perturbed state restores exact antisymmetry and idempotency $\mathbf{Y}^2 \to -\mathbf{I}$ in only $3\text{--}4$ matrix multiplications instead of $22\text{--}24$.
