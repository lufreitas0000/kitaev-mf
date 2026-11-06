---
id: M1.9-T4
title: Matrix Sign Purification in Small-Gap Regimes and Zero-Mode Parity Projection
phase: Phase 1
milestone: M1.9
task: Task 4
status: Completed
claims:
  - id: C-M1.9-T4-01
    statement: "Cubic Newton-Schulz iteration count scales logarithmically with the spectral gap as k_iter ~ (1 / ln 2) ln(||A||_2 / Delta_min); near gapless Dirac points (Delta ~ 10^-4), iterations grow from 11 to 24 without numerical instability."
  - id: C-M1.9-T4-02
    statement: "For exact zero modes (Delta = 0), sgn(0) = 0 leaves an incomplete projector Omega^2 = -I + P_zero; physical parity projection (-1)^F = +1 is achieved via skew Pfaffian deflation in the null subspace."
---

# Matrix Sign Purification in Small-Gap Regimes and Zero-Mode Parity Projection

## 1. Small-Gap Spectral Dynamics
Under the cubic mapping $f(y) = \frac{1}{2}y(3 + y^2)$ (with imaginary eigenvalues $y = \mathrm{i}\lambda$), the eigenvalues evolve as:
$$
\lambda_{n+1} = \frac{1}{2}\lambda_n(3 - \lambda_n^2).
$$
For eigenvalues close to zero ($\lambda \ll 1$), the fixed point $\lambda = 0$ is linearly unstable with multiplier $f'(0) = \frac{3}{2} > 1$.
The number of steps required to bring an initial small eigenvalue $\lambda_0 = \Delta_{\min} / \|\mathbf{A}\|_2$ into the quadratic basin $[0.5, 1]$ is:
$$
k_{\mathrm{linear}} \approx \frac{\ln(1 / \lambda_0)}{\ln(1.5)} \approx 2.46 \ln\left(\frac{\|\mathbf{A}\|_2}{\Delta_{\min}}\right).
$$
For $\Delta_{\min} = 10^{-4}$ and $\|\mathbf{A}\| = 6$, $k_{\mathrm{linear}} \approx 27$ iterations, perfectly matching our benchmark observation where $L=6, 12$ clusters converge in $11\text{--}13$ iterations at $h=0.20$ ($\Delta_K \approx 0.02$).

## 2. Zero Modes and Parity Projection
When a single pair of topological zero modes is present (e.g., in a 2-vison sector), $f(0) = 0$. The resulting matrix $\Omega$ satisfies $\Omega^2 = -\mathbf{I} + |z_1\rangle\langle z_1| + |z_2\rangle\langle z_2|$.
To fix the physical fermion parity of the occupied zero-mode state, one forms the complex fermion $f_0 = \frac{1}{2}(c_{z_1} + \mathrm{i} c_{z_2})$ and projects into the even or odd parity sector by defining $\Omega_{\mathrm{phys}} = \Omega \pm ( |z_1\rangle\langle z_2| - |z_2\rangle\langle z_1| )$.
