---
id: M1.9-T4
title: "Matrix Sign Purification in Small-Gap Regimes and Zero-Mode Parity Projection"
status: Completed
date: "2026-11-06"
questions:
  - "How does the spectral gap Delta_min govern the cubic Newton-Schulz iteration depth in real-space torus clusters?"
  - "How is physical parity projection (-1)^F = +1 preserved when exact zero modes arise in the single-particle skew Hamiltonian?"
depends_on:
  - M1.8-T2
  - M1.9-T1
contracts:
  - K-018
  - H10
artifacts:
  - src/kitaev_mf/services/purification.py
  - doc/labbook/chapters/ch_m1_9_real_space_and_purification.tex
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
