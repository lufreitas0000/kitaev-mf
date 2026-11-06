---
id: RND-P1-M1.6-T6
title: "Formal Statistical Scorecards, Basin Metrics, and Hypothesis Contracts for Nonlinear Solvers"
phase: 1
milestone: M1.6
task: T6
status: completed
claims:
  - CLM-M1-STAT-001
  - CLM-M1-STAT-002
  - CLM-M1-STAT-003
  - CLM-M1-STAT-004
---

# Formal Statistical Scorecards, Basin Metrics, and Hypothesis Contracts for Nonlinear Solvers

## 1. Context and Problem Statement

Previous evaluations of nonlinear solvers in Milestone M1.6 relied on single deterministic perturbation scans (e.g., perturbing with a fixed seed across $\delta \in \{0.001, 0.01, 0.05, 0.20\}$). As mandated by Rule 02 (Section 5) and scientific rigor, a single perturbation trajectory cannot certify the basin of attraction radius, stability, or convergence reliability of a nonlinear Grassmannian solver.

In this note, we establish the formal statistical framework, quantitative evaluation metrics, and pre-registered **Hypothesis Contracts (`H-01` to `H-12`)** to benchmark nonlinear solvers across statistical ensembles of initial states.

---

## 2. Statistical Metrics and Definitions

Let $x_0(\delta, \xi)$ denote an initial state perturbed from a target ground state $x^* = (\Omega^*, \boldsymbol{\lambda}^*)$ with amplitude $\delta \ge 0$ along a random direction $\xi$, where $\xi$ satisfies all physical symmetries ($C_3$ rotation, sublattice equivalence $A \leftrightarrow B$, and $\mathfrak{so}(4)$ tangent space constraints).

For an ensemble of $N_{\mathrm{trials}}$ independent realizations $\{\xi_j\}_{j=1}^{N_{\mathrm{trials}}}$ at a fixed perturbation amplitude $\delta$:

### 2.1 Convergence Probability $P_{\mathrm{conv}}(\delta)$ and Wilson Confidence Interval
Let $I_{\mathrm{conv}}^{(j)} \in \{0, 1\}$ indicate whether solver $S$ converged within $K_{\max}$ iterations ($r_K < \tau = 10^{-10}$).
The empirical convergence probability is:
$$
\hat{P}_{\mathrm{conv}}(\delta) = \frac{1}{N_{\mathrm{trials}}} \sum_{j=1}^{N_{\mathrm{trials}}} I_{\mathrm{conv}}^{(j)}.
$$
Because $\hat{P}$ can be close to 1 or 0, Gaussian Wald intervals fail. We compute the exact **Wilson score interval** at significance level $\alpha = 0.05$ ($z = 1.96$):
$$
w_\pm = \frac{\hat{P} + \frac{z^2}{2 N} \pm z \sqrt{\frac{\hat{P}(1-\hat{P})}{N} + \frac{z^2}{4 N^2}}}{1 + \frac{z^2}{N}}.
$$

### 2.2 Ground-State Fidelity $P_{\mathrm{same}}(\delta)$
Even when a solver converges, it may fall into an excited metastable saddle, nematic state, or competing flux sector.
Let $I_{\mathrm{same}}^{(j)} = 1$ if the solver converged to the target branch ($\|\Omega_K - \Omega^*\|_F < 10^{-6}$ and $|E_K - E^*| < 10^{-8}$), and 0 otherwise. The conditional target fidelity is:
$$
\hat{P}_{\mathrm{same}}(\delta) = \frac{\sum_{j} I_{\mathrm{same}}^{(j)}}{\sum_{j} I_{\mathrm{conv}}^{(j)}}.
$$

### 2.3 Basin Radii ($r_{50}, r_{90}$)
We fit the ensemble survival curve $\hat{P}_{\mathrm{conv}}(\delta)$ to a sigmoid logistic model:
$$
P_{\mathrm{fit}}(\delta) = \frac{1}{1 + \exp\left( \frac{\delta - r_{50}}{\sigma_{\mathrm{basin}}} \right)},
$$
from which we extract:
- $r_{50}$: The perturbation amplitude where $50\%$ of symmetric random seeds converge.
- $r_{90}$: The certified high-reliability radius where $90\%$ of random seeds converge ($P = 0.90$).

### 2.4 Asymptotic Contraction Factor $\rho_{\mathrm{eff}}$
For a residual sequence $\{r_k\}_{k=1}^K$, the instantaneous step contraction ratio is:
$$
q_k \equiv \frac{\|r_{k+1}\|_2}{\|r_k\|_2}.
$$
In the asymptotic linear tail ($\|r_k\| \in [10^{-11}, 10^{-5}]$), the effective contraction factor is the geometric mean:
$$
\rho_{\mathrm{eff}} = \exp\left( \frac{1}{|K_{\mathrm{tail}}|} \sum_{k \in K_{\mathrm{tail}}} \ln q_k \right).
$$

---

## 3. Master Pre-Registered Hypothesis Contracts (`H-01` to `H-12`)

The statistical benchmark script `scripts/m1_6/statistical_solver_benchmark.py` will test the following 12 pre-registered contracts:

1. **`H-01` (Asymptotic Contraction Regimes):** Fixed-Point and Damped Picard contract strictly linearly ($R^2 > 0.99$ on $\ln \|r_k\|$), whereas quasi-Newton and accelerated methods (DIIS, Broyden, L-Broyden, JFNK) achieve superlinear bursts with $\rho_{\mathrm{eff}} < 0.65\,\rho_{\mathrm{Picard}}$.
2. **`H-02` (Norm Equivalence Invariance):** The asymptotic contraction factor $\rho_{\mathrm{eff}}$ is invariant within $3\%$ across $\ell_1$, $\ell_2$, $\ell_\infty$, and normalized Frobenius norms.
3. **`H-03` (Critical Slowing Down vs. Saddle Trapping):** At $h \to h_{c1}^+$, slowing down is accompanied by positive semi-definite directional curvature $v^T \nabla^2 E v \to 0^+$, whereas intermediate-field stagnation ($h = 0.20$, $\delta \ge 0.20$) corresponds to negative directional curvature $v^T \nabla^2 E v < -10^{-4}$ (saddle trapping).
4. **`H-04` (Speed–Radius Trade-off):** Quasi-Newton methods (DIIS, Broyden) have $3\times$ faster asymptotic speed ($K_{50}$) than Damped Picard at $\delta \le 0.02$, but smaller high-reliability radius $r_{90}$ unless damped.
5. **`H-05` (Basin Fidelity $P_{\mathrm{same}}$):** At $h = 0.40$ (deep polarized), $P_{\mathrm{same}}(\delta) = 1.0$ for all converged runs. At $h = 0.20$, $P_{\mathrm{same}}(\delta)$ drops for $\delta \ge 0.10$ due to basin jumping into competing nematic/flux states.
6. **`H-06` (Iteration Survival ECDF Dominance):** Joint DIIS ($m=5$) and L-Broyden ($m=10$) stochastically dominate Plain Fixed-Point and JFNK in total $\Phi$-evaluations (Kolmogorov–Smirnov $p < 0.01$).
7. **`H-07` (Curvature-Checked Saddle Escape):** Augmenting quasi-Newton solvers with a directional curvature diagnostic eliminates false convergence to nematic saddles at $h=0.20, \delta=0.20$.
8. **`H-08` (Analytical vs. Profiled FLOP Agreement):** Operation profiler matches exact BLAS/LAPACK formulas within $<1\%$.
9. **`H-09` (Host CPU Roofline Calibration):** Empirical sustained GFLOP/s on the Intel i5-13500 predicts wall-clock times within $\pm 25\%$ for $M \ge 256$.
10. **`H-10` ($\Phi$-Evaluation Cost Dominance):** The spectral projection map $\Phi$ accounts for $>95\%$ of total FLOPs and $>90\%$ of wall time for all solvers.
11. **`H-11` (Full Broyden $\mathcal{O}(D^2)$ vs. L-Broyden $\mathcal{O}(mD)$):** Full Broyden scales quadratically in state dimension $D = 15 n_c$, creating a memory wall for large supercells, while L-Broyden scales strictly as $\mathcal{O}(mD)$.
12. **`H-12` (JFNK Krylov Work Ratio):** JFNK requires $>2.5\times$ more total $\Phi$-evaluations than DIIS due to inner GMRES matrix-free finite-difference evaluations.
