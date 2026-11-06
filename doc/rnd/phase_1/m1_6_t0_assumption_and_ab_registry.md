---
id: RND-P1-M1.6-T0
title: "Assumption, Hypothesis, and A/B Testing Registry: Solvers, Seeds, and Weak-Field Physics"
phase: 1
milestone: M1.6
task: T0
status: in_progress
claims:
  - CLM-M1-HYP-001
  - CLM-M1-HYP-002
  - CLM-M1-HYP-003
  - CLM-M1-HYP-004
  - CLM-M1-HYP-005
  - CLM-M1-ALG-001
  - CLM-M1-ALG-002
  - CLM-M1-ALG-003
  - CLM-M1-ALG-004
  - CLM-M1-ALG-005
---

# Assumption, Hypothesis, and A/B Testing Registry: Solvers, Seeds, and Weak-Field Physics

## 1. Context and Methodological Protocol

In Milestone M1.6, we establish an explicit hypothesis-driven verification framework. We do not adopt algorithmic rules or physical conclusions based on intuition alone. Every algorithmic design choice (damping, history depth, secant updates, Krylov tolerance, skew symmetrization) and every physical property (gauge multiplier scaling, basin separation, $C_3$ symmetry preservation or breaking) is formulated as a falsifiable hypothesis and subjected to quantitative testing against empirical simulation data.

---

## 2. Hypothesis Registry for Weak-Field Exploration ($h \in \{0.01, 0.05\}$)

### 2.1 Physics Hypotheses (PH)

| ID | Claim | Hypothesis Statement | Test & Metric | Status |
| :--- | :--- | :--- | :--- | :---: |
| **PH-01** | `CLM-M1-HYP-001` | **Weak-Field Gauge Multiplier Scaling:** Converged Lagrange multipliers $\boldsymbol{\lambda}^*(h)$ scale quadratically $\sim \mathcal{O}(h^2)$ rather than linearly $\mathcal{O}(h)$, because the linear Zeeman term couples to spin generators $N^a$, which commute with gauge generators $G^\gamma$ at lowest order. | Compare ratio $\|\boldsymbol{\lambda}^*(0.05)\| / \|\boldsymbol{\lambda}^*(0.01)\|$ against expected $(0.05/0.01)^2 = 25$. | `untested` |
| **PH-02** | `CLM-M1-HYP-002` | **Weak-Field Basin Separation:** At $h \in \{0.01, 0.05\}$, the zero-flux QSL fixed point is the global ground state with $W_p \approx 1 - \mathcal{O}(h^2)$, while the polarized seed converges to a distinct metastable local minimum with $W_p \ll 1$ and higher energy density ($E_{\mathrm{POL}} > E_{\mathrm{QSL}}$). | Compare converged $E_{\mathrm{site}}$ and $W_p$ between QSL and Polarized initial seeds. | `untested` |
| **PH-03** | `CLM-M1-HYP-003` | **Isotropic Sublattice Magnetization:** Under isotropic field $\mathbf{h} = (h, h, h)/\sqrt{3}$, the induced spin magnetization satisfies $m_x = m_y = m_z$ to machine precision ($< 10^{-10}$) and preserves sublattice symmetry $m_A = m_B$. | Evaluate $\|m_x - m_y\|_\infty$ and $\|m_A - m_B\|_\infty$. | `untested` |
| **PH-04** | `CLM-M1-HYP-004` | **Non-Vanishing Multipliers at Infinitesimal Field:** While $\boldsymbol{\lambda}^* = 0$ is exact at $h = 0$, any finite field $h > 0$ breaks time-reversal and dynamically requires $\boldsymbol{\lambda}^* \neq 0$ to satisfy $\operatorname{tr}(G^\gamma V_s) = 0$. | Verify $\|\boldsymbol{\lambda}^*\| > 10^{-6}$ for $h=0.01$. | `untested` |
| **PH-05** | `CLM-M1-HYP-005` | **$C_3$ Symmetry vs. Spontaneous Nematicity:** For field along $[111]$, the single-particle Hamiltonian commutes with $C_3$. In the weak-field QSL ground state, bond correlators preserve $C_3$ ($U_x = U_y = U_z$), but candidate metastable or perturbed seeds may exhibit spontaneous nematic symmetry breaking ($U_x \neq U_y$). | Measure bond anisotropy $\max_{a \neq b} \|U_a - U_b\|_F$. | `untested` |

### 2.2 Algorithm Hypotheses (AH)

| ID | Claim | Hypothesis Statement | Test & Metric | Status |
| :--- | :--- | :--- | :--- | :---: |
| **AH-01** | `CLM-M1-ALG-001` | **Baseline Fixed-Point (FP) Contraction Rate:** Plain FP iteration ($x_{k+1} = \Phi(x_k)$, $\alpha=1$) converges for $h = 0.01$ with asymptotic contraction factor $\mu < 0.9$, but shows noticeable slowdown or oscillations at $h = 0.05$ due to cross-coupling between $V$ and $\boldsymbol{\lambda}$. | Measure $\mu_{\mathrm{tail}}$ and iteration count $n_\Phi$ across $h \in \{0.01, 0.05\}$. | `untested` |
| **AH-02** | `CLM-M1-ALG-002` | **Joint Acceleration Efficiency in $n_\Phi$:** Joint DIIS ($m=6$) and Good Broyden converge in significantly fewer total map evaluations ($n_\Phi$) than plain FP and damped Picard on both $h = 0.01$ and $h = 0.05$. | Compare $n_\Phi$ to convergence tolerance $10^{-10}$. | `untested` |
| **AH-03** | `CLM-M1-ALG-003` | **JFNK Work Trade-Off:** JFNK requires fewer outer iterations than FP/DIIS, but its total map evaluations $n_\Phi = n_{\mathrm{outer}} \times (1 + n_{\mathrm{GMRES}})$ make it computationally more expensive than DIIS in the weak-field regime. | Benchmark total $n_\Phi$ and wall time of JFNK vs DIIS. | `untested` |
| **AH-04** | `CLM-M1-ALG-004` | **Skew Symmetrization Invariance & Leakage:** For skew seeds, Arm B (output symmetrization) and Arm C (output + post-accelerator retraction) produce numerically identical trajectories ($|x_B - x_C| < 10^{-14}$). For non-skew seeds, Arm A (no symmetrization) fails or leaks, while Arms B and C extinguish symmetric errors in 1 step. Arm D (intrinsic $\mathfrak{so}(4)$) eliminates skewness by construction. | Run 4-arm A/B benchmark on $h \in \{0.01, 0.05\}$. | `untested` |
| **AH-05** | `CLM-M1-ALG-005` | **Seed Basin Robustness:** Perturbations $\delta \le 0.05$ added to the QSL seed always converge back to the QSL basin at $h \in \{0.01, 0.05\}$, confirming a stable finite basin of attraction. | Sweep $\delta \in \{10^{-3}, 10^{-2}, 10^{-1}\}$ across 5 random perturbation seeds. | `untested` |

---

## 3. Algorithmic Assumption & Design Choice Registry

| Assumption ID | Statement | Origin | Rationale | Test Plan |
| :--- | :--- | :--- | :--- | :--- |
| **A-01** | Linear Picard damping $\alpha \in (0, 1]$ stabilizes oscillatory modes at the cost of slower asymptotic convergence for monotone modes. | `agent-proposed` | Classical fixed-point relaxation theory. | Test $\alpha \in \{0.2, 0.5, 0.8, 1.0\}$. |
| **A-02** | Joint DIIS should extrapolate on the damped target rather than the raw output. | `inherited-M1.3` | Smoothes high-frequency residual spikes before subspace projection. | A/B test DIIS on raw vs damped target. |
| **A-03** | DIIS subspace depth $m = 6$ is optimal for small unit cells. | `inherited-M1.3` | Balances memory and Gram matrix conditioning. | Test $m \in \{3, 6, 10\}$. |
| **A-04** | Tikhonov regularization $\epsilon_{\mathrm{reg}} = 10^{-12}$ prevents singular KKT systems in DIIS. | `inherited-M1.3` | Standard numerical linear algebra safeguard. | Test $\epsilon \in \{0, 10^{-12}, 10^{-8}\}$. |
| **A-05** | Dual multiplier step $\eta_\lambda = 0.2$ provides stable dual gradient ascent. | `inherited-M1.3` | Matches scale of gauge residual to single-particle band gap. | Sweep $\eta_\lambda \in \{0.05, 0.2, 0.5\}$. |
| **A-06** | Good Broyden with inverse Jacobian update $H_k$ does not require explicit restarts for $h \le 0.05$. | `agent-proposed` | Weak-field perturbations remain within local quadratic basin. | Monitor Broyden residual monotonicity. |
| **A-07** | Finite-difference step $\epsilon = 10^{-7}(1 + \|x\|)$ provides an accurate directional derivative for matrix-free GMRES in JFNK. | `agent-proposed` | Balances truncation and floating-point roundoff error. | Verify JFNK quadratic convergence ratio. |
| **A-08** | Monkhorst-Pack $24 \times 24$ grid containing Dirac points does not cause discontinuous jumps during weak-field updates. | `agent-proposed` | Finite field $h > 0$ opens a topological gap at Dirac points. | Compare $24 \times 24$ vs $25 \times 25$ grid. |
