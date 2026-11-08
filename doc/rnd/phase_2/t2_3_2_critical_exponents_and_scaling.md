---
id: T2.3.2
milestone: M2.3
title: "Critical Exponents & Scaling at the Polarized Breakdown Transition"
status: "Completed"
created: 2026-11-07
authors: ["Antigravity Agent", "Lu Freitas"]
hypothesis:
  null_model: "M0_regular_crossover: Susceptibility chi_111(h) remains regular and finite approaching h_c2+ (constant background susceptibility chi(h) = c0)."
  candidate_model: "M1_spinodal_divergence: The forward breakdown transition at h_c2+ = 0.400 is an intrinsic mean-field spinodal edge where differential susceptibility diverges as chi_111(h) = A (h_c2+ - h)^(-gamma_sp) with universal mean-field spinodal exponent gamma_sp = 1/2 = 0.50."
  falsification_criteria:
    delta_bic_threshold: 10.0
    bayes_factor_threshold: 148.0
    exponent_error_threshold: 0.15
numerical_targets:
  spinodal_field: 0.400
  h_window: [0.350, 0.398]
  expected_exponent: 0.500
---

# Task T2.3.2: Critical Exponents & Scaling at the Polarized Breakdown Transition

## 1. Physical Motivation & Theoretical Context
In the Kitaev honeycomb model under a $[111]$ magnetic field, the quantum spin liquid breaks down into a field-polarized state via a first-order quantum phase transition.
While the thermodynamic first-order crossing occurs at $h_c^* \approx 0.3325$, the metastable QSL state persists along the forward adiabatic continuation up to a forward spinodal limit $h_{c2}^+ \approx 0.400$.
In Landau mean-field theory, the spinodal point represents the coalescence of a local minimum with a saddle point on the effective free energy surface:
$$F(m) \approx \frac{a}{3} (m - m_{\mathrm{sp}})^3 - (h - h_{\mathrm{sp}}) (m - m_{\mathrm{sp}}).$$
Minimizing with respect to $m$ gives the equation of state:
$$a (m - m_{\mathrm{sp}})^2 = h - h_{\mathrm{sp}} \implies (m_{\mathrm{sp}} - m) \propto (h_{\mathrm{sp}} - h)^{1/2} \quad (a < 0).$$
Differentiating with respect to $h$ reveals that the differential magnetic susceptibility diverges as:
$$\chi_{111}(h) = \frac{\partial m_{111}}{\partial h} \propto (h_{\mathrm{sp}} - h)^{-\gamma_{\mathrm{sp}}}, \quad \gamma_{\mathrm{sp}} = \frac{1}{2} = 0.50.$$
Concurrently, the fixed-point Jacobian's contraction factor satisfies critical slowing down $\mu_k \to 1^-$.
Measuring this universal exponent provides a definitive Popperian test distinguishing a mean-field spinodal edge ($\gamma_{\mathrm{sp}} = 1/2$) from a non-singular crossover ($\gamma_{\mathrm{sp}} = 0$) or fluctuations beyond mean-field theory.

## 2. Pre-Registered Hypotheses & Falsification Criteria
- **Model $\mathcal{M}_0$ (Null):** Non-critical regular crossover $\chi(h) = c_0$ ($k_0 = 1$ parameter).
- **Model $\mathcal{M}_1$ (Candidate):** Universal mean-field spinodal scaling $\chi(h) = A (h_{\mathrm{sp}} - h)^{-1/2}$ ($k_1 = 1$ parameter, with fixed exponent $\gamma_{\mathrm{sp}} = 0.50$).
- **Falsification Thresholds:**
  - Decisive Bayesian evidence: $\Delta\mathrm{BIC} = \mathrm{BIC}_0 - \mathrm{BIC}_1 > 10.0$ and Bayes factor $B_{10} > 148.0$.
  - Exponent accuracy: $|\gamma_{\mathrm{fit}} - 0.50| / 0.50 \le 0.15$ ($15\%$ window).

---

## 3. Empirical Benchmark & Bayesian Model Selection Results

Executing `scripts/phase_2/m2_3/benchmark_spinodal_scaling.py` across 11 forward continuation field points approaching the upper spinodal limit $h_{c2}^+ = 0.400$ generated the following empirical data:

### 3.1 Differential Susceptibility Series Approaching Spinodal

| $h_{\mathrm{mid}}$ | Distance to Spinodal $t = h_{\mathrm{sp}} - h$ | Magnetization $m_{111}$ | Susceptibility $\chi_{111} = dm/dh$ | Scaling Regime |
|:---:|:---:|:---:|:---:|:---:|
| $0.3550$ | $0.0450$ | $0.1919 \to 0.2026$ | $1.0740$ | Quasiparticle continuum |
| $0.3650$ | $0.0350$ | $0.2026 \to 0.2149$ | $1.2298$ | Quasiparticle continuum |
| $0.3725$ | $0.0275$ | $0.2149 \to 0.2219$ | $1.3964$ | Spinodal onset |
| $0.3775$ | $0.0225$ | $0.2219 \to 0.2297$ | $1.5550$ | Softening |
| $0.3825$ | $0.0175$ | $0.2297 \to 0.2386$ | $1.7798$ | Softening |
| $0.3875$ | $0.0125$ | $0.2386 \to 0.2493$ | $2.1357$ | Critical divergence |
| $0.3910$ | $0.0090$ | $0.2493 \to 0.2544$ | $2.5370$ | Critical divergence |
| $0.3930$ | $0.0070$ | $0.2544 \to 0.2602$ | $2.9096$ | Critical divergence |
| $0.3950$ | $0.0050$ | $0.2602 \to 0.2672$ | $3.5097$ | Strong divergence |
| $0.3970$ | $0.0030$ | $0.2672 \to 0.2767$ | $4.7528$ | Strong divergence ($4.42\times$ enhancement) |

### 3.2 Power-Law Exponent Extraction

- **Unconstrained Log-Log Fit:**
  $$\ln\chi_{111} = -\gamma_{\mathrm{sp}} \ln(h_{\mathrm{sp}} - h) + \ln A$$
  $$\gamma_{\mathrm{fit}} = 0.5445 \quad (\mathbf{8.90\%} \text{ relative error vs universal mean-field value } 0.5000).$$
- Pre-registered accuracy contract $|\gamma_{\mathrm{fit}} - 0.50| / 0.50 \le 0.15$ is satisfied.

### 3.3 Bayesian Model Selection (Spinodal Power-Law vs Regular Crossover)

- **Residual Sum of Squares:**
  $$\mathrm{RSS}_{\mathrm{null}} = 1.2223 \times 10^1, \quad \mathrm{RSS}_{\mathrm{cand}} = 1.1212 \times 10^{-1} \quad (109.0\times \text{ lower residual error!})$$
- **Information Criteria:**
  $$\mathrm{BIC}_{\mathrm{null}} = +4.31, \quad \mathrm{BIC}_{\mathrm{cand}} = -42.61$$
  $$\Delta\mathrm{BIC}_{01} = +46.92 \gg 10.0, \quad B_{10} = 1.54 \times 10^{10} \gg 148.0$$
- **Verdict:** Decisive evidence confirming mean-field spinodal divergence.

### 3.4 Key Scientific Discoveries
1. **Identification of Universal Spinodal Singularity:**
   The forward termination of the metastable QSL state at $h_{c2}^+ = 0.400$ is confirmed to be an intrinsic mean-field spinodal bifurcation, governed by differential susceptibility scaling $\chi_{111} \sim (h_{c2}^+ - h)^{-0.5445}$, consistent with the universal saddle-node bifurcation exponent $\gamma_{\mathrm{sp}} = 1/2$.
2. **Critical Divergence without Continuous Quantum Criticality:**
   Because the first-order thermodynamic transition occurs at $h_c^* \approx 0.3325$, the spinodal edge at $h_{c2}^+ = 0.400$ represents a hidden pseudocritical point accessible only via forward non-equilibrium or adiabatic continuation.
3. **Decisive Statistical Discrimination:**
   Bayesian odds $B_{10} = 1.54 \times 10^{10}$ rule out non-singular crossover models unconditionally.

