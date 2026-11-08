---
id: T2.3.3
milestone: M2.3
title: "Global Phase Diagram Mapping & First-Order Boundary Anisotropy"
status: "Completed"
created: 2026-11-07
authors: ["Antigravity Agent", "Lu Freitas"]
hypothesis:
  null_model: "M0_isotropic_boundary: The first-order thermodynamic critical field h_c*(theta) is isotropic (independent of field angle theta in the (1, 1, 0) plane)."
  candidate_model: "M1_anisotropic_boundary: Point-group exchange symmetry induces an angular modulation h_c*(theta) = h_base + Delta_h sin^4(theta) cos^2(theta) / (4/27), with maximum stability along the high-symmetry [111] axis (theta_max = arccos(1/sqrt(3)) ~ 54.74 deg)."
  falsification_criteria:
    delta_bic_threshold: 10.0
    bayes_factor_threshold: 148.0
    symmetry_maximum_tolerance_deg: 3.0
numerical_targets:
  angular_resolution: 13
  h_window: [0.330, 0.340]
  expected_h_c_111: 0.3326
---

# Task T2.3.3: Global Phase Diagram Mapping & First-Order Boundary Anisotropy

## 1. Physical Motivation & Theoretical Context
In the Kitaev honeycomb model under a magnetic field $\mathbf{h} = (h_x, h_y, h_z)$, the ground state undergoes a first-order quantum phase transition between the non-Abelian quantum spin liquid (QSL) and the high-field polarized phase.
A key open question in the theory of Kitaev magnets is how the thermodynamic transition field $h_c^*(\hat{\mathbf{h}})$ depends on field orientation:
- Is the transition isotropic across all angles?
- Or does the Kitaev exchange anisotropy create an orientation-dependent stability surface $h_c^*(\theta)$?
Along the trajectory connecting the three high-symmetry axes $[001] \to [111] \to [110]$ in the $(1, 1, 0)$ reflection plane:
$$\mathbf{h}(\theta) = h_{\mathrm{mag}} \left( \frac{\sin\theta}{\sqrt{2}}, \frac{\sin\theta}{\sqrt{2}}, \cos\theta \right),$$
the $[111]$ direction corresponds to equal field projection onto all three orthogonal Kitaev exchange links ($x, y, z$). By contrast, fields along $[001]$ or $[110]$ selectively polarize a subset of bonds, inducing bond frustration imbalance.
Determining the exact profile $h_c^*(\theta)$ and evaluating Bayesian model selection tests whether the QSL phase is maximally robust along $[111]$.

## 2. Pre-Registered Hypotheses & Falsification Criteria
- **Model $\mathcal{M}_0$ (Null):** Isotropic first-order boundary $h_c^*(\theta) = h_0$ ($k_0 = 1$ parameter).
- **Model $\mathcal{M}_1$ (Candidate):** Anisotropic boundary $h_c^*(\theta) = h_{\mathrm{base}} + \Delta h \cdot \sin^4\theta \cos^2\theta / (4/27)$ ($k_1 = 2$ parameters, peaking at $\theta = \arccos(1/\sqrt{3}) \approx 54.74^\circ$).
- **Falsification Thresholds:**
  - Decisive Bayesian evidence: $\Delta\mathrm{BIC} = \mathrm{BIC}_0 - \mathrm{BIC}_1 > 10.0$ and Bayes factor $B_{10} > 148.0$.
  - Peak location: $|\theta_{\mathrm{max}} - 54.74^\circ| \le 3.0^\circ$.

---

## 3. Empirical Benchmark & Bayesian Model Comparison Results

Executing `scripts/phase_2/m2_3/benchmark_phase_boundary.py` across 13 orientations connecting $[001] \to [111] \to [110]$ in the $(1, 1, 0)$ reflection plane produced the following data:

### 3.1 First-Order Thermodynamic Critical Field Profile $h_c^*(\theta)$

| $\theta$ ($^\circ$) | Field Unit Vector $\hat{\mathbf{h}}$ | Critical Field $h_c^*$ ($|K|$) | Ground Energy $E_{\mathrm{cross}}$ ($|K|$) | Physical Symmetry |
|:---:|:---:|:---:|:---:|:---:|
| $0.00^\circ$ | $[0, 0, 1]$ | $0.33162$ | $-0.83162$ | $[001]$ high-symmetry axis |
| $7.50^\circ$ | $(0.092, 0.092, 0.991)$ | $0.33163$ | $-0.83163$ | Anisotropic transition |
| $15.00^\circ$ | $(0.183, 0.183, 0.966)$ | $0.33169$ | $-0.83169$ | Anisotropic transition |
| $22.50^\circ$ | $(0.271, 0.271, 0.924)$ | $0.33182$ | $-0.83182$ | Anisotropic transition |
| $30.00^\circ$ | $(0.354, 0.354, 0.866)$ | $0.33202$ | $-0.83202$ | Intermediate angle |
| $37.50^\circ$ | $(0.430, 0.430, 0.793)$ | $0.33226$ | $-0.83226$ | Approaching $[111]$ |
| $45.00^\circ$ | $(0.500, 0.500, 0.707)$ | $0.33247$ | $-0.83247$ | Approaching $[111]$ |
| $52.50^\circ$ | $(0.561, 0.561, 0.609)$ | $0.33260$ | $-0.83260$ | **Global Maximum (near $[111]$)** |
| $60.00^\circ$ | $(0.612, 0.612, 0.500)$ | $0.33256$ | $-0.83256$ | Decreasing |
| $67.50^\circ$ | $(0.653, 0.653, 0.383)$ | $0.33236$ | $-0.83236$ | Decreasing |
| $75.00^\circ$ | $(0.683, 0.683, 0.259)$ | $0.33206$ | $-0.83206$ | Decreasing |
| $82.50^\circ$ | $(0.701, 0.701, 0.131)$ | $0.33176$ | $-0.83176$ | Approaching $[110]$ |
| $90.00^\circ$ | $[1/\sqrt{2}, 1/\sqrt{2}, 0]$ | $0.33162$ | $-0.83162$ | $[110]$ high-symmetry axis |

### 3.2 Peak Location & Maximum QSL Stability

- **Fitted Baseline:** $h_{\mathrm{base}} = 0.33166 \pm 0.00001$.
- **Modulation Amplitude:** $\Delta h = +0.00098 \pm 0.00001$ ($+0.98 \times 10^{-3} |K|$ enhancement at $[111]$).
- **Peak Location:** $\theta_{\mathrm{max}} = 52.50^\circ$ (discrete grid maximum at $2.24^\circ$ from the exact $[111]$ angle $\arccos(1/\sqrt{3}) \approx 54.74^\circ$).
- Satisfies pre-registered contract tolerance $|\theta_{\mathrm{max}} - 54.74^\circ| \le 3.0^\circ$.

### 3.3 Bayesian Model Selection (Anisotropic Boundary vs Isotropic Null)

- **Residual Sum of Squares:**
  $$\mathrm{RSS}_{\mathrm{null}} = 1.6901 \times 10^{-6}, \quad \mathrm{RSS}_{\mathrm{cand}} = 1.0975 \times 10^{-8} \quad (154\times \text{ lower variance})$$
- **Information Criteria:**
  $$\mathrm{BIC}_{\mathrm{null}} = -198.81, \quad \mathrm{BIC}_{\mathrm{cand}} = -261.73$$
  $$\Delta\mathrm{BIC}_{01} = +62.92 \gg 10.0, \quad B_{10} = 4.59 \times 10^{13} \gg 148.0$$
- **Verdict:** Decisive evidence confirming anisotropic modulation of the first-order thermodynamic boundary.

### 3.4 Key Scientific Discoveries
1. **Maximal Topological Stability Along $[111]$:**
   The non-Abelian quantum spin liquid is maximally protected against field-induced breakdown along the $[111]$ axis ($h_c^* = 0.33260$), compared to $[001]$ ($0.33162$) and $[110]$ ($0.33162$). Symmetrical distribution of field across all three Kitaev bond types suppresses preferential bond-polarization channels.
2. **Phase Boundary Point-Group Symmetry in 3D Spin Space:**
   The Hamiltonian possesses cubic point-group symmetry ($O_h$) permuting the three orthogonal spin axes $(x, y, z)$. Within the $(1, 1, 0)$ sweep plane where $h_x = h_y$, $[111]$ is the 3-fold rotation axis ($C_3$). The leading cubic invariant under $h_x = h_y = h\sin\theta/\sqrt{2}$ and $h_z = h\cos\theta$ is $(h_x h_y h_z)^2 \propto \sin^4\theta \cos^2\theta$, attaining its global maximum precisely at $\theta = \arccos(1/\sqrt{3}) \approx 54.74^\circ$.
3. **Decisive Statistical Discrimination:**
   Bayesian odds $B_{10} = 4.59 \times 10^{13} \gg 148.0$ decisively falsify isotropic boundary models in favor of the $\sin^4\theta \cos^2\theta$ harmonic modulation.

