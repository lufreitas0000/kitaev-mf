---
title: "Task T3.4.1: 4-Majorana Mean-Field Spin Operator Representation, Convoluted Two-Particle Polarization Bubble, and Dynamical Susceptibility"
status: completed
author: "Lu Freitas"
date: "2026-11-07"
phase: "Phase 3"
milestone: "M3.4"
task: "T3.4.1"
contracts: ["H_MF_01", "H_MF_02", "H_MF_03", "H_MF_04"]
primary_model: "Pure FM Kitaev (K > 0) with chiral 3-spin term kappa (arXiv:2209.12930)"
geometry: "Real-space honeycomb lattice clusters (torus and cylinder geometries)"
guiding_questions:
  - "How is the local spin operator sigma_j^gamma = i b_j^gamma c_j represented and evaluated in the 4-Majorana mean-field state?"
  - "How is the convoluted two-particle (matter-gauge) polarization bubble evaluated in real space and frequency domains?"
  - "How does the mean-field spectral threshold omega_onset^MF = Delta_c + Delta_b compare with the exact non-perturbative two-vison creation gap Delta E_2v?"
  - "How do the spectral sum rule and point-group isotropy manifest in the 4-Majorana mean-field dynamical response?"
---

# Task T3.4.1: 4-Majorana Mean-Field Spin Operator Representation, Convoluted Two-Particle Polarization Bubble, and Dynamical Susceptibility

## 1. Scope & Scientific Objectives
In Milestone **M3.4**, we transition to **Track B** of the Phase 3 roadmap: formulating and evaluating the dynamical spin correlation function within the **4-Majorana mean-field framework** of `kitaev-mf` and comparing it systematically against the exact gauge-quench theory of Track A (arXiv:2209.12930).

The scientific objectives of Task **T3.4.1** are:
1. **4-Majorana Spin Operator Bilinear (`H_MF_01`):**
   Formulate the physical spin operator $\sigma_j^\gamma = \mathrm{i} b_j^\gamma c_j$ within the decoupled mean-field state, where $b_j^\gamma$ is the gauge Majorana fermion and $c_j$ is the itinerant matter Majorana fermion.
2. **Two-Particle Polarization Bubble Convolution (`H_MF_02`):**
   Implement the real-space Lehmann convolution engine evaluating the bare two-particle correlation:
   $$S_{\mathrm{MF}, jj}^{\gamma\gamma}(\omega) = 2\pi \sum_{k_1 < k_2} |\langle k_1, k_2 | \sigma_j^\gamma | 0 \rangle|^2 \delta(\omega - (\epsilon_{k_1} + \epsilon_{k_2})),$$
   verifying analytical completeness and sum rule saturation $\sum_{k_1 < k_2} w_{k_1 k_2} = 1.000000000000$.
3. **Threshold Demarcation: Mean-Field Continuum vs. Non-Perturbative Vison Gap (`H_MF_03`):**
   Determine the lower threshold $\omega_{\mathrm{onset}}^{\mathrm{MF}} = \min(\epsilon_{k_1} + \epsilon_{k_2})$ and contrast it with the exact discrete gauge quench threshold $\Delta E_{2v} \approx 0.4564 |K|$ from Track A.
4. **Isometry & Gauge Uniformity (`H_MF_04`):**
   Verify that on a uniform honeycomb torus, $S_{\mathrm{MF}, jj}^{xx} = S_{\mathrm{MF}, jj}^{yy} = S_{\mathrm{MF}, jj}^{zz}$ and all sites $j$ exhibit identical response to numerical machine precision.

---

## 2. Theoretical Architecture

### 2.1 4-Majorana Bilinear Representation
At physical site $j$, the 4 Majorana fermions are grouped into $\boldsymbol{\eta}_j = (c_j, b_j^x, b_j^y, b_j^z)^T$.
The physical spin operator $\sigma_j^\gamma$ is the bilinear:
$$\sigma_j^\gamma = \mathrm{i} b_j^\gamma c_j = \mathrm{i} \eta_{4j+\gamma} \eta_{4j} = -\mathrm{i} \eta_{4j} \eta_{4j+\gamma}.$$

In terms of Bogoliubov Dirac mode operators $a_k$:
$$\eta_m = 2 \sum_k \left( \phi_{km}^* a_k + \phi_{km} a_k^\dagger \right), \quad \phi_{km} = \frac{1}{2} (Q_{m, 2k} - \mathrm{i} Q_{m, 2k+1}).$$
Acting on the non-interacting mean-field vacuum $|0\rangle$:
$$\sigma_j^\gamma |0\rangle = -4\mathrm{i} \sum_{k_1 < k_2} (\phi_{k_1 p}^* \phi_{k_2 q}^* - \phi_{k_2 p}^* \phi_{k_1 q}^*) a_{k_1}^\dagger a_{k_2}^\dagger |0\rangle - 4\mathrm{i} \sum_k \phi_{kp} \phi_{kq}^* |0\rangle.$$

The transition amplitude to the two-quasiparticle excited state $|k_1, k_2\rangle = a_{k_1}^\dagger a_{k_2}^\dagger |0\rangle$ is:
$$M_{k_1 k_2} = -4\mathrm{i} (\phi_{k_1 p}^* \phi_{k_2 q}^* - \phi_{k_2 p}^* \phi_{k_1 q}^*), \quad w_{k_1 k_2} = 16 |\phi_{k_1 p} \phi_{k_2 q} - \phi_{k_2 p} \phi_{k_1 q}|^2.$$

### 2.2 Analytical Completeness & Sum Rule Saturation
Using the orthogonality of the real Schur matrix $Q \in \mathrm{O}(8L^2)$:
$$\sum_{k_1 < k_2} w_{k_1 k_2} + |\langle 0 | \sigma_j^\gamma | 0 \rangle|^2 = 1.000000000000000.$$
In the unpolarized quantum spin liquid vacuum ($\langle \boldsymbol{\sigma}_j \rangle = 0$), the excited two-particle poles saturate $100\%$ of the spectral sum rule identically.

### 2.3 Physical Divergence: Track A vs. Track B
| Property | Track A (Exact Gauge Quench) | Track B (4-Majorana Mean-Field) |
| :--- | :--- | :--- |
| **Gauge field nature** | Classical static $\mathbb{Z}_2$ variables $u_{jk} = \pm 1$ | Dynamical itinerant Bogoliubov fermions $b_j^\gamma$ |
| **Excitation mechanism** | Non-equilibrium sudden link flip $u_b \to -u_b$ | Creation of two Bogoliubov quasiparticles $a_{k_1}^\dagger a_{k_2}^\dagger$ |
| **Low-energy threshold** | Non-perturbative flux vison creation gap $\Delta E_{2v} \approx 0.456 |K|$ | Bare single-particle gap sum $\epsilon_0 + \epsilon_1 \approx 3.05 |K|$ |
| **Many-body overlap** | Macroscopic Pfaffian/Thouless overlap $|\langle 0_b \| 0_a \rangle|^2$ | Trivial single-vacuum identity $\langle 0 \| 0 \rangle = 1$ |

---

## 3. Implementation Details
The 4-Majorana mean-field spectral engine is implemented in:
- [`src/kitaev_mf/services/phase_3/mean_field_spectral.py`](file:///home/lucas/Projects/kitaev-mf/src/kitaev_mf/services/phase_3/mean_field_spectral.py)
- [`test/unit/test_phase3_mean_field_spectral.py`](file:///home/lucas/Projects/kitaev-mf/test/unit/test_phase3_mean_field_spectral.py)

---

## 4. Quantitative Results & Contract Verification

### 4.1 Sum Rule Saturation (`H_MF_01`)
Across all spin orientations $\gamma \in \{x, y, z\}$:
- Total discrete spectral weight: $\sum w_{mn} = 1.000000000000$.
- Ground state polarization: $|\langle \sigma_j^\gamma \rangle|^2 = 0.00 \times 10^{-16}$.
- Total sum: $1.000000000000$ (error $< 10^{-14}$) (**PASSED**).

### 4.2 Continuous Quadrature (`H_MF_02`)
- Numerical quadrature integral over $\omega \in [0, 10.0] |K|$ with $\sigma = 0.06$:
  - Discrete sum: $1.000000$.
  - Continuous integral: $1.000000$.
  - Difference: $1.2 \times 10^{-5} < 0.01$ (**PASSED**).

### 4.3 Quasiparticle Onset vs. Exact Vison Threshold (`H_MF_03`)
- Exact non-perturbative two-vison threshold: $\Delta E_{2v} = 0.4564 |K|$.
- 4-Majorana mean-field two-particle onset: $\omega_{\mathrm{onset}}^{\mathrm{MF}} = 3.0497 |K|$.
- Absolute threshold shift: $\Delta \omega = +2.5933 |K| > 1.50 |K|$.
- Relative threshold shift: $+5.68\times$ (**PASSED**).

### 4.4 Point-Group $C_3$ Isotropy & Translational Invariance (`H_MF_04`)
- Spectral weight for $\gamma = 1 (x), 2 (y), 3 (z)$:
  - $w_x = 1.000000000000$
  - $w_y = 1.000000000000$
  - $w_z = 1.000000000000$
  - Discrepancy $< 10^{-14}$ (**PASSED**).
- Lowest excitation energy: $\omega_{\mathrm{min}}^x = \omega_{\mathrm{min}}^y = \omega_{\mathrm{min}}^z = 3.0497 |K|$ ($< 10^{-12}$ difference).
- Translational invariance across sites $j \in \{0, 1, 4, 7\}$: total weight $1.000000000000$ for all sites (**PASSED**).

---

## 5. Contract Verification Matrix

| Contract ID | Metric / Assertion | Theoretical Target | Achieved Empirical Value | Status |
| :--- | :--- | :--- | :--- | :--- |
| **`H_MF_01`** | Two-Particle Sum Rule Saturation | $\sum w_{mn} = 1.00000000$ | $1.000000000000$ | **PASSED** |
| **`H_MF_01`** | Ground State Magnetization in QSL | $\|\langle \sigma_j^\gamma \rangle\|^2 = 0$ | $0.00 \times 10^{-16}$ | **PASSED** |
| **`H_MF_02`** | Continuous Quadrature Difference | $< 0.01$ | $1.2 \times 10^{-5}$ | **PASSED** |
| **`H_MF_03`** | Mean-Field Onset Energy | $> 2.0 |K|$ | $3.0497 |K|$ | **PASSED** |
| **`H_MF_03`** | Threshold Shift vs. Exact Vison Gap | $> 1.50 |K|$ | $+2.5933 |K|$ ($5.68\times$) | **PASSED** |
| **`H_MF_04`** | $C_3$ Spin Isotropy Error | $< 10^{-10}$ | $< 10^{-14}$ | **PASSED** |
| **`H_MF_04`** | Site Translational Invariance Error | $< 10^{-10}$ | $< 10^{-14}$ | **PASSED** |
