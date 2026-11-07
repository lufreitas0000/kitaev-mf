---
id: T2.3.1
milestone: M2.3
title: "Global Field-Directional Magnetization Surface & Torque Discontinuities"
status: "In Progress"
created: 2026-11-07
authors: ["Antigravity Agent", "Lu Freitas"]
hypothesis:
  null_model: "M0_isotropic: Magnetic torque tau = m x h is identically zero across all field angles (isotropic magnetization response)."
  candidate_model: "M1_anisotropic_torque: Kitaev exchange anisotropy induces an angular magnetic torque tau_perp(theta) with nodes at high-symmetry axes [001], [111], and [110], exhibiting discontinuous jumps across the first-order spinodal boundary."
  falsification_criteria:
    delta_bic_threshold: 10.0
    bayes_factor_threshold: 148.0
    symmetry_node_tolerance: 1.0e-5
numerical_targets:
  field_magnitudes: [0.10, 0.20, 0.30, 0.35, 0.40]
  angular_resolution: 13
  symmetry_directions: ["[001]", "[111]", "[110]"]
---

# Task T2.3.1: Global Field-Directional Magnetization Surface & Torque Discontinuities

## 1. Physical Motivation & Theoretical Context
In the pure ferromagnetic Kitaev honeycomb model ($K = -1.0$), application of a magnetic field $\mathbf{h}$ in arbitrary spatial directions breaks the exchange symmetry.
While along the high-symmetry axes:
- $[001] \quad (\theta = 0^\circ)$
- $[111] \quad (\theta \approx 54.74^\circ)$
- $[110] \quad (\theta = 90^\circ)$
the induced magnetization vector $\mathbf{m} = \frac{1}{N} \sum_s \langle \boldsymbol{\sigma}_s \rangle$ aligns with the applied field ($\mathbf{m} \parallel \mathbf{h}$), off-axis field orientations induce magnetic torque:
$$\boldsymbol{\tau} = \mathbf{m} \times \mathbf{h} \ne \mathbf{0}.$$
In magnetic torque magnetometry experiments (such as cantilever torque magnetometry on $\alpha\text{-RuCl}_3$), this torque directly probes the anisotropic free energy gradient $\boldsymbol{\tau} = -\partial F / \partial \hat{\mathbf{h}}$. Furthermore, across the first-order phase boundary separating the quantum spin liquid (QSL) from the field-polarized phase, the magnetization magnitude undergoes a first-order jump $\Delta m > 0$, causing discontinuous jumps in the magnetic torque profile.

## 2. Pre-Registered Hypothesis & Falsification Criteria
- **Model $\mathcal{M}_0$ (Null):** Isotropic linear response $\boldsymbol{\tau} = \mathbf{0}$ ($\mathrm{RSS}_{\mathrm{null}} = \sum_i \|\boldsymbol{\tau}_i\|^2$).
- **Model $\mathcal{M}_1$ (Candidate):** Anisotropic cubic-symmetric torque $\tau_\perp(\theta) = C \sin(2\theta)(3\cos^2\theta - 1)$ preserving high-symmetry nodes at $\theta = 0, \arccos(1/\sqrt{3}), \pi/2$.
- **Decisive Falsification Thresholds:**
  - $\Delta\mathrm{BIC} = \mathrm{BIC}_0 - \mathrm{BIC}_1 > 10.0$
  - Bayes factor $B_{10} = \exp(0.5 \Delta\mathrm{BIC}) > 148.0$
  - Node alignment: $\|\boldsymbol{\tau}(\theta_{\mathrm{node}})\| < 10^{-5}$ at all 3 high-symmetry axes.

## 3. Mathematical & Algorithmic Implementation
1. Parameterize field tilt in the $(1, 1, 0)$ plane connecting $[001]$ and $[110]$:
   $$\mathbf{h}(\theta) = h_{\mathrm{mag}} \left( \frac{\sin\theta}{\sqrt{2}}, \frac{\sin\theta}{\sqrt{2}}, \cos\theta \right)$$
2. Solve self-consistent mean-field state with Grassmannian DIIS accelerator (`solve_self_consistent_mean_field`).
3. Compute sublattice magnetization vectors $\langle \boldsymbol{\sigma}_s \rangle_G$ from on-site Grassmann matrices $V_s$.
4. Calculate magnetic torque $\boldsymbol{\tau} = \mathbf{m} \times \mathbf{h}$ and its projection along $(1, -1, 0)/\sqrt{2}$.
5. Perform Bayesian model comparison between isotropic null $\mathcal{M}_0$ and anisotropic candidate $\mathcal{M}_1$.

---

## 4. Empirical Benchmark & Bayesian Model Comparison Results

Executing `scripts/phase_2/m2_3/benchmark_magnetization_torque.py` across 13 orientations in the $(1, 1, 0)$ plane and 3 representative fields ($h \in \{0.10, 0.20, 0.35\}$) generated 39 benchmark data points:

### 4.1 Magnetic Torque Summary

| $h$ | $\theta$ ($[001] \to [110]$) | $\|\mathbf{m}\|$ | $\tau_{\perp}$ ($\times 10^{-4}$) | $\|\boldsymbol{\tau}\|$ ($\times 10^{-4}$) | Physical Regime |
|:---:|:---:|:---:|:---:|:---:|:---:|
| $0.10$ | $0.0^\circ$ ($[001]$) | $0.0760$ | $+0.0000$ | $0.0000$ | High-symmetry node |
| $0.10$ | $30.0^\circ$ | $0.0761$ | $+0.0415$ | $0.0415$ | Weak anisotropic response |
| $0.10$ | $54.7^\circ$ ($[111]$) | $0.0761$ | $+0.0000$ | $0.0000$ | High-symmetry node |
| $0.10$ | $75.0^\circ$ | $0.0761$ | $-0.0212$ | $0.0212$ | Inverted torque |
| $0.10$ | $90.0^\circ$ ($[110]$) | $0.0761$ | $+0.0000$ | $0.0000$ | High-symmetry node |
| $0.20$ | $30.0^\circ$ | $0.1576$ | $+0.8214$ | $0.8214$ | Intermediate field |
| $0.20$ | $54.7^\circ$ ($[111]$) | $0.1582$ | $+0.0000$ | $0.0000$ | High-symmetry node |
| $0.20$ | $75.0^\circ$ | $0.1577$ | $-0.7932$ | $0.7932$ | Inverted torque |
| $0.35$ | $0.0^\circ$ ($[001]$) | $0.3122$ | $+0.0000$ | $0.0000$ | High-symmetry node |
| $0.35$ | $30.0^\circ$ | $0.3214$ | $+17.341$ | $17.341$ | Critical enhancement ($\approx 418\times$ vs $h=0.10$) |
| $0.35$ | $54.7^\circ$ ($[111]$) | $0.3318$ | $+0.0000$ | $0.0000$ | High-symmetry node |
| $0.35$ | $75.0^\circ$ | $0.3228$ | $-24.234$ | $24.234$ | Critical enhancement |
| $0.35$ | $90.0^\circ$ ($[110]$) | $0.3114$ | $+0.0000$ | $0.0000$ | High-symmetry node |

### 4.2 Bayesian Model Selection

- **Residual Sum of Squares:**
  $$\mathrm{RSS}_{\mathrm{null}} = 2.4945 \times 10^{-5}, \quad \mathrm{RSS}_{\mathrm{cand}} = 9.1490 \times 10^{-6}$$
- **Bayesian Information Criterion:**
  $$\mathrm{BIC}_{\mathrm{null}} = -556.23, \quad \mathrm{BIC}_{\mathrm{cand}} = -584.36$$
  $$\Delta\mathrm{BIC}_{01} = +28.13 \gg 10.0, \quad B_{10} = 1.28 \times 10^6 \gg 148.0$$
- **Verdict:** Decisive evidence favoring anisotropic Kitaev magnetic torque over isotropic linear response.

### 4.3 Key Scientific Discoveries
1. **Exact High-Symmetry Nodes:**
   Along $[001]$, $[111]$, and $[110]$, the magnetic torque $\|\boldsymbol{\tau}\|$ vanishes to numerical precision ($< 10^{-11}$), confirming exact field alignment $\mathbf{m} \parallel \mathbf{h}$ enforced by the honeycomb point-group symmetries.
2. **Torque Sign Inversion and Easy-Axis Stability of $[111]$:**
   In the $[001] \to [111]$ quadrant, $\tau_\perp > 0$; in the $[111] \to [110]$ quadrant, $\tau_\perp < 0$. Because $\tau_\perp = -\partial F/\partial\theta$, the negative derivative $d\tau_\perp/d\theta < 0$ implies positive free-energy curvature $\partial^2 F/\partial\theta^2 > 0$. The torque is therefore a restoring torque driving magnetization toward $[111]$, establishing that $[111]$ is a **stable local energy minimum (easy axis)** along the $(1, 1, 0)$ sweep trajectory.
3. **Critical Non-linear Enhancement near Spinodal Breakdown:**
   As $h$ approaches the first-order spinodal field $h_{c2}^+ \approx 0.400$, the torque amplitude grows non-linearly by over $\approx 418\times$ between $h=0.10$ and $h=0.35$ (from $4.15 \times 10^{-6}$ to $1.734 \times 10^{-3}$ at $\theta=30^\circ$), reflecting the critical softening and divergence of the magnetic susceptibility tensor $\chi_{\alpha\beta}$ near the spinodal boundary.

