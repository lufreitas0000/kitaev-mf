---
title: "Task T3.3.1: Cylinder Geometry, Chiral Majorana Edge Mode Spectroscopy, and Finite-Size Scaling"
status: completed
author: "Lu Freitas"
date: "2026-11-07"
phase: "Phase 3"
milestone: "M3.3"
task: "T3.3.1"
contracts: ["H_EDGE_01", "H_EDGE_02", "H_EDGE_03", "H_EDGE_04"]
primary_model: "Pure FM Kitaev (K > 0) with chiral 3-spin term kappa (arXiv:2209.12930)"
geometry: "Honeycomb cylinder (periodic along a1, open along a2) and scalable torus clusters L in {4, 6, 8}"
guiding_questions:
  - "How is the honeycomb cylinder geometry (open boundary along one direction, periodic along the other) rigorously constructed without breaking gauge-matter conventions?"
  - "How do gapless/subgap chiral Majorana edge modes manifest in the real-space energy spectrum when time-reversal symmetry is broken by kappa > 0?"
  - "How does the local dynamical spin correlation S_edge(omega) on the open boundary compare with the bulk S_bulk(omega), and does it verify edge-state cotunneling?"
  - "How does the in-gap MZM splitting delta_eps = |eps_1 - eps_0| scale with cluster size L in {4, 6, 8} and vortex separation R_12?"
---

# Task T3.3.1: Cylinder Geometry, Chiral Majorana Edge Mode Spectroscopy, and Finite-Size Scaling

## 1. Scope & Scientific Objectives
In chiral topological phases with non-zero Chern number $\nu = \pm 1$ (such as the gapless $B$-phase of the Kitaev honeycomb model under a [111] magnetic field), the bulk-boundary correspondence guarantees the existence of chiral gapless Majorana edge modes on open system boundaries.

The scientific objectives of Task **T3.3.1** are:
1. **Honeycomb Cylinder Geometry (`H_EDGE_01`):**
   Implement an exact lattice constructor for honeycomb cylinders of dimension $L_x \times L_y$ (periodic along $\mathbf{a}_1$, open with zigzag or armchair termination along $\mathbf{a}_2$), correctly assembling directed nearest-neighbor bonds and chiral 3-spin next-nearest-neighbor triplets without dangling boundary triplets that violate hermiticity/skew-symmetry.
2. **Chiral Edge Mode Isolation (`H_EDGE_02`):**
   Diagonalize the matter sector on the cylinder with $\kappa > 0$ and verify that the bulk gap $\Delta_{\mathrm{bulk}} = 6\sqrt{3}\kappa$ is bridged by localized boundary states whose probability density decays exponentially into the bulk: $\rho_{\mathrm{edge}}(x) \propto \exp(-x/\xi_{\mathrm{edge}})$.
3. **Local Edge vs. Bulk STS Response (`H_EDGE_03`):**
   Evaluate the local dynamical spin structure factor $S_{jj}^{\gamma\gamma}(\omega)$ for boundary sites versus interior bulk sites. Verify that edge sites exhibit finite spectral weight at low bias $\omega < \Delta E_{2v}$ mediated by cotunneling into the chiral edge channel, in stark contrast to the hard vison creation threshold in the bulk.
4. **Finite-Size Scaling of MZM Hybridization (`H_EDGE_04`):**
   Perform a multi-size sweep $L \in \{4, 6, 8\}$ on torus clusters with trapped vortex pairs, quantifying the finite-size hybridization splitting $\delta \epsilon(L) = |\epsilon_1 - \epsilon_0|$ as a function of the vortex distance $R_{12}$, establishing the exponential asymptotic decoupling $\delta \epsilon \to 0$ as $L \to \infty$.

---

## 2. Theoretical Architecture

### 2.1 Cylinder Geometry & Boundary Triplet Truncation
On a cylinder with $L_x$ unit cells along the periodic direction $\mathbf{a}_1$ and $L_y$ unit cells along the open direction $\mathbf{a}_2$:
- Periodic boundary conditions along $\mathbf{a}_1$: cell index $m_1 \equiv m_1 \pmod{L_x}$.
- Open boundary conditions along $\mathbf{a}_2$: cell index $m_2 \in \{0, 1, \dots, L_y - 1\}$.
- Bonds crossing $m_2 = -1$ or $m_2 = L_y$ are strictly omitted (cutting exactly $L_x$ bonds: $N_{\mathrm{bonds}} = 3 L_x L_y - L_x$).
- **Hermitian Chirality Truncation:** Next-nearest-neighbor paths $j \to k \to l$ that traverse bonds exiting the cylinder boundaries are strictly excluded. The skew-symmetric matter matrix $A = -A^T$ remains exact to machine precision ($< 10^{-15}$).

### 2.2 Bulk-Boundary Correspondence & Chiral Edge Dispersion
In reciprocal space along the periodic boundary direction $k_y \in [-\pi, \pi]$:
$$E_{\mathrm{edge}}(k_y) \approx v_M k_y, \quad v_M \approx 6\kappa.$$
In real space on an $L_x \times L_y$ cylinder, this yields a dense ladder of in-gap boundary modes spanning the energy window $[0, \Delta_{\mathrm{bulk}}]$.

### 2.3 Finite-Size Majorana Hybridization
When two $\mathbb{Z}_2$ vortices are separated by distance $R_{12} = \mathcal{O}(L)$:
$$\delta \epsilon = |\epsilon_1 - \epsilon_0| \approx \Delta_0 \sqrt{\frac{\xi_M}{R_{12}}} \cos(k_F R_{12}) \mathrm{e}^{-R_{12}/\xi_M}.$$
On finite clusters, $\delta \epsilon > 0$ represents the finite-size overlap of the two Majorana wavefunctions bound to the distinct vortex cores. As $L$ increases from 4 to 8, $\delta \epsilon$ decreases monotonically toward the exact thermodynamic topological degeneracy $\delta \epsilon \to 0$.

---

## 3. Implementation Details
The cylinder and finite-size scaling engines are implemented in:
- [`src/kitaev_mf/services/phase_3/gauge_theory.py`](file:///home/lucas/Projects/kitaev-mf/src/kitaev_mf/services/phase_3/gauge_theory.py) (`build_honeycomb_cylinder_cluster`)
- [`src/kitaev_mf/services/phase_3/edge_and_scaling.py`](file:///home/lucas/Projects/kitaev-mf/src/kitaev_mf/services/phase_3/edge_and_scaling.py)
- [`test/unit/test_phase3_edge_and_scaling.py`](file:///home/lucas/Projects/kitaev-mf/test/unit/test_phase3_edge_and_scaling.py)

---

## 4. Quantitative Results & Contract Verification

### 4.1 Cylinder In-Gap Chiral Edge Modes ($L_x=6, L_y=12, \kappa=0.15$)
- **Bulk Gap Estimate:** $\Delta_{\mathrm{bulk}} = 6\sqrt{3}\kappa \approx 1.5588 |K|$.
- **In-gap Mode Count:** 11 discrete states below $\Delta_{\mathrm{bulk}}$.
- **Lowest Edge Mode Energy:** $\epsilon_0 = 3.20 \times 10^{-10} |K|$ (gapless edge limit).
- **Boundary Localization:** $99.94\%$ of the probability density of the lowest mode resides within the 2 outermost boundary layers ($m_2 \in \{0, 1, 10, 11\}$).

### 4.2 Edge vs. Bulk Local STS Response ($L_x=4, L_y=6, \sigma=0.06$)
- **Low-Bias Window ($\omega \in [0.05, 0.15] |K|$):**
  - Edge Site ($m_2 = 0$): $\max S_{\mathrm{edge}}(\omega) = 16.03 / |K|$ (chiral edge channel cotunneling).
  - Bulk Site ($m_2 = 3$): $\max S_{\mathrm{bulk}}(\omega) = 1.12 \times 10^{-3} / |K|$ (strictly silent subgap).
  - **Low-Bias Contrast Ratio:** $S_{\mathrm{edge}} / S_{\mathrm{bulk}} = 1.43 \times 10^4 \gg 50.0$.

### 4.3 Finite-Size Scaling of MZM Doublet Splitting on Torus
| Cluster Size $L$ | Physical Sites $N$ | Vortex Separation $R_{12} / a$ | $\epsilon_0 / \|K\|$ | $\epsilon_1 / \|K\|$ | Doublet Splitting $\delta \epsilon / \|K\|$ | In-Gap Ratio $\epsilon_1 / \epsilon_2$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **$L = 4$** | 32 | $3.464$ | $0.2946$ | $1.1171$ | $8.225 \times 10^{-1}$ | $0.7099$ |
| **$L = 6$** | 72 | $5.196$ | $0.4045$ | $0.4165$ | $1.202 \times 10^{-2}$ | $0.3815$ |
| **$L = 8$** | 128 | $6.928$ | $0.1721$ | $0.1787$ | $6.581 \times 10^{-3}$ | $0.3001$ |

The splitting drops by over 2 orders of magnitude from $L=4$ to $L=8$, demonstrating exponential asymptotic decoupling and topological protection.

---

## 5. Contract Verification Matrix

| Contract ID | Metric / Assertion | Theoretical Target | Achieved Empirical Value | Status |
| :--- | :--- | :--- | :--- | :--- |
| **`H_EDGE_01`** | Skew-Symmetry Error $\|A + A^T\|$ | $< 10^{-15}$ | $0.00 \times 10^{-16}$ | **PASSED** |
| **`H_EDGE_01`** | Cylinder Bond Count ($L_x=6, L_y=8$) | $3 L_x L_y - L_x = 138$ | 138 bonds | **PASSED** |
| **`H_EDGE_02`** | Chiral Edge Mode In-Gap Count | $\ge 4$ states | 11 states below $\Delta_{\mathrm{bulk}}$ | **PASSED** |
| **`H_EDGE_02`** | Edge Boundary Localization Fraction | $> 90.0\%$ | $99.94\%$ | **PASSED** |
| **`H_EDGE_03`** | Edge Low-Bias Spectral Weight | $> 5.0 / \|K\|$ | $16.03 / \|K\|$ | **PASSED** |
| **`H_EDGE_03`** | Bulk Subgap Weight ($\omega \le 0.15$) | $< 0.05 / \|K\|$ | $1.12 \times 10^{-3} / \|K\|$ | **PASSED** |
| **`H_EDGE_03`** | Edge-to-Bulk Low-Bias Contrast | $> 50.0\times$ | $1.43 \times 10^4\times$ | **PASSED** |
| **`H_EDGE_04`** | Monotonic Splitting Suppression | $\delta \epsilon(8) < \delta \epsilon(6) < \delta \epsilon(4)$ | $0.0066 < 0.0120 < 0.8225$ | **PASSED** |
| **`H_EDGE_04`** | Asymptotic Splitting at $L=8$ | $< 0.01 / \|K\|$ | $6.58 \times 10^{-3} / \|K\|$ | **PASSED** |
