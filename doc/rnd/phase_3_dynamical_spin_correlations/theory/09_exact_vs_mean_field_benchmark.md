---
title: "Task T3.5.1: Systematic Exact (Track A) vs. 4-Majorana Mean-Field (Track B) Cross-Validation Benchmark"
status: completed
author: "Lu Freitas"
date: "2026-11-07"
phase: "Phase 3"
milestone: "M3.5"
task: "T3.5.1"
contracts: ["H_COMP_01", "H_COMP_02", "H_COMP_03", "H_COMP_04"]
primary_model: "Pure FM Kitaev (K > 0) with chiral 3-spin term kappa (arXiv:2209.12930)"
geometry: "Real-space honeycomb lattice clusters L in {4, 6, 8}"
guiding_questions:
  - "How do the exact non-equilibrium gauge quench (Track A) and 4-Majorana mean-field bubble (Track B) compare quantitatively across excitation thresholds, peak structures, and sum rules?"
  - "How does the spectral threshold shift Delta omega = omega_onset^MF - Delta E_2v scale with chiral 3-spin coupling kappa in [0.05, 0.25] |K|?"
  - "How is the vortex-core MZM peak selectivity in Track A contrasted against the bare mean-field vacuum of Track B?"
  - "What are the physical implications for gauge projection and vertex corrections in 4-Majorana mean-field theories of quantum spin liquids?"
---

# Task T3.5.1: Systematic Exact (Track A) vs. 4-Majorana Mean-Field (Track B) Cross-Validation Benchmark

## 1. Scope & Scientific Objectives
Milestone **M3.5** delivers the final synthesis of the **Phase 3** program: a systematic, quantitative cross-validation comparing the **exact non-equilibrium gauge quench theory** (Track A, reproducing arXiv:2209.12930) and the **4-Majorana mean-field theory** (Track B, implemented in `kitaev-mf`).

The scientific objectives of Task **T3.5.1** are:
1. **Unified Cross-Validation Benchmark Engine (`H_COMP_01`):**
   Implement [`scripts/phase_3/benchmark_exact_vs_mean_field.py`](file:///home/lucas/Projects/kitaev-mf/scripts/phase_3/benchmark_exact_vs_mean_field.py) executing synchronized comparisons between Track A and Track B across cluster sizes $L \in \{4, 6, 8\}$ and chirality parameters $\kappa \in [0.05, 0.25] |K|$.
2. **Threshold Shift & Scaling Analysis (`H_COMP_02`):**
   Evaluate the energy threshold difference $\Delta \omega(\kappa) = \omega_{\mathrm{onset}}^{\mathrm{MF}}(\kappa) - \Delta E_{2v}(\kappa)$ across the parameter sweep, demonstrating that the upward shift is an intrinsic consequence of treating gauge Majoranas as deconfined Bogoliubov quasiparticles.
3. **MZM Core Selectivity Contrast (`H_COMP_03`):**
   Showcase the selective emergence of the MZM zero-bias peak in Track A ($S_{\mathrm{core}}(0) > 10 / |K|$) versus the complete absence of subgap weight in Track B ($S_{\mathrm{MF}}(\omega \le 1.0) = 0$), explaining the necessity of gauge vortex pinning.
4. **Master Comparative Table & Structured Data (`H_COMP_04`):**
   Export comprehensive numerical datasets to `doc/labbook/figures/phase_3/benchmark_exact_vs_mf.dat` and record certified findings in Front 1 and Front 2 documentation.

---

## 2. Theoretical & Methodological Comparison

| Feature / Dimension | Track A: Exact Non-Equilibrium Gauge Quench | Track B: 4-Majorana Mean-Field Polarization Bubble |
| :--- | :--- | :--- |
| **Foundational Reference** | Knolle et al. (2014), arXiv:2209.12930 | `kitaev-mf` HFB mean-field theory |
| **Gauge link operator $u_{jk}^\gamma$** | Classical static $\mathbb{Z}_2$ variable ($u_{jk} = \pm 1$) | Decoupled dynamical Majorana pair $\mathrm{i} b_j^\gamma b_k^\gamma$ |
| **Spin operator action** | Sudden link quench $u_b \to -u_b$ | Creation of two-quasiparticle pair $a_{k_1}^\dagger a_{k_2}^\dagger$ |
| **Intermediate states** | Non-equilibrium sector $G' = G \cdot \sigma_j$ | Single self-consistent quasiparticle vacuum $|0\rangle$ |
| **Transition matrix element** | Macroscopic Pfaffian overlap $\langle 0_b | 0_a \rangle$ & Thouless coefficients | Bilinear Wick determinant $4(\phi_{k_1 p} \phi_{k_2 q} - \phi_{k_2 p} \phi_{k_1 q})$ |
| **Lowest excitation energy** | Two-vison gap $\Delta E_{2v} \approx 0.26\text{--}0.48 |K|$ | Two-quasiparticle onset $\epsilon_0 + \epsilon_1 \approx 3.05 |K|$ |
| **Majorana Zero Mode (MZM)** | Sharp zero-bias peak ($S \approx 12.8$) bound to vortex core | Subgap silent; requires explicit localized flux ansatz |
| **Sum rule $\int \frac{\mathrm{d}\omega}{2\pi} S(\omega)$** | Saturated ($1.000 \pm 0.05$) via Lehmann cutoff | Saturated analytically ($1.000000000000$) to machine precision |

---

## 3. Implementation Details
The master benchmark suite is implemented in:
- [`scripts/phase_3/benchmark_exact_vs_mean_field.py`](file:///home/lucas/Projects/kitaev-mf/scripts/phase_3/benchmark_exact_vs_mean_field.py)
- [`test/unit/test_phase3_cross_validation.py`](file:///home/lucas/Projects/kitaev-mf/test/unit/test_phase3_cross_validation.py)
- Graphical assets:
  - `doc/labbook/figures/phase_3/fig5_exact_vs_mean_field.png` (300 DPI, 571 kB)
  - `doc/labbook/figures/phase_3/fig5_exact_vs_mean_field.pdf` (vector, 45 kB)
  - `doc/labbook/figures/phase_3/benchmark_exact_vs_mf.dat` (ASCII table, 0.8 kB)

---

## 4. Quantitative Results & Master Sweep Table

| Chirality $\kappa / \|K\|$ | Track A: Core Onset $\omega_{\mathrm{core}} / \|K\|$ | Track A: Vison Gap $\Delta E_{2v} / \|K\|$ | Track B: MF Threshold $\omega_{\mathrm{onset}}^{\mathrm{MF}} / \|K\|$ | Threshold Shift $\Delta \omega / \|K\|$ | Relative Shift $\Delta \omega / \Delta E_{2v}$ | MZM Splitting $\delta \epsilon / \|K\|$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **$0.05$** | $0.2140$ | $0.4805$ | $3.0497$ | $+2.5692$ | $5.35\times$ | $0.0542$ |
| **$0.10$** | $0.1872$ | $0.4714$ | $3.0497$ | $+2.5783$ | $5.47\times$ | $0.0245$ |
| **$0.15$** | $0.1564$ | $0.4564$ | $3.0497$ | $+2.5934$ | $5.68\times$ | $0.0120$ |
| **$0.20$** | $0.1250$ | $0.4361$ | $3.0497$ | $+2.6136$ | $5.99\times$ | $0.0058$ |
| **$0.25$** | $0.0911$ | $0.4117$ | $3.0497$ | $+2.6381$ | $6.41\times$ | $0.0027$ |

### 4.1 Key Findings:
1. **Intrinsic Quasiparticle Continuum Shift:** The mean-field threshold $\omega_{\mathrm{onset}}^{\mathrm{MF}} \approx 3.05 |K|$ strictly exceeds the exact two-vison threshold $\Delta E_{2v} \approx 0.41\text{--}0.48 |K|$ by over $2.56 |K|$ across all tested parameters ($5.35\times\text{--}6.41\times$ larger). This discrepancy is an inescapable feature of the bare 4-Majorana mean-field approximation, which neglects gauge confinement and projects gauge link fluctuations into high-energy itinerant quasiparticle pairs.
2. **Selective MZM Binding:** Track A clearly separates the vortex translation channel ($\omega \approx 0.09\text{--}0.21 |K|$) from the bulk vison creation threshold ($\Delta E_{2v} \ge 0.41 |K|$), yielding a sharp zero-bias STS signature. Track B on a translationally invariant state has strictly zero subgap weight ($S_{\mathrm{MF}}(\omega \le 1.0) = 0$).
3. **MZM Hybridization Scaling:** As chiral coupling $\kappa$ increases from $0.05$ to $0.25 |K|$, the bulk band gap $\Delta_{\mathrm{bulk}} = 6\sqrt{3}\kappa$ widens, sharpening the Majorana localization length $\xi_M \sim 1/\kappa$ and suppressing the finite-size hybridization splitting by over a factor of 20 (from $0.0542 |K|$ down to $0.0027 |K|$).

---

## 5. Contract Verification Matrix

| Contract ID | Metric / Assertion | Expected Value | Observed Empirical Value | Status |
| :--- | :--- | :--- | :--- | :--- |
| **`H_COMP_01`** | Benchmark Engine Execution | Valid outputs across $\kappa$ sweep | 5 records, zero errors or NaNs | **PASSED** |
| **`H_COMP_02`** | Threshold Gap Discrepancy | $\Delta \omega > 1.50 |K|$ | $\Delta \omega \ge 2.5692 |K|$ ($5.35\times\text{--}6.41\times$) | **PASSED** |
| **`H_COMP_02`** | Exact Vison Gap Upper Bound | $\Delta E_{2v} < 0.50 |K|$ | $\max \Delta E_{2v} = 0.4805 |K|$ | **PASSED** |
| **`H_COMP_03`** | Core MZM Low-Bias Onset | $\omega_{\mathrm{core}} \le 0.25 |K|$ | $\omega_{\mathrm{core}} \in [0.091, 0.214] |K|$ | **PASSED** |
| **`H_COMP_03`** | Bulk Vison Lower Bound | $\Delta E_{2v} \ge 0.40 |K|$ | $\Delta E_{2v} \in [0.412, 0.481] |K|$ | **PASSED** |
| **`H_COMP_04`** | Asset File Sizes | PNG $> 25\text{ kB}$, PDF $> 15\text{ kB}$ | PNG $= 571\text{ kB}$, PDF $= 45\text{ kB}$ | **PASSED** |
| **`H_COMP_04`** | Master Tabular Formatting | ASCII table with 7 columns | 5 rows, 7 columns, valid headers | **PASSED** |
