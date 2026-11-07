---
id: T2.4.1
milestone: M2.4
title: "Two-Peak Specific Heat C_v(T) & Fractionalization Entropy Stages (S = 1/2 ln 2)"
status: "Completed"
created: 2026-11-07
authors: ["Antigravity Agent", "Lu Freitas"]
hypothesis:
  null_model: "M0_single_peak: Thermal excitations in the Kitaev magnet are governed by a single characteristic energy scale, yielding a single Schottky/Debye specific heat peak C_v(T) without an intermediate entropy plateau."
  candidate_model: "M1_fractionalized_two_peak: Spin fractionalization into itinerant Majorana fermions (scale T_H ~ 0.5 |K|) and localized Z2 gauge fluxes/visons (scale T_L ~ 0.05 |K|) produces two well-separated specific heat peaks separated by a wide entropy plateau at S(T_mid) = 1/2 ln 2."
  falsification_criteria:
    peak_separation_ratio_min: 5.0
    intermediate_plateau_entropy_rel_tol: 0.05
    high_temp_asymptotic_entropy_rel_tol: 0.02
    delta_bic_threshold: 10.0
    bayes_factor_threshold: 148.0
numerical_targets:
  temperature_range: [0.003, 6.0]
  temperature_points: 40
  vison_gap: 0.137
  plateau_entropy_target: 0.34657
  high_t_entropy_target: 0.69315
---

# Task T2.4.1: Two-Peak Specific Heat $C_v(T)$ & Fractionalization Entropy Stages ($S = \frac{1}{2}\ln 2$)

## 1. Physical Motivation & Theoretical Context

In conventional quantum magnets, magnetic degrees of freedom freeze out simultaneously at a single phase transition or crossover temperature $T_c \sim \mathcal{O}(J)$, releasing the full spin entropy $S = \ln(2S+1) = \ln 2$ (for $S=1/2$) in a single specific heat anomaly.

In the pure Kitaev quantum spin liquid ($K = -1.0, \mathbf{h} = 0$), the physical spin operators fractionalize into two distinct emergent species with hierarchically separated energy scales:
1. **Localized $\mathbb{Z}_2$ Gauge Fluxes (Visons):**
   Associated with hexagonal plaquette loop operators $\hat{W}_p = \pm 1$. The ground state lies in the vortex-free sector ($\hat{W}_p = +1$). Creating a pair of visons costs a finite flux excitation gap $\Delta_v \approx 0.137\,|K|$ (determined from the long-range asymptotic limit of the two-vortex potential $V(\infty) = 2\Delta_v$ in Task T2.2.1). Thermal proliferation of visons occurs at a low crossover scale:
   $$T_L \sim 0.4 \Delta_v \approx 0.05 - 0.06\,|K|.$$
   Because there is 1 independent flux per hexagonal plaquette and $N_p = N_{\mathrm{spins}}/2$ plaquettes, freezing of the vison gauge field releases exactly half of the total spin entropy:
   $$\Delta S_{\mathrm{vison}} = \frac{N_p}{N_{\mathrm{spins}}} \ln 2 = \frac{1}{2}\ln 2 \approx 0.34657 \quad (\text{per spin}).$$

2. **Itinerant Matter Majorana Fermions ($c$-Majoranas):**
   Disperse across the honeycomb lattice with a gapless Dirac cone and bandwidth $W \approx 1.5 - 3.0\,|K|$. Thermal excitation of the itinerant Majorana Fermi sea occurs at a much higher characteristic scale:
   $$T_H \sim 0.3 - 0.6\,|K|.$$
   There is 1 $c$-Majorana per site, corresponding to $N/2$ complex fermions for $N$ spins. Thermal excitation across this band releases the remaining half of the spin entropy:
   $$\Delta S_{\mathrm{Majorana}} = \frac{N/2}{N} \ln 2 = \frac{1}{2}\ln 2 \approx 0.34657 \quad (\text{per spin}).$$

Together, the thermal specific heat exhibits the hallmark **two-peak signature** (Nasu, Udagawa, Motome, PRL 113, 197205 (2014)):
$$C_v(T) = C_v^{\mathrm{vison}}(T) + C_v^{\mathrm{Majorana}}(T),$$
characterized by an intermediate temperature window $T_L \ll T \ll T_H$ where the entropy plateaus at:
$$S(T) \approx \frac{1}{2}\ln 2,$$
and an asymptotic limit $\lim_{T \to \infty} S(T) = \ln 2$.

---

## 2. Pre-Registered Hypotheses & Falsification Criteria

- **Model $\mathcal{M}_0$ (Null Hypothesis — Un-fractionalized Single-Scale Magnet):**
  A single Schottky peak model governed by one excitation gap $\Delta$:
  $$C_v^{(0)}(T) = A \cdot \left(\frac{\Delta}{T}\right)^2 \frac{e^{\Delta/T}}{\left(1 + e^{\Delta/T}\right)^2}$$
  with $k_0 = 2$ free parameters ($A, \Delta$).

- **Model $\mathcal{M}_1$ (Candidate Hypothesis — Fractionalized Two-Peak QSL):**
  A two-component model with separated vison and Majorana scales:
  $$C_v^{(1)}(T) = A_1 \left(\frac{\Delta_1}{T}\right)^2 \frac{e^{\Delta_1/T}}{(1 + e^{\Delta_1/T})^2} + A_2 \left(\frac{\Delta_2}{T}\right)^2 \frac{e^{\Delta_2/T}}{(1 + e^{\Delta_2/T})^2}$$
  with $k_1 = 4$ free parameters ($A_1, \Delta_1, A_2, \Delta_2$).

- **Popperian Falsification Criteria:**
  1. **Peak Separation:** Peak temperature ratio $T_H / T_L \ge 5.0$.
  2. **Intermediate Entropy Plateau:** At the specific heat minimum $T_{\mathrm{mid}} = \mathrm{argmin}_{T \in (T_L, T_H)} C_v(T)$:
     $$\frac{|S(T_{\mathrm{mid}}) - \frac{1}{2}\ln 2|}{\frac{1}{2}\ln 2} \le 5.0\%.$$
  3. **High-Temperature Asymptotic Entropy:** At $T = 6.0\,|K|$:
     $$\frac{|S(T_{\mathrm{max}}) - \ln 2|}{\ln 2} \le 2.0\%.$$
  4. **Bayesian Model Comparison:**
     $$\Delta\mathrm{BIC} = \mathrm{BIC}_0 - \mathrm{BIC}_1 > 10.0 \quad \text{and} \quad B_{10} = \exp\left(\frac{1}{2}\Delta\mathrm{BIC}\right) > 148.0.$$

---

## 3. Empirical Benchmark & Bayesian Model Comparison Results

Executing `scripts/phase_2/m2_4/benchmark_specific_heat.py` on a $32 \times 32$ momentum grid with $\Delta_v = 0.137\,|K|$ yields:

### 3.1 Characteristic Peak Temperatures & Extrema
- **Low-Temperature Vison Peak:** $T_L = 0.0588\,|K|$, with peak specific heat $C_v(T_L) = 0.2229$.
- **High-Temperature Majorana Peak:** $T_H = 0.6090\,|K|$, with peak specific heat $C_v(T_H) = 0.1711$.
- **Peak Temperature Ratio:**
  $$\frac{T_H}{T_L} = \frac{0.6090}{0.0588} = 10.36 \ge 5.0 \quad (\text{PASS}).$$
- **Intermediate Specific Heat Minimum:** $T_{\mathrm{mid}} = 0.1892\,|K|$, with $C_v(T_{\mathrm{mid}}) = 0.1022$.

### 3.2 Entropy Stages & Plateau Exactness
- **Intermediate Plateau Entropy:**
  $$S(T_{\mathrm{mid}}) = 0.33646 \quad \text{vs} \quad \frac{1}{2}\ln 2 = 0.34657 \quad (\Delta_{\mathrm{rel}} = 2.92\% \le 5.0\% \text{ PASS}).$$
- **High-Temperature Entropy:**
  $$S(T_{\mathrm{max}} = 6.31\,|K|) = 0.68848 \quad \text{vs} \quad \ln 2 = 0.69315 \quad (\Delta_{\mathrm{rel}} = 0.67\% \le 2.0\% \text{ PASS}).$$

### 3.3 Bayesian Model Selection (Two-Peak vs Single-Peak)
- **Residual Sum of Squares:**
  $$\mathrm{RSS}_0 = 1.9496 \times 10^{-1}, \quad \mathrm{RSS}_1 = 6.2511 \times 10^{-3} \quad (31.2\times \text{ variance reduction}).$$
- **Bayesian Information Criteria:**
  $$\mathrm{BIC}_0 = -205.51, \quad \mathrm{BIC}_1 = -335.73$$
  $$\Delta\mathrm{BIC}_{01} = +130.22 \gg 10.0, \quad B_{10} = 1.90 \times 10^{28} \gg 148.0.$$
- **Verdict:** Decisive, overwhelming empirical refutation of the un-fractionalized single-scale hypothesis.

---

## 4. Key Scientific Discoveries

1. **Hierarchy of Fractionalized Excitations:**
   The Kitaev quantum spin liquid exhibits an order-of-magnitude separation of thermal scales ($T_H / T_L = 10.36$). Vison excitations freeze out at $T_L \approx 0.06\,|K|$, while Majorana fermions remain itinerant and active up to $T_H \approx 0.61\,|K|$.
2. **Quantized Intermediate Entropy Plateau:**
   The thermodynamic entropy displays a pronounced plateau at $S(T) \approx \frac{1}{2}\ln 2$ around $T \approx 0.19\,|K|$, providing direct thermodynamic proof of $S=1/2$ spin fractionalization into discrete topological sectors.
3. **Overwhelming Bayesian Discrimination:**
   Bayesian odds $B_{10} = 1.90 \times 10^{28} \gg 148.0$ decisively falsify any single-peak description, establishing the two-peak thermal profile as a defining signature of fractionalization in Kitaev materials.
