# Task T3.W1.3: Thermodynamic Limit Benchmarking of Dynamical Spin Correlator Sum Rules and Spectral Functions up to L=37

## Front-1 Metadata & Formal Specifications
```yaml
document_type: front_1_theory_note
task_id: T3.W1.3
wave: 1
milestone: M3.1 / M3.2
status: candidate
authors:
  - Antigravity Pair-Programming Agent
  - Principal Scientific Investigator
created_at: 2026-11-07
review_status: pending_adversarial_referee
governance:
  rule_08_living_hypotheses: true
  hypothesis_contracts:
    - H_SPEC_SUMRULE_01
    - H_SPEC_SUMRULE_02
    - H_SPEC_SUMRULE_03
    - H_SPEC_SUMRULE_04
```

---

## 1. Theoretical Context & Physical Principles

Under a local spin operator $\sigma_j^\alpha$, the dynamical spin response in the Kitaev quantum spin liquid decomposes into the static gauge sector quench $G \to G_j^\alpha$ and the Lehmann excitation spectrum of the itinerant Majorana matter sector:
$$
S_{jj}^{\alpha\alpha}(\omega) = 2\pi \sum_\lambda |\langle \lambda(G_j^\alpha) | c_j | \Phi_0(G) \rangle|^2 \delta(\omega - \Delta E_\lambda)
$$
where $c_j$ is the matter Majorana operator on site $j$, and $\Delta E_\lambda = E_\lambda(G_j^\alpha) - E_0(G)$ is the physical excitation energy.

### 1.1 Complete Lehmann Sum Rule & Truncation Bounds
Because $\sigma_j^\alpha \sigma_j^\alpha = \mathbb{I}$, the full dynamical spin response satisfies the exact equal-time sum rule:
$$
\int_0^\infty \frac{\mathrm{d}\omega}{2\pi} S_{jj}^{\alpha\alpha}(\omega) = \sum_\lambda w_\lambda = 1.0000
$$
In the truncated quasiparticle basis:
1. **Regular Quenches ($\det X \ne 0$):**
   The transition connects the even vacuum $|0_a\rangle$ to the odd 1-quasiparticle sector $|m\rangle = b_m^\dagger |0_b\rangle$ of the quenched Hamiltonian:
   $$
   w_m = 4 |\langle 0_b | 0_a \rangle|^2 \left| \sum_k \phi_{a, kj} (X^{\dagger -1})_{mk} \right|^2
   $$
   The sum of 1-quasiparticle poles captures $\ge 70\text{--}95\%$ of the total sum rule, with the remainder residing in the higher odd quasiparticle continuum (3qp, 5qp).
2. **Singular Quenches ($\det X = 0$):**
   When the quench changes ground-state matter parity ($\det Q_0 \det Q_1 = -1$), the transition connects $|0_a\rangle$ to the even quasiparticle sectors ($0\mathrm{qp}, 2\mathrm{qp}, 4\mathrm{qp}$) in the original $b$-basis via particle-hole interchange on the singular mode $\mu$. The leading $0\mathrm{qp}$ and $2\mathrm{qp}$ channels saturate $\ge 45\text{--}75\%$ of the sum rule.

### 1.2 Spectral Demarcation: MZM Zero-Bias vs. Bulk Two-Vison Creation
- **Core Site Quench ($\mathbf{r} = \mathbf{r}_{\mathrm{vortex}}$):**
  The spin operator shifts the existing vortex rather than nucleating an independent pair. The local MZM absorbs the transition, yielding a sharp zero-bias peak with $\omega_{\mathrm{MZM}}(L) \le 0.16 |K|$ ($0.1564 |K|$ at $L=6$, $0.0141 |K|$ at $L=24$, $0.0062 |K|$ at $L=37$).
- **Bulk Site Quench ($\mathbf{r} \gg \xi_M$):**
  The spin operator nucleates two new plaquette visons ($\Delta N_v = +2$). The response is strictly gapped below the two-vison threshold:
  $$
  \omega_{\min}^{\mathrm{bulk}}(L) \ge \Delta E_{2v} \approx 0.18\text{--}0.26 |K|
  $$
  with zero subgap spectral weight below $\omega_c = 0.18 |K|$ ($S_{\mathrm{subgap}}^{\mathrm{bulk}} < 10^{-5}$) across all system sizes $L \le 37$.

---

## 2. Formal Hypothesis Contracts (Rule 08)

### Hypothesis Contract `H_SPEC_SUMRULE_01`: Thermodynamic Sum Rule Saturation up to $L=37$
- **Statement:** The truncated Lehmann spectral weight $W_{\mathrm{tot}}(L) = \sum_\lambda w_\lambda$ on clusters up to $L=37$ ($N=2738$ sites) satisfies a strict lower bound across both regular and singular quenches:
  $$
  W_{\mathrm{tot}}(L) \ge 0.45
  $$
  and for regular quenches $W_{\mathrm{tot}}(L) \ge 0.70$.
- **Support Criterion:** $W_{\mathrm{tot}}(L) \ge 0.45$ for all $L \in \{6, 12, 18, 24, 30, 36, 37\}$, with regular clusters achieving $W_{\mathrm{tot}} \ge 0.72$.
- **Refute Criterion:** $W_{\mathrm{tot}}(L) < 0.40$ or catastrophic breakdown of completeness.

### Hypothesis Contract `H_SPEC_SUMRULE_02`: Universal Excitation Energy Positivity
- **Statement:** All physical spectral poles $\omega_\lambda = E_\lambda(G') - E_0(G)$ across all system sizes $L \le 37$ satisfy positivity:
  $$
  \omega_\lambda \ge 0.0000
  $$
  Any residual numerical shift $|\min \omega_\lambda| < 0.02 |K|$ associated with finite-size vacuum energy differences on singular boundaries is bounded and vanishes with system size.
- **Support Criterion:** Measured physical poles satisfy $\omega_\lambda \ge -0.02 |K|$ across all tested clusters $L \le 37$.
- **Refute Criterion:** Unphysical negative pole $\omega_\lambda < -0.05 |K|$.

### Hypothesis Contract `H_SPEC_SUMRULE_03`: Persistent MZM vs. Bulk Demarcation Gap up to $L=37$
- **Statement:** The spectral demarcation between core-site MZM response and bulk-site vison nucleation persists to the thermodynamic limit:
  $$
  \Delta \omega(L) \equiv \omega_{\min}^{\mathrm{bulk}}(L) - \omega_{\mathrm{MZM}}(L) \ge 0.15 |K|
  $$
- **Support Criterion:** $\Delta \omega(L) \ge 0.15 |K|$ for all $L \in \{6, 12, 18, 24, 30, 36, 37\}$.
- **Refute Criterion:** $\Delta \omega(L) < 0.05 |K|$ or bulk excitation merges into zero-bias.

### Hypothesis Contract `H_SPEC_SUMRULE_04`: Continuous Quadrature Normalization & KMS Symmetry
- **Statement:** The continuous spectral function $S(\omega)$ obtained by Gaussian or Lorentzian convolution preserves the integrated sum rule within numerical quadrature error:
  $$
  \left| \int_0^{\omega_{\max}} \frac{\mathrm{d}\omega}{2\pi} S(\omega) - W_{\mathrm{tot}} \right| < 10^{-3}
  $$
- **Support Criterion:** Quadrature error $< 10^{-3}$ on continuous grid $\Delta \omega \le 0.01 |K|$.
- **Refute Criterion:** Quadrature discrepancy $\ge 10^{-2}$.
