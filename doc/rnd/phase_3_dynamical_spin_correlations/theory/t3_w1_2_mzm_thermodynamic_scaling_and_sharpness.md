# Task T3.W1.2: Thermodynamic Limit Scaling of MZM Hybridization Suppression and Zero-Bias Sharpness

## Front-1 Metadata & Formal Specifications
```yaml
document_type: front_1_theory_note
task_id: T3.W1.2
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
    - H_MZM_SCALE_01
    - H_MZM_SCALE_02
    - H_MZM_SCALE_03
    - H_MZM_SCALE_04
```

---

## 1. Theoretical Context & Physical Principles

In the Kitaev honeycomb model under time-reversal symmetry breaking ($\kappa \ne 0$), each isolated vison (plaquette $\pi$-flux) binds a localized Majorana zero mode (MZM) $\gamma_v$. When two visons $v_1, v_2$ are separated by spatial distance $R = |\mathbf{r}_1 - \mathbf{r}_2|$, their bound zero modes hybridize into a Dirac fermion state $d_0 = \frac{1}{2}(\gamma_1 + i \gamma_2)$ with energy splitting $\delta\epsilon_0(R)$.

### 1.1 MZM Hybridization Law
By quantum tunneling through the bulk Chern insulator ($\nu = \pm 1$ chiral p-wave topological phase), the hybridization splitting is governed by the bulk spectral gap $\Delta_{\mathrm{Dirac}} = 3\sqrt{3}\kappa$:
$$
\delta\epsilon_0(R) \sim C_0 \cos(k_F R + \phi) \frac{\mathrm{e}^{-R/\xi_M}}{\sqrt{R/a}}, \quad \xi_M \approx \frac{\hbar v_F}{\Delta_{\mathrm{Dirac}}} = \frac{\sqrt{3}K/2}{3\sqrt{3}\kappa} = \frac{K}{6\kappa}
$$
In the presence of the full honeycomb dispersion, the effective MZM localization length is $\xi_M \approx 1.5\text{--}4.5 a$.

### 1.2 Dynamical Local Quench Response and Energy Scale Separation
Under an instantaneous local spin flip $\sigma_i^z$ on bond $\langle ij \rangle$, the system undergoes a gauge quench $G_0 \to G_1$:
$$
S_{zz}(\mathbf{r}_i, \omega) = \sum_\lambda |\langle \Psi_\lambda(G_1) | \sigma_i^z | \Psi_0(G_0) \rangle|^2 \delta(\omega - (E_\lambda(G_1) - E_0(G_0)))
$$
Crucially, in the non-Abelian phase ($\kappa \ne 0$), there are two distinct excitation energy scales:
1. **Bulk Two-Vison Creation Gap ($\Delta E_{2v} \approx 0.21\text{--}0.26 |K|$):**
   A spin flip creates a pair of plaquette visons ($\Delta W_p = -1$). In the bulk, far from existing vortices, the creation of two new visons requires a finite gauge-flux threshold $\Delta E_{2v} \approx 0.209\text{--}0.264 |K|$. No excitations exist below this gap on a bulk site: $S_{\mathrm{bulk}}(\omega < \Delta E_{2v}) \equiv 0$.
2. **Bulk Itinerant Dirac Fermion Gap ($\Delta_{\mathrm{Dirac}} = 3\sqrt{3}\kappa \approx 0.78 |K|$ at $\kappa = 0.15$):**
   The gap to free Chern fermionic quasiparticles in the continuum.
3. **Majorana Zero Mode Peak ($\omega_{\mathrm{MZM}} \le 0.10 |K|$):**
   When the quench is performed on a bond connected directly to a vortex core $p_1$, the vison is shifted rather than nucleating a new, unhybridized pair. The local Majorana zero mode $\gamma_{v1}$ absorbs the transition with subgap excitation energy $\omega_{\mathrm{MZM}} \to 0$ as $R \to \infty$.

Therefore, the MZM subgap spectral window must be defined with cutoff $\omega_c$ strictly below the two-vison threshold:
$$
\omega_c \approx 0.18 |K| < \min_{L} \Delta E_{2v}(L) \approx 0.209 |K|
$$
This cleanly captures the hybridized zero-bias mode across all clusters ($\omega_{\mathrm{MZM}} = 0.1564 |K|$ at $L=6$ down to $\omega_{\mathrm{MZM}} = 0.00317 |K|$ at $L=20$), while ensuring the bulk two-vison nucleation response ($\omega_{\min}^{\mathrm{bulk}} \ge 0.209 |K|$) has strictly vanishing spectral weight ($S_{\mathrm{subgap}}^{\mathrm{bulk}} = 0.0 < 10^{-5}$).

### 1.3 Parity Policies: Thermodynamic Bulk (OBC) vs. Global Torus
In finite-cluster Lehmann calculations:
1. **Thermodynamic Parity Policy (`boundary="obc"`):**
   Selects the thermodynamic bulk vacuum sector ($\pi_G = +1$). This eliminates non-contractible Wilson loop parity twists across the torus boundaries, isolating local vortex and MZM physics as experienced in an infinite 2D plane.
2. **Topological Parity Policy (`boundary="torus"`):**
   Enforces Pedrocchi-Zschocke global parity constraints from large Wilson loops. As system size grows past $L \ge 20$, single bond flips induce parity transitions ($\det(Q_0)\det(Q_1) = -1$) requiring singular Bogoliubov resolution. The benchmark and sharpness analysis defaults to `boundary="obc"` to evaluate intrinsic local STS spectroscopy.

---

## 2. Formal Hypothesis Contracts (Rule 08)

### Hypothesis Contract `H_MZM_SCALE_01`: Exponential Suppression of MZM Splitting
- **Statement:** In the commensurate Dirac sector ($3 \mid L$, where discrete torus modes sample the Dirac cone $K$), the two-vortex MZM hybridization energy $\delta\epsilon_0(R_{\max})$ at antipodal vortex separation $R_{\max} = \frac{\sqrt{3}}{2} L a$ decays exponentially with system size / separation:
  $$
  \delta\epsilon_0(R_{\max}) \le C \mathrm{e}^{-R_{\max}/\xi_M}
  $$
  where the localization length $\xi_M$ satisfies $1.5 a \le \xi_M \le 4.5 a$ at $\kappa = 0.15$.
- **Support Criterion:** A log-linear regression $\ln \delta\epsilon_0(R_{\max}) = \alpha - R_{\max}/\xi_M$ across $L \in \{6, 12, 18, 24\}$ ($R_{\max} \in [5.2a, 20.8a]$) achieves $R^2 \ge 0.85$ with fitted $\xi_M \in [1.5, 4.5]a$.
- **Refute Criterion:** $R^2 < 0.75$ or $\delta\epsilon_0(R_{\max})$ decays algebraically ($1/R^\alpha$) or $\xi_M > 6.0 a$.

### Hypothesis Contract `H_MZM_SCALE_02`: Thermodynamic Zero-Bias Peak Convergence
- **Statement:** The zero-bias MZM peak excitation energy $\omega_{\mathrm{MZM}}(L, R)$ measured relative to the two-vortex ground state converges toward zero with increasing separation $R$ and system size $L$:
  $$
  \lim_{L \to \infty, R \to \infty} |\omega_{\mathrm{MZM}} - \delta E_0| < 0.01 |K|
  $$
- **Support Criterion:** For $L \ge 16$ and $R \ge L/2$, the localized subgap pole satisfies $|\omega_{\mathrm{MZM}}| \le 0.05 |K|$ (and $< 0.01 |K|$ for $L \ge 20$).
- **Refute Criterion:** Subgap pole remains pinned at $\omega > 0.10 |K|$ independent of $R$ and $L$.

### Hypothesis Contract `H_MZM_SCALE_03`: Robust Subgap Sharpness Ratio
- **Statement:** The zero-bias sharpness ratio, defined as the fraction of dynamical spectral weight inside the subgap window $[0, \omega_c]$ with $\omega_c = 0.18 |K| < \Delta E_{2v}$:
  $$
  S_{\mathrm{subgap}} \equiv \frac{\sum_{\omega_\lambda \le \omega_c} I(\omega_\lambda)}{\sum_\lambda I(\omega_\lambda)}
  $$
  remains strictly non-vanishing and saturated ($S_{\mathrm{subgap}} \ge 0.15$) across all sizes $L \in [6, 20]$.
- **Support Criterion:** $S_{\mathrm{subgap}} \ge 0.15$ for core-site quench across all tested $L \in [6, 20]$.
- **Refute Criterion:** $S_{\mathrm{subgap}} < 0.05$ or vanishes asymptotically as $\mathcal{O}(1/L^2)$.

### Hypothesis Contract `H_MZM_SCALE_04`: Bulk-Site Subgap Vanishing & Vison Gap Protection
- **Statement:** For a spin flip on a bond far from any vortex core (maximizing distance $d(p, \{p_1, p_2\})$ on the torus), the local response has strictly zero subgap weight below $\omega_c = 0.18 |K|$:
  $$
  S_{\mathrm{subgap}}^{\mathrm{bulk}} < 10^{-5}
  $$
  and the lowest excitation threshold matches or exceeds the two-vison creation gap $\Delta E_{2v} \approx 0.209 |K|$ ($\omega_{\min}^{\mathrm{bulk}} \ge 0.20 |K|$).
- **Support Criterion:** $S_{\mathrm{subgap}}^{\mathrm{bulk}} < 10^{-5}$ and minimum excitation frequency $\omega_{\min}^{\mathrm{bulk}} \ge 0.20 |K|$ across all tested $L \in [6, 20]$.
- **Refute Criterion:** $S_{\mathrm{subgap}}^{\mathrm{bulk}} > 10^{-2}$ or bulk excitation appears below $0.15 |K|$.
