# Author Response and Revision Record (Round 1 -> Round 2)

We thank Referee 1, Referee 2, and the Handling Editor for their constructive, rigorous, and searching reviews of manuscript `PHYS-2026-MF-0412`. We have incorporated all four mandatory directives:

## 1. Response to Directive A: Epistemic Demarcation
- **Theorem 1 Demotion:** Demoted from "Theorem" to `Proposition 1: First-Order Critical Boundary and Hysteresis`, explicitly framed as the result of forward and backward adiabatic numerical continuation on finite meshes.
- **Lagrange Multiplier Demotion:** Demoted from "Theorem" to `Proposition 3: Gauge Susceptibility and Linear Counter-Torque Response`. The linear response relation is derived analytically via the Implicit Function Theorem on $\mathbf{g}_s = \mathbf{0}$, while the numerical relation $\|\boldsymbol{\lambda}^*\| \approx \sqrt{5}h$ is explicitly demarcated as an empirical regression (<0.8% error).
- **Proof vs Regression Clarification:** In Section 3.2 and Section 1, we made the strict separation between the algebraic symmetry theorem (which proves that $C_3$ and $\mathcal{T}$ selection rules enforce an $h^3$ exponent) and the numerical prefactor $2.592$, which is an empirical regression on finite grids.
- **Pearson Walk Clarification:** Section 2.3 now cleanly separates Borwein et al.'s continuous thermodynamic 3-step random walk integral $\bar{U} = W_3(1)/3 \approx 0.5248657$ from the discrete finite-$N_k$ mesh residual quadrature scaling $\epsilon_{\mathrm{gauge}} \approx \sqrt{3\pi / (8N_k)} \propto N_k^{-1/2}$.
- **Fock vs. Hilbert Space:** Section 2.1 explicitly states that mean-field states are unprojected quasi-free states in the $4^N$-dimensional Majorana Fock space $\mathcal{F}$ satisfying the local gauge constraint only on average ($\langle \mathbf{G}_i \rangle = \mathbf{0}$).

## 2. Response to Directive B: Algebraic and Operator Normalization
- **Spin/Gauge Bilinear Normalization:** Eqs. (13)–(14) have been rescaled with exact prefactors $S_i^\alpha = -\frac{\ii}{8}\mathbf{c}_i^T M^\alpha \mathbf{c}_i$ and $G_i^\gamma = -\frac{\ii}{8}\mathbf{c}_i^T G^\gamma \mathbf{c}_i$, ensuring that the spin algebra satisfies $[S^\alpha, S^\beta] = \ii \epsilon^{\alpha\beta\gamma} S^\gamma$ with Casimir $\mathbf{S}^2 = 3/4$.
- **Correction of Parity Projector Sign:** Eq. (21) now defines $\mathcal{D}_i = + c_i^0 c_i^1 c_i^2 c_i^3 = 2\mathbf{S}_i^2 - 2\mathbf{G}_i^2 = +1$, correctly projecting onto the physical spin-$1/2$ odd-parity doublet ($n_1+n_2=1$) and excluding the even-parity unphysical gauge doublet.
- **Cross-Susceptibility Derivation:** Proposition 3 now provides the exact matrix derivative $\left.\frac{\mathrm{d}\boldsymbol{\lambda}^*}{\mathrm{d}\mathbf{h}}\right|_{\mathbf{h}=\mathbf{0}} = -\boldsymbol{\chi}_{\mathrm{gauge}}^{-1} \boldsymbol{\chi}_{G,S}$, explicitly incorporating the driving cross-susceptibility $\boldsymbol{\chi}_{G,S} = \frac{\partial \langle \mathbf{g}_s \rangle}{\partial \mathbf{h}}$.

## 3. Response to Directive C: Perturbation Theory & Decoupling Mechanisms
- **Vison Gap vs Fock Gap:** Remark 1 in Section 3.2 now provides a detailed physical discussion of the prefactor ratio $\Delta_{\mathrm{pert}} / \Delta_K^{\mathrm{MF}} \approx 11.4$. We explain that Kitaev's exact spin perturbation creates virtual pairs of static $\mathbb{Z}_2$ visons with threshold $2\Delta_{\mathrm{flux}} \approx 0.26|K|$ giving energy denominators $(2\Delta_{\mathrm{flux}})^2 \approx 0.068|K|^2$, whereas Majorana mean-field theory decouples the 3-spin operator into bilinear fermions with a single-particle gap set by the flat gauge bands ($\sim |K|$), explaining the suppression of the mean-field prefactor.

## 4. Response to Directive D: Computational Complexity, Roofline, and Resonances
- **Table 1 FLOP Overhauls:**
  - Good Broyden overhead corrected to exact $8 D^2 \text{ FLOPs}$.
  - Joint DIIS overhead corrected to $4 m D + \mathcal{O}(m^3) \text{ FLOPs}$ (strictly linear in $D$).
  - JFNK separates Krylov overhead ($2 k_{\mathrm{kry}} D$) from map evaluation cost ($k_{\mathrm{kry}} \Phi$).
  - Added explicit 95% Wilson confidence intervals $[w^-, w^+]$ to Table 1.
- **Track A Real Eigensolve FLOP Fallacy:** Section 4.2 explicitly refutes the "4.0x FLOP reduction" claim, noting that dense matrix multiplications to assemble $\mathbf{S} = -\mathbf{A}^2$ and $-\mathbf{A}\mathbf{S}^{-1/2}$ bring the total cost to $\frac{26}{3}M^3 \approx 8.67M^3 \text{ FLOPs}$ (higher than complex `zheevd` at $\frac{16}{3}M^3 \approx 5.33M^3$). Real runtime improvements are attributed to Level-3 BLAS SIMD cache locality.
- **Dirac Resonance at $3 \mid L$:** Section 4.2 documents the Dirac resonance where condition number explodes to $\kappa(\mathbf{A}) \approx 6 \times 10^{15}$. In Track A, $\kappa(\mathbf{S}) \approx 3.6 \times 10^{31}$ causes catastrophic mantissa loss, spiking Frobenius error to $9.5 \times 10^{-9}$, while Track B purification is immune ($\sim 4 \times 10^{-15}$).
- **Newton-Schulz Roofline & Deflation:** Section 4.3 clarifies that purification consumes $48 M^3 \text{ FLOPs}$ ($9\times$ more than eigensolvers), and its speedup is a hardware roofline effect ($114.6\text{ GFLOP/s}$ cache-blocked `dgemm` vs $20\text{ GFLOP/s}$ memory-bound eigensolvers). The secular deflation stalling of `zheevd` under $C_{6v}$ symmetry multiplets is explicitly explained.
- **Table Formatting:** Both tables resized using `\resizebox{\columnwidth}{!}{...}`, eliminating all overfull horizontal boxes (0 errors, 0 warnings).
