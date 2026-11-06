# SciPost Physics Handling Editor: Final Acceptance Decision

**Manuscript Reference:** `PHYS-2026-MF-0412`
**Manuscript Title:** Majorana Parton Mean-Field Theory of the Extended Kitaev Honeycomb Model: Algebraic Symmetries, Grassmannian Solvers, and Kinematic Invariants
**Author:** Lu Freitas
**Editorial Decision:** **ACCEPT FOR PUBLICATION**

---

### Handling Editor Synthesis & Justification

Following two rigorous, uncompromising rounds of adversarial peer review by Referee 1 (Senior Mathematical Physicist) and Referee 2 (Expert Computational Physicist), all editorial directives and technical requirements have been comprehensively resolved:

1. **Epistemic Rigor & Provenance:**
   - The distinction between exact algebraic proofs (such as the cubic selection rule $\Delta_K \propto h^3$ and flat-band gauge neutrality $\boldsymbol{\lambda}^* \equiv \mathbf{0}$) and empirical regressions (the prefactor $2.592$ and the linear relation $\|\boldsymbol{\lambda}^*\| \approx \sqrt{5}h$) is rigorously enforced across the introduction, main text, and conclusions.
   - Numerical continuation findings are properly classified as Propositions rather than Theorems.
   - The Pearson random walk continuous integral is cleanly decoupled from discrete finite-mesh $\mathcal{O}(N_k^{-1/2})$ quadrature noise.

2. **Algebraic Consistency:**
   - The parity projector $\mathcal{D}_i = + c_i^0 c_i^1 c_i^2 c_i^3 = \frac{4}{3}(\mathbf{S}_i^2 - \mathbf{G}_i^2) = +1$ correctly selects the physical spin-$1/2$ doublet and eliminates the unphysical gauge sector.
   - The spin and gauge bilinear normalizations accurately satisfy the canonical Lie algebra $[S^\alpha, S^\beta] = \ii \epsilon^{\alpha\beta\gamma} S^\gamma$ with quadratic Casimir $\mathbf{S}^2 = 3/4$.

3. **Physical Perturbation Mechanisms:**
   - The physical origin of the prefactor discrepancy between exact spin-space perturbation theory ($\Delta_{\mathrm{pert}} \approx 29.5 h^3$) and mean-field theory ($\Delta_K^{\mathrm{MF}} \approx 2.592 h^3$) is clearly explained by the mismatch in energy denominators: non-local two-vison threshold $(2\Delta_{\mathrm{flux}})^2 \approx 0.068 |K|^2$ versus local single-particle gauge Majorana bandwidth $\sim |K|^2$.

4. **Computational Complexity & Roofline Analysis:**
   - All FLOP formulas for Full Broyden ($8D^2$), Joint DIIS ($4mD + \mathcal{O}(m^3)$), and JFNK are mathematically exact.
   - The fallacy of "4.0x FLOP reduction" in real-space eigensolves has been fully corrected by accounting for the $\mathcal{O}(M^3)$ matrix multiplications in $S = -A^2$ and $A S^{-1/2}$.
   - Matrix sign purification is accurately presented as a hardware roofline victory ($114.64\text{ GFLOP/s}$ sustained throughput in cache-blocked `dgemm`), and the catastrophic mantissa loss under $3 \mid L$ Dirac resonances in eigensolvers is thoroughly documented.

5. **Editorial Formatting:**
   - Compiled PDF `scipost_kitaev_mf.pdf` has **0 errors, 0 warnings, and 0 overfull horizontal boxes**.
   - Tables and mathematical typography conform to the highest SciPost Physics standards.

The paper is hereby **accepted for publication** in SciPost Physics.
