# SciPost Physics Editorial Decision Letter

**Manuscript Reference:** `PHYS-2026-MF-0412`
**Manuscript Title:** Majorana Parton Mean-Field Theory of the Extended Kitaev Honeycomb Model: Algebraic Symmetries, Grassmannian Solvers, and Kinematic Invariants
**Author:** Lu Freitas
**Handling Editor:** SciPost Physics Editorial College
**Date:** November 6, 2026
**Editorial Decision:** **Major Revision**

---

Dear Dr. Freitas,

The peer review process for your manuscript submitted to *SciPost Physics* is now complete. We have received detailed, independent evaluations from two expert referees:
- **Referee 1** (Senior Mathematical Physicist and Quantum Many-Body Theorist), who evaluated the operator algebras, physical Hilbert space projections, symmetry arguments, and analytical perturbation theory.
- **Referee 2** (Expert Computational Physicist and Numerical Analyst), who audited the algorithmic complexity, fixed-point accelerators, FLOP counters, hardware roofline mechanics, and finite-lattice Dirac resonances.

Both referees recognize the ambition and potential value of your work. In particular, they commend the geometric formulation of the self-consistent Majorana mean-field equations on the orthogonal Grassmannian manifold bundle and the practical demonstration of matrix sign purification for high-throughput real-space cluster calculations.

However, both reports independently conclude that the manuscript in its current form suffers from **critical mathematical errors, misleading nomenclature, severe arithmetic undercounting, and the omission of fundamental numerical failure modes**. The mathematical foundation of the operator representation is compromised by a fatal sign error in the local parity projector and factor-of-two discrepancies in the spin algebra. Furthermore, empirical regressions are erroneously designated as mathematical theorems; FLOP counts in the numerical section suppress leading-order matrix multiplications by factors between 1.6× and 4×; and the discussion of matrix sign purification obscures the true hardware roofline mechanism while ignoring catastrophic mantissa loss at the $L \equiv 0 \pmod 3$ Dirac point resonance.

In accordance with *SciPost Physics* criteria—which mandate absolute scientific clarity, rigorous self-contained mathematical exposition, complete reproducibility, and zero tolerance for unsubstantiated claims or arithmetic fallacies—the Editorial College has decided that this manuscript requires **Major Revision**.

To guide your revision, the required modifications have been synthesized into **four mandatory thematic directives**. A revised submission will only be considered if each directive is thoroughly and rigorously addressed in both the manuscript and an exhaustive, point-by-point Response to Referees.

---

## Thematic Editorial Directives

### Directive A: Epistemic Demarcation (Theorems vs. Empirical Regressions vs. Conjectures)

*SciPost Physics* enforces strict epistemic demarcation between analytical proofs, exact algebraic statements, empirical numerical regressions, and physical conjectures. The manuscript repeatedly blurs these boundaries.

1. **Demotion of Theorem 1 (Critical Boundaries and Hysteresis):**
   Under Section 3.1, "Theorem 1" lists numerical continuation outputs: $h_c^* / |K| = 0.3325 \pm 0.0005$, spinodal boundaries $h_{s1} \approx 0.28 |K|$ and $h_{s2} \approx 0.40 |K|$, and finite magnetization/energy jumps. Numerical estimates obtained from iterative solvers on finite discretizations accompanied by error bars are **empirical simulation results**, not mathematical theorems. You must demote Theorem 1 to **"Computational Result 1: First-Order Transition and Hysteresis"** (or **"Numerical Observation 1"**). The text must explicitly distinguish exact thermodynamic phase boundaries from finite-grid continuation estimates.
2. **Restructuring of Theorem 3 (Linear Gauge Susceptibility Law):**
   Theorem 3 improperly conflates a formal linear response relation with an empirical regression formula $\|\boldsymbol{\lambda}^*\| \approx \sqrt{5} h$ (bounded by $< 0.8\%$ deviation across $h \le 0.10$). You must separate this section:
   - Present the formal vector linear response relation governing the saddle-point Lagrange multipliers as a **Proposition** derived via the Implicit Function Theorem (incorporating the cross-susceptibility tensor; see Directive B).
   - Designate the proportionality $\|\boldsymbol{\lambda}^*\| \approx \sqrt{5} h$ strictly as an **Empirical Scaling Relation** or numerical fit, clearly stating that $\sqrt{5}$ is an empirical regression parameter unless an exact algebraic derivation is provided.
3. **Retraction of Proof Claims for Numerical Prefactors:**
   In Section 1 (Line 16) and Section 3.2, retract all assertions of having "proven" the cubic prefactor $\Delta_K = (2.592 \pm 0.003) h^3$ or the ratio $\|\boldsymbol{\lambda}^*\| = \sqrt{5} h$. While symmetry selection rules rigorously prove that the leading-order Dirac gap opening must scale as $h^3$ (Theorem 2), the numerical prefactor $2.592 \pm 0.003$ is an empirical fit. Explicitly qualify that symmetry proves the exponent, whereas the prefactor is determined numerically (or derived analytically per Directive C).
4. **Resolution of the Pearson Walk Discretization Law Confusion:**
   In Section 1 and Section 2.3 (Eq. 52), the manuscript claims to have "proven the Pearson Walk Discretization Law ($\mathcal{O}(N_k^{-1/2})$)" governing Brillouin zone quadrature errors. You must completely untangle two distinct concepts:
   - The classical 3-step Pearson planar random walk integral ($W_3(1) \approx 1.574597$), evaluated by Borwein et al. (2012), provides the exact continuous thermodynamic average of the nearest-neighbor bond correlator $\bar{U} \approx 0.5248657$.
   - The residual scaling $\epsilon_{\mathrm{gauge}}(N_k) \propto N_k^{-1/2}$ on regular deterministic momentum meshes is a numerical quadrature property over conical Dirac singularities. You must remove all claims of having "proven" this as an exact mathematical law unless you supply a rigorous analytical proof demonstrating that deterministic midpoint sampling of the gauge residual is mathematically isomorphic to a 2D random flight with Rayleigh-distributed phases. Otherwise, designate it as an empirical quadrature error model.
5. **Physical Hilbert Space vs. Unprojected Fock Space:**
   Explicitly clarify in Section 1 and Section 2.1 the dimensional hierarchy between the physical spin Hilbert space ($\dim \mathcal{H}_{\mathrm{phys}} = 2^N$) and the unprojected 4-flavor Majorana Fock space ($\dim \mathcal{F} = 4^N = 2^{2N}$). Emphasize that quadratic Majorana mean-field theory operates strictly on unprojected quasi-free states in $\mathcal{F}$, enforcing the local gauge constraint $\langle \mathbf{G}_i \rangle = \mathbf{0}$ only on average via Lagrange multipliers, which does not eliminate gauge fluctuations or enforce exact physical projection ($P_G |\Psi_{\mathrm{MF}}\rangle$).

---

### Directive B: Algebraic and Operator Normalization

The algebraic foundation in Section 2 must be made fully consistent, rigorous, and free of sign or normalization errors.

1. **Normalization of Spin and Gauge Bilinears ($\mathfrak{so}(4) \cong \mathfrak{su}(2)_{\mathrm{spin}} \oplus \mathfrak{su}(2)_{\mathrm{gauge}}$):**
   In Section 2.1 (Eqs. 13–14), the bilinear operators are defined as $S_i^\alpha = -\frac{\ii}{4} \mathbf{c}_i^T M^\alpha \mathbf{c}_i$ and $G_i^\gamma = -\frac{\ii}{4} \mathbf{c}_i^T G^\gamma \mathbf{c}_i$, while claiming to satisfy the standard $\mathfrak{su}(2)$ Lie algebra $[S^\alpha, S^\beta] = \ii \epsilon^{\alpha\beta\gamma} S^\gamma$. As proven by Referee 1, this definition evaluates to:
   $$[S_i^\alpha, S_i^\beta] = 2\ii \epsilon^{\alpha\beta\gamma} S_i^\gamma, \qquad \mathbf{S}_i^2 = 3,$$
   which corresponds to **Pauli operators** $\boldsymbol{\sigma}_i$, not spin-$1/2$ operators $\mathbf{S}_i = \frac{1}{2}\boldsymbol{\sigma}_i$. You must systematically resolve this across the entire manuscript:
   - Either normalize the spin-$1/2$ generators as $S_i^\alpha = \frac{\ii}{8} \mathbf{c}_i^T M^\alpha \mathbf{c}_i = \frac{\ii}{4}(c_i^0 c_i^\alpha - \frac{1}{2}\epsilon^{\alpha\beta\gamma} c_i^\beta c_i^\gamma)$,
   - Or explicitly denote them as Pauli operators $\sigma_i^\alpha$ with commutation relations $[\sigma^\alpha, \sigma^\beta] = 2\ii \epsilon^{\alpha\beta\gamma} \sigma^\gamma$, adjusting all Hamiltonian coupling prefactors accordingly.
   - Clean up index notation to strictly separate 4-flavor Majorana indices ($\mu, \nu \in \{0, 1, 2, 3\}$) from 3-space Cartesian/spin indices ($\alpha, \beta, \gamma \in \{1, 2, 3\}$).
2. **Fatal Sign Inversion in the Local Parity Projector $D_i$:**
   In Eq. (21), the physical parity projector is defined as $D_i = - c_i^0 c_i^1 c_i^2 c_i^3 = 1 - 2\mathbf{G}_i^2 = +1$. Under the canonical Dirac mapping $d_1 = \frac{1}{2}(c^0 + \ii c^1)$, $d_2 = \frac{1}{2}(c^2 + \ii c^3)$, the product of four Majoranas is:
   $$c^0 c^1 c^2 c^3 = -(-1)^{n_1 + n_2}.$$
   The physical spin-$1/2$ sector $(1/2, 0)$ is the odd-parity subspace ($n_1 + n_2 = 1$), where $c^0 c^1 c^2 c^3 = +1$. Your definition $D_i = - c_i^0 c_i^1 c_i^2 c_i^3 = +1$ forces $c^0 c^1 c^2 c^3 = -1$, which selects the **even-parity** subspace ($n_1 + n_2 \in \{0, 2\}$) where $\mathbf{S}^2 = 0$ and $\mathbf{G}^2 = 3/4$ (the unphysical gauge doublet). You must:
   - Correct the definition to $\mathcal{D}_i \equiv + c_i^0 c_i^1 c_i^2 c_i^3 = +1$.
   - Explicitly include the Dirac occupation mapping in Section 2.1 to prove that $\mathcal{D}_i = +1$ correctly isolates the physical spin doublet $\mathcal{F}_+ \cong (1/2, 0)$.
3. **Rigorous Formulation of the Gauge Cross-Susceptibility Tensor:**
   In Section 3.3 (Eq. 75), the linear response relation $\boldsymbol{\lambda}^*(\mathbf{h}) = -\boldsymbol{\chi}_{\mathrm{gauge}}^{-1} \mathbf{h}$ silently sets the driving cross-susceptibility tensor to identity. Applying the Implicit Function Theorem to the stationarity constraint $\mathbf{g}_s(\boldsymbol{\lambda}^*(\mathbf{h}), \mathbf{h}) = \langle \mathbf{G}_s \rangle = \mathbf{0}$ yields:
   $$\frac{\mathrm{d}\boldsymbol{\lambda}_s^*}{\mathrm{d}\mathbf{h}} = - \boldsymbol{\chi}_{\mathrm{gauge}}^{-1} \, \boldsymbol{\chi}_{G, S}, \qquad \text{where } \boldsymbol{\chi}_{G, S} \equiv \frac{\partial \langle \mathbf{G}_s \rangle}{\partial \mathbf{h}}.$$
   You must formally derive $\boldsymbol{\chi}_{G, S}$ from the unprojected Majorana band eigenstates and clarify whether the empirical proportionality constant $\sqrt{5}$ originates from the contraction $\boldsymbol{\chi}_{\mathrm{gauge}}^{-1} \boldsymbol{\chi}_{G, S}$ or is an accidental numerical regression.

---

### Directive C: Perturbation Theory & Physical Decoupling Mechanisms

The comparison between mean-field gap opening and Kitaev's exact spin perturbation theory must be elevated from qualitative remarks to quantitative mathematical physics.

1. **Physical Decoupling: Vison Flux Gap vs. Hartree-Fock Fock Gap:**
   Expand Remark 1 in Section 3.2 into a detailed physical and mathematical exposition contrasting the perturbation mechanisms:
   - In Kitaev's exact spin perturbation theory, the ground state is in the vortex-free sector of the physical spin Hilbert space. A single Zeeman operator $\sigma_j^\alpha$ anticommutes with the adjacent plaquette flux operators $W_p$, creating a pair of localized $\mathbb{Z}_2$ vison excitations. The 3rd-order perturbation pathway $\sigma_j^x \sigma_l^z \sigma_k^y$ around a vertex proceeds via intermediate states containing **two static visons** at energy $2\Delta_{\mathrm{flux}} \approx 0.26 |K|$, yielding an energy denominator $(2\Delta_{\mathrm{flux}})^2 \approx 0.068 |K|^2$.
   - In the 4-flavor Majorana mean-field theory, the $\mathbb{Z}_2$ gauge fields are replaced by static bond expectation values $\chi, \eta$. The intermediate virtual states in the single-particle Rayleigh-Schrödinger expansion of $\mathcal{H}_{\mathrm{MF}}(\mathbf{k})$ are the **flat bands of the gauge Majoranas** ($c^1, c^2, c^3$) situated at energy $\pm |K|$. The energy denominator is governed by $|K|^2$, which is larger than the two-vison denominator by a factor of $|K|^2 / (2\Delta_{\mathrm{flux}})^2 \approx 1 / 0.068 \approx 14.7$.
   - Show that this ratio quantitatively explains why the mean-field Dirac gap prefactor ($\Delta_K^{\mathrm{MF}} \approx 2.59 h^3$) is suppressed by an order of magnitude relative to the exact result ($\Delta_{\mathrm{pert}} \approx 29.5 h^3$).
2. **Analytical Single-Particle Perturbation Expansion:**
   Because the mean-field Hamiltonian $\mathcal{H}_{\mathrm{MF}}(\mathbf{k})$ is an explicit $8 \times 8$ matrix in momentum space, you must provide the **exact analytical 3rd-order perturbation expansion** at the Dirac point $\mathbf{K} = (\frac{1}{3}, \frac{2}{3})$. Derive the closed-form algebraic expression for the mean-field coefficient $C_{\mathrm{MF}}$ such that $\Delta_K^{\mathrm{MF}}(h) = C_{\mathrm{MF}} h^3$, and compare it directly with $C_{\mathrm{pert}} = 6\sqrt{3}/\Delta E^2$, demonstrating exact algebraic agreement with your numerical continuation value ($2.592$).

---

### Directive D: Computational Complexity, FLOP Formulas, Roofline, and Dirac Resonances

Section 4 and Table 1 require a complete overhaul to eliminate arithmetic errors, retract fallacious complexity claims, and incorporate omitted numerical failure modes.

1. **Correction of Linear Algebra FLOP Formulas in Table 1:**
   Correct the leading-order arithmetic operation counts for all non-linear accelerators:
   - **Good Broyden:** Update the step overhead from $2 D^2$ to **$8 D^2 \text{ FLOPs}$**. Explicitly detail the operational breakdown: computing $H_k y_k$ ($2 D^2$), $s_k^T H_k$ ($2 D^2$), inner product ($2 D$), rank-1 outer product update ($2 D^2$), and step application $H_k R_k$ ($2 D^2$).
   - **Joint DIIS (Anderson Acceleration):** Update the step overhead from $\mathcal{O}(m^2 D)$ to **$4 m D + \mathcal{O}(m^3) \text{ FLOPs}$**. Explain that updating the Gram matrix requires $m$ vector inner products ($2 m D$), solving the $(m+1) \times (m+1)$ KKT system requires $\mathcal{O}(m^3)$ operations, and forming the new state requires linear vector combination ($2 m D$). Eliminate the erroneous $m^2 D$ term.
   - **JFNK:** In Table 1, strictly separate the internal linear algebra overhead of Modified Gram-Schmidt GMRES ($\mathcal{O}(k_{\mathrm{kry}}^2 D)$) from the external physics operator evaluation calls ($k_{\mathrm{kry}}$ evaluations of $\Phi$).
2. **Retraction and Rectification of the $4.0\times$ Real Eigensolve FLOP Fallacy (Track A):**
   In Section 4.3 (Lines 40–41), the claim of a "theoretical $4.0\times$ FLOP reduction" for real skew-symmetric eigensolves (`dsyevd` vs. `zheevd`) is mathematically false and must be formally retracted.
   - You must publish the full, honest arithmetic complexity of Track A:
     $$\mathbf{S} = -\mathbf{A}^2 \,(2 M^3) + \text{`dsyevd`} \,(\tfrac{8}{3} M^3) + \mathbf{S}^{-1/2} \,(2 M^3) + \Omega = -\mathbf{A}\mathbf{S}^{-1/2} \,(2 M^3) = \frac{26}{3} M^3 \approx \mathbf{8.67 M^3 \text{ FLOPs}}.$$
   - Compare this directly against standard complex Hermitian eigensolves (`zheevd`), which consume $\frac{16}{3} M^3 \approx \mathbf{5.33 M^3 \text{ FLOPs}}$.
   - Explicitly state that Track A requires **$1.63\times$ MORE floating-point operations** than complex diagonalization. Any observed wall-clock speedup must be truthfully attributed to cache locality and single-precision memory streaming, not to an arithmetic FLOP advantage.
3. **Reconciling Matrix Sign Purification: Roofline Saturation vs. Cuppen Deflation Stalling (Track B):**
   - Retract the claim in Section 4.4 and Section 5.1 that cubic Newton-Schulz purification "breaks the $\mathcal{O}(M^3)$ scaling bottleneck." Purification is strictly an $\mathcal{O}(M^3)$ algorithm consuming $4 M^3 \text{ FLOPs}$ per iteration, which totals **$48 M^3 \text{ FLOPs}$** across 12 iterations—nearly **an order of magnitude more arithmetic work** than eigensolvers.
   - Reframe the discussion entirely around the **hardware roofline model**: explain that Newton-Schulz speedup ($63.53\times$) stems from replacing memory-bound, scalar secular equation searches in eigensolvers with compute-bound, cache-blocked BLAS-3 `dgemm` matrix multiplications operating at peak theoretical throughput ($>104\text{ GFLOP/s}$ vs. $20.2\text{ GFLOP/s}$).
   - Address the Cuppen divide-and-conquer deflation stalling pathology: explain that the 44.6-second runtime of `zheevd` at $M=1152$ is driven by thread roundoff preventing deflation tolerance $|d_i - d_j| \le \text{tol}$ on high-symmetry $C_{6v}$ degenerate eigenvalue multiplets, causing `dlaed4` secular equation stagnation, rather than normal cubic scaling.
   - Document small-gap linear slowdown: provide the iteration scaling formula $k_{\mathrm{linear}} \approx 2.466 \ln(\|\mathbf{A}\|_2 / \Delta_{\min})$ for gapless or near-gapless regimes.
   - Document projector incompleteness at exact zero modes: state that for $\Delta = 0$, $\Omega^2 = -\mathbf{I} + \mathbf{P}_{\mathrm{zero}}$, and specify the skew Pfaffian deflation procedure required to fix physical fermion parity $(-1)^F = \pm 1$.
4. **Incorporation of the $L \equiv 0 \pmod 3$ Dirac Resonance and Mantissa Precision Degradation:**
   Add a dedicated subsection and include the laboratory benchmark data (Figure F18) analyzing the Dirac point resonance on $L \times L$ periodic tori:
   - When $3 \mid L$, the discrete lattice exactly samples the physical Dirac nodes $\mathbf{K}, \mathbf{K}'$, causing the smallest singular value to collapse to $s_{\min} \approx 10^{-16}$ and the Hamiltonian condition number to explode to $\kappa(\mathbf{A}) \approx \mathbf{6 \times 10^{15}}$.
   - Prove that Track A fails catastrophically in this regime because forming $\mathbf{S} = -\mathbf{A}^2$ squares the condition number to $\kappa(\mathbf{S}) \approx \mathbf{3.6 \times 10^{31}} \gg \epsilon_{\mathrm{mach}}^{-1}$, losing the lower half of the 53-bit IEEE 754 mantissa (8 decimal digits of precision) and causing the covariance error to explode from $1.3 \times 10^{-14}$ to $9.5 \times 10^{-9}$.
   - Contrast this with Track B (Newton-Schulz), which avoids explicit matrix inversion and maintains stable convergence.
5. **Reporting Wilson Score Confidence Intervals in Table 1:**
   Fulfill the promise made in Section 4.2 by replacing point values in Table 1 with exact 95% Wilson score confidence intervals $[w^-, w^+]$ for all solver basin survival rates. Remove the unsubstantiated "$\ge 0.500$" lower-bound claims unless explicit numerical scans beyond $\delta = 0.500$ are reported.

---

## Submission of Revised Manuscript

When preparing your revised manuscript, you must provide:
1. A **clean version** of the revised manuscript adhering strictly to SciPost formatting guidelines. Ensure that all LaTeX tables are properly formatted without overfull `\hbox` warnings.
2. A **tracked-changes version** (or detailed diff) highlighting all textual, mathematical, and tabular modifications.
3. A comprehensive, point-by-point **Response to Referees and Editorial Directives**, detailing the exact changes made in response to each numbered objection of Referee 1, Referee 2, and the four Directives above.

The revised manuscript will be returned to both referees for re-evaluation. Given the depth of the required algebraic, analytical, and computational overhaul, a failure to fully resolve any of the directives above will result in the immediate rejection of the manuscript.

Please submit your revision within **60 days**. If you require an extension, please contact the Editorial College well in advance.

Yours sincerely,

**Handling Editor, SciPost Physics**
*On behalf of the SciPost Physics Editorial College*
`https://scipost.org/SciPostPhys`
