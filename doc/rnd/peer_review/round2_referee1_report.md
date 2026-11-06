# SciPost Physics Referee Report 2

**Manuscript Title:** Majorana Parton Mean-Field Theory of the Extended Kitaev Honeycomb Model: Algebraic Symmetries, Grassmannian Solvers, and Kinematic Invariants
**Author:** Lu Freitas
**Referee:** Referee 1 (Senior Mathematical Physicist and Quantum Many-Body Theorist)
**Date:** November 6, 2026
**Recommendation:** **Minor Revision (Acceptance Pending Synchronization of Introduction and Minor Pre-Print Polish)**

---

## 1. Executive Summary of Round 2 Assessment

In this second round of review, I evaluated the revised manuscript (`doc/manuscript/sections/`) alongside the author's detailed response record (`doc/rnd/peer_review/round2_author_response.md`).

The author has undertaken a serious, rigorous, and highly commendable revision that directly resolves the primary mathematical and physical objections raised in Round 1:
1. **Operator-Level Algebra & Parity Projector Sign:** The spin and gauge bilinears in Eqs. (13)–(14) are now normalized with the exact factor of $1/8$, establishing the true $\mathfrak{su}(2)$ algebra $[S_i^\alpha, S_i^\beta] = \ii \epsilon^{\alpha\beta\gamma} S_i^\gamma$ with quadratic Casimir $\mathbf{S}_i^2 = 3/4$. Crucially, the fatal sign inversion in the local parity projector (Eq. 24) has been corrected to $\mathcal{D}_i = + c_i^0 c_i^1 c_i^2 c_i^3 = +1$, correctly selecting the physical odd-parity Dirac doublet ($n_1+n_2=1$) and excluding the gauge doublet.
2. **Epistemic Demarcation (Theorems vs. Propositions):** In Section 3, former Theorem 1 (continuation data, spinodal fields, critical hysteresis) has been appropriately demoted to **Proposition 1**, and former Theorem 3 has been demoted to **Proposition 3**, cleanly separating the exact linear response derivative from the empirical regression $\|\boldsymbol{\lambda}^*\| \approx \sqrt{5}h$. In Section 2.3, the thermodynamic bond integral of Borwein et al. has been cleanly decoupled from finite-mesh quadrature scaling.
3. **Gauge Response Theory:** Proposition 3 now correctly incorporates the driving cross-susceptibility tensor $\boldsymbol{\chi}_{G, S} \equiv \frac{\partial \langle \mathbf{g}_s \rangle}{\partial \mathbf{h}}$ derived from the Implicit Function Theorem on $\mathbf{g}_s = \mathbf{0}$.
4. **Perturbative Decoupling Mechanism:** Remark 1 in Section 3.2 now provides a physically sound explanation for the prefactor ratio $\Delta_{\mathrm{pert}} / \Delta_K^{\mathrm{MF}} \approx 11.4$, contrasting non-local virtual two-vison intermediate states with single-particle localized Majorana bandwidths.

However, a close comparative reading reveals **one critical editorial omission**: the author updated Sections 2, 3, and 4, but **failed to synchronize Section 1 (Introduction)**. The Introduction still contains several unrevised, obsolete claims (such as claiming to "prove" empirical regressions and quoting outdated FLOP/speedup numbers) that directly contradict the revised body of the paper. Additionally, a minor numerical prefactor in the Casimir identity in Eq. (24) requires correction.

---

## 2. Detailed Technical Verification Against Round 1 Directives

### 1. Parity Projector Sign and Generator Normalization (Section 2.1)
* **Status:** **Resolved (with one minor scalar typo).**
* In Eqs. (13)–(14), the definitions:
  $$S_i^\alpha = -\frac{\ii}{8} \mathbf{c}_i^T M^\alpha \mathbf{c}_i = \frac{\ii}{4} c_i^0 c_i^\alpha - \frac{\ii}{8} \epsilon^{\alpha\beta\gamma} c_i^\beta c_i^\gamma,$$
  $$G_i^\gamma = -\frac{\ii}{8} \mathbf{c}_i^T G^\gamma \mathbf{c}_i = \frac{\ii}{4} c_i^0 c_i^\gamma + \frac{\ii}{8} \epsilon^{\gamma\alpha\beta} c_i^\alpha c_i^\beta,$$
  now reproduce the exact $\mathfrak{su}(2)$ commutation relations $[S_i^\alpha, S_i^\beta] = \ii \epsilon^{\alpha\beta\gamma} S_i^\gamma$ and $[G_i^\gamma, G_i^\delta] = \ii \epsilon^{\gamma\delta\rho} G_i^\rho$, with $[S_i^\alpha, G_i^\gamma] = 0$.
* In Eq. (24), defining $\mathcal{D}_i \equiv + c_i^0 c_i^1 c_i^2 c_i^3 = +1$ correctly identifies the physical spin-$1/2$ sector.
* **Minor Typo in Eq. (24):** The identity is written as $\mathcal{D}_i = 2 \mathbf{S}_i^2 - 2 \mathbf{G}_i^2 = +1$. Since $\mathbf{S}_i^2 = 3/4$ on the physical space and $\mathbf{G}_i^2 = 0$, evaluating $2(3/4) - 2(0) = 3/2 \neq 1$. To equal $\mathcal{D}_i = +1$, the scalar factor must be $\frac{4}{3}$, i.e.:
  $$\mathcal{D}_i = \frac{4}{3}\left(\mathbf{S}_i^2 - \mathbf{G}_i^2\right) = 1 - \frac{8}{3}\mathbf{G}_i^2 = +1.$$
  This scalar factor must be corrected in the final text.
* The unprojected Fock space status of MMFT is now explicitly and accurately stated in lines 26–30.

### 2. Epistemic Demotions and Nomenclature (Section 3)
* **Status:** **Resolved in Sections 2 and 3; Failed Synchronization in Section 1.**
* The demotions of Theorem 1 to Proposition 1 and Theorem 3 to Proposition 3 in Section 3 successfully resolve the category errors. The manuscript now properly distinguishes analytical proofs (Theorem 1 on flat-band neutrality and Theorem 2 on symmetry selection rules) from computational continuation and empirical regressions.
* In Section 2.3, Item 4, the Borwein continuous 3-step walk integral $W_3(1)/3$ is now cited properly and separated from finite-$N_k$ quadrature scaling.
* **Critical Issue in Section 1 (Introduction):**
  While Section 3 was corrected, Section 1 (`sec01_introduction.tex`) was left unrevised:
  - *Line 15:* Still asserts: *"Furthermore, we prove the Pearson Walk Discretization Law ($\mathcal{O}(N_k^{-1/2})$) governing numerical Brillouin zone quadrature errors."* (Should state: *We establish the Pearson walk scaling...*)
  - *Line 16:* Still asserts: *"we prove that the topological Dirac mass opens with a strictly cubic power law $\Delta_K = (2.592 \pm 0.003) h^3$ ... while the physical gauge Lagrange multiplier satisfies the universal asymptotic law $\|\boldsymbol{\lambda}^*\| = \sqrt{5} h$."* (Should state: *we prove that the gap opening is cubic, with numerical prefactor $2.592$, while the Lagrange multiplier exhibits linear scaling $\approx \sqrt{5}h$...*)
  - *Line 18:* Still asserts: *"halving memory ($2.0\times$) and cutting FLOPs by $4.0\times$ ... achieving a 63.53x speedup"*, which directly contradicts Section 4.2 (which refutes the $4.0\times$ FLOP claim as a cache artifact) and Table 2 (which lists measured speedups of $\approx 2.31\times$).
  Section 1 must be edited immediately to match the rigor of Sections 2–4.

### 3. Physical Rigor of Perturbation Theory (Section 3.2, Remark 1)
* **Status:** **Resolved.**
* The explanation in Remark 1 accurately identifies the fundamental physical difference between spin perturbation theory and Majorana mean-field theory:
  - In exact spin perturbation theory, the Zeeman operator flips spins and creates pairs of static $\mathbb{Z}_2$ visons on adjacent plaquettes, giving an intermediate-state energy denominator $(2\Delta_{\mathrm{flux}})^2 \approx (0.26|K|)^2 \approx 0.068|K|^2$.
  - In 4-flavor Majorana mean-field theory, the $\mathbb{Z}_2$ fluxes are replaced by static quadratic expectations, and intermediate states in single-particle perturbation theory correspond to the flat gauge Majorana bands at energy $\pm |K|$, giving a denominator of order $|K|^2 \approx 14.7 \times (2\Delta_{\mathrm{flux}})^2$.
  This physically elucidates the prefactor ratio $\Delta_{\mathrm{pert}} / \Delta_K^{\mathrm{MF}} \approx 11.4$.

---

## 3. Actionable Demands for Final Acceptance

Before the manuscript is formally recommended for production in SciPost Physics, the author must implement two simple, localized fixes:

1. **Synchronize Section 1 (`sec01_introduction.tex`):**
   - Update Line 15 to remove the claim of having "proven" the Pearson Walk Discretization Law; phrase it as an empirical scaling law for finite-mesh quadrature.
   - Update Line 16 to state that symmetry selection rules prove the *cubic exponent* $\mathcal{O}(h^3)$, while the prefactor $2.592$ and the Lagrange multiplier relation $\|\boldsymbol{\lambda}^*\| \approx \sqrt{5}h$ are empirical regressions.
   - Update Line 18 to replace the obsolete "cutting FLOPs by $4.0\times$" and "63.53x speedup" with the verified cache-locality explanation and the sustained throughput/speedups reported in Section 4 ($114.6\text{ GFLOP/s}$ and measured $1.4\times\text{--}2.3\times$ speedups).
2. **Correct the Casimir Scalar Prefactor in Eq. (24) (`sec02`):**
   - Change $\mathcal{D}_i = 2 \mathbf{S}_i^2 - 2 \mathbf{G}_i^2$ to $\mathcal{D}_i = \frac{4}{3}(\mathbf{S}_i^2 - \mathbf{G}_i^2) = 1 - \frac{8}{3}\mathbf{G}_i^2 = +1$.

---

## 4. Final Recommendation

**Recommendation:** **Minor Revision (Acceptance Pending Synchronization)**

The author has successfully transformed the manuscript into a rigorous, well-grounded work of mathematical physics and computational quantum many-body theory. Once the localized synchronization fixes in Section 1 and Eq. (24) are applied, the manuscript will meet all criteria for publication in **SciPost Physics**.
