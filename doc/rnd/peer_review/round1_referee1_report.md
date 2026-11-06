# SciPost Physics Referee Report 1

**Manuscript Title:** Majorana Parton Mean-Field Theory of the Extended Kitaev Honeycomb Model: Algebraic Symmetries, Grassmannian Solvers, and Kinematic Invariants
**Author:** Lu Freitas
**Referee:** Referee 1 (Senior Mathematical Physicist and Quantum Many-Body Theorist)
**Date:** November 6, 2026
**Recommendation:** **Major Revision (Substantial Mathematical and Formal Overhaul Required)**

---

## 1. Executive Summary & General Scientific Assessment

The manuscript presents an extensive Hartree-Fock-Bogoliubov mean-field investigation of the extended Kitaev honeycomb model formulated in terms of a four-flavor Majorana parton representation per physical lattice site. The author constructs a geometric framework on the orthogonal Grassmannian manifold bundle, examines non-linear iterative fixed-point solvers (Picard, DIIS, Broyden, JFNK), details real-space skew-symmetric reformulations, and benchmarks matrix sign purification methods.

From a mathematical physics and quantum many-body perspective, the manuscript contains interesting computational results, particularly regarding the Grassmannian formulation of the saddle-point equations and the high-throughput matrix sign purification on real-space clusters.

However, in its present form, the manuscript suffers from **severe formal, algebraic, and conceptual defects**. Most critically:
1. **Operator Normalization and Parity Projection Errors:** Equations (13), (14), and (21) in Section 2 contain conflicting normalization prefactors and an erroneous sign in the local parity projector. If taken literally, Eq. (21) projects onto the unphysical gauge-doublet sector rather than the physical spin-1/2 doublet.
2. **Category Errors in Mathematical Nomenclature:** Multiple empirical numerical regressions and continuation data (e.g., critical fields with floating-point error bars, spinodal boundaries, and an empirical multiplier scaling $\approx \sqrt{5}h$) are improperly designated as **Theorems** (Theorem 1 and Theorem 3 in Section 3). In mathematical physics, empirical fits cannot be labeled as theorems.
3. **Incomplete Gauge Response Theory:** The linear response relation in Theorem 3 omits the driving cross-susceptibility tensor that couples the external magnetic field to the gauge constraint.
4. **Superficial Perturbative Comparison:** The discussion comparing mean-field Dirac mass opening ($\Delta_K^{\mathrm{MF}} \approx 2.592 h^3$) with Kitaev's exact third-order chiral perturbation theory ($\Delta_{\mathrm{pert}} \approx 29.5 h^3$) overlooks the fundamental contrast between non-local virtual two-vison intermediate states in spin space and single-particle band hybridization in Majorana Fock space, while missing an obvious opportunity to provide an exact analytical perturbation expansion of the single-particle Hamiltonian.
5. **Conflation of Thermodynamic Integrals with Discretization Noise:** The claim in Section 1 and Section 2.3 of having "proven the Pearson Walk Discretization Law" mischaracterizes Borwein et al.'s classical random walk evaluation of the continuous bond correlator integral by confusing it with numerical quadrature noise on discrete meshes.

Below, I provide a detailed, numbered exposition of these objections, complete with specific line and equation citations, followed by mandatory actionable revision demands.

---

## 2. Major Conceptual, Algebraic, and Mathematical Objections

### Objection 1: Operator Inconsistencies in the $\mathfrak{so}(4) \cong \mathfrak{su}(2)_{\mathrm{spin}} \oplus \mathfrak{su}(2)_{\mathrm{gauge}}$ Algebra and Fatal Sign Inversion in Local Parity Projector $D_i$

**Citations:**
* `doc/manuscript/sections/sec02_algebraic_foundations_and_invariants.tex`, Lines 6–16 (Eqs. 13–14) and Lines 19–22 (Eq. 21).

**Mathematical Scrutiny:**
1. **Prefactor Factor-of-Two Discrepancy in Spin Bilinears:**
   In Section 2.1, the author introduces four Majorana operators $\{c_i^\mu, c_j^\nu\} = 2\delta_{ij}\delta^{\mu\nu}\mathbf{1}$ and writes the spin and gauge bilinears in Eq. (13) and Eq. (14) as:
   $$S_i^\alpha = -\frac{\ii}{4} \mathbf{c}_i^T M^\alpha \mathbf{c}_i = \frac{\ii}{2} c_i^0 c_i^\alpha - \frac{\ii}{4} \epsilon^{\alpha\beta\gamma} c_i^\beta c_i^\gamma,$$
   $$G_i^\gamma = -\frac{\ii}{4} \mathbf{c}_i^T G^\gamma \mathbf{c}_i = \frac{\ii}{2} c_i^0 c_i^\gamma + \frac{\ii}{4} \epsilon^{\gamma\alpha\beta} c_i^\alpha c_i^\beta.$$
   The text immediately claims (Line 16) that these satisfy the standard $\mathfrak{su}(2)$ Lie algebra:
   $$[S_i^\alpha, S_i^\beta] = \ii \epsilon^{\alpha\beta\gamma} S_i^\gamma, \qquad [G_i^\gamma, G_i^\delta] = \ii \epsilon^{\gamma\delta\rho} G_i^\rho.$$
   **This is mathematically false.**
   Let us compute the commutator directly from the Clifford CAR algebra. For $\alpha=1, \beta=2$, the expansion in Eq. (13) gives:
   $$S^1 = \frac{\ii}{2}(c^0 c^1 - c^2 c^3), \qquad S^2 = \frac{\ii}{2}(c^0 c^2 - c^3 c^1).$$
   Evaluating the commutator yields:
   $$[S^1, S^2] = -\frac{1}{4} [c^0 c^1 - c^2 c^3, \, c^0 c^2 - c^3 c^1] = -(c^0 c^3 - c^1 c^2) = 2 \ii S^3.$$
   Notice the factor of **$2$**. The operators defined in Eq. (13) and Eq. (14) satisfy the algebra of **Pauli operators** $\boldsymbol{\sigma}$, not spin-$1/2$ operators $\mathbf{S} = \frac{1}{2}\boldsymbol{\sigma}$. Indeed, the quadratic Casimir of the operators in Eq. (13) evaluates to $\mathbf{S}^2 = 3$ on the physical space, rather than $s(s+1) = 3/4$.
   To represent physical spin-$1/2$ generators satisfying $[S_i^\alpha, S_i^\beta] = \ii \epsilon^{\alpha\beta\gamma} S_i^\gamma$, the bilinear must be normalized as:
   $$S_i^\alpha = \frac{\ii}{8} \mathbf{c}_i^T M^\alpha \mathbf{c}_i = \frac{\ii}{4}(c_i^0 c_i^\alpha - \epsilon^{\alpha\beta\gamma} c_i^\beta c_i^\gamma), \qquad G_i^\alpha = \frac{\ii}{8} \mathbf{c}_i^T G^\alpha \mathbf{c}_i.$$
   The manuscript has conflated the Pauli matrix $\sigma_i^\alpha \equiv \ii c_i^0 c_i^\alpha$ with the spin operator $S_i^\alpha$.

2. **Fatal Sign Inversion in Local Parity Projector $D_i$:**
   In Eq. (21), the author defines the local physical projector via:
   $$D_i = - c_i^0 c_i^1 c_i^2 c_i^3 = 1 - 2 \mathbf{G}_i^2 = +1.$$
   This sign is completely erroneous under the standard CAR representation.
   Introduce the canonical Dirac fermionic modes:
   $$d_1 = \frac{1}{2}(c^0 + \ii c^1), \qquad d_2 = \frac{1}{2}(c^2 + \ii c^3),$$
   with occupation number operators $n_a = d_a^\dagger d_a \in \{0, 1\}$. Direct substitution shows:
   $$c^0 c^1 = \ii (1 - 2n_1), \qquad c^2 c^3 = \ii (1 - 2n_2).$$
   Consequently:
   $$c^0 c^1 c^2 c^3 = -(1 - 2n_1)(1 - 2n_2) = -(-1)^{n_1 + n_2}.$$
   The physical spin-$1/2$ Hilbert space $\mathcal{H}_{\mathrm{phys}, i} \cong (1/2, 0)_{\mathrm{spin, gauge}}$ corresponds to the **odd-parity** Dirac subspace $\operatorname{span}\{|10\rangle, |01\rangle\}$ where $n_1 + n_2 = 1$. In this subspace:
   $$c^0 c^1 c^2 c^3 = -(-1)^1 = +1.$$
   Therefore, the physical projector requires $\mathcal{D}_i \equiv + c_i^0 c_i^1 c_i^2 c_i^3 = +1$.
   By inserting an extraneous minus sign and defining $D_i = - c_i^0 c_i^1 c_i^2 c_i^3 = +1$, the author forces $c_i^0 c_i^1 c_i^2 c_i^3 = -1$, which selects the **even-parity** Dirac subspace $\operatorname{span}\{|00\rangle, |11\rangle\}$!
   In that subspace, the gauge Casimir is $\mathbf{G}^2 = 3/4$ (gauge doublet) and the physical spin Casimir vanishes identically ($\mathbf{S}^2 = 0$). The definition in Eq. (21) literally projects onto the **completely unphysical gauge-variant sector**. This error must be rectified immediately.

---

### Objection 2: Mislabeling of Empirical Numerical Fits and Continuation Data as "Theorems"

**Citations:**
* `doc/manuscript/sections/sec01_introduction.tex`, Line 16;
* `doc/manuscript/sections/sec03_field_response_and_multistability.tex`, Lines 13–29 (Theorem 1) and Lines 70–83 (Theorem 3).

**Mathematical Scrutiny:**
SciPost Physics and the broader mathematical physics community maintain strict criteria separating analytical theorems from empirical or numerical simulations:
1. **Theorem 1 is not a Theorem:**
   Under Section 3.1, the heading "Theorem 1 (First-Order Critical Boundary and Hysteresis)" is attached to a list of numerical continuation outputs:
   - "$h_c^* / |K| = 0.3325 \pm 0.0005$"
   - "$h_{s2} \approx 0.40 |K|$, $h_{s1} \approx 0.28 |K|$"
   - "$\Delta M = 0.412 \pm 0.005$, $\Delta E = 0.057 |K|$"
   A numerical data point obtained from numerical iterations on finite meshes with statistical/discretization uncertainty intervals ($\pm 0.0005$) is an **empirical simulation result**, not a mathematical theorem. There is no proof of the existence of the thermodynamic phase transition, no analytical proof of the discontinuity of the order parameter, and no proof of the spinodal boundaries. Labeling this a "Theorem" is an unacceptable breach of editorial and scientific rigor.
2. **Theorem 3 Contains Empirical Regressions:**
   Theorem 3 ("Linear Gauge Susceptibility Law") asserts:
   $$\|\boldsymbol{\lambda}^*\| \approx \sqrt{5} \cdot h \approx 2.236 \cdot h, \quad \text{with deviations bounded by } < 0.8\% \text{ across } h \le 0.10.$$
   An approximate empirical proportionality containing an unproven numerical factor ($\sqrt{5}$) and an empirical error margin ($< 0.8\%$) cannot be stated as a Theorem.
3. **Introduction Claim of Proof for Numerical Regressions:**
   In Section 1, Line 16, the manuscript claims: *"we prove that the topological Dirac mass opens with a strictly cubic power law $\Delta_K = (2.592 \pm 0.003) h^3$ ... while the physical gauge Lagrange multiplier satisfies the universal asymptotic law $\|\boldsymbol{\lambda}^*\| = \sqrt{5} h$."*
   One can prove from symmetry selection rules that the leading term in an analytic expansion must be proportional to $h^3$ (Theorem 2). One **cannot** "prove" an empirical coefficient with numerical regression error bars $(2.592 \pm 0.003)$. The prefactor is a numerical fit, not a proved mathematical constant.

---

### Objection 3: Incomplete Linear Response Formulation of the Gauge Lagrange Multipliers

**Citations:**
* `doc/manuscript/sections/sec03_field_response_and_multistability.tex`, Lines 70–83 (Theorem 3, Eq. 75).

**Mathematical Scrutiny:**
In Section 3.3, Eq. (75), the author writes:
$$\boldsymbol{\lambda}_s^*(\mathbf{h}) = -\boldsymbol{\chi}_{\mathrm{gauge}}^{-1} \mathbf{h} + \mathcal{O}(h^3),$$
where $\boldsymbol{\chi}_{\mathrm{gauge}} = \frac{\partial \langle \mathbf{g}_s \rangle}{\partial \boldsymbol{\lambda}_s}$ is the static local gauge susceptibility tensor.
This formula is dimensionally and physically incomplete.
In the Lagrangian / KKT formulation:
$$\mathcal{L}(\Omega, \boldsymbol{\lambda}; \mathbf{h}) = E(\Omega) - \sum_{s} \mathbf{h} \cdot \langle \mathbf{S}_s \rangle - \sum_s \boldsymbol{\lambda}_s \cdot \langle \mathbf{G}_s \rangle.$$
The saddle-point condition enforcing gauge neutrality is the vector constraint equation:
$$\mathbf{g}_s(\boldsymbol{\lambda}^*(\mathbf{h}), \mathbf{h}) \equiv \langle \mathbf{G}_s \rangle = \mathbf{0}.$$
Differentiating this identity with respect to the external field $\mathbf{h}$ via the Implicit Function Theorem yields:
$$\frac{\mathrm{d}\mathbf{g}_s}{\mathrm{d}\mathbf{h}} = \frac{\partial \langle \mathbf{G}_s \rangle}{\partial \boldsymbol{\lambda}_s} \frac{\mathrm{d}\boldsymbol{\lambda}_s^*}{\mathrm{d}\mathbf{h}} + \frac{\partial \langle \mathbf{G}_s \rangle}{\partial \mathbf{h}} = \mathbf{0},$$
which implies:
$$\frac{\mathrm{d}\boldsymbol{\lambda}_s^*}{\mathrm{d}\mathbf{h}} = - \left( \frac{\partial \langle \mathbf{G}_s \rangle}{\partial \boldsymbol{\lambda}_s} \right)^{-1} \left( \frac{\partial \langle \mathbf{G}_s \rangle}{\partial \mathbf{h}} \right) = - \boldsymbol{\chi}_{\mathrm{gauge}}^{-1} \, \boldsymbol{\chi}_{G, S},$$
where $\boldsymbol{\chi}_{G, S} \equiv \frac{\partial \langle \mathbf{G}_s \rangle}{\partial \mathbf{h}}$ is the **cross-susceptibility tensor** measuring the linear response of the unprojected gauge bilinears to the physical Zeeman field.
By writing simply $-\boldsymbol{\chi}_{\mathrm{gauge}}^{-1}\mathbf{h}$ in Eq. (75), the author has implicitly and silently assumed that $\boldsymbol{\chi}_{G, S} = \mathbf{I}$, equating the spin and gauge operator matrix elements without justification. The cross-susceptibility $\boldsymbol{\chi}_{G, S}$ must be derived explicitly from the Majorana band eigenvectors, and the origin of the numerical constant $\sqrt{5}$ must be clarified: is it an exact algebraic property of the cross-susceptibility product $\boldsymbol{\chi}_{\mathrm{gauge}}^{-1} \boldsymbol{\chi}_{G, S}$, or merely an accidental numerical coincidence of the fit?

---

### Objection 4: Perturbation Theory Claims: Spin Space Visons vs. Fock Space Majorana Hybridization

**Citations:**
* `doc/manuscript/sections/sec03_field_response_and_multistability.tex`, Lines 51–65 (Remark 1).

**Mathematical Scrutiny:**
In Remark 1, the author compares the mean-field Dirac gap opening $\Delta_K^{\mathrm{MF}} = 2.592 (h/|K|)^3$ with Kitaev's exact third-order perturbation result $\Delta_{\mathrm{pert}} \approx 29.5 (h/|K|)^3$, noting a ratio of $\sim 11.4$.
While the remark correctly attributes the discrepancy to the difference between virtual vison intermediate states and single-particle Majorana bandwidths, the discussion remains qualitatively superficial and mathematically incomplete:
1. **The Nature of Intermediate States:**
   In Kitaev's exact spin perturbation theory (Kitaev 2006, Section 8), the unperturbed ground state lies in the vortex-free sector of the physical spin Hilbert space. A single Zeeman operator $\sigma_j^\alpha$ anticommutes with the two plaquette flux operators $W_p$ sharing the $\alpha$-bond, creating a pair of localized $\mathbb{Z}_2$ visons. The third-order perturbation process $\sigma_j^x \sigma_l^z \sigma_k^y$ around a vertex involves intermediate states with **two static visons** at energy $2\Delta_{\mathrm{flux}} \approx 0.26 |K|$. The perturbation denominator is therefore $(2\Delta_{\mathrm{flux}})^2 \approx 0.068 |K|^2$.
   In sharp contrast, in the 4-flavor Majorana mean-field theory, the $\mathbb{Z}_2$ gauge field is replaced by static mean-field bond expectation values $\chi, \eta$. The intermediate states in the single-particle Rayleigh-Schrödinger perturbation expansion of the $8 \times 8$ Bloch Hamiltonian $\mathcal{H}_{\mathrm{MF}}(\mathbf{k})$ are the **flat bands of the gauge Majoranas** $c^1, c^2, c^3$ situated at energy $\pm |K|$. The single-particle energy denominator is therefore governed by $|K|^2$, which is larger than the two-vison denominator by a factor of $|K|^2 / (2\Delta_{\mathrm{flux}})^2 \approx 1 / 0.068 \approx 14.7$. This explains why the mean-field gap is smaller by an order of magnitude.
2. **Missing Analytical Single-Particle Expansion:**
   Because the mean-field Hamiltonian is a transparent $8 \times 8$ matrix in momentum space, why did the author not calculate the 3rd-order degenerate perturbation expansion of $\mathcal{H}_{\mathrm{MF}}(\mathbf{K})$ **analytically**? Such a calculation would yield the exact closed-form algebraic prefactor $C_{\mathrm{MF}}$ such that $\Delta_K^{\mathrm{MF}}(h) = C_{\mathrm{MF}} h^3$, allowing an exact analytical comparison with $C_{\mathrm{pert}} = 6\sqrt{3}/\Delta E^2$. Instead, the manuscript relies entirely on a numerical power-law regression $(2.592 \pm 0.003)$. A mathematical physics paper should provide the exact algebraic calculation of $C_{\mathrm{MF}}$.

---

### Objection 5: Conflation of Pearson Walk Integrals (Borwein et al.) with Finite-Mesh Quadrature Noise

**Citations:**
* `doc/manuscript/sections/sec01_introduction.tex`, Line 15;
* `doc/manuscript/sections/sec02_algebraic_foundations_and_invariants.tex`, Lines 49–54 (Item 4, Eq. 52).

**Mathematical Scrutiny:**
In Section 1, Line 15, the manuscript asserts: *"we prove the Pearson Walk Discretization Law ($\mathcal{O}(N_k^{-1/2})$) governing numerical Brillouin zone quadrature errors."*
In Section 2.3, Item 4, Eq. (52), this is presented as:
$$\epsilon_{\mathrm{gauge}}(N_k) = \sqrt{\frac{3\pi}{8 N_k}} \propto N_k^{-1/2}.$$
This text presents a severe confusion:
1. **Borwein et al. Pearson Walk Integral:**
   The classical 3-step Pearson planar random walk integral (Borwein, Straub, Wan, Zudilin, *Canad. J. Math.* 2012) evaluates the continuous thermodynamic integral of the honeycomb band dispersion:
   $$3\bar{U} = \int_{[0, 1)^2} |1 + e^{2\pi \ii k_1} + e^{2\pi \ii k_2}| \, \mathrm{d}k_1 \mathrm{d}k_2 = W_3(1) = \frac{3 \cdot 2^{1/3}}{16 \pi^4} \Gamma^6(1/3) + \frac{16 \cdot 2^{2/3}\pi^2}{\Gamma^6(1/3)} \approx 1.574597237.$$
   This is an exact evaluation of the **thermodynamic matter bond correlator** $\bar{U} \approx 0.52486574585$.
2. **Discretization Noise on a Deterministic Mesh:**
   Evaluating an integral over the 2D Brillouin zone on a regular $N_k$-point grid using the midpoint rule is a deterministic quadrature procedure. If an integrand has conical Dirac singularities, Euler-Maclaurin expansions fail, producing algebraic convergence. However, calling the finite-grid residual a "Pearson Walk Law" and claiming to have "proved" it with prefactor $\sqrt{3\pi / 8N_k}$ is mathematically unsubstantiated in the text. There is no proof anywhere in the manuscript establishing that regular grid sampling of the gauge residual behaves as a 2D random flight with Rayleigh-distributed step phases. If this is an empirical hypothesis, it must be labeled as such; if it is a mathematical theorem, the full proof must be provided.

---

### Objection 6: Hilbert Space vs. Fock Space Dimensionality and Physical Projection in Mean-Field Theory

**Citations:**
* `doc/manuscript/sections/sec01_introduction.tex`, Line 9;
* `doc/manuscript/sections/sec02_algebraic_foundations_and_invariants.tex`, Lines 18–23.

**Mathematical Scrutiny:**
In Section 1, Line 9, the manuscript states that four Majoranas per site *"doubles the local state space from $\dim \mathcal{H}_{\mathrm{phys}, i} = 2$ to $\dim \mathcal{F}_i = 4$."*
While this is true for a single isolated site, the author fails to state the consequences for the many-body system:
- For a lattice of $N$ sites, the physical spin Hilbert space dimension is $\dim \mathcal{H}_{\mathrm{phys}} = 2^N$.
- The unprojected Majorana Fock space dimension is $\dim \mathcal{F} = 4^N = 2^{2N}$.
In exact quantum many-body theory, physical states must satisfy $P_G |\Psi_{\mathrm{phys}}\rangle = |\Psi_{\mathrm{phys}}\rangle$, where $P_G = \prod_i \frac{1 + \mathcal{D}_i}{2}$ retains only a fraction $2^{-N}$ of the Fock space.
In quadratic mean-field theory, Wick's theorem is applied to an unprojected quasi-free state $|\Psi_{\mathrm{MF}}\rangle \in \mathcal{F}$. The gauge constraint is enforced **only on average** via Lagrange multipliers: $\langle \Psi_{\mathrm{MF}} | \mathbf{G}_i | \Psi_{\mathrm{MF}} \rangle = \mathbf{0}$.
The manuscript must explicitly clarify that:
1. The mean-field state $|\Psi_{\mathrm{MF}}\rangle$ is **not** an element of the physical Hilbert space $\mathcal{H}_{\mathrm{phys}}$, but an unprojected quasi-free state on $\mathcal{F}$.
2. Enforcing $\langle \mathbf{G}_i \rangle = \mathbf{0}$ does not eliminate gauge fluctuations; it merely ensures that the state is a stationary saddle point of the KKT Lagrangian. Gauge-projected variational Monte Carlo or Gutzwiller projection $P_G |\Psi_{\mathrm{MF}}\rangle$ would be required to restore strict physical state representability.

---

## 3. Minor Editorial and Notational Deficiencies

1. **Section 1, Line 18:** The text claims a *"63.53x speedup over standard divide-and-conquer eigensolvers at machine-precision accuracy ($4.22 \times 10^{-15}$)."* The error metric should be explicitly stated as the relative Frobenius norm error $\|\mathbf{Y}^2 + \mathbf{I}\|_F / \sqrt{M}$, rather than unqualified "accuracy".
2. **Section 2.1, Eq. (10):** The index condition $(\mu < \nu)$ is specified, but in Eq. (13) and (14) the summation conventions over Greek indices are mixed between $0, 1, 2, 3$ and $1, 2, 3$. The manuscript should use distinct index sets (e.g., $\mu, \nu \in \{0, 1, 2, 3\}$ for Majorana flavors and $\alpha, \beta, \gamma \in \{1, 2, 3\}$ for spatial/spin coordinates).
3. **Section 4.4, Eq. (46):** The Newton-Schulz iteration is written as:
   $$\mathbf{Y}_{k+1} = \frac{1}{2} \mathbf{Y}_k (3\mathbf{I} + \mathbf{Y}_k^2).$$
   The author should explicitly note for the reader that the plus sign arises because $\mathbf{Y}_k$ is real skew-symmetric ($\mathbf{Y}_k^2 \le 0$), so that $\ii \mathbf{Y}_k$ is Hermitian and the standard iteration $X_{k+1} = \frac{1}{2} X_k (3\mathbf{I} - X_k^2)$ maps directly to this formula under $X_k = \ii \mathbf{Y}_k$.

---

## 4. Required Actionable Revisions

To meet the rigorous standards of SciPost Physics, the author must execute the following revisions:

1. **Correct the Majorana Lie Algebra and Parity Projector:**
   - Correct Eq. (13) and Eq. (14) to either normalize $S_i^\alpha = \frac{\ii}{8} \mathbf{c}_i^T M^\alpha \mathbf{c}_i$ (if defining spin-$1/2$ generators with $[S^\alpha, S^\beta] = \ii \epsilon^{\alpha\beta\gamma} S^\gamma$) or explicitly denote them as Pauli operators $\sigma_i^\alpha$ with commutation relations $[\sigma^\alpha, \sigma^\beta] = 2\ii \epsilon^{\alpha\beta\gamma} \sigma^\gamma$.
   - Fix Eq. (21): remove the incorrect minus sign and define $\mathcal{D}_i = + c_i^0 c_i^1 c_i^2 c_i^3 = +1$. Explicitly demonstrate the mapping to the Dirac occupation basis to prove that this selects the physical doublet $\mathcal{F}_+ \cong (1/2, 0)$.
2. **Eliminate All Misnomers of "Theorems" for Empirical Fits:**
   - Demote **Theorem 1** to **"Numerical Observation 1"** or **"Computational Result: First-Order Transition and Hysteresis"**. Clearly state that the critical and spinodal fields are numerical estimates with finite-mesh uncertainties.
   - Demote **Theorem 3** to **"Proposition: Linear Gauge Susceptibility and Empirical Response"**. Separate the formal linear response relation from the empirical regression $\|\boldsymbol{\lambda}^*\| \approx \sqrt{5} h$.
   - In Section 1 and Section 3.2, rephrase all claims of having "proven" the empirical prefactor $\Delta_K = (2.592 \pm 0.003) h^3$ or the numerical ratio $\sqrt{5}$. State clearly that symmetry proves the cubic exponent $h^3$, while the prefactor is determined numerically.
3. **Rigorous Formulation of the Gauge Response Law:**
   - In Section 3.3, provide the complete vector linear response derivation including the cross-susceptibility tensor $\boldsymbol{\chi}_{G, S} = \frac{\partial \langle \mathbf{G}_s \rangle}{\partial \mathbf{h}}$. Explain whether the prefactor $\sqrt{5}$ can be derived analytically from $\boldsymbol{\chi}_{\mathrm{gauge}}^{-1} \boldsymbol{\chi}_{G, S}$ or if it is purely an empirical fit.
4. **Deepen the Perturbation Theory Section:**
   - Expand Remark 1 to clearly contrast the non-local 2-vison intermediate states in Kitaev's spin perturbation theory with the local single-particle gauge Majorana band hybridization in mean-field theory.
   - Provide the analytical 3rd-order perturbation expansion of the $8 \times 8$ single-particle mean-field Hamiltonian $\mathcal{H}_{\mathrm{MF}}(\mathbf{K})$ to determine the exact analytical prefactor $C_{\mathrm{MF}}$.
5. **Clarify the Pearson Walk Discussion:**
   - Disentangle the exact thermodynamic bond correlator evaluation (Borwein et al. 2012) from the discrete finite-mesh residual scaling. If Eq. (52) is an empirical model for the quadrature error, label it as an empirical scaling law and remove claims of having mathematically "proven" it.
6. **Clarify Hilbert vs. Fock Space Constraints:**
   - Add an explicit discussion in Section 1 or Section 2.1 clarifying that MMFT operates on the $4^N$-dimensional unprojected Fock space $\mathcal{F}$, enforcing the gauge constraint only at the expectation-value level $\langle \mathbf{G}_i \rangle = \mathbf{0}$.

---

## 5. Recommendation

**Overall Recommendation:** **Major Revision**

The manuscript addresses a timely and computationally demanding problem in frustrated quantum magnetism and offers valuable algorithmic insights. However, the presence of algebraic errors in the foundational operator definitions, the conflation of empirical fits with mathematical theorems, and the oversimplified perturbation theory comparisons prevent publication in SciPost Physics in its current form.

The manuscript can be recommended for publication only after a rigorous, thorough revision addressing all the objections detailed above.
