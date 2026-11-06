---
id: T01
title: "Epistemic Demarcation of Theorems vs. Empirical Propositions and Structural Reconciliation of the Cubic Dirac Mass Gap"
status: Completed
date: "2026-11-06"
questions:
  - "Do any \\begin{theorem} or \\begin{proof} environments in the Living Lab Book contain empirical regressions, finite-grid error bars, or numerical script benchmarks?"
  - "How does the 4-Majorana self-consistent mean-field Dirac mass gap Delta_K^MF(h) = c_MF (h/|K|)^3 structurally differ from Kitaev's 3rd-order spin-perturbation gap without invoking ad-hoc hand-waved arithmetic?"
depends_on:
  - M1.6-T7
  - M1.7-T1
  - M1.10-T1
contracts:
  - K-010
  - K-019
  - Contract H6
artifacts:
  - doc/labbook/chapters/ch_m1_6_solvers_and_basins.tex
  - doc/labbook/chapters/ch_m1_8_supercell_bloch_scaling.tex
  - doc/rnd/phase_1/m1_10_t1_exact_perturbation_and_normalization_reconciliation.md
  - doc/manuscript/sections/sec03_field_response_and_multistability.tex
---

# T01: Epistemic Demarcation of Theorems vs. Empirical Propositions and Structural Reconciliation of the Cubic Dirac Mass Gap

## 1. Scope & Adversarial Audit Findings (`[B3]`, `[B5]`)

During the three-agent adversarial peer review of Waves 1--4 (`Referee_1_Theory`, `Referee_2_Numerics`, `Editor_SciPost`), two blocking epistemic and theoretical defects were identified across Front 1 (`doc/rnd/`), Front 2 (`doc/labbook/`), and Front 3 (`doc/manuscript/`):

1. **Defect `[B3]` (Epistemic Conflation in Living Lab Book Theorems):**
   Under Rule 03 (`03_editorial_rigor.md`) and `AGENTS.md` Section 4, `\begin{theorem}` and `\begin{proof}` environments are reserved strictly for exact algebraic identities and analytical proofs; finite-grid numerical benchmarks, power-law regressions with $\pm$ standard errors, and $R^2$ goodness-of-fit statistics must reside in `\begin{proposition}` or `\begin{remark}` environments. In the Living Lab Book (`doc/labbook/`), four environments violated this demarcation:
   - `ch_m1_6_solvers_and_basins.tex`:
     - `Theorem 7.2` (`thm:m1_6_ab_failure`): Reports numerical A/B benchmark outcomes (`Arm A` vs. `Arms B, C`) and roundoff thresholds.
     - `Theorem 7.3` (`thm:m1_6_critical_transition`): Reports finite-step ($\Delta h = 0.002$) numerical sweep values $h_c^* = 0.3325 \pm 0.0002$ and $h_{c2}^+ = 0.400 \pm 0.002$.
     - `Theorem 7.4` (`thm:m1_6_cubic_gap`): Combines the exact symmetry proof of the cubic exponent ($\mathcal{I}^{(3)}(\mathbf{h}) = h_x h_y h_z \propto h^3$) with the numerical regression $(2.5921 \pm 0.005)(h/|K|)^{2.9936 \pm 0.0010}$ and cites `scripts/m1_6/benchmark_majorana_gap.py` inside `\begin{proof}`.
   - `ch_m1_8_supercell_bloch_scaling.tex`:
     - `Theorem 8.2` (`thm:m1_7_stability`): Reports 10-trial empirical random-seed stability metrics ($P_{\mathrm{conv}} = 1.000$, $N_{\mathrm{nem}} \le 7.32 \times 10^{-8}$).

2. **Defect `[B5]` (Hand-Waved Arithmetic in Cubic Dirac Mass Gap Reconciliation):**
   In `doc/rnd/phase_1/m1_10_t1_exact_perturbation_and_normalization_reconciliation.md` (line 99), `ch_m1_6_solvers_and_basins.tex` (lines 289--295), and `sec03_field_response_and_multistability.tex` (line 73), the comparison between Kitaev's third-order spin-perturbation gap formula and the 4-Majorana self-consistent mean-field gap $\Delta_{\mathrm{canonical}}^{\mathrm{MF}}(h) = 4 \Delta_K^{\mathrm{MF}}(h) = 10.3684 h^3/|K|^2$ substituted a naive gauge flat-band gap $2\Delta_{\mathrm{gauge}} \approx 1.05 |K|$ into the spin-perturbation formula $\frac{4}{\Delta E^2} h^3$ to get $3.63 h^3$, and then bridged $3.63 h^3 \to 10.37 h^3$ via an un-derived "screening factor of $2.8510$".

---

## 2. Guiding Questions & Falsifiable Hypotheses

- **Guiding Question Q1 (Epistemic Demarcation):** Can every statement in Chapters 7 and 8 of the Living Lab Book be cleanly partitioned so that `\begin{theorem}` contains exclusively exact algebraic/symmetry propositions with analytical proofs, while all finite-size numerical fits and solver benchmarks are housed in `\begin{proposition}` environments?
- **Guiding Question Q2 (Operator-Level Mechanism of the Mean-Field Cubic Gap):** Why is a naive substitution $\Delta E \to 2\Delta_{\mathrm{gauge}}$ into the spin-perturbation formula $\frac{4}{\Delta E^2}h^3$ structurally incomplete in 4-Majorana Hartree-Fock-Bogoliubov (HFB) theory?
  - *Structural Analysis:*
    1. **Factor of 4 (Canonical vs. Bloch Generator Normalization):** In `src/kitaev_mf/services/hamiltonian.py`, the single-particle Bloch matrix is $H_{\mathrm{Bloch}}(\mathbf{k}) = \frac{\mathrm{i}}{4} A(\mathbf{k})$, whereas canonical Majorana normal-mode excitation energies are eigenvalues of $\mathrm{i} A(\mathbf{k})$. This gives the exact algebraic identity $\Delta_{\mathrm{canonical}}^{\mathrm{MF}}(h) = 4 \Delta_K^{\mathrm{MF}}(h) = 10.3684 \, h^3/|K|^2$.
    2. **Structural Difference Between Exact Spin Perturbation Theory and HFB Mean-Field Perturbation Theory:**
       - In **exact spin perturbation theory** on the physical subspace $\mathcal{F}_{\mathrm{phys}}$, $\sigma_i^\gamma = \mathrm{i} c_i^0 c_i^\gamma$ flips the $\mathbb{Z}_2$ bond gauge variable $u_{\langle ij\rangle_\gamma} \mapsto -u_{\langle ij\rangle_\gamma}$, exciting a non-local two-vison pair plus arbitrary gapless matter fermion pairs across the Brillouin zone. The exact third-order Feshbach resolvent $Q(E_0 - H_K)^{-1}Q$ sums over the full two-vison + matter continuum (whose isolated two-vison threshold is $\Delta_{\mathrm{2v}} \approx 0.264 |K|$, while Kitaev's static single-mode approximation $\Delta_{\mathrm{flux}} \approx 0.264 |K|$ neglects the dispersive matter continuum; computing the exact continuum-integrated coupling $\kappa_{\mathrm{exact}}$ is the dedicated target of Track 4 / Tasks `T12`--`T15`).
       - In **self-consistent 4-Majorana HFB mean-field theory**, there are no non-local $\mathbb{Z}_2$ flux flips. Instead, the bare Zeeman operator $-\sum_\gamma h_\gamma N^\gamma$ only couples the matter flavor $c^0$ to the gauge flavors $c^\gamma$ ($N^\gamma_{0\gamma} = -1$). Because the unperturbed $h=0$ QSL fixed point is block-diagonal between the dispersive matter band ($c^0$) and the three localized flat gauge bands ($c^x, c^y, c^z$ with two-particle energy $2\Delta_{\mathrm{gauge}} = 2|K|\bar{U} \approx 1.04973 |K|$), a bare perturbation $-\sum_\gamma h_\gamma N^\gamma$ alone **cannot** generate a $\mathcal{T}$-odd chiral mass at $\mathcal{O}(h^3)$ without breaking the gauge-sector degeneracy! Rather, the chiral mass at $\mathbf{K}$ arises from the **self-consistent KKT linear response**: at $\mathcal{O}(h)$, the induced on-site expectation $V_s$ triggers both the linear KKT Lagrange multiplier counter-field $\lambda_s^\gamma = \mathcal{O}(h)$ (which couples directly to the dual gauge generators $G^\gamma$ having matrix elements $G^\gamma_{\alpha\beta} = -\epsilon_{\gamma\alpha\beta}$ that mix different gauge flavors $c^\alpha \leftrightarrow c^\beta$) and the off-diagonal Fock bond self-energy $\delta U_b^{0\gamma}, \delta U_b^{\alpha\beta} = \mathcal{O}(h)$.
       - Consequently, in the third-order path $c_A^0 \xrightarrow{\mathcal{O}(h)} c^\alpha \xrightarrow{\mathcal{O}(h)} c^\beta \xrightarrow{\mathcal{O}(h)} c_A^0$ around the Dirac point $\mathbf{K}$, the middle vertex is dressed by the self-consistent KKT multiplier $\lambda^\gamma G^\gamma$ and Fock bond self-energy $\delta U_b$, rather than a bare Zeeman vertex. Expanding the KKT stationarity condition $[\mathcal{H}(\mathbf{k}), \Pi^-(\mathbf{k})] = 0$ analytically to $\mathcal{O}(h^3)$ is the formal target of `Contract H6` (Task `T08`).

---

## 3. Resolution Summary

1. **`doc/labbook/chapters/ch_m1_6_solvers_and_basins.tex`:**
   - Demoted `Theorem 7.2` (`thm:m1_6_ab_failure`) to `\begin{proposition}[Failure of Unprojected Iteration under Krylov Subspaces]`.
   - Demoted `Theorem 7.3` (`thm:m1_6_critical_transition`) to `\begin{proposition}[First-Order Thermodynamic Transition and Metastability Window]`.
   - Split `Theorem 7.4` (`thm:m1_6_cubic_gap`) into:
     - `\begin{theorem}[Algebraic Exponent Selection Rule for Dirac Mass Opening]` (`\label{thm:m1_6_cubic_gap}`) containing **only** the exact symmetry selection rules ($\mathcal{O}(h) = 0$, $\mathcal{O}(h^2) = 0$, and $\mathcal{I}^{(3)}(\mathbf{h}) = h_x h_y h_z = \frac{h^3}{3\sqrt{3}}$) and its pure symmetry `\begin{proof}`, and
     - `\begin{proposition}[Empirical Cubic Dirac Mass Scaling in 4-Flavor Mean-Field Theory]` (`\label{prop:m1_6_cubic_gap_fit}`) housing the numerical regression $(2.5921 \pm 0.0034)(h/|K|)^{2.9936 \pm 0.0006}$, $R^2 > 0.99999$, and global Brillouin-zone minimum check.
   - Rewrote Section `\ref{sec:m1_6_perturbative_reconciliation}` to state the exact operator-level distinction between the Feshbach two-vison + matter continuum resolvent and the HFB KKT self-consistent perturbation expansion (`Contract H6` / `T08`), eliminating the hand-waved `3.63 -> 10.37` arithmetic.
2. **`doc/labbook/chapters/ch_m1_8_supercell_bloch_scaling.tex`:**
   - Demoted `Theorem 8.2` (`thm:m1_7_stability`) to `\begin{proposition}[Empirical Ground-State Isotropy and Stability Against Modulation]`.
3. **`doc/rnd/phase_1/m1_10_t1_exact_perturbation_and_normalization_reconciliation.md` & `doc/manuscript/sections/sec03_field_response_and_multistability.tex`:**
   - Added complete Rule 03 YAML frontmatter to `m1_10_t1_exact_perturbation_and_normalization_reconciliation.md`.
   - Replaced line 99 (`3.63 h^3 ... 10.37 h^3`) and the corresponding remark in `sec03_field_response_and_multistability.tex` with the rigorous operator-level comparison between the exact Feshbach two-vison + matter continuum resolvent (`T12`--`T15`) and the HFB KKT self-consistent perturbation pathway (`T08` / `Contract H6`).
