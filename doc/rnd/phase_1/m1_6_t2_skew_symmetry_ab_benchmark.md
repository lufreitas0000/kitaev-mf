---
id: RND-P1-M1.6-T2
title: "Skew-Symmetry Handling A/B Testing: Round-off Leakage, JFNK Stability, and Retraction Necessity"
phase: 1
milestone: M1.6
task: T2
status: completed
claims:
  - CLM-M1-ALG-004
---

# Skew-Symmetry Handling A/B Testing: Round-off Leakage, JFNK Stability, and Retraction Necessity

## 1. Context and Motivation

In the 4-flavor Majorana mean-field theory, the on-site correlation matrix $V_s \in \mathfrak{so}(4)$ represents expectations of skew-symmetric fermion bilinears:
$$V_s^{\mu\nu} = \langle i c_{s}^\mu c_s^\nu \rangle - i \delta^{\mu\nu}.$$
Physical validity and the `MeanFieldState` invariant require strict real skew-symmetry:
$$V_s^T = -V_s \iff \|V_s + V_s^T\|_\infty = 0.$$
In Milestone M1.6, we tested whether explicit post-accelerator retraction ($V \leftarrow \frac{1}{2}(V - V^T)$) is an indispensable mathematical safeguard or redundant overhead, as formulated in **Hypothesis AH-04** (`CLM-M1-ALG-004`).

---

## 2. Experimental Design (A/B Testing Matrix)

We evaluated three skew-symmetry handling arms across all four solvers (Plain FP, Joint DIIS, Good Broyden, JFNK):
1. **Arm A (`raw`):** No output symmetrization ($V_{\mathrm{out}}$ taken directly from $\operatorname{Im} P_{ss}$) and no post-accelerator step retraction.
2. **Arm B (`symmetrize_output`):** Output symmetrization enforced at BZ integration ($V_s \leftarrow \frac{1}{2}(V_s - V_s^T)$), but no step retraction.
3. **Arm C (`retract_step`):** Both output symmetrization and post-accelerator step retraction.

Each arm was evaluated on:
- Fields: $h \in \{0.01, 0.05\}$ along $[111]$.
- Seeds: Exact QSL, perturbed QSL ($\delta = 0.01$), exact Polarized, and perturbed Polarized ($\delta = 0.01$).

---

## 3. Empirical Results and Benchmark Evidence

The benchmark data was exported to `data/processed/benchmarks/m1_6/ab_skew_retraction.dat` with run logs in `build/rnd/m1_6/runs_skew_ab.jsonl`.

### 3.1 Key Observations

1. **Failure of Arm A (`raw`) in Krylov Subspaces:**
   On perturbed states under $h = 0.05$, Arm A crashed inside JFNK:
   ```
   Failed raw jfnk 0.05 qsl_pert: On-site correlation matrix V must be skew-symmetric, max error: 3.5487e-12
   Failed raw jfnk 0.05 pol_pert: On-site correlation matrix V must be skew-symmetric, max error: 1.2580e-12
   ```
   During matrix-free GMRES iterations, finite-difference directional perturbations $J(x) v \approx \frac{G(x + \epsilon v) - G(x)}{\epsilon}$ accumulate floating-point asymmetry in $V$. Without symmetrization, this round-off error leaks into the single-particle Hamiltonian, destroying self-adjointness and violating the $10^{-12}$ domain assertion.

2. **Equivalence of Arm B and Arm C:**
   For all convergent runs across FP, DIIS, Broyden, and JFNK, Arm B (`symmetrize_output`) and Arm C (`retract_step`) produced **identical trajectories** to machine precision ($< 10^{-14}$):
   - Map evaluations ($n_\Phi$) were identical in every test.
   - Commutator stationarity residuals $\epsilon_{\mathrm{comm}}$ matched to all recorded digits.
   - Gauge residuals $\epsilon_{\mathrm{gauge}}$ matched to all recorded digits.

3. **Status of Hypothesis AH-04 (`CLM-M1-ALG-004`):**
   - **Supported:** Arm A is mathematically unstable under Krylov methods due to round-off leakage.
   - **Supported:** Arm B and Arm C are numerically equivalent. Enforcing skew-symmetrization at the correlator integration boundary ($V \leftarrow \frac{1}{2}(V - V^T)$) is necessary and sufficient to guarantee strict $\mathfrak{so}(4)$ invariance across all fixed-point algorithms.
