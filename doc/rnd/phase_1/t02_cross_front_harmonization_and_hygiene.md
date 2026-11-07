---
id: T02
title: "Cross-Front Notation, FLOP Complexity, Statistical Speedup, and Hypothesis Registry Harmonization"
status: Completed
date: "2026-11-06"
questions:
  - "Are local fermion parity and spin/gauge operator conventions completely consistent between the Living Lab Book (Chapter 1) and the SciPost Monograph (Section 2)?"
  - "Do the Living Lab Book (Chapters 7, 9, 11), SciPost Monograph (Section 4), PROJECT_STATE.yaml (K-017, K-018, K-019), HYPO_CONTRACTS.md, and CHANGELOG.md report identical FLOP counts and N_repeats=5 statistical speedups?"
  - "Does every \\section{} across all 11 Lab Book chapters and 5 Manuscript sections carry an explicit % EDITORIAL QUESTION: tag, and are orphaned build files removed?"
depends_on:
  - T01
  - M1.8-T1
  - M1.8-T2
  - M1.11
contracts:
  - K-002
  - K-008
  - K-017
  - K-018
  - K-019
artifacts:
  - doc/manuscript/sections/sec02_algebraic_foundations_and_invariants.tex
  - doc/labbook/chapters/ch_m1_1_foundations_algebra.tex
  - doc/labbook/chapters/ch_m1_6_solvers_and_basins.tex
  - doc/labbook/chapters/ch_m1_9_real_space_and_purification.tex
  - doc/labbook/chapters/ch_m1_11_hypothesis_contracts.tex
  - doc/HYPOTHESIS_CONTRACTS.md
  - PROJECT_STATE.yaml
  - CHANGELOG.md
  - Makefile
---

# T02: Cross-Front Notation, FLOP Complexity, Statistical Speedup, and Hypothesis Registry Harmonization

## 1. Scope & Adversarial Audit Findings (`[W1]`--`[W5]`)

Following Atomic Task `T01`, the three-agent adversarial review (`Referee_1_Theory`, `Referee_2_Numerics`, `Editor_SciPost`) identified cross-front inconsistencies (`[W1]`--`[W5]`) between the Living Lab Book (`doc/labbook/`), the SciPost Executive Monograph (`doc/manuscript/`), and the governance records (`PROJECT_STATE.yaml`, `doc/HYPOTHESIS_CONTRACTS.md`, `CHANGELOG.md`):

1. **Defect `[W1]` (`PH-01` & `H6` Cross-Front Synchronization):**
   - In `ch_m1_11_hypothesis_contracts.tex` (line 38) and `doc/HYPOTHESIS_CONTRACTS.md` (line 15), `PH-01` referenced `K-019` instead of `K-008` (the dedicated invariant for $\|\boldsymbol{\lambda}^*(h)\|_2 \approx \sqrt{5} h$), and `doc/HYPOTHESIS_CONTRACTS.md` (lines 26, 48, 76--78) and `ch_m1_11_hypothesis_contracts.tex` (line 262) still contained remnants of the `2.85x` screening arithmetic eliminated in `T01`.
2. **Defect `[W2]` (Track A FLOP Accounting & Track B Statistical Speedup Harmonization in `K-017`, `K-018`):**
   - **Track A (`K-017`):** In `ch_m1_9_real_space_and_purification.tex` (`thm:m1_8_real_vs_complex_algebra`, lines 31--37), `PROJECT_STATE.yaml` (`K-017`), and `CHANGELOG.md` (line 23), Track A claimed a "$4.0\times$ FLOP reduction ($\frac{4}{3}M^3$ vs. $\frac{16}{3}M^3$)" by counting only the bare tridiagonal reduction of a real symmetric matrix and ignoring the $6M^3$ real matrix multiplications required to form $\mathbf{S} = -\mathbf{A}^2$, assemble $\mathbf{S}^{-1/2} = \mathbf{V}\boldsymbol{\Lambda}^{-1/2}\mathbf{V}^T$, and reconstruct $\Omega = -\mathbf{A}\mathbf{S}^{-1/2}$. As documented in `sec04_grassmannian_solvers_and_real_space_scaling.tex` (line 43), total Track A arithmetic is $\frac{8}{3}M^3 + 6M^3 = \frac{26}{3}M^3 \approx 8.67 M^3$ real FLOPs (or $\frac{4}{3}M^3 + 6M^3 = \frac{22}{3}M^3$ if eigenvectors are omitted), so Track A's $2.0\times$ advantage is strictly in **memory footprint** ($8M^2$ vs. $16M^2$ bytes) and real-valued Level-3 BLAS cache streaming rather than total FLOP count.
   - **Track B (`K-018`):** `sec04_grassmannian_solvers_and_real_space_scaling.tex` (Table 2, line 68) reports the $N_{\mathrm{repeats}}=5$ statistical benchmark at $M=1152$: `zheevd` $1418.6 \pm 420\text{ ms}$ vs. Newton-Schulz purification $613.0 \pm 42\text{ ms}$ ($\mathbf{2.31\times}$ mean speedup, $109.73\text{ GFLOP/s}$ at $M=1152$ and $114.64\text{ GFLOP/s}$ at $M=1568$, Frobenius error $4.01 \times 10^{-15}$), whereas `ch_m1_9_real_space_and_purification.tex` (line 118), `PROJECT_STATE.yaml` (`K-018`), `doc/HYPOTHESIS_CONTRACTS.md` (line 28), and `CHANGELOG.md` (line 24) still foregrounded a single-run $44596\text{ ms}$ secular-equation deflation stall ($63.53\times$).
   - **Full Broyden Step FLOPs:** `sec04` (Table 1, line 27) wrote `8 D^2 FLOPs` while `ch_m1_6` (Table 7.1, line 169 & Table 7.2, line 395) wrote `2 D^2 + 4D` / `2 D^2`. For the Sherman-Morrison update $H_{k+1} = H_k + \frac{s_k - H_k y_k}{s_k^T H_k y_k} s_k^T H_k$ plus step $H_{k+1} R_{k+1}$, the exact leading-order cost is three matrix-vector products ($H_k y_k$, $H_k^T s_k$, $H_{k+1} R_{k+1} \implies 3 \times 2D^2 = 6D^2$) plus one rank-1 outer product update ($2D^2$), totaling $6D^2 + \mathcal{O}(D)$ (or $4D^2 + \mathcal{O}(D)$ if $H_k R_k$ is reused). We harmonize both tables to $4D^2 + \mathcal{O}(D)$ (one matvec $H_k y_k$ + one transposed matvec $s_k^T H_k$ + one rank-1 outer update $2D^2$).
3. **Defect `[W3]` (Phase 2 Scope Alignment):**
   - Verified that `sec05_conclusion_and_outlook.tex` (lines 17--23), `AGENTS.md`, `PROJECT_STATE.yaml`, and `CHANGELOG.md` all consistently define Phase 2 as **Physics Production in the Pure Ferromagnetic Kitaev $[111]$ Model ($K = -1.0$)** (`T=0` bulk ground state, `T=0` isolated vortex pairs, `T>0` thermal covariance, and dynamical spin correlations) and Phase 3 as Antiferromagnetic Kitaev ($K > 0$).
4. **Defect `[W4]` (Parity Operator Sign Convention $\mathcal{D}_i = -c_i^0 c_i^1 c_i^2 c_i^3$):**
   - In `ch_m1_1_foundations_algebra.tex` (Theorem 1.4, line 346) and `src/kitaev_mf/domain/majorana.py` (`K-002`), the local fermion parity operator is $\mathcal{D}_i \equiv -c_i^0 c_i^1 c_i^2 c_i^3$, which satisfies $\mathcal{D}_i = +1$ on the physical odd-Dirac-fermion gauge-singlet subspace ($\bm{\tau}_i^2 = \frac{3}{2}(\Id - \mathcal{D}_i) = 0$, $\bm{\sigma}_i^2 = \frac{3}{2}(\Id + \mathcal{D}_i) = 3\Id$).
   - However, `sec02_algebraic_foundations_and_invariants.tex` (Eq. 4, line 24) erroneously wrote $\mathcal{D}_i \equiv +c_i^0 c_i^1 c_i^2 c_i^3$ and used $1/8$ prefactors in Eqs. (2)--(3) (`S_i^\alpha = -\frac{\ii}{8}\mathbf{c}_i^T M^\alpha \mathbf{c}_i`) that contradicted Eq. (4). We harmonize `sec02` to match `ch_m1_1` (`\mathcal{D}_i \equiv -c_i^0 c_i^1 c_i^2 c_i^3`, $\sigma_i^\alpha = \frac{\ii}{4}\mathbf{c}_i^T M^\alpha \mathbf{c}_i = \frac{\ii}{2}(c_i^0 c_i^\alpha - c_i^{\alpha^+} c_i^{\alpha^-})$, $\tau_i^\gamma = \frac{\ii}{4}\mathbf{c}_i^T G^\gamma \mathbf{c}_i = \frac{\ii}{2}(c_i^0 c_i^\gamma + c_i^{\gamma^+} c_i^{\gamma^-})$).
5. **Defect `[W5]` (`% EDITORIAL QUESTION:` Coverage and Orphaned `doc/manuscript/main.*` Cleanup):**
   - Added the missing `% EDITORIAL QUESTION:` header before `\section{Summary and Theoretical Transition}` in `ch_m1_1_foundations_algebra.tex` (line 494) so that 100% of `\section{}` headers across all 11 Lab Book chapters and 5 Manuscript sections have an explicit guiding question.
   - Removed the 7 orphaned build files `doc/manuscript/main.{aux,bbl,blg,log,out,pdf,toc}` and updated `Makefile` line 21 help text to reference `doc/manuscript/scipost_kitaev_mf.pdf`.

---

## 2. Guiding Questions & Verification Plan

- **Q1:** Does `sec02_algebraic_foundations_and_invariants.tex` reproduce the exact algebraic sign identities of `ch_m1_1_foundations_algebra.tex` and `test_local_majorana_representations.py` (`\mathcal{D}_i = -c_i^0 c_i^1 c_i^2 c_i^3 = +1$ on $\mathcal{H}_{\mathrm{phys}, i}$)?
- **Q2:** Are `K-017`, `K-018`, and `K-019` completely synchronized across `PROJECT_STATE.yaml`, `ch_m1_9`, `ch_m1_11`, `sec04`, `doc/HYPOTHESIS_CONTRACTS.md`, and `CHANGELOG.md`?
- **Verification:** Run `make verify && make labbook && make manuscript` and audit with `Editor_SciPost` and `Referee_2_Numerics`.
