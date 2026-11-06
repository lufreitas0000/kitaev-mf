# SciPost Physics Referee Report 2 (Round 2)

**Manuscript Title:** Majorana Parton Mean-Field Theory of the Extended Kitaev Honeycomb Model: Algebraic Symmetries, Grassmannian Solvers, and Kinematic Invariants
**Author:** Lu Freitas
**Reviewer:** Referee 2 (Expert Computational Physicist / Numerical Analyst)
**Date:** November 6, 2026
**Recommendation:** **Accept with Minor Typographical Reconciliations** (Clear for publication upon harmonizing introductory/concluding summary sentences with Section 4)

---

## 1. Summary of Revisions and Referee Evaluation

In Round 1, I raised severe objections regarding:
1. Significant arithmetic undercounting of step FLOPs in the solver scorecard (Good Broyden, Joint DIIS, JFNK).
2. The erroneous assertion of a "theoretical $4.0\times$ FLOP reduction" in Track A real skew-symmetric eigensolves, which completely ignored three dense $M \times M$ matrix multiplications.
3. Misleading claims that matrix sign purification "breaks the $\mathcal{O}(M^3)$ scaling bottleneck," while running an algorithm that executes $9\times$ *more* FLOPs than standard eigensolvers.
4. Total omission of the $L \equiv 0 \pmod 3$ Dirac point resonance, matrix condition number explosion ($\kappa(\mathbf{A}) \approx 6 \times 10^{15}$), condition squaring ($\kappa(\mathbf{S}) \approx 3.6 \times 10^{31}$), and mantissa bit loss in Track A.
5. Lack of explicit 95% Wilson score confidence intervals and unjustified upper-bound extrapolations ($R_{90} \ge 0.500$).
6. Attribution of the 44.6-second `zheevd` runtime at $M=1152$ to normal eigensolver scaling rather than Cuppen divide-and-conquer deflation stalling.
7. Table horizontal overflow (overfull `\hbox`).

The author has undertaken an exemplary, technically uncompromising overhaul of Section 4 and provided a detailed point-by-point response in `doc/rnd/peer_review/round2_author_response.md`. I have re-audited the manuscript source files and checked the TeX compilation logs. My evaluation of the seven points is detailed below.

---

## 2. Point-by-Point Verification of Required Numerical Corrections

### 1. Correction of Solver Arithmetic Complexity (Table 1 & Section 4.1)
- **Good Broyden:** The step overhead has been corrected from the erroneous $2 D^2$ to the mathematically exact **$8 D^2 \text{ FLOPs}$**, properly accounting for $w_k = H_k y_k$ ($2 D^2$), $v_k^T = s_k^T H_k$ ($2 D^2$), the rank-1 outer product update ($2 D^2$), and step application ($2 D^2$).
- **Joint DIIS (Anderson Acceleration):** The step overhead was corrected from $\mathcal{O}(m^2 D)$ to **$4 m D + \mathcal{O}(m^3) \text{ FLOPs}$**, strictly reflecting the $2 m D$ dot-product update of the Gram matrix, $2 m D$ trial vector combination, and small $(m+1) \times (m+1)$ linear solve.
- **JFNK:** The linear algebra Krylov step overhead ($2 k_{\mathrm{kry}} D$) has been disentangled from the physical mean-field map evaluation calls ($k_{\mathrm{kry}} \Phi$).
- **Statistical Intervals:** Table 1 now reports explicit 95% Wilson confidence intervals $[w^-, w^+]$ at $\delta = 0.35$ (e.g., $[0.886, 1.000]$ for Broyden, DIIS, and JFNK; $[0.192, 0.512]$ for L-Broyden; $[0.006, 0.167]$ for Picard). The basin radius is honestly bounded at the tested radius $R_{90} = 0.350$.

### 2. Retraction of the Track A Real Eigensolve FLOP Fallacy (Section 4.3)
The author has fully retracted the false "4.0x FLOP reduction" claim. In Section 4.3 (item 2), the arithmetic cost is derived with exact prefactors:
$$\mathcal{W}_{\mathrm{Track\,A}} = \frac{8}{3} M^3 + 6 M^3 = \frac{\mathbf{26}}{\mathbf{3}} M^3 \approx \mathbf{8.67 M^3 \text{ FLOPs}}.$$
The author explicitly notes that Track A requires $1.63\times$ *more* arithmetic FLOPs than direct complex diagonalization via `zheevd` ($\frac{16}{3} M^3 \approx 5.33 M^3$). The observed wall-clock speedups ($1.36\times\text{--}6.5\times$) are now correctly attributed to Level-3 BLAS cache-line streaming and lack of complex branch stalls, rather than an arithmetic complexity advantage.

### 3. Matrix Sign Purification: Algorithmic Complexity vs. Hardware Roofline
The author has completely eliminated claims of "breaking the $\mathcal{O}(M^3)$ scaling bottleneck." Section 4.4 explicitly states:
- Newton-Schulz purification consumes $4 k_{\mathrm{iter}} M^3 \text{ FLOPs}$, which at $k_{\mathrm{iter}} = 12$ equals **$48 M^3 \text{ FLOPs}$** ($9.0\times$ more arithmetic than `zheevd`).
- The speedup is properly explained through arithmetic intensity and hardware roofline dynamics: dense matrix multiplication (`dgemm`) achieves sustained throughput of **$114.64\text{ GFLOP/s}$** (67% of theoretical hardware peak), completely avoiding memory stalls.

### 4. Divide-and-Conquer Eigensolver Secular Deflation Stalling
The anomalous single-run 44.6 s runtime at $M=1152$ in Round 1 has been replaced by genuine multi-repeat benchmark statistics ($N_{\mathrm{repeats}} = 5$). In Table 2:
- At $M=1152$, `zheevd` takes $1418.6 \pm 420\text{ ms}$, while purification takes $613.0 \pm 42\text{ ms}$, representing a realistic and credible speedup of **$2.31\times$**.
- At $M=1568$, `zheevd` takes $3581.9 \pm 187\text{ ms}$ vs $1748.7 \pm 35\text{ ms}$ for purification (**$2.05\times$ speedup**).
- The text explicitly notes that the high run-to-run variance in `zheevd` ($\pm 420\text{ ms}$) is driven by secular equation deflation checks failing across multi-threaded workers on high-symmetry $C_{6v}$ degenerate eigenvalue multiplets.

### 5. Dirac Resonance ($3 \mid L$), Condition Number Squaring, and Mantissa Loss
Section 4.3 (item 3) now comprehensively documents the Dirac point resonance:
- At $3 \mid L$, momentum discretization samples the Dirac node $\mathbf{K} = (\frac{1}{3}, \frac{2}{3})$, causing $s_{\min}(\mathbf{A}) \approx 10^{-16}$ and condition number explosion to $\kappa(\mathbf{A}) \approx 6 \times 10^{15}$.
- Track A squares the condition number to $\kappa(\mathbf{S}) \approx 3.6 \times 10^{31} \gg \epsilon_{\mathrm{mach}}^{-1}$, losing 8 mantissa digits and spiking Frobenius discrepancy to $9.5 \times 10^{-9}$.
- Track B purification acts as a smooth spectral sign filter rather than an inversion, maintaining machine precision ($\approx 4 \times 10^{-15}$) across all geometries.

### 6. Table Layout and Overfull Boxes
Both Table 1 and Table 2 have been reformatted using `\resizebox{\columnwidth}{!}{...}` with `\scriptsize`. Independent verification of the LaTeX compilation transcript reveals:
- **0 LaTeX compilation errors**
- **0 Overfull `\hbox` warnings**
- **0 Underfull `\hbox` warnings**

---

## 3. Remaining Minor Editorial Reconciliation (Mandatory Before Publication)

While Section 4 has been completely corrected and now represents a gold standard of numerical rigor, the author appears to have overlooked synchronizing the summary text in the Introduction and Conclusion:

1. **Section 1 (Introduction, line 18):**
   The text still reads:
   *"...halving memory ($2.0\times$) and cutting FLOPs by $4.0\times$. Finally, we implement cubic Newton-Schulz matrix sign purification... achieving a $63.53\times$ speedup..."*
   **Required Edit:** Update this sentence to match Section 4: remove the debunked "cutting FLOPs by $4.0\times$" and replace the obsolete "$63.53\times$ speedup" with the verified sustained benchmark speedup ($2.05\times\text{--}2.31\times$, reaching $114.6\text{ GFLOP/s}$ at 67% CPU hardware roofline peak).

2. **Section 5 (Conclusion, line 13):**
   Item 5 still reads:
   *"...reduced FLOPs by $4.0\times$. Newton-Schulz matrix sign purification delivered a $63.53\times$ speedup at $M=1152$..."*
   **Required Edit:** Harmonize with Section 4: state that memory is halved ($2.0\times$), while purification achieves sustained throughput of $114.64\text{ GFLOP/s}$ and outpaces divide-and-conquer eigensolvers by up to $2.31\times$ by bypassing secular deflation stalling.

3. **Section 4 Editorial Comment (line 1):**
   Remove the LaTeX comment `% EDITORIAL QUESTION: ... how does eigenvector-free purification break the O(M^3) scaling bottleneck?` to maintain consistency with the main text.

---

## 4. Final Recommendation

**Accept with Minor Typographical Reconciliations.**
The revised manuscript is numerically sound, mathematically transparent, and computationally rigorous. The author has addressed all core numerical analyst concerns in an exemplary manner. Once the minor text synchronizations in Section 1 and Section 5 are applied, the manuscript will be fully suitable for publication in *SciPost Physics*.
