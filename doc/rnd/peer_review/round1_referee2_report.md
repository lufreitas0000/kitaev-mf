# SciPost Physics Referee Report 2

**Manuscript Title:** Majorana Parton Mean-Field Theory of the Extended Kitaev Honeycomb Model: Algebraic Symmetries, Grassmannian Solvers, and Kinematic Invariants
**Author:** Lu Freitas
**Reviewer:** Referee 2 (Expert Computational Physicist / Numerical Analyst)
**Date:** November 6, 2026
**Recommendation:** Major Revision (Reject unless all computational and numerical discrepancies are rigorously reconciled)

---

## 1. Overview and Executive Summary

In this manuscript, the author formulates the self-consistent Majorana mean-field equations of the extended Kitaev honeycomb model as a geometric fixed-point problem on the orthogonal Grassmannian bundle. The manuscript presents an empirical comparison of non-linear fixed-point solvers (Picard, DIIS, Full Broyden, Limited-Memory Broyden, and JFNK) and introduces real skew-symmetric arithmetic and cubic Newton-Schulz matrix sign purification for real-space torus clusters.

While the mathematical formulation on the Grassmannian bundle and the application of matrix sign purification to Majorana mean-field theory are timely and of clear interest to the computational quantum materials community, the manuscript in its current form contains **severe arithmetic errors, undercounted FLOP formulas, unsubstantiated scaling claims, and omissions of critical numerical failure modes**. Most egregiously:
1. The author claims that real skew-symmetric eigensolves achieve a "theoretical $4.0\times$ FLOP reduction" over complex Hermitian eigensolves, completely ignoring three dense $M \times M$ matrix multiplications that render the real algorithm *arithmetically more expensive* than LAPACK's complex `zheevd`.
2. The manuscript repeatedly asserts that cubic Newton-Schulz purification "breaks the $\mathcal{O}(M^3)$ scaling bottleneck," which is mathematically false: the algorithm is strictly $\mathcal{O}(M^3)$ and executes nearly an order of magnitude *more* floating-point operations than standard eigensolvers, succeeding purely due to BLAS-3 cache locality and memory bandwidth roofline saturation.
3. The manuscript promises 95% Wilson score confidence intervals for solver basin reliability, yet Table 1 contains no confidence intervals whatsoever.
4. The manuscript completely ignores the catastrophic $L \equiv 0 \pmod 3$ Dirac point resonance, matrix condition number explosion ($\kappa(\mathbf{A}) \approx 6 \times 10^{15}$), mantissa precision degradation, zero-mode projector incompleteness ($\Omega^2 = -\mathbf{I} + \mathbf{P}_0$), and LAPACK Cuppen divide-and-conquer deflation stalling.

The manuscript cannot be recommended for publication in *SciPost Physics* until the author addresses the following formal, numbered critiques and incorporates all requested benchmarks.

---

## 2. Detailed Technical Critiques

### Critique 1: Erroneous Arithmetic Complexity and Under-Counted Solver FLOP Formulas

In Section 4.1, 4.3, and Table 1, the arithmetic complexity formulas for the non-linear fixed-point accelerators are riddled with errors and inconsistencies:

1. **Good Broyden Step Overhead Undercounting:**
   In Table 1, the step overhead for Full Broyden is stated as $2 D^2 \text{ FLOPs}$. This is an egregious undercount.
   The Good Broyden inverse Jacobian update is:
   $$H_{k+1} = H_k + \frac{(s_k - H_k y_k) s_k^T H_k}{s_k^T H_k y_k}.$$
   Evaluating the denominator and numerator requires:
   - Computing $w_k = H_k y_k$: matrix-vector product requiring $2 D^2$ FLOPs.
   - Computing $v_k^T = s_k^T H_k$: vector-matrix product requiring $2 D^2$ FLOPs.
   - Inner product $s_k^T w_k = v_k^T y_k$: $2 D$ FLOPs.
   - Rank-1 outer product update $H_{k+1} = H_k + u_k v_k^T$: $D^2$ multiplications and $D^2$ additions $\implies 2 D^2$ FLOPs.
   - Applying $H_k R_k$ for the quasi-Newton step: $2 D^2$ FLOPs.
   The actual leading-order operational overhead is **$8 D^2 \text{ FLOPs}$**, not $2 D^2$. Table 1 undercounts the linear algebra overhead by a factor of 4.

2. **Incorrect DIIS (Anderson Acceleration) Complexity Scaling:**
   In Table 1, the step overhead for Joint DIIS is listed as $\mathcal{O}(m^2 D)$. This is mathematically incorrect.
   In DIIS with history depth $m$, the Gram matrix $B_{pq} = \langle G_p, G_q \rangle$ is updated at step $k$ by computing $m$ vector inner products between the new residual $G_k$ and the $m$ stored residuals, requiring exactly $2 m D \text{ FLOPs}$.
   The $(m+1) \times (m+1)$ KKT linear system is solved via Cholesky or LU factorization in $\frac{1}{3} m^3$ or $\frac{2}{3} m^3 \text{ FLOPs}$.
   The new trial state is formed as a linear combination of $m$ history vectors: $\sum_{j=0}^{m-1} c_j X_{k-j}$, requiring $2 m D \text{ FLOPs}$.
   Therefore, the total DIIS operational overhead is strictly:
   $$\mathcal{W}_{\mathrm{DIIS}} = 4 m D + \mathcal{O}(m^3) \text{ FLOPs}.$$
   There is **no** $m^2 D$ leading-order term in standard DIIS. Listing $\mathcal{O}(m^2 D)$ misrepresents the complexity of Anderson mixing.

3. **Conflation of Solver Overhead and Evaluation Cost in JFNK:**
   In Table 1, under JFNK, the author enters $\mathcal{O}(k_{\mathrm{kry}} \Phi)$ in the "Step Overhead" column. This conflates the algorithmic overhead of the Krylov accelerator with the physical evaluation of the mean-field map $\Phi$. The actual linear algebra step overhead of GMRES is the modified Gram-Schmidt Arnoldi orthogonalization, which scales as $\mathcal{O}(k_{\mathrm{kry}}^2 D)$, whereas $k_{\mathrm{kry}}$ evaluations of $\Phi$ represent the function evaluation workload. Table 1 must separate internal solver overhead from external operator calls.

---

### Critique 2: Mathematical Invalidation of the Claimed $4.0\times$ Real Eigensolve FLOP Reduction (Track A)

In Section 4.3 (lines 40–41), the author writes:
> *"Real symmetric eigensolves (`dsyevd`) consume $\frac{4}{3} M^3 \text{ FLOPs}$ compared to $\frac{16}{3} M^3 \text{ FLOPs}$ for complex `zheevd`, yielding a theoretical $4.0\times$ FLOP reduction and measured wall-clock speedups between $1.36\times$ and $6.5\times$."*

This statement is mathematically false and constitutes an unphysical idealization. The author compares only the isolated tridiagonalization step of a real matrix against a complex matrix, completely omitting the operations required to execute the mean-field map in Track A:
1. To obtain a real symmetric matrix from the real skew-symmetric Hamiltonian $\mathbf{A} \in \mathbb{R}^{M \times M}$, one must compute $\mathbf{S} = -\mathbf{A}^2 = \mathbf{A} \mathbf{A}^T$. A dense matrix multiplication of two $M \times M$ matrices consumes **$2 M^3 \text{ FLOPs}$**.
2. Computing the eigenvalues and eigenvectors of $\mathbf{S}$ via LAPACK `dsyevd` (divide-and-conquer) consumes $\frac{4}{3} M^3 \text{ FLOPs}$ for tridiagonal reduction plus $\frac{4}{3} M^3 \text{ FLOPs}$ to accumulate the orthogonal eigenvectors $\mathbf{V}$, totaling **$\frac{8}{3} M^3 \text{ FLOPs}$**.
3. Reconstructing the inverse square root $\mathbf{S}^{-1/2} = \mathbf{V} \boldsymbol{\Lambda}^{-1/2} \mathbf{V}^T$ requires multiplying $\mathbf{V}$ by $\boldsymbol{\Lambda}^{-1/2} \mathbf{V}^T$, consuming another **$2 M^3 \text{ FLOPs}$**.
4. Finally, forming the physical covariance matrix $\Omega = -\mathbf{A} \mathbf{S}^{-1/2}$ requires a third dense matrix multiplication, consuming another **$2 M^3 \text{ FLOPs}$**.

Summing these leading-order contributions yields:
$$\mathcal{W}_{\mathrm{Track\,A}} = 2 M^3 + \frac{8}{3} M^3 + 2 M^3 + 2 M^3 = \frac{26}{3} M^3 \approx \mathbf{8.67 M^3 \text{ FLOPs}}.$$
In contrast, computing the spectral decomposition of $\ii \mathbf{A}$ directly using complex LAPACK `zheevd` consumes:
$$\mathcal{W}_{\mathbb{C}} = \frac{16}{3} M^3 \approx \mathbf{5.33 M^3 \text{ FLOPs}}.$$
**Conclusion:** Track A real skew arithmetic requires **$1.63\times$ MORE arithmetic FLOPs** than direct complex diagonalization!
Indeed, the author's own laboratory record (`doc/labbook/chapters/ch_m1_6_solvers_and_basins.tex`, line 205) documents that the real skew eigensolve costs $(\frac{4}{3} M^3 + 6 M^3) = \frac{22}{3} M^3 \approx 7.33 M^3 \text{ FLOPs}$. Suppressing the $6 M^3$ matrix multiplication overhead in the submitted manuscript while claiming a "theoretical $4.0\times$ FLOP reduction" is an unscientific misrepresentation of numerical reality.

---

### Critique 3: Misleading Scaling Claims for Matrix Sign Purification (Track B)

In Section 4.4 and Section 5.1 (item 5), the author claims that cubic Newton-Schulz matrix sign purification "breaks the $\mathcal{O}(M^3)$ scaling bottleneck."

1. **Newton-Schulz is an $\mathcal{O}(M^3)$ Algorithm:**
   The iteration $\mathbf{Y}_{k+1} = \frac{1}{2} \mathbf{Y}_k (3\mathbf{I} + \mathbf{Y}_k^2)$ requires:
   - Computing $\mathbf{Y}_k^2$: for skew-symmetric $\mathbf{Y}_k$, this is a matrix multiplication costing $2 M^3 \text{ FLOPs}$.
   - Computing $\mathbf{Y}_k (3\mathbf{I} + \mathbf{Y}_k^2)$: a second matrix multiplication costing $2 M^3 \text{ FLOPs}$.
   Thus, each iteration consumes **$4 M^3 \text{ FLOPs}$**.
   For the benchmarked $M=1152$ system, the algorithm requires $k_{\mathrm{iter}} = 12$ iterations, consuming:
   $$\mathcal{W}_{\mathrm{purif}} = 12 \times 4 M^3 = \mathbf{48 M^3 \text{ FLOPs}}.$$
   Purification performs **$9.0\times$ more arithmetic FLOPs** than complex `zheevd` ($\frac{16}{3} M^3$)!

2. **Algorithmic Complexity vs. Hardware Roofline Throughput:**
   The observed $63.53\times$ wall-clock speedup does *not* arise from breaking asymptotic $\mathcal{O}(M^3)$ complexity. It arises exclusively from **arithmetic intensity and hardware roofline saturation**:
   - `dgemm` operations operate on regular cache-blocked memory layouts, saturating CPU SIMD vector units at $> 104\text{ GFLOP/s}$ (near peak theoretical floating-point performance).
   - Divide-and-conquer eigensolvers (`zheevd` / `dlaed0`) are bottlenecked by scalar secular equation root-finding, irregular memory access, and severe multi-threaded lock contention, dropping effective throughput to $20\text{--}23\text{ GFLOP/s}$.
   The manuscript must replace all false assertions about "breaking $\mathcal{O}(M^3)$ scaling" with a precise discussion of BLAS-3 cache locality versus BLAS-1/2 memory stalls.

---

### Critique 4: Total Omission of the $L \equiv 0 \pmod 3$ Dirac Resonance, Condition Number Explosion, and Mantissa Loss

Section 4 completely omits one of the most critical numerical phenomena governing real-space honeycomb simulations: the **$L \equiv 0 \pmod 3$ Dirac point resonance**.

1. **Condition Number Explosion:**
   On a periodic $L \times L$ torus, the fractional Brillouin zone Dirac points $\mathbf{K} = (\frac{1}{3}, \frac{2}{3})$ and $\mathbf{K}' = (\frac{2}{3}, \frac{1}{3})$ are sampled if and only if $3 \mid L$. When $3 \mid L$ ($L \in \{6, 9, 12, \dots\}$), the smallest singular value of the skew-symmetric Hamiltonian plunges to $s_{\min} \approx 10^{-16}$, causing the condition number to explode to:
   $$\kappa(\mathbf{A}) = \frac{s_{\max}(\mathbf{A})}{s_{\min}(\mathbf{A})} \approx \mathbf{6 \times 10^{15}}.$$

2. **Catastrophic Precision Loss in Track A:**
   In Track A, forming $\mathbf{S} = -\mathbf{A}^2$ squares the condition number:
   $$\kappa(\mathbf{S}) = (\kappa(\mathbf{A}))^2 \approx \mathbf{3.6 \times 10^{31}} \gg \epsilon_{\mathrm{mach}}^{-1} \approx 10^{16}.$$
   Inverting $\mathbf{S}^{1/2}$ near this singular boundary causes catastrophic cancellation, losing 8 digits of precision (the entire lower half of the 53-bit IEEE 754 double-precision mantissa). As recorded in the author's internal labbook (Figure F18), the Frobenius error of Track A explodes from $1.3 \times 10^{-14}$ (for $3 \nmid L$) to $9.5 \times 10^{-9}$ (for $3 \mid L$).

3. **Suppression in the Manuscript:**
   None of this appears in Section 4 of the manuscript! The author presents Track A as a viable method without warning the reader that it fails with catastrophic precision degradation whenever the cluster samples the physical Dirac nodes.

---

### Critique 5: Small-Gap Slowdown and Incompleteness of Matrix Sign Purification at Zero Modes

While Track B (Newton-Schulz) avoids the condition number squaring of Track A, the author's depiction of purification as an unconditional machine-precision panacea is misleading:

1. **Linear Slowing Down Near Gapless Nodes:**
   The scalar iteration $f(x) = \frac{1}{2} x (3 - x^2)$ has an unstable fixed point at $x=0$ with derivative $f'(0) = \frac{3}{2}$. In the small-gap regime ($\lambda_0 = \Delta_{\min} / \|\mathbf{A}\|_2 \ll 1$), cubic convergence does not begin until the iterate enters the basin $[0.5, 1.0]$. The initial phase is governed by slow linear expansion:
   $$k_{\mathrm{linear}} \approx \frac{\ln(1/\lambda_0)}{\ln(1.5)} \approx 2.466 \ln\left( \frac{\|\mathbf{A}\|_2}{\Delta_{\min}} \right).$$
   Near a gapless Dirac resonance where $\Delta_{\min} \sim 10^{-4}$ or $10^{-6}$, the required iterations double or triple ($k \approx 25\text{--}40$), consuming $100\text{--}160 M^3 \text{ FLOPs}$ and substantially eroding the wall-clock speedup. The manuscript reports benchmark results only at the safe, gapped point $h=0.20$ ($\Delta_K \approx 0.02, k=11\text{--}12$) without documenting this gap-dependent iteration penalty.

2. **Projector Incompleteness at Exact Zero Modes ($\Delta = 0$):**
   At zero magnetic field ($h=0$) or in topological defect sectors (such as pairs of $\mathbb{Z}_2$ vison vortices hosting Majorana zero modes), the Hamiltonian possesses exact zero eigenvalues $\epsilon_0 = 0$. Because $\operatorname{sgn}(0) = 0$, the cubic iteration yields:
   $$\lim_{k \to \infty} \mathbf{Y}_k = \Omega, \qquad \Omega^2 = -\mathbf{I} + \mathbf{P}_{\mathrm{zero}},$$
   where $\mathbf{P}_{\mathrm{zero}}$ is the orthogonal projector onto the zero-mode null subspace. The purified covariance matrix is **not** an orthogonal complex structure and fails to define a pure Gaussian state!
   The manuscript does not describe how physical fermion parity projection $(-1)^F = \pm 1$ is resolved (e.g., via skew Pfaffian deflation in the null space). Asserting that eigensolvers can be "completely bypassed" without addressing zero-mode projection is unacceptable.

---

### Critique 6: Pathology of Cuppen Eigensolver Deflation Stalling

In Section 4.4, Table 2 lists an enormous wall-clock runtime for the divide-and-conquer eigensolver at $M=1152$: **$44,596\text{ ms}$** (44.6 seconds), compared to $702\text{ ms}$ for purification.

The author treats this 44.6-second runtime as a baseline for eigensolver scaling, but fails to explain the underlying numerical pathology:
- On an $L=12$ honeycomb torus, the Hamiltonian possesses non-abelian crystallographic point-group symmetries ($C_{6v}$ and translations), leading to massive eigenvalue degeneracies (up to 12-fold degenerate multiplets).
- In Cuppen's divide-and-conquer algorithm (`dlaed0`/`zlaed0`), degenerate eigenvalues should theoretically deflate. However, due to parallel floating-point roundoff errors across worker threads, clustered eigenvalues differ by $\delta \sim 10^{-14}$.
- The deflation tolerance test $|d_i - d_j| \le \text{tol}$ fails, forcing the secular equation root-finder (`dlaed4`) into pathological stagnation while resolving nearly colliding poles. Thread barriers and lock contention cause wall-clock runtimes to spike erratically (from $\sim 400\text{ ms}$ at $M=800$ to $>44,000\text{ ms}$ at $M=1152$).
- Presenting this unmitigated deflation stalling as standard eigensolver scaling without explanation or comparison to robust algorithms (such as MRRR / `zheevr` or randomized block Lanczos) is scientifically misleading.

---

### Critique 7: Missing Wilson Confidence Intervals and Stagnation in Solver Comparisons

1. **Unfulfilled Promise of Confidence Intervals:**
   In Section 4.2 (lines 16–17), the text explicitly promises:
   > *"Point estimates $\hat{p}$ are reported with 95\% Wilson score confidence intervals:"*
   Yet in Table 1, **not a single confidence interval is reported!** The table lists only point values for $R_{90}$ and $R_{50}$ (e.g., "$0.350$", "$0.010$"). The Wilson bounds $[w^-, w^+]$ for each solver across the tested perturbation radii $\delta$ are completely missing.

2. **Unverified Upper Bound Extrapolations:**
   Table 1 lists $R_{90} \ge 0.500$ for Full Broyden and JFNK. However, $0.500$ was simply the maximum perturbation radius tested ($\delta \in [0.005, 0.500]$). In numerical benchmarking, one cannot claim a lower bound "$\ge 0.500$" without testing beyond the search boundary to identify where the basin actually breaks down.

3. **Misattribution of Picard Failure:**
   The author reports that Picard takes $K_{50} = 115$ iterations and stagnates at $\sim 10^{-6}$. In standard self-consistent field literature, damped Picard operates on a dual manifold where the gauge Lagrange multipliers $\boldsymbol{\lambda}$ are updated via dual gradient ascent ($\delta \boldsymbol{\lambda} = \eta_\lambda \mathbf{g}$). The manuscript fails to disclose the value of $\eta_\lambda$ (internally fixed to $\eta_\lambda = 0.2$). Picard's slow convergence is primarily caused by an ill-conditioned, unoptimized step size $\eta_\lambda$ on the dual constraint, rather than an intrinsic failure of fixed-point contraction.

---

## 3. Required Revisions and Mandatory Benchmarks

To meet the technical standards of *SciPost Physics*, the author must execute the following revisions before this manuscript can be considered:

1. **Correct all FLOP formulas and Table 1 overhead entries:**
   - Update Good Broyden step overhead to **$8 D^2 \text{ FLOPs}$** (accounting for $H_k y_k$, $s_k^T H_k$, rank-1 update, and step application).
   - Update Joint DIIS step overhead to **$4 m D + \mathcal{O}(m^3) \text{ FLOPs}$**, eliminating the incorrect $\mathcal{O}(m^2 D)$ term.
   - Separate solver linear algebra overhead from operator evaluation calls ($\Phi$) in Table 1.

2. **Correct the Track A FLOP narrative:**
   - Retract the claim of a "theoretical $4.0\times$ FLOP reduction" in Section 4.3.
   - State the full, honest arithmetic complexity of Track A: $\mathcal{W}_{\mathrm{Track\,A}} = \frac{26}{3} M^3 \text{ FLOPs}$, proving that Track A is $1.63\times$ more expensive in arithmetic FLOPs than complex `zheevd` ($\frac{16}{3} M^3$).

3. **Reconcile Matrix Sign Purification Complexity vs Throughput:**
   - Retract the claim that purification "breaks the $\mathcal{O}(M^3)$ scaling bottleneck."
   - Explicitly report that purification consumes $48 M^3 \text{ FLOPs}$ ($9\times$ more than eigensolvers) and attribute the $63.53\times$ speedup to hardware cache blocking, SIMD vectorization, and roofline throughput ($104.5\text{ GFLOP/s}$ vs. $20.2\text{ GFLOP/s}$).

4. **Incorporate the $L \equiv 0 \pmod 3$ Resonance and Track A Failure:**
   - Add a dedicated subsection in Section 4 analyzing the Dirac point resonance at $3 \mid L$.
   - Include Figure F18 from the laboratory records showing the condition number explosion $\kappa(\mathbf{A}) \approx 6 \times 10^{15}$, condition number squaring $\kappa(\mathbf{S}) \approx 3.6 \times 10^{31}$, and the loss of 8 mantissa digits in Track A.

5. **Document Small-Gap Scaling and Zero-Mode Parity Projection:**
   - Provide the logarithmic iteration scaling formula $k \approx 2.466 \ln(\|\mathbf{A}\|_2 / \Delta_{\min})$ for small-gap regimes.
   - Explicitly state the incompleteness of the sign projector at zero modes ($\Omega^2 = -\mathbf{I} + \mathbf{P}_{\mathrm{zero}}$) and define the skew Pfaffian deflation procedure required to fix physical fermion parity $(-1)^F = \pm 1$.

6. **Explain Cuppen Deflation Stalling:**
   - Add a technical note explaining that the 44.6-second runtime of `zheevd` at $M=1152$ is driven by Cuppen deflation failure on symmetry-degenerate honeycomb stars, rather than standard cubic scaling.

7. **Report Exact Wilson Score Intervals in Table 1:**
   - Provide the explicit 95% Wilson confidence intervals $[w^-, w^+]$ for all solver basin survival rates across the tested $\delta$ grid.

---

## 4. Overall Recommendation

**Major Revision.**
The numerical techniques and physical insights presented in this manuscript are valuable, but the presentation in Section 4 falls far short of the rigorous standards expected in numerical physics. The arithmetic undercounting, misleading complexity claims, and suppression of known numerical failure modes must be thoroughly rectified before publication.
