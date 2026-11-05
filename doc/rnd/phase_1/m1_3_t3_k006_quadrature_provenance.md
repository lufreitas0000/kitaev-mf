---
id: RND-P1-M1.3-T3
title: "Reproducible K-006 Brillouin-Zone Quadrature"
phase: 1
milestone: M1.3
task: T3
status: completed
claims:
  - CLM-M1-COR-006
artifacts:
  - data/processed/benchmarks/k006_quadrature.dat
---

# Note 04d: Reproducible K-006 Brillouin-Zone Quadrature

**Status:** M1.4.3b numerical provenance
**Implementation:** services/quadrature.py; scripts/k006_quadrature.py
**Tests:** test/unit/test_k006_quadrature.py

## 1. Integral and numerical method

The matter bond correlator magnitude is
\[
\bar U=\frac13\int_0^1\!dk_1\int_0^1\!dk_2\,
\left|1+e^{-2\pi i k_1}+e^{-2\pi i k_2}\right|.
\]
This is the continuous-integrand form of K-006. It is preferable numerically to the equivalent normalized expression containing the form-factor magnitude in the denominator: the absolute-value integrand is continuous, whereas the ratio is undefined at the Dirac points.

Use a midpoint grid \(k_j=(j+1/2)/N\). For \(N\) divisible by three, the half-shifted coordinates cannot equal 1/3 or 2/3, so no grid point samples the Dirac points. The cusp remains and controls the nonanalytic convergence behavior. The implementation evaluates rows in chunks, bounding temporary memory by \(O(N n_{\rm chunk})\).

## 2. Convergence and empirical error

The deterministic midpoint implementation gives:

| N | \(\bar U_N\) |
|---:|---:|
| 96 | 0.524865911187376 |
| 192 | 0.524865766524707 |
| 384 | 0.524865748435109 |
| 768 | 0.524865746173698 |
| 1536 | 0.524865745891015 |
| 3072 | 0.524865745855679 |
| 6144 | 0.524865745851262 |
| 12288 | 0.524865745850710 |

For three successive meshes, define \(\Delta_N=\bar U_{2N}-\bar U_N\) and infer
\[
p_N=\log_2\left|\frac{\Delta_{N/2}}{\Delta_N}\right|.
\]
If the asymptotic error is \(cN^{-p}\), the estimated error of the finer value \(\bar U_N\) is
\[
\epsilon_N\simeq\frac{|\Delta_{N/2}|}{2^{p_N}-1},
\qquad
\bar U_{\rm R}\simeq\bar U_N+
\frac{\Delta_{N/2}}{2^{p_N}-1}.
\]
Using the triple (3072, 6144, 12288) gives \(p\simeq2.99967\), estimated error \(6.31\times10^{-13}\) for the 6144 mesh, and Richardson value
\[
\bar U_{\rm R}\simeq0.524865745850631.
\]
The independent 24576 mesh gives 0.524865745850642, consistent with the extrapolation at the displayed precision.

This is numerical convergence evidence. The error estimate is conditional on an asymptotic power law and is not a rigorous quadrature bound. The stable reported value is approximately \(\bar U=0.5248657458506\); do not claim more certified digits.

## 3. Physics convention and use

The directed zero-flux gauge link remains \(U_b^{\gamma\gamma}=+1\), and the matter covariance is \(U_b^{00}=-\bar U\). The quadrature therefore updates the magnitude used in K-006 without changing its sign convention or the pure-Kitaev energy normalization established in Note 04c.

Run from the repository root with the project environment active:

    python scripts/k006_quadrature.py

The CLI accepts an explicit doubling sequence with --meshes and a memory-control parameter --row-chunk-size. Automated tests check direct small-mesh agreement, chunk invariance, invalid inputs, and observed convergence order.
