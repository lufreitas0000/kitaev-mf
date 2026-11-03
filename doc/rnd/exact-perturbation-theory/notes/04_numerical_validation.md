# 04 — Numerical validation plan

## Stage A: fixed-gauge Majorana solver

For an $L_1\times L_2$ periodic torus, there are $N=L_1L_2$ unit cells.
Construct the $N\times N$ matrix $T$ from the link variables and compute
its singular values. This gives all nonnegative matter-mode energies
without diagonalizing a $2N\times2N$ Hermitian matrix.

Analytical reference targets:

- Flux-free maximum positive matter energy: exactly $6K$ for periodic
  boundary conditions; the universal band bound is $6K$ for any twist.
- Flux-free thermodynamic Gaussian vacuum energy: approximately
  $-1.5746K$ per two-site unit cell, or $-0.7873K$ per spin, in the
  Pauli convention used here.
- Equivalent spin-operator convention:
  $-0.3936J$ per two-site unit cell for
  $H=-J\sum S_i^\alpha S_j^\alpha$ with $S^\alpha=\sigma^\alpha/2$.
- Unprojected matter Fock-space width: approximately $3.1492K$
  per unit cell in the thermodynamic limit.
- One flipped link changes exactly two adjacent plaquette fluxes.
- Pair energies approach the approximate targets in note 02 after
  parity and finite-size effects are handled.


- Unprojected finite-torus Fock width:
  $W_{\mathrm{Fock}}=\sum_n\epsilon_n$, exactly.
- Fixed-parity width for a gauge representative with required parity
  $p\in\{0,1\}$:
  $W_p=\sum_n\epsilon_n-
  [p+\mathbf 1_{N\bmod2\ne p}]\epsilon_{\min}$.
- Thermodynamic width density:
  $W_{\mathrm{Fock}}/N\to3.14919448K$ per two-site unit cell,
  equal to minus twice the vacuum energy density.
- Confirm numerically that the parity correction is controlled by
  $\epsilon_{\min}$ and vanishes per unit cell in the gapless limit.

The thermodynamic vacuum-energy quadrature to implement independently is
$$
-\frac{K}{4\pi^2}
\int_{-\pi}^{\pi}\!\!d\theta_1
\int_{-\pi}^{\pi}\!\!d\theta_2\,
\sqrt{3+2\cos\theta_1+2\cos\theta_2
+2\cos(\theta_1-\theta_2)}.
$$
The quadrature should be converged with respect to grid refinement or a
controlled integration tolerance; its output is compared with the
literature value, not used as its own reference.

The Dirac points cause pronounced finite-size dependence. Compare several
$L_1,L_2$, especially different values modulo three, and vary boundary
twists/global Wilson loops.

## Stage B: physical projection and parity

Derive the relation between matter-fermion parity and the gauge
configuration from the product of local constraints $D_i=+1$. The result
depends on site ordering and boundary conditions. Compute the Gaussian
vacuum parity independently, for example from the Pfaffian of the
covariance matrix, and combine it with the gauge-dependent constraint.
If the vacuum parity is incompatible, the lowest physical state requires
occupying a matter mode.

Do not infer physical-sector energies by comparing Gaussian-vacuum
energies alone.

## Stage C: independent spin-space exact diagonalization

For small clusters, construct the spin Hamiltonian directly and use sparse
exact diagonalization. Check the plaquette flux labels of eigenstates.
When degeneracies mix sectors, resolve the commuting flux operators within
the degenerate eigenspace. Compare sector-resolved energies against the
Majorana calculation.

## Output policy

Generated arrays, caches, logs, and plots belong under
build/rnd/exact-perturbation-theory/, not beside source notes. Keep only
small curated reference tables and scripts in version control. Do not
commit raw parameter sweeps or large binary outputs.
