# 04 — Numerical validation plan

## Stage A: fixed-gauge Majorana solver

For an \(L_1\times L_2\) periodic torus, there are \(N=L_1L_2\) unit cells.
Construct the \(N\times N\) matrix \(T\) from the link variables and compute
its singular values. This gives all positive matter energies without
diagonalizing a \(2N\times2N\) Hermitian matrix.

Initial validation targets:

- Flux-free maximum positive matter energy: \(6K\).
- Flux-free thermodynamic matter vacuum energy: approximately
  \(-1.5746K\) per unit cell.
- Unprojected matter Fock-space width: approximately \(3.1492K\)
  per unit cell.
- One flipped link changes exactly two adjacent plaquette fluxes.
- Pair energies approach the approximate targets in note 02 after
  parity and finite-size effects are handled.

The Dirac points cause pronounced finite-size dependence. Compare several
\(L_1,L_2\), especially different values modulo three, and vary boundary
twists/global Wilson loops.

## Stage B: physical projection and parity

Derive the relation between matter-fermion parity and the gauge
configuration from the product of local constraints \(D_i=+1\). The result
depends on site ordering and boundary conditions. Compute the Gaussian
vacuum parity independently, for example from the Pfaffian of the
covariance matrix, and combine it with the gauge-dependent constraint.
If the vacuum parity is incompatible, the lowest physical state requires
occupying a matter mode.

Do not infer physical-sector energies by comparing Gaussian vacuum
energies alone.

## Stage C: independent spin-space exact diagonalization

For small clusters, construct the spin Hamiltonian directly and use sparse
exact diagonalization. Check the plaquette flux labels of eigenstates.
When degeneracies mix sectors, resolve the commuting flux operators within
the degenerate eigenspace. Compare sector-resolved energies against the
Majorana calculation.

## Output policy

Generated arrays, caches, logs, and plots belong under
`build/rnd/exact-perturbation-theory/`, not beside source notes. Keep only
small curated reference tables and scripts in version control. Do not
commit raw parameter sweeps or large binary outputs.
