# Research status

## Current stage

M0 — establish the unperturbed spectrum, flux-sector thresholds, and
admissible low-energy resolvent domain. No perturbative coefficients have
been implemented.

## Starting conventions

- Spin Hamiltonian:
  $H_K=-K\sum_{\langle ij\rangle_\alpha}
  \sigma_i^\alpha\sigma_j^\alpha$, with $K>0$.
- Perturbation:
  $V=-\sum_{i,\alpha}h_\alpha\sigma_i^\alpha$.
- Energies are reported in units of $K$.
- Each oriented $A\to B$ bond is counted once.
- Physical states satisfy $D_i=+1$ at every site.

## Milestones

| ID | Task | Status |
|---|---|---|
| M0.1 | Derive and verify the flux-free matter dispersion and one-particle upper edge | pending |
| M0.2 | Verify the thermodynamic matter vacuum energy per unit cell | pending |
| M0.3 | Distinguish one-particle band range from extensive many-body width | pending |
| M0.4 | Compute adjacent-pair energy from one link flip | pending |
| M0.5 | Compute separated-pair energies and finite-size convergence | pending |
| M0.6 | Derive physical matter parity and check against spin-space ED | pending |
| M0.7 | Establish the Feshbach resolvent domain and perturbative conditions | pending |

## Important distinctions

1. The maximum positive matter quasiparticle energy is $6K$ in the
   isotropic convention. This is a one-particle energy scale.
2. The unprojected matter Fock-space spectral width is extensive. Its
   thermodynamic value is approximately $3.1492K$ per unit cell.
3. A vortex-pair creation energy is an $O(K)$ defect energy, independent
   of system size in the thermodynamic limit.
4. A matter-energy cutoff $P_\Lambda$ does not isolate a spectral block
   from its full orthogonal complement: that complement contains flux-free
   matter states immediately above the cutoff.
5. The flux threshold must be defined relative to the same physical
   Hilbert space and boundary conditions as the reference energy.
6. The Feshbach operator is defined only for energies in the resolvent set
   of the eliminated block. The chosen external-energy interval must lie
   strictly below that block's spectrum for the simple real-energy bound.

## Numerical status

No numerical results have been generated in this research track. Values
quoted in the notes are reference targets from Kitaev's Appendix A.

## Open issues

- Derive the finite-torus physical matter-parity condition, including gauge
  configuration, site ordering, and boundary conditions.
- Compare global Wilson-loop sectors and boundary twists.
- Quantify finite-size effects near the Dirac points.
- Establish a useful perturbative bound for the field-dependent eliminated
  block. The norm of an extensive field perturbation is generally extensive,
  so a global operator-norm estimate may be too strong in the thermodynamic
  limit.
