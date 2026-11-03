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
  $V=-\sum_i\sum_\alpha h_\alpha\sigma_i^\alpha$.
- Energies are reported in units of $K$.
- Each oriented $A\to B$ bond is counted once.
- Physical states satisfy $D_i=+1$ at every site.

## Milestones

| ID | Task | Status |
|---|---|---|
| M0.1 | Derive and verify the flux-free matter dispersion and one-particle upper edge | analytically complete; numerical implementation pending |
| M0.2 | Derive the thermodynamic matter-vacuum energy per unit cell | analytically complete; independent quadrature pending |
| M0.3 | Distinguish one-particle band range from extensive many-body width | analytically complete; finite-size numerical check pending |
| M0.4 | Compute adjacent-pair energy from one link flip | pending |
| M0.5 | Compute separated-pair energies and finite-size convergence | pending |
| M0.6 | Derive physical matter parity and check against spin-space ED | pending |
| M0.7 | Establish the Feshbach resolvent domain and perturbative conditions | pending |

## M0.1 result

In the vortex-free translational gauge, the Bloch matrix has eigenvalues
$\pm\epsilon(\mathbf q)$, with
$$
\epsilon(\mathbf q)=2K\left|
1+e^{i\theta_1}+e^{i\theta_2}
\right|,
\qquad \theta_a=\mathbf q\cdot\mathbf n_a.
$$
The triangle inequality proves $0\leq\epsilon(\mathbf q)\leq6K$,
with equality at $\theta_1=\theta_2=0$. The periodic finite-torus
maximum is exactly $6K$ because the zero-momentum point is always
allowed. This is the one-particle matter-band edge, not a many-body
bandwidth. The two Dirac points occur at
$(\theta_1,\theta_2)=(\pm2\pi/3,\mp2\pi/3)$ modulo $2\pi$.

## M0.2 result

The thermodynamic Gaussian vacuum energy per two-site unit cell is
$$
e_{\mathrm{vac}}
=-\frac{K}{4\pi^2}
\int_{-\pi}^{\pi}\!\!d\theta_1
\int_{-\pi}^{\pi}\!\!d\theta_2\,
\sqrt{3+2\cos\theta_1+2\cos\theta_2
+2\cos(\theta_1-\theta_2)}
\simeq -1.5746K.
$$
This is the value reported in Kitaev (2006), Appendix A. It is
$-0.7873K$ per spin in the present Pauli convention. For
$H=-J\sum S_i^\alpha S_j^\alpha$ with $S^\alpha=\sigma^\alpha/2$,
the equivalent value is $-0.3936J$ per two-site unit cell, since
$K=J/4$. The analytic integral and normalization are established;
an independent numerical quadrature remains part of the solver tests.

## M0.3 result

The unprojected matter Fock-space width is exactly
$W_{\mathrm{Fock}}(N)=\sum_n\epsilon_n$. Its thermodynamic density is
$-2e_{\mathrm{vac}}\simeq3.14919448K$ per two-site unit cell. In a fixed
gauge representative, physical projection retains one matter-parity
sector; the finite-size width differs from the unprojected width by at
most two times the smallest matter-mode energy. Since the isotropic model
is gapless and $\epsilon_{\min}\to0$ along regular thermodynamic
sequences, this does not change the width per unit cell.

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

No numerical results have been generated in this research track.
The values quoted above are analytical consequences or reference values
reported by Kitaev (2006), Appendix A; they are not outputs of the
project's numerical implementation.

## Open issues

- Derive the finite-torus physical matter-parity condition, including gauge
  configuration, site ordering, and boundary conditions.
- Compare global Wilson-loop sectors and boundary twists.
- Quantify finite-size effects near the Dirac points.
- Establish a useful perturbative bound for the field-dependent eliminated
  block. The norm of an extensive field perturbation is generally extensive,
  so a global operator-norm estimate may be too strong in the thermodynamic
  limit.
