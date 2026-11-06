# 00 — Model and conventions

## Spin Hamiltonian

We use the ferromagnetic isotropic convention
$$
H_K=-K\sum_{\alpha=x,y,z}\sum_{\langle ij\rangle_\alpha}
\sigma_i^\alpha\sigma_j^\alpha,\qquad K>0,
$$
with every bond counted once. The field perturbation is
$$
V=-\sum_i\sum_{\alpha=x,y,z}h_\alpha\sigma_i^\alpha.
$$
All energies are in units of $K$, unless stated otherwise.

This convention must be distinguished from conventions using
$S_i^\alpha=\sigma_i^\alpha/2$, or absorbing factors into the exchange
couplings.

## Majorana representation

At each site introduce four Majoranas $b_i^x,b_i^y,b_i^z,c_i$, with
$$
\sigma_i^\alpha=i b_i^\alpha c_i,\qquad
D_i=b_i^x b_i^y b_i^z c_i=+1
$$
on physical states. For an oriented $A\to B$ bond of type $\alpha$,
define
$$
u_{ij}=i b_i^\alpha b_j^\alpha=\pm1.
$$
Then
$$
H_K[u]=iK\sum_{\langle ij\rangle_\alpha}u_{ij}c_i c_j.
$$
Equivalently, in $H=(i/4)\sum_{jk}A_{jk}c_jc_k$, the oriented matrix
elements are $A_{ij}=2K u_{ij}$ and $A_{ji}=-2K u_{ij}$.

The link variables are gauge dependent. Plaquette fluxes and
noncontractible Wilson loops are gauge invariant.

## Honeycomb unit cell

Use one $A$ and one $B$ site per unit cell $\mathbf r$, with primitive
translations $\mathbf n_1,\mathbf n_2$. Choose
$$
A_{\mathbf r}\leftrightarrow B_{\mathbf r}\;(z),\quad
A_{\mathbf r}\leftrightarrow B_{\mathbf r+\mathbf n_1}\;(x),\quad
A_{\mathbf r}\leftrightarrow B_{\mathbf r+\mathbf n_2}\;(y).
$$
In the vortex-free translational gauge, $u^x=u^y=u^z=+1$.
Different orientation conventions can change Bloch phases but not singular
values.

## Reference

Kitaev (2006), especially the model definition, Majorana representation,
matter dispersion, and Appendix A. Equation numbering and conventions
must be checked against the selected PDF before final publication.
