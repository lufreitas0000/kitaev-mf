# 01 — Flux sectors and matter spectrum

## Quadratic matter problem

For a fixed link configuration $u$, the Majorana Hamiltonian is quadratic.
With bipartite ordering $(A_1,\ldots,A_N,B_1,\ldots,B_N)$, define
$$
T_{\mathbf r,\mathbf r'}=2K\left[
u^z_{\mathbf r}\delta_{\mathbf r',\mathbf r}
+u^x_{\mathbf r}\delta_{\mathbf r',\mathbf r+\mathbf n_1}
+u^y_{\mathbf r}\delta_{\mathbf r',\mathbf r+\mathbf n_2}
\right],
$$
with periodic translations on a torus. The antisymmetric Majorana matrix is
$$
A=\begin{pmatrix}0&T\\-T^\mathsf T&0\end{pmatrix}.
$$
The $N$ nonnegative matter-mode energies are the singular values $s_n(T)$.
The unprojected Gaussian matter-vacuum energy is
$$
E_{\mathrm{vac}}(u)=-\frac12\sum_{n=1}^N s_n(T).
$$
This is not necessarily the lowest physical energy in that gauge sector:
the Gaussian vacuum must satisfy the physical matter-parity constraint.

## M0.1 — Flux-free dispersion and its upper edge

In the translationally invariant flux-free gauge, $u^x=u^y=u^z=+1$.
The Fourier transform in the two-site unit cell gives
$$
i\widetilde A(\mathbf q)=
\begin{pmatrix}
0&if(\mathbf q)\\
-if(\mathbf q)^*&0
\end{pmatrix},
\qquad
f(\mathbf q)=2K\left(
1+e^{i\mathbf q\cdot\mathbf n_1}
+e^{i\mathbf q\cdot\mathbf n_2}
\right).
$$
Its eigenvalues are $\pm |f(\mathbf q)|$. Thus the positive matter-mode
energy is
$$
\epsilon(\mathbf q)=2K\left|
1+e^{i\theta_1}+e^{i\theta_2}
\right|,
\qquad \theta_a=\mathbf q\cdot\mathbf n_a.
$$
The triangle inequality gives
$$
0\leq\epsilon(\mathbf q)
\leq 2K(1+1+1)=6K.
$$
The upper bound is attained at $\theta_1=\theta_2=0$ modulo $2\pi$.
Consequently, $6K$ is the exact upper edge of the positive one-particle
matter band. On a periodic finite torus the $\mathbf q=0$ mode is allowed,
so the finite-size maximum is also exactly $6K$. With twisted boundary
conditions the allowed momenta are shifted; the universal bound remains
$6K$, but a finite-size sample need not attain it.

The band touches zero at the two inequivalent Dirac momenta
$(\theta_1,\theta_2)=(2\pi/3,-2\pi/3)$ and
$(-2\pi/3,2\pi/3)$ modulo $2\pi$. A finite torus samples these points
only for compatible sizes and twists. This accounts for the pronounced
modulo-three finite-size effects in the ground-state energy.

The number $6K$ is a one-quasiparticle excitation energy, not the
many-body bandwidth. For $N$ unit cells, the unprojected matter Fock
spectrum is
$$
E(\{n_{\mathbf q}\})=
\sum_{\mathbf q}\epsilon(\mathbf q)
\left(n_{\mathbf q}-\frac12\right),
\qquad n_{\mathbf q}\in\{0,1\},
$$
and its full width is $\sum_{\mathbf q}\epsilon(\mathbf q)$, which is
extensive. The physical parity constraint can remove one of the two
Fock-space endpoints at finite size; it does not change the
thermodynamic energy density.

## M0.2 — Thermodynamic vortex-free vacuum energy

For $N$ unit cells, the unprojected vacuum energy is
$$
E_{\mathrm{vac}}=-\frac12\sum_{\mathbf q}\epsilon(\mathbf q).
$$
In the thermodynamic limit, the momentum sum per unit cell becomes the
normalized Brillouin-zone integral. Using $\theta_a=\mathbf q\cdot
\mathbf n_a$, each angle ranges over a $2\pi$-periodic coordinate, and
$$
e_{\mathrm{vac}}
=\lim_{N\to\infty}\frac{E_{\mathrm{vac}}}{N}
=-\frac{K}{4\pi^2}
\int_{-\pi}^{\pi}\!\!d\theta_1
\int_{-\pi}^{\pi}\!\!d\theta_2\,
\sqrt{3+2\cos\theta_1+2\cos\theta_2
+2\cos(\theta_1-\theta_2)}.
$$
Numerical evaluation gives
$$
e_{\mathrm{vac}}\simeq -1.5746K
\quad\text{per two-site unit cell}.
$$
This value is reported by Kitaev (2006), Appendix A, for the
thermodynamic vortex-free phase. The integral above fixes the normalization
and provides the direct numerical target for an independent quadrature
implementation.

There are two sites per unit cell, so the same result is
$-0.7873K$ per spin. Be careful when comparing with papers using
$S_i^\alpha=\sigma_i^\alpha/2$: a Hamiltonian written as
$-J\sum_{\langle ij\rangle_\alpha}S_i^\alpha S_j^\alpha$ corresponds
to $K=J/4$ here and has energy $-0.3936J$ per two-site unit cell.
These numbers differ only by the spin-operator and coupling convention.

The integral computes the Gaussian matter-vacuum energy density. At finite
size, the physical parity constraint may require occupying a matter mode,
so the physical ground-state energy need not equal $E_{\mathrm{vac}}$.
That correction is subextensive and does not alter the thermodynamic
energy per unit cell. The vortex-free sector is the ground-state flux
sector; see Lieb's theorem as discussed in Kitaev (2006), Section 6.

## Flux sectors

A single link flip changes the fluxes of the two plaquettes sharing that
link, creating an adjacent vortex pair. A longer string of flipped links
separates the endpoints; the string itself is gauge dependent.

On a torus the total number of vortices is even. An isolated vortex
therefore requires an open boundary or an infinite-plane limit, with the
compensating flux at the boundary or infinity.

The matter Fock-space width and vortex-pair creation energy are different
quantities: the former grows linearly with system size, while the latter
approaches an $O(K)$ defect energy.

## Reference

A. Kitaev, “Anyons in an exactly solved model and beyond,” *Annals of
Physics* **321**, 2–111 (2006), Sec. 6, Eqs. (31)–(34), and Appendix A,
pp. 44–45, doi:10.1016/j.aop.2005.10.005.
