# 01 — Flux sectors and matter spectrum

## Quadratic matter problem

For a fixed link configuration \(u\), the Majorana Hamiltonian is quadratic.
With bipartite ordering \((A_1,\ldots,A_N,B_1,\ldots,B_N)\), define
\[
T_{\mathbf r,\mathbf r'}=2K\left[
u^z_{\mathbf r}\delta_{\mathbf r',\mathbf r}
+u^x_{\mathbf r}\delta_{\mathbf r',\mathbf r+\mathbf n_1}
+u^y_{\mathbf r}\delta_{\mathbf r',\mathbf r+\mathbf n_2}
\right],
\]
with periodic translations on a torus. The antisymmetric Majorana matrix is
\[
A=\begin{pmatrix}0&T\\-T^\mathsf T&0\end{pmatrix}.
\]
The \(N\) positive matter energies are the singular values \(s_n(T)\).
The unprojected Gaussian matter-vacuum energy is
\[
E_{\mathrm{vac}}(u)=-\frac12\sum_{n=1}^N s_n(T).
\]
This is not necessarily the lowest physical energy in that gauge sector:
the Gaussian vacuum must satisfy the physical matter-parity constraint.

## Flux-free dispersion

For the translationally invariant flux-free gauge,
\[
f(\mathbf q)=2K\left(1+e^{i\mathbf q\cdot\mathbf n_1}
+e^{i\mathbf q\cdot\mathbf n_2}\right),\qquad
\epsilon(\mathbf q)=|f(\mathbf q)|.
\]
Therefore
\[
0\leq\epsilon(\mathbf q)\leq6K.
\]
The upper edge \(6K\) is the maximum energy of one positive matter
quasiparticle, not the many-body bandwidth.

For \(N\) unit cells, the unprojected matter Fock spectrum is
\[
E(\{n_{\mathbf q}\})=
\sum_{\mathbf q}\epsilon(\mathbf q)
\left(n_{\mathbf q}-\frac12\right),
\qquad n_{\mathbf q}\in\{0,1\}.
\]
Its full spectral width is \(\sum_{\mathbf q}\epsilon(\mathbf q)\), which
is extensive. In the thermodynamic isotropic model, the vortex-free
ground-state energy is approximately \(-1.5746K\) per unit cell (two
spins per cell). Thus the unprojected Fock-space width is approximately
\(3.1492K\) per unit cell. Physical parity can change finite-size endpoints,
but not this thermodynamic width per cell.

## Flux sectors

A single link flip changes the fluxes of the two plaquettes sharing that
link, creating an adjacent vortex pair. A longer string of flipped links
separates the endpoints; the string itself is gauge dependent.

On a torus the total number of vortices is even. An isolated vortex
therefore requires an open boundary or an infinite-plane limit, with the
compensating flux at the boundary or infinity.

The matter Fock-space width and vortex-pair creation energy are different
quantities: the former grows linearly with system size, while the latter
approaches an \(O(K)\) defect energy.
