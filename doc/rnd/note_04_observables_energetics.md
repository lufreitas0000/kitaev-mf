# Note 04: Gaussian Observables and Variational Energetics

## 1. Scope and conventions

This note fixes the observable and energy conventions for the translationally invariant Majorana-Gaussian mean-field ansatz used by the finite-cell implementation. It distinguishes exact operator identities in the enlarged Majorana Fock space, Wick expectations in an unprojected Gaussian state, physical spin expectation values in the exactly projected spin Hilbert space, and finite Brillouin-zone sums versus thermodynamic-limit integrals.

These objects must not be identified without an additional argument. Satisfying local gauge constraints only on average does not make an unprojected Gaussian state an exactly projected physical state.

At each site,
\[
\sigma_i^a=i c_i^0c_i^a=\frac{i}{2}c_i^T N^a c_i,
\qquad
S_i^a=\frac12\sigma_i^a.
\]
The on-site and directed-bond covariance matrices are
\[
V_i^{\mu\nu}=\langle i c_i^\mu c_i^\nu\rangle-i\delta^{\mu\nu},
\qquad
U_b^{\mu\nu}=\langle i c_i^\mu c_j^\nu\rangle,
\]
where bond $b=(i,j)$ has the unit-cell orientation. Thus $V_i$ is real and skew-symmetric, while $U_b$ is real but need not be symmetric.

Define
\[
m_i^a:=\operatorname{tr}(V_i^T N^a),
\qquad
g_i^\gamma:=\operatorname{tr}(G^\gamma V_i).
\]
The first is the expectation of $i c_i^T N^a c_i=2\sigma_i^a$; the second is the local gauge residual in the convention used by the correlator service.

## 2. One-point spin observables

Because $N^a$ is real and skew-symmetric and has zero trace,
\[
\begin{aligned}
\langle \sigma_i^a\rangle_{\mathrm G}
&=\frac12\sum_{\mu,\nu}N^a_{\mu\nu}
  \langle i c_i^\mu c_i^\nu\rangle\\
&=\frac12\operatorname{tr}(V_i^T N^a)
=\frac12 m_i^a,
\\
\langle S_i^a\rangle_{\mathrm G}
&=\frac14 m_i^a.
\end{aligned}
\]
The subscript $\mathrm G$ denotes expectation in the unprojected Gaussian state. In a translation-invariant ansatz these moments are constant within each sublattice, but need not be equal between sublattices.

They are not automatically expectation values in the exactly projected state
\[
|\Psi_{\mathrm{phys}}\rangle
=\frac{P_G|\Psi_{\mathrm G}\rangle}
{\sqrt{\langle\Psi_{\mathrm G}|P_G|\Psi_{\mathrm G}\rangle}}.
\]
Replacing projected expectations by Gaussian contractions is an approximation unless a separate projection identity or controlled limit is established. The fact that local spin operators commute with the local gauge constraint does not by itself equate expectations in the projected and unprojected states.

## 3. Two-site spin correlators from Wick's theorem

For distinct sites $i,j$ connected by directed bond $b$, define
\[
C_b^{ab}:=
\operatorname{tr}\!\left[(N^a)^T U_b N^b U_b^T\right]
=\sum_{\mu,\nu,\rho,\sigma}
N^a_{\mu\nu}N^b_{\rho\sigma}
U_b^{\mu\rho}U_b^{\nu\sigma}.
\]
For a fermionic Gaussian state, Wick's theorem gives
\[
\boxed{
\langle \sigma_i^a\sigma_j^b\rangle_{\mathrm G}
=\frac14 m_i^a m_j^b-\frac12 C_b^{ab}.
}
\]
Consequently,
\[
\boxed{
\langle S_i^a S_j^b\rangle_{\mathrm G}
=\frac1{16}m_i^a m_j^b-\frac18 C_b^{ab}.
}
\]
The sign follows from the ordering
$c_i^\mu c_i^\nu c_j^\rho c_j^\sigma$ and the convention
$U_b^{\mu\nu}=\langle i c_i^\mu c_j^\nu\rangle$.
This formula assumes distinct sites; it must not be extended to an on-site product by simply setting $i=j$.

For exchange tensor $J_b^{ab}$, the Gaussian bond energy is
\[
E_{b,\mathrm{ex}}^{\mathrm G}
=\sum_{a,b'=1}^3 J_b^{ab'}
\left(\frac14 m_i^a m_j^{b'}-\frac12 C_b^{ab'}\right).
\]
It includes both the Hartree contribution and the intersite Wick contraction. Omitting either term changes the variational energy and its derivatives.

## 4. Gaussian variational energy

Use the spin Hamiltonian and directed-edge convention of Note 00d: each physical edge appears once in the directed bond set $\mathcal E$, equivalently as one half of the two-orientation ordered sum. For a Gaussian state with covariances $(U,V)$, Wick's theorem gives
\[
\boxed{
\begin{aligned}
E_{\mathrm{spin}}^{\mathrm G}[U,V]
={}&\sum_{b=(i,j)\in\mathcal E}
\sum_{a,b'=1}^3 J_b^{ab'}
\left[
\frac14 m_i^a m_j^{b'}
-\frac12
\operatorname{tr}\!\left((N^a)^T U_b N^{b'}U_b^T\right)
\right]\\
&-\frac12\sum_i\sum_{a=1}^3 h_i^a m_i^a.
\end{aligned}
}
\]
Here $b'$ is a spin-component index, not a bond label. The Zeeman term follows from $\langle\sigma_i^a\rangle_{\mathrm G}=m_i^a/2$.

If the quadratic constrained mean-field Hamiltonian includes the multiplier convention $H_\lambda^{\mathrm{constr}}=2\sum_{i,\gamma}\lambda_i^\gamma T_i^\gamma$, its Gaussian expectation is
\[
E_\lambda^{\mathrm G}
=\sum_{i,\gamma}2\lambda_i^\gamma\langle T_i^\gamma\rangle_{\mathrm G}
=-\frac14\sum_{i,\gamma}\lambda_i^\gamma g_i^\gamma.
\]
The sign follows from
$\langle i c_i^T G^\gamma c_i\rangle
=\operatorname{tr}((G^\gamma)^T V_i)
=-\operatorname{tr}(G^\gamma V_i)$.
Thus
\[
E_{\mathrm{MF}}^{\mathrm G}[U,V,\lambda]
=E_{\mathrm{spin}}^{\mathrm G}[U,V]
-\frac14\sum_{i,\gamma}\lambda_i^\gamma g_i^\gamma.
\]
The multiplier term vanishes when the constraints hold exactly. At a saddle point that satisfies them only approximately, its residual contribution must be reported rather than discarded.

This is the expectation of the original spin Hamiltonian in an unprojected Gaussian state, evaluated by Wick factorization. It is a variational mean-field functional, not the exact ground-state energy of the interacting spin Hamiltonian and not, by itself, the energy of the exactly projected state.

## 5. Quadratic spectrum and zero-point energy convention

For the real-space quadratic Majorana operator
\[
H_{\mathrm{quad}}=\frac{i}{4}\mathbf c^T M_{\mathrm{real}}\mathbf c,
\]
let $\lambda_m>0$ denote the positive eigenvalues of the Hermitian matrix $iM_{\mathrm{real}}$. Its ground-state contribution, with no additional scalar constant, is
\[
E_{\mathrm{quad},0}=-\frac12\sum_{\lambda_m>0}\lambda_m.
\]
In the Bloch convention of Note 02,
$\mathcal H(\mathbf k)=\frac{i}{4}\mathcal M(\mathbf k)$ and
$\varepsilon_n(\mathbf k)=4\mathcal E_n(\mathbf k)$.
For a full Brillouin-zone sum, count only positive eigenvalues of each
$\mathcal H(\mathbf k)$:
\[
E_{\mathrm{quad},0}
=-2\sum_{\mathbf k\in\mathrm{BZ}}
\sum_{n:\,\mathcal E_n(\mathbf k)>0}\mathcal E_n(\mathbf k).
\]
The sum is over all allowed momenta of the finite periodic system. There is no additional factor of one half for the $\mathbf k,-\mathbf k$ pair after restricting the inner sum to positive eigenvalues. This convention must be checked against the number of positive real-space Majorana frequencies, $2N_{\mathrm{sub}}N_{\mathrm{cell}}$.

If the many-body quadratic Hamiltonian is written with a scalar term as in Note 02,
\[
H_{\mathrm{MF}}=
\frac12\sum_{\mathbf k}\Psi_{\mathbf k}^\dagger
\mathcal H(\mathbf k)\Psi_{\mathbf k}-\mathrm C,
\]
its ground-state energy is the quadratic spectral contribution together with that scalar term, $E_{\mathrm{quad},0}-\mathrm C$.
The scalar $\mathrm C$ must be derived from the chosen decoupling and normal-ordering convention and compared with the Gaussian functional in Section 4. It must not be set to zero or inferred from band eigenvalues alone. Note 02 introduces $\mathrm C$ but does not give its explicit expression; deriving and testing it is an outstanding M1.4.3 task.

## 6. Brillouin-zone and finite-size normalization

For a periodic system with $N_{\mathrm{cell}}$ unit cells, the allowed fractional momenta form a discrete mesh in $[0,1)^2$ and
\[
\frac1{N_{\mathrm{cell}}}\sum_{\mathbf k}
\longrightarrow
\int_0^1 dk_1\int_0^1 dk_2
\]
in the thermodynamic limit. The fractional coordinates are dimensionless; reciprocal vectors obey
$\mathbf a_m\cdot\mathbf b_n=2\pi\delta_{mn}$.
The normalized BZ average is independent of the Cartesian area of the chosen unit cell.

The correlator implementation uses
\[
U_b=\frac1{N_k}\sum_{\mathbf k}(\cdots),
\qquad
V_s=\frac1{N_k}\sum_{\mathbf k}(\cdots),
\]
with $N_k=N_1N_2$ for its regular mesh. These are per-bond and per-sublattice correlators, not extensive quantities. Total energy is extensive; energy per unit cell is $E/N_{\mathrm{cell}}$, and energy per site is $E/(N_{\mathrm{sub}}N_{\mathrm{cell}})$.

For a supercell, the BZ volume is reduced by the index of the supercell and the number of sublattices increases by the same factor. Comparisons between primitive and supercell calculations must compare total energies per physical site (or consistently normalized energy densities), not raw per-cell energies.

A finite regular mesh approximates the normalized BZ integral by a Riemann sum. Mesh convergence is numerical evidence. At a gapless point, non-smooth integrands or a zero-mode prescription can affect convergence; a single mesh value is not a thermodynamic-limit error estimate.

## 7. Analytic checks and limiting cases

The following checks are independent targets for M1.4.2 and M1.4.3.

1. **Zero covariance:** $U_b=0$, $V_i=0$ gives $m_i^a=g_i^\gamma=C_b^{ab}=0$. The Gaussian spin and multiplier energies vanish, leaving only any separately specified scalar convention.
2. **Pure Kitaev reference state:** for $K<0$, zero field, zero other couplings, $V=0$, $U_b^{00}=u_{\mathrm{matter}}$, and $U_b^{\gamma\gamma}=+1$, the local moments vanish. The exchange functional reduces to the matter-bond contribution and must be consistent with K-006 and the pure-Kitaev energy per site.
3. **Unpolarized gauge-neutral state:** $V_i=0$ implies all one-point spin moments and gauge residuals vanish. It does not imply all two-site spin correlators vanish: the intersite Wick contraction can remain nonzero.
4. **Zero modes:** the half-weight prescription in the spectral projector gives a particle-hole-symmetric covariance. It is a specified Gaussian covariance for the zero-mode subspace; when exact zero modes are present, it need not represent a unique pure many-body ground state.
5. **Field-only response:** with exchange couplings zero, the Gaussian energy contains only Zeeman and constraint terms. A test should verify the sign and factor using a state with nonzero $V_i$ and independently constructed Majorana matrices.
6. **Cell normalization:** primitive-cell and supercell calculations describing the same state must agree in energy per physical site after folding and consistent covariance mapping.

The pure-Kitaev energy normalization should be derived explicitly in the follow-up energetics batch rather than inferred from the correlator magnitude alone. A comparison with a spin-model ground-state energy must state whether it is per bond, per site, or per unit cell and use the same Hamiltonian sign convention.

## 8. Scope and unresolved items

This note specifies Gaussian observables and the variational energy functional. It does not claim that unprojected Gaussian expectations equal exactly projected spin expectations. It does not derive the scalar double-counting term $\mathrm C$, a finite-temperature free energy, response functions, or a thermodynamic-limit numerical error bound for K-006.

The explicit $\mathrm C$ expression and a reproducible high-precision K-006 quadrature/convergence artifact remain open energetics/provenance tasks. They should be completed and tested before the corresponding numerical energy benchmark is treated as a strict reference.
