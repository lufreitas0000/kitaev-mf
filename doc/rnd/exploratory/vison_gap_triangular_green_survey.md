---
id: RND-EXP-VISON-GAP
title: "Analytical Calculation of the Kitaev Vison Gap via Triangular Green Functions"
phase: 1
milestone: M1.3
status: active
questions:
  - "How can the two-vison excitation energy in the Kitaev honeycomb model be evaluated analytically via triangular lattice Green functions?"
references:
  - "Panigrahi, Coleman, Tsvelik, PRB 108, 045151 (2023)"
---

# Analytical rederivation of the Kitaev vison gap
## Research programme, mathematical framework, and development roadmap

**Primary reference:** A. Panigrahi, P. Coleman, and A. Tsvelik, *Analytic calculation of the vison gap in the Kitaev spin liquid*, Phys. Rev. B **108**, 045151 (2023), [arXiv:2303.00093](https://arxiv.org/abs/2303.00093), [published version](https://doi.org/10.1103/PhysRevB.108.045151).

**Objective:** Completely rederive the adjacent-vison gap calculation using lattice Green functions, their random-walk/moment representations, and the exact theory of triangular-lattice Green functions. Establish the mathematical results needed for the derivation, independently verify the published result, and develop a subsequent extension to separated visons and long gauge strings.

---

# 1. Research objectives and scope

The main objective is to obtain an analytical and independently verifiable calculation of the energy cost of creating two adjacent visons in the isotropic Kitaev honeycomb model.

The reference paper formulates the problem as scattering of free Majorana fermions from a local bond reversal. Its central result is

\[
\Delta_v=0.1311656(3)K=0.263313(6)J,
\qquad K=2J.
\]

The calculation uses a local Green function, a finite-dimensional scattering determinant, and a frequency integral over the associated scattering phase shift. The local Green-function integrals are reduced to one-dimensional expressions and evaluated numerically.

Our programme is to preserve the physical scattering formulation while replacing the intermediate Green-function treatment with a more systematic mathematical framework.

The programme has four objectives:

1. Reconstruct the complete calculation in the reference paper, including conventions, normalizations, analytic continuation, and the zero-temperature limit.
2. Derive exact representations for the local and real-space Green functions from triangular-lattice Green-function theory.
3. Establish the mathematical statements underlying the calculation, including resolvent identities, lattice Green-function recurrences, determinant identities, and spectral-shift formulae.
4. Extend the formalism to two separated visons, first at finite separation and eventually in the asymptotic regime.

The adjacent-vison calculation is the first required deliverable. The separated-vison problem is a subsequent research direction, not an assumption required for the first result.

# 2. Reference literature

## 2.1 Primary physical reference

**Panigrahi, Coleman, and Tsvelik (2023)**

A. Panigrahi, P. Coleman, and A. Tsvelik, “Analytic calculation of the vison gap in the Kitaev spin liquid,” *Physical Review B* **108**, 045151 (2023).

- [arXiv:2303.00093](https://arxiv.org/abs/2303.00093)
- [HTML version, v2](https://arxiv.org/html/2303.00093v2)
- [DOI: 10.1103/PhysRevB.108.045151](https://doi.org/10.1103/PhysRevB.108.045151)

The relevant structure is:

- Majorana representation and static \(\mathbb Z_2\) gauge fields.
- Flux-free honeycomb Hamiltonian.
- A local bond reversal producing two adjacent visons.
- A finite-rank scattering potential.
- Local Green functions \(g_0(z)\) and \(g_2(z)\).
- A scattering determinant and phase shift.
- A zero-temperature frequency integral for the vison gap.

The paper's isotropic convention is

\[
K=2J,\qquad V=4J=2K,
\]

and

\[
\gamma_{\mathbf k}=iK(1+e^{ik_1}+e^{ik_2}).
\]

The paper obtains

\[
\Delta_v
=
K\int_0^\infty\frac{dx}{2\pi}
\operatorname{Im}\log
\left[
(1+2g_2(x-i0^+))^2
-
(2g_0(x-i0^+))^2
\right].
\]

This expression, including its branch prescription and normalization, will be independently rederived rather than adopted without verification.

## 2.2 Triangular- and honeycomb-lattice Green functions

**Horiguchi (1972)**

T. Horiguchi, “Lattice Green's Functions for the Triangular and Honeycomb Lattices,” *Journal of Mathematical Physics* **13**, 1411–1419 (1972).

[DOI: 10.1063/1.1666155](https://doi.org/10.1063/1.1666155)

This is the principal mathematical reference for the project. It establishes that triangular-lattice Green functions at arbitrary lattice sites can be expressed using complete elliptic integrals of the first and second kinds, and relates honeycomb-lattice Green functions to triangular-lattice Green functions.

Its central relevance is the reduction of our two-sublattice honeycomb resolvent to a scalar triangular-lattice resolvent.

**Morita and Horiguchi (1972)**

T. Morita and T. Horiguchi, “Analytic properties of the lattice Green function,” *Journal of Physics A* **5**, 67 (1972).

[DOI: 10.1088/0305-4470/5/1/011](https://doi.org/10.1088/0305-4470/5/1/011)

Relevant for analytic continuation, singularities associated with critical points of the dispersion, and the local behaviour of lattice Green functions near their spectral singularities.

**Horiguchi (1991)**

T. Horiguchi, “Lattice Green's function for anisotropic triangular lattice,” *Physica A* **178**, 351–363 (1991).

[DOI: 10.1016/0378-4371(91)90024-7](https://doi.org/10.1016/0378-4371(91)90024-7)

Relevant for extensions beyond the isotropic model. The arbitrary-site Green function is represented using a finite set of fundamental lattice Green functions and their parameter derivatives, expressed through complete elliptic integrals of the first, second, and third kinds.

## 2.3 Further mathematical references

**K. Uchiyama (1998)**

K. Uchiyama, “Green's Functions for Random Walks on \(\mathbb Z^N\),” *Proceedings of the London Mathematical Society* **77**, 215–240 (1998).

Relevant background for connecting lattice Green functions, random walks, and potential kernels.

**D. Kapanadze (2021)**

D. Kapanadze, “The far-field behaviour of Green's function for a triangular lattice and radiation conditions,” *Mathematical Methods in the Applied Sciences* **44**, 12746–12759 (2021).

[DOI: 10.1002/mma.7575](https://doi.org/10.1002/mma.7575)

Potentially useful for the large-distance asymptotics of real-space triangular-lattice Green functions, especially after retarded boundary conditions are imposed.

These references serve different purposes. Horiguchi supplies the exact lattice Green-function machinery; Morita and Horiguchi supply analytic-structure results; random-walk literature provides the combinatorial interpretation of resolvent moments; far-field lattice Green-function work is relevant to the separated-vison extension.

# 3. Analytical findings established so far

## 3.1 Free Majorana Hamiltonian

In the flux-free sector, the isotropic honeycomb Hamiltonian is represented by

\[
H_0(\mathbf k)=\vec\gamma_{\mathbf k}\cdot\vec\tau.
\]

Using the reference-paper convention,

\[
\gamma_{\mathbf k}=iKf(\mathbf k),
\qquad
f(\mathbf k)=1+e^{ik_1}+e^{ik_2}.
\]

The spectrum is

\[
\epsilon_{\mathbf k}=\pm K|f(\mathbf k)|.
\]

The free resolvent is

\[
G_0(z,\mathbf k)
=
[z-H_0(\mathbf k)]^{-1}
=
\frac{z+H_0(\mathbf k)}
{z^2-K^2|f(\mathbf k)|^2}.
\]

Here \(z\) is a complex spectral parameter outside the spectrum. Retarded and advanced Green functions are obtained as boundary values from the upper and lower half-planes, respectively.

The spectral support is contained in

\[
-3K\leq \epsilon\leq 3K.
\]

The actual density of states has additional structure, including the Dirac points and van Hove singularities. The resolvent is analytic off the spectrum.

## 3.2 The triangular-lattice reduction

The central algebraic identity is

\[
|f(\mathbf k)|^2
=
3+2\cos k_1+2\cos k_2
+2\cos(k_1-k_2).
\]

This is the dispersion polynomial of a triangular Bravais lattice, up to an affine rescaling.

Define the dimensionless triangular-lattice Green function

\[
\mathcal T_{\mathbf R}(s)
=
\int_{[0,2\pi)^2}
\frac{dk_1\,dk_2}{(2\pi)^2}
\frac{e^{i\mathbf k\cdot\mathbf R}}
{s-|f(\mathbf k)|^2}.
\]

The diagonal honeycomb Green function becomes

\[
G_{AA}(\mathbf R,z)
=
\frac{z}{K^2}\mathcal T_{\mathbf R}
\left(\frac{z^2}{K^2}\right).
\]

In particular,

\[
g_0(z)
=
\frac{z}{K^2}\mathcal T_{\mathbf 0}
\left(\frac{z^2}{K^2}\right).
\]

Equivalently, defining the dimensionless frequency \(\zeta=z/K\),

\[
g_0(z)=\frac{1}{K}\zeta\mathcal T_{\mathbf 0}(\zeta^2).
\]

This identifies the local honeycomb Green function with a triangular-lattice Green function at the origin.

The off-diagonal Green function also reduces to a finite combination of triangular-lattice Green functions. With the Fourier convention \(e^{i\mathbf k\cdot\mathbf R}\) and the corresponding convention for \(H_{AB}(\mathbf k)=Kf(\mathbf k)\),

\[
G_{AB}(\mathbf R,z)
=
\frac{1}{K}
\left[
\mathcal T_{\mathbf R}(\zeta^2)
+
\mathcal T_{\mathbf R+\mathbf a_1}(\zeta^2)
+
\mathcal T_{\mathbf R+\mathbf a_2}(\zeta^2)
\right].
\]

The shifts must be checked against the precise Fourier and sublattice conventions used in each derivation. Reversing the Fourier convention reverses the shifts.

Thus, all components of the honeycomb resolvent can be constructed from scalar triangular-lattice Green functions.

## 3.3 Local Green functions in the reference paper

The paper defines

\[
g_0(z)
=
\int_{\mathrm{BZ}}\frac{d^2k}{(2\pi)^2}
\frac{z}{z^2-K^2|f(\mathbf k)|^2}
\]

and

\[
g_2(z)
=
\int_{\mathrm{BZ}}\frac{d^2k}{(2\pi)^2}
\frac{K(1+\cos k_1+\cos k_2)}
{z^2-K^2|f(\mathbf k)|^2}.
\]

The second expression is not an independent scalar lattice Green function. In the isotropic case, the symmetry among the three cosine terms in \(|f|^2\) gives an exact relation.

Writing \(\zeta=z/K\), one obtains

\[
\boxed{
g_2(z)=\frac{\zeta g_0(z)-1}{3}.
}
\]

This relation follows from the triangular-lattice symmetry and the identity

\[
|f|^2=3+2q,\qquad
q=\cos k_1+\cos k_2+\cos(k_1-k_2).
\]

Indeed, symmetry implies that the Brillouin-zone averages of the three cosine terms divided by any function of \(|f|^2\) are equal. Therefore,

\[
\left\langle
\frac{1+\cos k_1+\cos k_2}
{\zeta^2-|f|^2}
\right\rangle
=
\frac{\zeta^2\mathcal T_{\mathbf 0}(\zeta^2)-1}{3}.
\]

Since \(g_0= \zeta\mathcal T_{\mathbf 0}(\zeta^2)/K\), the stated identity follows.

This is an important simplification: in the isotropic problem, the determinant can ultimately be expressed entirely in terms of the local triangular-lattice Green function.

The relation must be rederived with the paper's exact normalization and Fourier convention in the formal derivation.

## 3.4 One-dimensional representation

Integrating one momentum variable gives

\[
\mathcal T_{\mathbf 0}(s)
=
\frac{1}{\pi}
\int_0^\pi
\frac{d\theta}
{
\sqrt{
(s-1-4\cos^2\theta)^2
-
16\cos^2\theta
}
}.
\]

This expression is valid with a branch prescription defined by analytic continuation from the resolvent domain. The square-root branch cannot be chosen independently at each real frequency.

The integral is algebraic and is expected to reduce to complete elliptic integrals. The exact reduction, including its branch structure, is a primary analytical task.

## 3.5 Random-walk interpretation

The expansion

\[
\mathcal T_{\mathbf 0}(s)
=
\frac{1}{s}
\sum_{n=0}^{\infty}
\frac{\langle |f(\mathbf k)|^{2n}\rangle}{s^n}
\]

is valid for sufficiently large \(|s|\), specifically outside the convergence region determined by the spectrum.

The coefficients

\[
\langle |f(\mathbf k)|^{2n}\rangle
\]

are weighted closed-walk counts on the associated triangular lattice. This follows because Brillouin-zone integration annihilates every nonzero lattice translation.

More generally, for a hopping function

\[
h(\mathbf k)=\sum_{j=1}^{q}t_j e^{i\mathbf k\cdot\boldsymbol\delta_j},
\]

the moments of \(|h(\mathbf k)|^2\) count weighted closed walks generated by differences of the hopping vectors.

The random-walk interpretation is therefore general for bipartite hopping problems of this form. It is not restricted to the Kitaev model.

A separate quantity considered earlier was

\[
\frac13
\left\langle
|1+e^{i\theta_1}+e^{i\theta_2}|
\right\rangle
\approx 0.5248657458506.
\]

This is a three-step planar random-walk endpoint-distance average. It is mathematically related to the same form factor, but it is not itself the local Green function or the vison gap. The moment expansion of the resolvent is the more directly useful random-walk connection for the present project.

# 4. Mathematical framework to be developed

The project should state general definitions, lemmas, and theorems before specializing to the honeycomb model. The following is the proposed mathematical structure.

## 4.1 Lattice Hamiltonians and resolvents

**Definition 1 — Periodic bipartite hopping Hamiltonian**

Let \(\Lambda\) be a Bravais lattice with two sublattices. Consider

\[
H(\mathbf k)=
\begin{pmatrix}
0 & h(\mathbf k)\\
h(\mathbf k)^* & 0
\end{pmatrix},
\qquad
h(\mathbf k)=\sum_j t_j e^{i\mathbf k\cdot\boldsymbol\delta_j}.
\]

The spectrum is

\[
E_\pm(\mathbf k)=\pm |h(\mathbf k)|.
\]

**Definition 2 — Lattice resolvent**

For \(z\) outside the spectrum, define

\[
G_0(z)=(z-H_0)^{-1}.
\]

For a periodic system, its real-space matrix elements are Fourier coefficients of \(G_0(z,\mathbf k)\).

**Lemma 1 — Bipartite resolvent identity**

For the above two-sublattice Hamiltonian,

\[
G_0(z,\mathbf k)
=
\frac{1}{z^2-|h(\mathbf k)|^2}
\begin{pmatrix}
z & h(\mathbf k)\\
h(\mathbf k)^* & z
\end{pmatrix}.
\]

This identity is the algebraic basis for the honeycomb-to-triangular reduction.

## 4.2 Moment expansion and closed walks

**Definition 3 — Spectral moments**

For a lattice Hamiltonian \(H\) and a site state \(|\mathbf R,\alpha\rangle\), define

\[
\mu_n(\mathbf R,\alpha)
=
\langle\mathbf R,\alpha|H^n|\mathbf R,\alpha\rangle.
\]

**Lemma 2 — Walk interpretation**

For finite-range hopping Hamiltonians, \(\mu_n\) is a sum over all length-\(n\) closed paths beginning and ending at the specified site, with each path weighted by the product of its hopping amplitudes.

**Theorem 1 — Resolvent moment expansion**

For \(|z|>\|H\|\),

\[
G_{\mathbf R\mathbf R}(z)
=
\sum_{n=0}^{\infty}
\frac{\mu_n(\mathbf R)}{z^{n+1}}.
\]

The expansion is exact in its domain of convergence. Its analytic continuation is the resolvent itself, not necessarily the termwise continuation of the series.

This theorem gives a systematic way to derive high-frequency expansions and validate the elliptic Green-function expressions.

## 4.3 Triangular-lattice Green functions

**Definition 4 — Triangular-lattice Green function**

For

\[
\varepsilon_\triangle(\mathbf k)
=
3+2\cos k_1+2\cos k_2+2\cos(k_1-k_2),
\]

define

\[
\mathcal T_{\mathbf R}(s)
=
\int_{\mathrm{BZ}}\frac{d^2k}{(2\pi)^2}
\frac{e^{i\mathbf k\cdot\mathbf R}}
{s-\varepsilon_\triangle(\mathbf k)}.
\]

**Lemma 3 — Lattice difference equation**

The Fourier definition implies a discrete difference equation relating the origin Green function to its nearest-neighbour values and, more generally, relating values at neighbouring lattice sites.

The exact coefficients and normalization must be derived explicitly for the chosen dispersion convention.

**Theorem 2 — Finite-generation representation**

For the isotropic triangular lattice, arbitrary-site Green functions can be reduced to a finite set of fundamental Green functions and their parameter derivatives, using the lattice difference equations and symmetry.

Horiguchi's 1972 result provides the principal reference and target for reproducing this construction.

The project must distinguish between:
- an exact representation in terms of elliptic integrals;
- a recurrence that generates arbitrary-site values;
- a numerical evaluation of those expressions.

These are related but not interchangeable claims.

## 4.4 Analytic continuation and spectral singularities

**Definition 5 — Resolvent boundary values**

For real \(E\), define

\[
G_0^{R/A}(E)
=
\lim_{\eta\downarrow0}G_0(E\pm i\eta).
\]

**Lemma 4 — Analyticity of the resolvent**

For a bounded self-adjoint lattice Hamiltonian, the resolvent is analytic outside the real spectrum and satisfies

\[
G_0(z)^\dagger=G_0(z^*).
\]

**Theorem 3 — Spectral representation**

For a local state,

\[
G_{\mathbf R\mathbf R}(z)
=
\int_{\mathbb R}
\frac{\rho_{\mathbf R}(E)}{z-E}\,dE,
\]

where \(\rho_{\mathbf R}\) is the corresponding spectral measure.

The boundary values determine the density of states through the discontinuity across the spectrum.

The project must identify the singular points of the triangular dispersion and derive the local behaviour of \(\mathcal T_{\mathbf 0}(s)\) around them. In particular, the Dirac-point and van Hove contributions should be treated separately.

## 4.5 Finite-rank perturbations and scattering determinants

**Definition 6 — Finite-rank perturbation**

Let

\[
H=H_0+V,
\]

where \(V\) has finite rank and can be written

\[
V=U\Lambda U^\dagger.
\]

**Theorem 4 — Finite-rank determinant identity**

For \(z\) outside the spectra of \(H_0\) and \(H\),

\[
\frac{\det(z-H)}{\det(z-H_0)}
=
\det\left[
I-\Lambda U^\dagger G_0(z)U
\right].
\]

The determinant on the right acts only in the finite-dimensional perturbation space.

This is the central mathematical reason that a local bond reversal can be treated using a small matrix rather than the full infinite lattice Hamiltonian.

**Lemma 5 — Resolvent identity**

The perturbed resolvent is

\[
G(z)
=
G_0(z)+
G_0(z)U
\left[
\Lambda^{-1}-U^\dagger G_0(z)U
\right]^{-1}
U^\dagger G_0(z),
\]

when the indicated inverse exists.

The equivalent \(T\)-matrix representation will also be derived.

## 4.6 Spectral shift and the energy difference

**Definition 7 — Scattering determinant**

Define

\[
D(z)=
\det\left[I-\Lambda U^\dagger G_0(z)U\right].
\]

The phase of \(D(E+i0^+)\) determines the spectral shift, subject to the convention used for the logarithm and the sign of the perturbation.

**Theorem 5 — Spectral-shift relation**

For a suitable trace-class perturbation, the spectral shift function \(\xi(E)\) is related to the boundary value of the perturbation determinant through

\[
\xi(E)
=
\frac{1}{\pi}
\operatorname{Im}\log D(E+i0^+),
\]

up to the sign convention adopted for the perturbation determinant.

The precise sign and normalization must be fixed by deriving the relation from the chosen definition of \(D\), not by copying a convention from another source.

**Theorem 6 — Ground-state energy shift**

For the Majorana quadratic Hamiltonian, the zero-temperature energy shift can be expressed as a frequency integral over the scattering phase. The factor of \(1/2\) associated with Majorana degrees of freedom must be retained throughout the derivation.

The project will derive the formula starting from the Matsubara determinant, perform the contour deformation explicitly, and only then take the zero-temperature limit.

This is essential because the Majorana spectrum is particle-hole symmetric and a careless use of a complex-fermion formula can introduce a factor-of-two error.

# 5. Detailed development roadmap

The work is organized into phases, milestones, and independently verifiable tasks. Each milestone should produce a self-contained mathematical note and a reproducible validation record.

## Phase I — Reconstruction of the reference calculation

### M1.1 — Conventions and physical setup

**Objective:** Establish an unambiguous definition of the unperturbed Hamiltonian, the gauge configuration, and the vison excitation.

Tasks:

- T1.1.1: Reproduce the Majorana representation and the static \(\mathbb Z_2\) gauge-field formulation.
- T1.1.2: Fix the normalization of Majorana operators and the quadratic Hamiltonian.
- T1.1.3: Derive the momentum-space Hamiltonian and identify the sublattice Pauli-matrix convention.
- T1.1.4: Reproduce the form factor \(\gamma_{\mathbf k}\) and the dispersion.
- T1.1.5: Specify the gauge transformation that moves a string without changing its endpoints.
- T1.1.6: Derive the bond-reversal perturbation and verify \(V=4J=2K\).
- T1.1.7: Distinguish the adjacent-vison pair energy from a single-vison energy.

**Deliverable:** A convention and model-definition note.

**Acceptance criteria:** Every factor of two, sign, sublattice convention, and Fourier normalization is independently checked.

### M1.2 — Scattering formulation

**Objective:** Reconstruct the finite-rank scattering problem.

Tasks:

- T1.2.1: Express the local bond perturbation as a finite-rank operator.
- T1.2.2: Derive the perturbed resolvent using the resolvent identity.
- T1.2.3: Derive the finite-dimensional determinant using the matrix determinant lemma.
- T1.2.4: Derive the Matsubara free-energy difference.
- T1.2.5: Reproduce the determinant in terms of the two local Green-function components.
- T1.2.6: Verify the sign of the determinant and its dependence on \(V\).

**Deliverable:** A scattering-theory note with all determinant manipulations shown.

**Acceptance criteria:** The resulting determinant agrees with the reference paper after translating conventions.

### M1.3 — Spectral-shift and zero-temperature formula

**Objective:** Derive the vison gap as a phase-shift integral.

Tasks:

- T1.3.1: Derive the Matsubara-frequency expression from the quadratic Majorana partition function.
- T1.3.2: Perform the contour deformation to the real-frequency axis.
- T1.3.3: Establish the correct retarded or advanced boundary prescription.
- T1.3.4: Fix the logarithm branch by continuity from large complex frequency.
- T1.3.5: Derive the zero-temperature occupation factor.
- T1.3.6: Verify the Majorana factor of \(1/2\).
- T1.3.7: Recover the paper's dimensionless integral.

**Deliverable:** A complete derivation of the phase-shift formula.

**Acceptance criteria:** The derivation fixes the sign, integration domain, normalization, and branch of the logarithm without relying on numerical fitting.

## Phase II — Exact triangular-lattice Green functions

### M2.1 — Honeycomb-to-triangular mapping

**Objective:** Establish the exact reduction of the honeycomb resolvent to triangular-lattice Green functions.

Tasks:

- T2.1.1: Derive the identity for \(|f(\mathbf k)|^2\).
- T2.1.2: Define the triangular-lattice Green function with explicit normalization.
- T2.1.3: Derive the diagonal honeycomb Green-function relation.
- T2.1.4: Derive the off-diagonal relation and all required lattice shifts.
- T2.1.5: Verify the relations using large-\(|z|\) expansions.
- T2.1.6: Check the symmetry properties under sublattice exchange and complex conjugation.

**Deliverable:** A mathematical note proving the honeycomb-to-triangular reduction.

**Acceptance criteria:** All matrix elements have consistent units and agree with the direct Fourier representation.

### M2.2 — Local triangular-lattice Green function

**Objective:** Obtain an exact expression for \(\mathcal T_{\mathbf 0}(s)\).

Tasks:

- T2.2.1: Derive the two-dimensional Brillouin-zone integral.
- T2.2.2: Perform the one-dimensional reduction independently.
- T2.2.3: Derive the polynomial under the square root and identify its roots.
- T2.2.4: Map the integral to a standard complete elliptic-integral form.
- T2.2.5: Derive the analytic continuation for complex \(s\).
- T2.2.6: Identify the branch cuts and the correct boundary values.
- T2.2.7: Derive asymptotic expansions for large \(|s|\).
- T2.2.8: Derive local expansions near the relevant spectral singularities.

**Deliverable:** An exact expression for \(\mathcal T_{\mathbf 0}(s)\), together with a branch and analytic-continuation specification.

**Acceptance criteria:** The expression reproduces the first several moments, the large-\(|s|\) asymptotics, and the direct two-dimensional integral.

### M2.3 — Exact local functions \(g_0\) and \(g_2\)

**Objective:** Obtain the local Green functions required by the scattering determinant.

Tasks:

- T2.3.1: Express \(g_0(z)\) in terms of \(\mathcal T_{\mathbf 0}\).
- T2.3.2: Prove the isotropic relation \(g_2=(\zeta g_0-1)/3\).
- T2.3.3: Substitute the exact elliptic expression into both functions.
- T2.3.4: Derive their real and imaginary parts on the physical frequency axis.
- T2.3.5: Analyze the low-frequency and band-edge behaviour.
- T2.3.6: Compare with the one-dimensional integral representation in the paper.

**Deliverable:** Exact local Green functions and a separate note proving their algebraic relation.

**Acceptance criteria:** Both functions reproduce the paper's integral representations and are numerically consistent with direct Brillouin-zone integration.

## Phase III — Exact adjacent-vison gap

### M3.1 — Scattering determinant in elliptic form

**Objective:** Express the complete scattering determinant using the exact local Green function.

Tasks:

- T3.1.1: Insert the exact \(g_0\) and \(g_2\) into the finite-rank determinant.
- T3.1.2: Simplify the determinant using the isotropic identity.
- T3.1.3: Identify zeros and poles in the complex frequency plane.
- T3.1.4: Determine the correct continuous phase of the determinant.
- T3.1.5: Identify possible resonance contributions and their relation to the local density of states.

**Deliverable:** An exact analytic representation of the scattering determinant.

### M3.2 — Phase shift and vison gap

**Objective:** Evaluate the phase-shift integral using the exact Green functions.

Tasks:

- T3.2.1: Derive the phase shift from the determinant.
- T3.2.2: Establish the phase normalization at large frequency.
- T3.2.3: Analyze the phase near zero frequency.
- T3.2.4: Analyze its behaviour near van Hove and band-edge singularities.
- T3.2.5: Evaluate the zero-temperature integral with controlled analytic continuation.
- T3.2.6: Establish numerical error estimates independently of the reference paper's interpolation procedure.
- T3.2.7: Recover \(\Delta_v/K\) and compare with \(0.1311656(3)\).

**Deliverable:** A complete independent derivation of the adjacent-vison gap.

**Acceptance criteria:** The result agrees with the reference value within controlled numerical uncertainty. The calculation must not rely on fitting the published result.

### M3.3 — Independent consistency checks

Tasks:

- T3.3.1: Compare the elliptic expression with direct two-dimensional integration.
- T3.3.2: Compare the one-dimensional quadrature with the elliptic expression.
- T3.3.3: Check high-frequency moments.
- T3.3.4: Check the determinant in the weak-perturbation limit.
- T3.3.5: Check the spectral-shift normalization.
- T3.3.6: Compare the resulting energy with finite-size diagonalization, if available.

**Deliverable:** A validation report.

The numerical comparison is a verification tool, not a substitute for the analytical derivation.

## Phase IV — Real-space Green functions

### M4.1 — Triangular-lattice recurrence relations

**Objective:** Construct arbitrary-site Green functions from a finite set of fundamental functions.

Tasks:

- T4.1.1: Derive the triangular-lattice difference equation.
- T4.1.2: Identify the independent symmetry classes of lattice sites.
- T4.1.3: Derive recurrence relations between neighbouring Green functions.
- T4.1.4: Reconstruct the Horiguchi representation.
- T4.1.5: Determine the minimal set of fundamental elliptic integrals.
- T4.1.6: Derive the parameter derivatives needed for arbitrary sites.

**Deliverable:** A recurrence-based representation of \(\mathcal T_{\mathbf R}(s)\).

### M4.2 — Honeycomb real-space Green functions

Tasks:

- T4.2.1: Obtain \(G_{AA}(\mathbf R,z)\).
- T4.2.2: Obtain \(G_{BB}(\mathbf R,z)\).
- T4.2.3: Obtain \(G_{AB}(\mathbf R,z)\) and \(G_{BA}(\mathbf R,z)\).
- T4.2.4: Check Hermiticity and sublattice symmetries.
- T4.2.5: Verify selected near-neighbour matrix elements directly.
- T4.2.6: Derive large-distance asymptotics in the relevant spectral regimes.

**Deliverable:** A complete real-space honeycomb Green-function formalism.

## Phase V — Two separated visons

This phase generalizes the local bond-reversal problem. It should begin only after the single-bond scattering theory and real-space Green functions have been established.

### M5.1 — Finite-separation gauge configurations

**Objective:** Construct a gauge configuration containing two visons at arbitrary finite separation.

Tasks:

- T5.1.1: Specify the two plaquette fluxes and their positions.
- T5.1.2: Construct a finite gauge string joining the two fluxes.
- T5.1.3: Derive the associated changes in the Majorana hopping Hamiltonian.
- T5.1.4: Identify the finite set of modified bonds.
- T5.1.5: Express the perturbation as a finite-rank operator.
- T5.1.6: Determine the rank and structure of the finite-dimensional scattering matrix.
- T5.1.7: Establish gauge equivalence between different strings with identical endpoints.

**Deliverable:** A general finite-separation perturbation operator and scattering determinant.

**Acceptance criteria:** Gauge-related strings produce identical physical spectra and energy differences.

### M5.2 — Finite-separation determinant

**Objective:** Express the two-vison energy in terms of real-space Green functions.

Tasks:

- T5.2.1: Construct the finite-dimensional matrix \(U^\dagger G_0(z)U\).
- T5.2.2: Express its entries using the real-space honeycomb Green functions.
- T5.2.3: Derive the determinant for the chosen string geometry.
- T5.2.4: Identify the dependence on vison separation.
- T5.2.5: Derive the associated spectral shift and energy integral.
- T5.2.6: Verify the adjacent-vison limit against M3.

**Deliverable:** An exact finite-separation determinant formulation.

### M5.3 — Large-separation asymptotics

**Objective:** Determine the asymptotic interaction between separated visons.

Tasks:

- T5.3.1: Derive the large-\(|\mathbf R|\) asymptotics of the required real-space Green functions.
- T5.3.2: Identify contributions from the Dirac points.
- T5.3.3: Determine whether van Hove regions contribute to the relevant energy integral.
- T5.3.4: Expand the finite-separation determinant in the large-distance regime.
- T5.3.5: Determine the asymptotic decay of the vison interaction.
- T5.3.6: Separate universal power-law contributions from lattice-scale corrections.
- T5.3.7: Compare the result with finite-size or direct numerical calculations.

**Deliverable:** A controlled large-separation expansion.

### M5.4 — Infinite-string limit

**Objective:** Formulate and analyze the limit in which the gauge string becomes arbitrarily long.

This is a distinct operator problem. A string of length \(L\) modifies \(O(L)\) bonds, so its perturbation rank generally grows with \(L\). The finite-rank determinant from the adjacent-vison problem is therefore not automatically a fixed-size determinant in this limit.

Tasks:

- T5.4.1: Define a sequence of finite strings and their corresponding Hamiltonians.
- T5.4.2: Determine which observables are invariant under changes of string path.
- T5.4.3: Establish whether the energy difference has a well-defined \(L\to\infty\) limit at fixed endpoints or fixed endpoint separation.
- T5.4.4: Identify whether the determinant sequence admits a Fredholm determinant representation.
- T5.4.5: Investigate whether a straight string gives a block Toeplitz structure.
- T5.4.6: Determine whether a Wiener–Hopf or related factorization is applicable.
- T5.4.7: Analyze the asymptotic determinant using an appropriate large-\(L\) theorem.
- T5.4.8: Separate endpoint contributions from possible extensive or boundary contributions.

**Deliverable:** A mathematically controlled infinite-string formulation, or a proof that the proposed limit must be reformulated.

# 6. Important conceptual issue: the gauge string

A gauge string connecting two visons is not itself a gauge-invariant physical object. Its endpoints determine the plaquette fluxes; different paths connecting the same endpoints are related by gauge transformations, provided the corresponding gauge transformation is implemented consistently on the matter Majoranas.

Consequently, the infinite-string problem should not be formulated simply as the energy of a string whose length tends to infinity.

There are at least two distinct limits:

1. **Fixed endpoint separation, changing string path or length.** This tests gauge invariance and the equivalence of different representatives of the same flux configuration.
2. **Increasing endpoint separation.** This probes the interaction and asymptotic energy of two physically separated visons.

These limits must not be conflated.

In particular, the infinite separation of two visons does not necessarily imply that a physical energy proportional to the string length survives. In the exactly solvable Kitaev model, the gauge string is not a conventional confining string with a physical tension. The energy is determined by the Majorana spectrum in the corresponding gauge sector.

The relevant asymptotic physical quantity is the energy of the two-vison configuration as their endpoints are separated, compared with the appropriate isolated-vison reference energy.

# 7. Toeplitz and asymptotic determinant strategy

The Toeplitz idea is plausible, but it is not yet an established result for this problem.

Suppose a straight string is oriented along a lattice direction. The perturbation matrix in the space of modified bonds may have entries of the form

\[
M_{mn}(z)
=
\delta_{mn}
-
\lambda\,G_0(\mathbf R_m-\mathbf R_n,z),
\qquad
m,n=1,\ldots,L,
\]

or a block-matrix generalization when each string segment contains several Majorana degrees of freedom.

If translation invariance along the string is preserved, then

\[
M_{mn}(z)=M_{m-n}(z),
\]

and the matrix becomes Toeplitz or block Toeplitz.

The determinant then takes the form

\[
D_L(z)=\det M_L(z).
\]

Potential methods include:

- Toeplitz determinant asymptotics.
- Block Toeplitz factorization.
- Wiener–Hopf factorization.
- Fredholm determinant limits.
- Fisher–Hartwig asymptotics, if the symbol has suitable singularities.
- Transfer-matrix or recurrence methods for special geometries.

However, several conditions must be checked before applying these results:

1. The string geometry must generate translation invariance along its interior.
2. The Green-function kernel must have the required decay or summability properties.
3. The matrix symbol must be correctly defined and its singularities characterized.
4. Endpoint corrections must be separated from the bulk string contribution.
5. The determinant must be defined with an appropriate normalization if the infinite-volume limit is not trace-class.

The gapless Majorana spectrum makes these conditions nontrivial. At real frequencies, the Green-function kernel can have algebraic long-distance behaviour. Even at complex frequency, the decay length can diverge as the spectral parameter approaches the real axis.

Therefore, the appropriate order is:

1. Derive finite-\(L\) determinants.
2. Establish their matrix structure.
3. Study the symbol in the complex-frequency domain.
4. Establish convergence or regularization.
5. Only then derive large-\(L\) asymptotics.

A Toeplitz determinant should be treated as a candidate mathematical structure, not an assumed property of the infinite-string problem.

# 8. Generalizations beyond the isotropic honeycomb model

Once the isotropic result is complete, the same architecture can be extended.

## 8.1 Anisotropic Kitaev honeycomb model

For

\[
\gamma_{\mathbf k}
=
2i(J_z+J_xe^{ik_1}+J_ye^{ik_2}),
\]

the denominator becomes

\[
|\gamma_{\mathbf k}|^2
=
4|J_z+J_xe^{ik_1}+J_ye^{ik_2}|^2.
\]

This corresponds to a weighted or anisotropic triangular-lattice Green function.

The relevant mathematical reference is Horiguchi (1991).

Tasks:

- Derive the anisotropic triangular-lattice mapping.
- Identify the fundamental elliptic integrals.
- Derive the anisotropic local scattering determinant.
- Determine the dependence of the vison energy on \(J_x,J_y,J_z\).
- Analyze the transition between gapless and gapped matter sectors.

## 8.2 Other bipartite lattices

For a general bipartite lattice with

\[
h(\mathbf k)=\sum_j t_j e^{i\mathbf k\cdot\boldsymbol\delta_j},
\]

the local resolvent depends on the distribution of \(|h(\mathbf k)|^2\).

The moments are weighted closed-walk counts. This provides a systematic route to high-frequency expansions, numerical moment methods, and potentially exact Green functions when the associated lattice Green functions are known.

## 8.3 More general finite-rank gauge defects

The finite-rank determinant framework applies whenever a gauge-sector change modifies a finite number of hopping terms.

This suggests applications to:
- multiple local bond reversals;
- more general vison configurations;
- other exactly solvable Majorana spin liquids;
- lattices with more than two sublattices;
- multi-component Majorana scattering problems.

The exact special-function reduction will depend on the lattice, but the scattering determinant formalism is much more general.

# 9. Numerical and analytical validation strategy

Numerical calculations should support the theory, not replace it.

Each major analytical result should be checked through independent routes.

**Green functions**
- Direct two-dimensional Brillouin-zone integration.
- One-dimensional integral reduction.
- Elliptic-integral representation.
- High-frequency moment expansion.
- Real-space recurrence relations.

**Scattering determinant**
- Direct finite-rank matrix evaluation.
- Analytic determinant formula.
- Weak-perturbation expansion.
- Spectral-shift consistency.

**Vison gap**
- Direct integration of the derived phase shift.
- Independent numerical evaluation of the original paper's representation.
- Finite-size diagonalization where practical.
- Convergence analysis with respect to frequency resolution and broadening.

All numerical calculations must record the frequency convention, branch prescription, integration method, and convergence parameters.

# 10. Proposed research-note structure

The theory should be documented in modular notes, rather than one oversized derivation.

| Note | Subject | Main output |
|---|---|---|
| N01 | Model and conventions | Majorana Hamiltonian, gauge sector, perturbation |
| N02 | Finite-rank scattering theory | Resolvent and determinant identities |
| N03 | Spectral shift and energy | Matsubara-to-real-frequency derivation |
| N04 | Honeycomb-to-triangular mapping | Exact Green-function reduction |
| N05 | Triangular local Green function | Elliptic expression and analytic continuation |
| N06 | Local scattering functions | Exact \(g_0\), \(g_2\), and determinant |
| N07 | Adjacent-vison gap | Complete independent derivation |
| N08 | Triangular real-space Green functions | Recurrences and arbitrary-site representation |
| N09 | Honeycomb real-space Green functions | All sublattice components |
| N10 | Finite-separation visons | Finite-rank determinant for two endpoints |
| N11 | Large-separation asymptotics | Interaction and decay laws |
| N12 | Infinite-string limit | Toeplitz/Fredholm/Wiener–Hopf analysis |

Each note should contain:
- definitions and conventions;
- assumptions and domains of validity;
- lemmas and theorems with proofs;
- explicit derivations;
- references to the original mathematical results;
- numerical checks where appropriate;
- unresolved questions and the next tasks.

# 11. Immediate next steps

The first development batch should be limited to Phase I and the beginning of Phase II.

**Batch A**
1. Reconstruct the reference paper's Majorana Hamiltonian and bond-flip perturbation.
2. Fix all Fourier, sublattice, and normalization conventions.
3. Derive the finite-rank scattering determinant.
4. Derive the Matsubara-to-real-frequency phase-shift formula.

**Batch B**
1. Prove the honeycomb-to-triangular Green-function mapping.
2. Derive the one-dimensional integral for \(\mathcal T_{\mathbf 0}(s)\).
3. Establish the exact relation between \(g_0\) and \(g_2\).
4. Identify the elliptic-integral reduction and its analytic branches.

**Batch C**
1. Insert the exact local Green function into the scattering determinant.
2. Derive the phase shift and its zero-temperature integral.
3. Recover the published adjacent-vison gap.
4. Produce independent numerical validation.

Only after Batch C is complete should we begin the arbitrary-site Green functions and separated-vison programme.

# 12. Completion criteria

The main project is complete when:

1. The reference-paper calculation has been rederived from the Majorana Hamiltonian without relying on intermediate formulae that have not been independently verified.
2. The local Green function has an exact elliptic-integral representation with a fully specified analytic continuation.
3. The isotropic relation between \(g_0\) and \(g_2\) has been proved.
4. The finite-rank scattering determinant has been derived with the correct Majorana normalization.
5. The zero-temperature phase-shift integral reproduces the published adjacent-vison gap.
6. Independent numerical calculations validate the analytical expressions.
7. The derivation is documented as a sequence of self-contained mathematical notes.

The extended project is complete when the finite-separation determinant is derived and the large-separation limit is controlled. An infinite-string Toeplitz or Fredholm formulation is a further theoretical result, not a prerequisite for the adjacent-vison calculation.

---

**Central methodological principle:** The physical scattering problem is finite-rank for a local bond reversal. The honeycomb lattice Green function reduces to a triangular-lattice Green function, whose exact mathematical structure is already well studied. The project should exploit this structure systematically, derive the adjacent-vison result independently, and only then address the more delicate large-separation and infinite-string limits.
