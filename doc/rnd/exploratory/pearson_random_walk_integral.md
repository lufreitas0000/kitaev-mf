---
id: RND-EXP-PEARSON-WALK
title: "Pearson Random Walk and Brillouin-Zone Correlator Integral"
phase: 1
milestone: M1.3
status: completed
questions:
  - "What is the closed-form evaluation of the pure Kitaev thermodynamic matter bond correlator?"
contracts:
  - A06
claims:
  - K-006
artifacts:
  - scripts/k006_quadrature.py
  - data/processed/benchmarks/k006_quadrature.dat
---

# Pearson Random Walk and Brillouin-Zone Correlator Integral

Yes. This integral is a classical object in the theory of random walks.

**Identification.** Write $\theta_j = 2\pi k_j$. Then

$$3\bar U = \int_{[0,1)^2}\!\left|1+e^{i\theta_1}+e^{i\theta_2}\right| dk_1\,dk_2 = \mathbb{E}\,|X_3|,$$

which is the **mean distance from the origin after a 3-step planar random walk with unit steps and uniformly random directions** (the Pearson random walk). The first step can be fixed to $1$ by rotation invariance, which leaves the other two angles uniform on $[0,2\pi)$. In the notation of Borwein–Straub–Wan–Zudilin, "Densities of short uniform random walks" (Canad. J. Math. 2012),

$$\bar U = \frac{W_3(1)}{3}, \qquad W_3(1) = 1.574597237\ldots,$$

and indeed $1.5745972/3 = 0.5248657$, matching your Richardson value.

**Closed form.** $W_3(1)$ is known exactly in terms of $\Gamma(1/3)$. Numerically, the following reproduces $1.5746$ (I checked it to about 4–5 digits, so confirm the exact coefficients against the paper before quoting it):

$$W_3(1) = \frac{3\cdot 2^{1/3}}{16\,\pi^4}\,\Gamma^6\!\left(\tfrac13\right) + \frac{16\cdot 2^{2/3}\,\pi^2}{\Gamma^6\!\left(\tfrac13\right)}.$$

**Further structure.**
- The probability density $p_3(x)$ of the walk's distance is expressed via the complete elliptic integral $K$. It has a logarithmic singularity at $x=1$, where $|1+e^{i\theta_1}+e^{i\theta_2}|=1$ is the saddle-type critical value. This corresponds to the **van Hove singularity** of the honeycomb band structure.
- $W_3(1)$ can also be written as a ${}_3F_2$ hypergeometric value and connects to Bessel moments, $W_3(s)=\int_0^\infty s\,x\,J_0(x)^3 (xs)^{\ldots}$-type integrals, so it is also a Bessel-function moment, $\int_0^\infty J_0(x)^3 J_1(x)\,dx$ style.
- The modulus $|1+e^{i\theta_1}+e^{i\theta_2}|$ is the magnitude of the honeycomb structure factor $f(\mathbf k)$, so $\bar U$ is the Brillouin-zone-averaged $|f|/3$, that is, the mean single-particle graphene/Majorana band energy normalized by its maximum $3$.

Searching for "$W_3(1)$" or "mean distance of a 3-step random walk" will lead you to the relevant literature (Pearson 1905, Kluyver 1906, Borwein et al.).


---


## The 3-step walk

Take a walker in the plane who makes $n$ steps of unit length, each in an independent uniformly random direction. This is the **Pearson random walk** (the "random flight" in 2D). Its endpoint is

$$X_n=\sum_{k=1}^n e^{i\theta_k},\qquad \theta_k\sim U[0,2\pi),$$

and the quantity of interest is the distance $|X_n|\in[0,n]$. For $n=3$ the endpoint satisfies $|X_3|\le 3$, with $|X_3|=3$ only when all three steps align. This is the $\Gamma$ point of the honeycomb Brillouin zone, where $f(\mathbf k)=3$. Because of rotation invariance, you can fix $\theta_1=0$, which gives exactly your integrand $|1+e^{i\theta_2}+e^{i\theta_3}|$.

## Definition of $W_n(s)$

$$W_n(s)=\int_{[0,1]^n}\Big|\sum_{k=1}^n e^{2\pi i x_k}\Big|^s\,dx_1\cdots dx_n=\mathbb E\,|X_n|^s ,$$

that is, the $s$-th moment of the distance. It is defined for complex $s$, and it is analytic for $\operatorname{Re}s>-2$ when $n\ge 3$ ($\operatorname{Re}s>-1$ for $n=2$). Equivalently, $W_n(s)=\int_0^n x^s p_n(x)\,dx$, where $p_n$ is the density of $|X_n|$.

**Bessel representation** (Kluyver 1906). The density is

$$p_n(x)=\int_0^\infty x\,t\,J_0(xt)\,J_0(t)^n\,dt ,$$

so $W_n(s)$ is a Mellin-type moment of $J_0^n$. A useful special case is $W_n(-1)=\int_0^\infty J_0(x)^n\,dx$. This replaces the loose "$J_0^3J_1$" remark in my previous message, which was not a correct statement.

## Known closed forms

- $n=1$: $W_1(s)=1$.
- $n=2$: $W_2(s)=\binom{s}{s/2}$, so $W_2(1)=4/\pi$.
- **Even moments** ($s=2k$) for any $n$:

$$W_n(2k)=\sum_{a_1+\dots+a_n=k}\binom{k}{a_1,\dots,a_n}^2 .$$

  This is an integer. For $n=3$ these are $1,3,15,93,639,\dots$ (Domb numbers), and they satisfy a linear recurrence.
- $n=3$, general $s$ (Borwein–Straub–Wan–Zudilin):

$$W_3(s)=\frac{3^{s+3/2}}{2\pi}\frac{\Gamma(1+\frac s2)^2}{\Gamma(s+2)}\;{}_3F_2\!\left(\tfrac s2+1,\tfrac s2+1,\tfrac s2+1;\,1,\tfrac{s+3}{2};\,\tfrac14\right).$$

  I am quoting this from memory, so check the parameters against the paper. It does pass the check $W_3(0)=1$.
- $n=3$, $s=1$: this is your constant,

$$W_3(1)=\frac{3\cdot2^{1/3}}{16\pi^4}\Gamma^6\!\left(\tfrac13\right)+\frac{16\cdot2^{2/3}\pi^2}{\Gamma^6(\frac13)}=1.5745972\ldots$$

  Numerically the first term is $0.8964\ldots$ and the second is $0.6781\ldots$. The first term is $W_3(-1)=\int_0^\infty J_0^3$.
- $n=4$: $W_4(s)$ also has a ${}_4F_3$ representation at argument $1/16$, and special values involve $\pi$ and $\zeta(3)$, but it is more intricate.

## Is there a general closed form?

No. For $n\ge5$ there is no known closed form, even at $s=1$. What does exist is:

1. **Hypergeometric or Meijer-$G$ series** in $s$ for each fixed $n$. Borwein et al. give a series of this type for $W_n(s)$ that converges for all $s$ not an odd integer below $-2$.
2. **Linear recurrences in $s$.** $W_3(s)$ and $W_4(s)$ satisfy them, relating $W(s+4)$, $W(s+2)$ and $W(s)$ with polynomial coefficients. This gives all odd moments from two seed values.
3. **Asymptotics.** For large $n$, $X_n$ is approximately a 2D Gaussian (Rayleigh distribution), so $W_n(s)\sim n^{s/2}\,\Gamma(1+s/2)$.

The standard references are Borwein–Straub–Wan–Zudilin (*Canad. J. Math.* 64, 2012) and Borwein–Straub–Wan, "Three-step and four-step random walk integrals" (*Experimental Math.* 2013).


---

Yes. This integral has a natural interpretation in several areas of mathematics, particularly planar random walks, lattice Green functions, and special-function theory. There is also an important distinction from the Mahler measure of the honeycomb-lattice polynomial.

Let us define

\\[ I_3 = \int_0^1\int_0^1 \left|1+e^{-2\pi i k_1}+e^{-2\pi i k_2}\right| \\,dk_1\\,dk_2, \\]

so that \\(\bar U=I_3/3\\).

The first important observation is that this is not an arbitrary Brillouin-zone integral.

## 1. Connection to a three-step planar random walk

Consider a planar random walk with three unit-length steps,

\\[ \mathbf R_3 = \mathbf e_1+\mathbf e^{i\theta_1} +\mathbf e^{i\theta_2}, \\]

where the first step is fixed by rotational invariance and the other two angles are independently uniform on \\([0,2\pi)\\).

Its mean end-to-end distance is

\\[ \mathbb E[|\mathbf R_3|] = \frac{1}{(2\pi)^2} \int_0^{2\pi}\int_0^{2\pi} |1+e^{i\theta_1}+e^{i\theta_2}| \\,d\theta_1d\theta_2. \\]

Changing variables to \\(\theta_j=2\pi k_j\\), we obtain exactly

\\[ \boxed{ \bar U=\frac13\mathbb E[|\mathbf R_3|]. } \\]

Thus, your integral is one-third of the mean distance from the origin after three independent unit steps in a two-dimensional isotropic random walk.

This gives it an established probabilistic interpretation.

## 2. An exact evaluation in terms of Gamma functions

There is a stronger result: the three-step random-walk problem has an exact closed-form solution.

Writing \\(W_3(1)=\mathbb E[|\mathbf R_3|]\\), the result is

\\[ \boxed{ W_3(1)= \frac{3\\,2^{1/3}}{16\pi^4} \Gamma\left(\frac13\right)^6 + \frac{27\\,2^{2/3}}{4\pi^4} \Gamma\left(\frac23\right)^6 . } \\]

Consequently, your integral is exactly

\\[ \boxed{ \begin{aligned} \bar U={}& \frac{2^{1/3}}{16\pi^4} \Gamma\left(\frac13\right)^6\\\ &+ \frac{9\\,2^{2/3}}{4\pi^4} \Gamma\left(\frac23\right)^6 . \end{aligned} } \\]

Numerically,

\\[ \boxed{\bar U=0.524865745850594\ldots} \\]

This is not merely a numerical coincidence.

The result is established in the literature on short planar random walks:

Jonathan M. Borwein, Dirk Nuyens, Armin Straub and James Wan, Some arithmetic properties of short random walk integrals, The Ramanujan Journal 26, 109–132 (2011).

DOI: [10.1007/s11139-011-9325-y](https://doi.org/10.1007/s11139-011-9325-y?utm_source=chatgpt.com).

The authors derive closed forms for the moments of three-step and four-step planar random walks, including the expected distance after three steps.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://arminstraub.com\&sz=32)

arminstraub.com+1



## 3. A Bessel-function representation

There is another useful representation, particularly interesting from a mathematical-physics perspective.

The characteristic function of a uniformly distributed unit vector in two dimensions is

\\[ \mathbb E[e^{i\mathbf q\cdot\mathbf e}] =J_0(|\mathbf q|), \\]

where \\(J_0\\) is the Bessel function of the first kind.

For three independent steps, independence gives

\\[ \mathbb E[J_0(t|\mathbf R_3|)] =J_0(t)^3. \\]

Using the identity

\\[ r=\int_0^\infty \frac{1-J_0(tr)}{t^2}\\,dt, \qquad r\geq 0, \\]

we obtain

\\[ \boxed{ W_3(1)= \int_0^\infty \frac{1-J_0(t)^3}{t^2}\\,dt. } \\]

Therefore,

\\[ \boxed{ \bar U= \frac13\int_0^\infty \frac{1-J_0(t)^3}{t^2}\\,dt. } \\]

This connects your Brillouin-zone integral to Bessel-function integrals, which are important in random-walk theory, Fourier analysis and lattice models.

## 4. Connection to the honeycomb lattice and Mahler measure

There is another interesting connection.

Define the Laurent polynomial

\\[ P(x,y)=1+x+y. \\]

Your integral is its \\(L^1\\) torus norm:

\\[ I_3=\int\_{\mathbb T^2}|P(x,y)|\\,d\mu, \\]

where \\(d\mu\\) is normalized Haar measure.

The polynomial \\(1+x+y\\) also appears in the theory of Mahler measures. However, the Mahler measure is defined using the logarithm:

\\[ m(P)= \int\_{\mathbb T^2} \log|1+x+y|\\,d\mu. \\]

Unlike your integral, this quantity has a particularly elegant arithmetic evaluation, due to Smyth:

\\[ \boxed{ m(1+x+y)= \frac{3\sqrt3}{4\pi} L(\chi\_{-3},2). } \\]

Here \\(L(\chi\_{-3},s)\\) is a Dirichlet \\(L\\)-function.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://academic.oup.com\&sz=32)

Oxford Academic+1



Thus, the two integrals involve the same polynomial but represent different moments:

- Your integral: first absolute moment, \\(L^1\\).
- Mahler measure: logarithmic moment.

They should not be identified with one another.

## 5. Why this is interesting for the Kitaev model

In the isotropic Kitaev honeycomb model, the flux-free Majorana dispersion has the form

\\[ \varepsilon(\mathbf k) \propto \left|1+e^{i k_1}+e^{i k_2}\right|, \\]

with the proportionality factor depending on the Hamiltonian and Fourier-transform conventions.

Your integral therefore measures the Brillouin-zone average of the dimensionless dispersion.

There is a useful conceptual distinction here:

- The Dirac points determine the low-energy behavior.
- The integral above measures a global property of the entire Brillouin zone.
- The random-walk interpretation provides an exact mathematical evaluation independent of the low-energy approximation.

My recommendation: for your Kitaev-model research notes, I would record this result under the name three-step planar random-walk integral, with the Gamma-function evaluation and Bessel representation.

The primary reference is the 2011 Borwein–Nuyens–Straub–Wan paper. It provides a rigorous mathematical foundation for the numerical value, rather than treating it as a special numerical constant obtained from Brillouin-zone integration.


---

Yes, and the reorganized argument is cleaner than the document's. The central step is that the **real-space** Green function reduces to one scalar function, the triangular-lattice Green function $T(\mathbf R;w)$, evaluated at $w=z^2-3$. I also have to correct something in my earlier message.

**Correction.** I called the numbers $1,3,15,93,639,\dots$ "Domb numbers". That was wrong. They are $a_n=\sum_k\binom nk^2\binom{2k}{k}$, which count closed walks on the honeycomb lattice. The Domb numbers $1,4,28,256,\dots$ are $W_4(2k)$.

## 1. Matrix structure (with $z=i\omega_n$, $K=1$ units)

Take $f(\mathbf k)=1+e^{ik_1}+e^{ik_2}$ and $H_{\mathbf k}=\begin{pmatrix}0& iKf\\ -iKf^*&0\end{pmatrix}$. Then $H_{\mathbf k}^2=K^2|f|^2$, so

$$G_0(z,\mathbf k)=\frac{z+H_{\mathbf k}}{z^2-K^2|f|^2}.$$

The Green function therefore has two scalar building blocks:
- the diagonal part, $z/(z^2-K^2|f|^2)$;
- the off-diagonal part, $\pm iK f^{(*)}/(z^2-K^2|f|^2)$.

## 2. Reduction to the triangular lattice

Expanding the modulus,

$$|f|^2=3+\varepsilon_{\mathbf k},\qquad \varepsilon_{\mathbf k}=2\cos k_1+2\cos k_2+2\cos(k_1-k_2).$$

This $\varepsilon_{\mathbf k}$ is the nearest-neighbor triangular-lattice dispersion with $t=1$ (six neighbors). Define

$$T(\mathbf R;w)=\int\!\frac{d^2k}{(2\pi)^2}\frac{e^{i\mathbf k\cdot\mathbf R}}{w-\varepsilon_{\mathbf k}},\qquad w=\frac{z^2}{K^2}-3 .$$

Then $\dfrac1{z^2-K^2|f|^2}=\dfrac1{K^2}\dfrac1{w-\varepsilon}$.

## 3. Real-space blocks

The factor $f\,e^{i\mathbf k\cdot\mathbf R}$ is a sum of three shifted phases, $e^{i\mathbf k\cdot\mathbf R}+e^{i\mathbf k\cdot(\mathbf R+\mathbf e_1)}+e^{i\mathbf k\cdot(\mathbf R+\mathbf e_2)}$. Hence, up to sign and phase conventions for $\mathbf R$ and the sublattice offset,

$$G_{AA}(\mathbf R)=G_{BB}(\mathbf R)=\frac{z}{K^2}\,T(\mathbf R;w),$$

$$G_{AB}(\mathbf R)=\frac{i}{K}\big[T(\mathbf R)+T(\mathbf R+\mathbf e_1)+T(\mathbf R+\mathbf e_2)\big],$$

$$G_{BA}(\mathbf R)=-\frac{i}{K}\big[T(\mathbf R)+T(\mathbf R-\mathbf e_1)+T(\mathbf R-\mathbf e_2)\big].$$

So the full real-space $G_0$ is a finite combination of the single function $T$ on the triangular lattice. This is the exact statement of the bipartite reduction.

## 4. The local identity, derived more simply

The document derived $g_2=(zg_0-1)/3$ from permutation symmetry of the three unit steps. The lattice equation $(w-\varepsilon)T=\delta_{\mathbf R,0}$ gives the same identity more directly:

$$wT(0)-6T(\mathrm{nn})=1\;\Rightarrow\; T(0)+2T(\mathrm{nn})=\frac{(w+3)T(0)-1}{3}=\frac{z^2T(0)-1}{3},$$

where $T(\mathrm{nn})=T(\mathbf e_1)=T(\mathbf e_2)$ by the $C_6$ symmetry. With $g_0=zT(0)$ (setting $K=1$) this is exactly $G_{AB}(0)/i=(zg_0-1)/3$. The identity is therefore just the lattice Green-function equation, and the random-walk symmetry is not needed.

At $z\to0$ this gives $\langle \mathrm{Re}f/|f|^2\rangle=1/3$, so the $\tau$-component has magnitude $1/(3K)$. I could not check the paper's Eq. (33), which the document says gives $-2/3$. The mismatch may come from a different normalization of $g_2$ or of $V$ (for example a factor of 2). I would check that convention before concluding the paper has an error.

## 5. What you get for Matsubara frequencies

For $z=i\omega_n$ we have $w=-\omega_n^2/K^2-3\le-3$. This lies **below** the triangular band $[-3,6]$, so:
- $T(\mathbf R;w)$ is real and analytic, with no branch-cut ambiguity.
- $T$ decays exponentially in $|\mathbf R|$.
- The correlation length diverges as $\omega_n\to0$, because $w\to-3$ is the band bottom of the triangular lattice, which corresponds to the Dirac points.

The analytic continuation $i\omega_n\to\omega+i0^+$ is needed only at the very end, for the phase shift.

## 6. How to compute $T(\mathbf R;w)$

- **On-site.** $T(0;w)$ is a complete elliptic integral (Horiguchi 1972; Morita–Horiguchi 1972). I have not verified the explicit modulus, so take it from those papers.
- **Off-site.** The lattice equation plus lattice symmetry gives a recursion in $\mathbf R$. This is the standard route by which the off-site values reduce to rational combinations of $K(\kappa)$ and $E(\kappa)$ with coefficients rational in $w$. I am citing this from the literature and have not worked it out here.
- **Series.** For large $|z|$,

$$T(0)=\frac{1}{z^2}\sum_{n\ge0}\frac{a_n}{z^{2n}},\qquad a_n=1,3,15,93,639,\dots$$

The $a_n$ satisfy $(n+1)^2a_{n+1}=(10n^2+10n+3)a_n-9n^2a_{n-1}$. I checked this for the first terms. Its growth rates $1$ and $9$ correspond to singularities at $z^2=1$ (the van Hove point) and $z^2=9$ (the band edge).

## 7. Role of the random-walk results

- The **even** moments $\langle|f|^{2n}\rangle=a_n$ count closed walks on the honeycomb lattice. The off-site analogue $\langle e^{i\mathbf k\cdot\mathbf R}|f|^{2n}\rangle$ counts open walks between two sites.
- The **odd** moment $\bar U=\langle|f|\rangle/3$ has no walk interpretation. It does not enter $G_0$ at all. It appears in energy-type quantities such as $\sum|\epsilon_{\mathbf k}|$.

## 8. Suggested next step

Choose between two routes:
1. Take $T(0;w)$ and $T(\mathrm{nn};w)$ in closed elliptic form from Horiguchi, which gives the local $g_0,g_2$ and hence $D(z)$ analytically.
2. Evaluate $T(\mathbf R;w)$ numerically for $w\le-3$ by a one-dimensional integral, which is stable and needs no continuation.

---
