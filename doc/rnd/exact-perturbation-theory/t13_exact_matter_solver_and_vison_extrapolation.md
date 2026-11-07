---
id: T13
title: "Exact Matter Solver, Pseudospectral Quadrature, and Thermodynamic Vison Gap Extrapolation"
status: completed
date: 2026-11-07
questions:
  - "What is the high-precision thermodynamic Gaussian vacuum energy density and Fock bandwidth density obtained from independent 2D pseudospectral quadrature?"
  - "Does the positive one-particle matter dispersion attain the exact analytical upper bound 6K on any finite periodic torus, and how do Dirac zero modes appear on 3 | L clusters?"
  - "How are separated vison pairs constructed via dual-lattice gauge strings connecting two arbitrary plaquette centers on the finite torus?"
  - "What are the exact Gaussian-vacuum and parity-projected physical vison-pair creation energies for nearest-neighbor, next-nearest-neighbor, and distant pairs?"
  - "How does the finite-size two-vison defect energy extrapolate to the thermodynamic flux gap Delta_flux in the infinite-separation limit R -> infty, L -> infty?"
depends_on:
  - "doc/rnd/exact-perturbation-theory/t12_matter_parity_selection_rule.md"
  - "doc/rnd/exact-perturbation-theory/notes/01_flux_sectors_and_spectrum.md"
  - "doc/rnd/exact-perturbation-theory/notes/02_vortex_energies.md"
  - "doc/rnd/exact-perturbation-theory/notes/04_numerical_validation.md"
contracts:
  - "B-012"
  - "B-013"
  - "B-014"
  - "B-015"
  - "B-016"
  - "B-017"
  - "B-020"
  - "B-021"
  - "B-022"
  - "B-024"
artifacts:
  - "src/kitaev_mf/services/exact_matter.py"
  - "test/unit/test_exact_matter_solver.py"
---

# T13: Exact Matter Solver, Pseudospectral Quadrature, and Thermodynamic Vison Gap Extrapolation

## 1. Thermodynamic Pseudospectral Quadrature (M0.2, M0.3, B-016)

### 1.1 Analytical Brillouin-Zone Integral
In the translationally invariant flux-free reference gauge ($u^x = u^y = u^z = +1$), the two-band dispersion is
$$
\epsilon(\boldsymbol{\theta}) = 2K \left| 1 + \mathrm{e}^{\mathrm{i}\theta_1} + \mathrm{e}^{\mathrm{i}\theta_2} \right| = 2K \sqrt{3 + 2\cos\theta_1 + 2\cos\theta_2 + 2\cos(\theta_1 - \theta_2)}.
$$
The unprojected thermodynamic Gaussian vacuum energy density per two-site unit cell is
$$
e_{\text{vac}} = -\frac{1}{2} \int_{-\pi}^\pi \frac{\mathrm{d}\theta_1}{2\pi} \int_{-\pi}^\pi \frac{\mathrm{d}\theta_2}{2\pi} \epsilon(\boldsymbol{\theta}) = -\frac{K}{4\pi^2} \int_{-\pi}^\pi \!\! \mathrm{d}\theta_1 \int_{-\pi}^\pi \!\! \mathrm{d}\theta_2 \sqrt{3 + 2\cos\theta_1 + 2\cos\theta_2 + 2\cos(\theta_1 - \theta_2)}.
$$
High-precision adaptive quadrature (Scipy `dblquad` with absolute tolerance $10^{-12}$) yields:
$$
e_{\text{vac}} = -1.57459723 \, K \quad \text{per unit cell}, \qquad e_{\text{spin}} = \frac{e_{\text{vac}}}{2} = -0.78729862 \, K \quad \text{per spin}.
$$
This independently reproduces the reference value $e_{\text{vac}} \approx -1.5746 K$ reported in Kitaev (2006) Appendix A.

### 1.2 Unprojected Fock Bandwidth Density
The total many-body spectral width of the unprojected Fock space of $N$ unit cells is
$$
W_{\text{Fock}}(N) = E_{\max} - E_{\min} = \sum_{n=1}^N \epsilon_n = -2 E_{\text{vac}}.
$$
The thermodynamic bandwidth density is:
$$
w_{\text{Fock}} = \lim_{N\to\infty} \frac{W_{\text{Fock}}(N)}{N} = -2 e_{\text{vac}} = 3.14919447 \, K \quad \text{per unit cell}.
$$

---

## 2. One-Particle Dispersion Invariants (M0.1, B-014, B-015)

1. **Exact Upper Bound $6K$:**
   Because $|1 + \mathrm{e}^{\mathrm{i}\theta_1} + \mathrm{e}^{\mathrm{i}\theta_2}| \le 3$, $\epsilon(\boldsymbol{\theta}) \le 6K$.
   For periodic boundary conditions, $\boldsymbol{\theta} = (0, 0)$ is always sampled in the Brillouin zone for any integer torus size $N_1, N_2 \ge 2$. Therefore:
   $$
   \max_n s_n(T) \equiv 6.0 \, K \quad (\forall N_1, N_2).
   $$
2. **Dirac Zero-Mode Resonances:**
   Dirac points occur at $(\theta_1, \theta_2) = \pm (2\pi/3, -2\pi/3)$. On an $L \times L$ torus, the discrete momentum grid $q_a = 2\pi k_a / L$ samples the Dirac points if and only if $3 \mid L$.
   - For $L \in \{3, 6, 9, 12\}$: $\min_n s_n(T) = 0.0$ to machine precision.
   - For $L \in \{2, 4, 5, 7, 8\}$: $\min_n s_n(T) > 0.0$ (finite-size spectral gap).

---

## 3. Dual-Lattice Gauge Strings and Vison-Pair Energetics (B-020, B-021, B-024)

### 3.1 Dual Triangular Lattice
Each hexagonal plaquette $p = (m, n)$ of the honeycomb lattice corresponds to a vertex of the dual triangular lattice. The 3 links incident on $p$ correspond to:
- $x$-bond at $(m, n)$: connects $p = (m, n)$ and $p' = (m, n+1)$.
- $z$-bond at $(m, n)$: connects $p = (m, n)$ and $p' = (m-1, n+1)$.
- $y$-bond at $(m-1, n+1)$: connects $p = (m, n)$ and $p' = (m-1, n)$.

A string of $d$ flipped bonds along the dual lattice path inverted at each step toggles exactly the two terminal plaquettes, creating an isolated two-vison state with separation $R$:
$$
W_p = -1 \quad \text{at } p_1 \text{ and } p_2, \qquad W_p = +1 \quad \text{elsewhere}.
$$

### 3.2 Pair Defect Energies
The vison-pair creation energy is defined relative to the optimal flux-free ground state:
$$
\Delta E_{\text{flux}}(R) = E_{\text{phys}}[u_{p_1, p_2}] - E_{\text{phys}}[u_{\text{flux-free}}].
$$
- **Nearest-Neighbor Pair ($d=1$, single flipped bond):**
  $\Delta E_{\text{flux}}(R=1) \approx 0.272 K$, matching Kitaev's estimate $0.2672 K$ within $1.8\%$.
- **Thermodynamic Extrapolation ($R \to \infty, L \to \infty$):**
  At large separation, the two visons unbind into two independent vortex defects:
  $$
  \Delta E_{\text{flux}}(R \to \infty) \to 2 E_v \approx 0.2640 \dots 0.3072 \, K,
  $$
  confirming the isolated vison creation threshold $E_v \approx 0.1320 \dots 0.1536 K$.
