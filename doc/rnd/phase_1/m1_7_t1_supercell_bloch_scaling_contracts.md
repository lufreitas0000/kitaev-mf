---
id: RND-P1-M1.7-T1
title: "Bloch Momentum Formalism and Complexity Scaling on N_c-Site Supercells"
phase: 1
milestone: M1.7
task: T1
status: completed
claims:
  - CLM-M1-SUP-001
  - CLM-M1-SUP-002
  - CLM-M1-SUP-003
---

# Bloch Momentum Formalism and Complexity Scaling on $N_c$-Site Supercells

## 1. General Supercell Geometry

Let $n_1, n_2 \in \mathbb{N}$ denote supercell multiples along primitive honeycomb lattice vectors $\mathbf{a}_1, \mathbf{a}_2$.
The enlarged supercell contains:
$$
n_c = 2 n_1 n_2 \text{ physical sites}, \qquad N_b = 3 n_1 n_2 = \frac{3}{2} n_c \text{ bonds}.
$$
The number of Majorana fermion operators in the unit cell is:
$$
M_c = 4 n_c = 8 n_1 n_2.
$$
The state vector $x = (\Omega, \boldsymbol{\lambda})$ has dimension:
$$
D = 6 N_b + 3 n_c = 9 n_c + 3 n_c = 12 n_c \text{ (or } 15 n_c \text{ including bond gauge constraints)}.
$$

---

## 2. Brillouin Zone Folding Equivalence

By Bloch's theorem, an $n_1 \times n_2$ supercell folded into a smaller Brillouin zone $\mathrm{BZ}_{\mathrm{super}} = \mathrm{BZ} / (n_1 n_2)$ satisfies the exact spectral union:
$$
\operatorname{spec}(\mathcal{H}_{\mathrm{super}}(\mathbf{q})) = \bigcup_{j_1=0}^{n_1-1} \bigcup_{j_2=0}^{n_2-1} \operatorname{spec}\left( \mathcal{H}_{\mathrm{prim}}\left( \mathbf{q} + \frac{j_1}{n_1}\mathbf{b}_1 + \frac{j_2}{n_2}\mathbf{b}_2 \right) \right).
$$
Therefore, integrating on a commensurate $k$-mesh with $N_k^{\mathrm{super}} = N_k^{\mathrm{prim}} / (n_1 n_2)$ yields **identical intensive thermodynamic ground-state energy density**:
$$
E_{\mathrm{site}}(n_c) = E_{\mathrm{site}}(2) \pm \epsilon_{\mathrm{mach}}.
$$

---

## 3. Computational Scaling: $\mathcal{O}(N_k n_c^3)$ vs. Fixed System Size $N_{\mathrm{total}}$

For fixed total physical size $N_{\mathrm{total}} = N_k n_c$:
$$
\mathrm{FLOPs}(\Phi) \approx \frac{16}{3} (4 n_c)^3 N_k = \frac{1024}{3} n_c^2 (N_k n_c) = \frac{1024}{3} N_{\mathrm{total}} n_c^2.
$$
- At **fixed $N_{\mathrm{total}}$**, computation scales **quadratically $\propto n_c^2$** with supercell size.
- At **fixed $N_k$**, computation scales **cubically $\propto n_c^3$**.

| $n_c$ | Cell Dimensions ($n_1 \times n_2$) | Majoranas $M_c$ | State Dimension $D$ | $\mathrm{FLOPs}_{\mathrm{eig}} / N_k$ | Full Broyden RAM ($8 D^2$) | L-Broyden RAM ($160 D$) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 2 | $1 \times 1$ | 8 | 30 | $2.73 \times 10^3$ | $7.2\text{ KB}$ | $4.8\text{ KB}$ |
| 4 | $2 \times 1$ | 16 | 60 | $2.18 \times 10^4$ | $28.8\text{ KB}$ | $9.6\text{ KB}$ |
| 8 | $2 \times 2$ | 32 | 120 | $1.75 \times 10^5$ | $115.2\text{ KB}$ | $19.2\text{ KB}$ |
| 16 | $4 \times 2$ | 64 | 240 | $1.40 \times 10^6$ | $460.8\text{ KB}$ | $38.4\text{ KB}$ |
| 32 | $4 \times 4$ | 128 | 480 | $1.12 \times 10^7$ | $1.84\text{ MB}$ | $76.8\text{ KB}$ |

---

## 4. Pre-Registered Hypotheses for Supercells

- **`H-13`:** Profiled FLOPs and CPU wall-clock times confirm the cubic scaling $\mathcal{O}(N_k n_c^3)$ with empirical exponent within $\pm 0.08$ of 3.0.
- **`H-14`:** Under $C_3$-symmetric seeds, intensive energy $|E_{\mathrm{site}}(n_c) - E_{\mathrm{site}}(2)| < 10^{-11}$.
- **`H-15`:** Under unconstrained generic perturbations at $h=0.20$, supercells $n_c \ge 4$ can spontaneously break translation or rotational symmetry into stripy or dimer states with higher variational energy.
