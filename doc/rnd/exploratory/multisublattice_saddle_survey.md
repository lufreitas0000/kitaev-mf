---
id: RND-EXP-SADDLE-SUB
title: "Saddle-Point Algorithms, Multi-Sublattice Order, and Parton Decoupling in the Extended Kitaev Model"
phase: 1
milestone: M1.6
status: active
questions:
  - "How do multi-sublattice magnetic orders and isospin Lagrange multiplier saddle points couple in the extended Kitaev model?"
contracts:
  - T1
  - T2
---

# Saddle-Point Algorithms, Multi-Sublattice Order, and Parton Decoupling in the Extended Kitaev Model

Mean-field treatments of quantum spin liquids and proximate magnetically ordered phases in bond-directional $J\text{-}K\text{-}\Gamma\text{-}\mathbf{h}$ systems, such as $\alpha\text{-RuCl}_3$ and honeycomb iridates, require careful management of auxiliary Hilbert space redundancies. In the foundational monograph of Lu Freitas (Universidade Federal do Rio Grande do Norte / Heinrich-Heine-Universität Düsseldorf) and the associated publication Freitas et al., Majorana mean-field theory (MMFT) provides an effective framework for computing spatial polarizations and local spin correlations around $\mathbb{Z}_2$ vortex cores in the extended Kitaev model. The quantitative accuracy of these parton methods depends on how the unphysical states introduced by fractionalization are eliminated. When local spin operators are expressed in terms of four Majorana fermions per lattice site, the Hilbert space dimension doubles from $2$ to $4$, necessitating an exact or average projection via local $SO(4)$ or isospin $SU(2)$ constraints. Benchmarking these frameworks requires analyzing the mathematical formulations of the underlying saddle-point equations, the multi-sublattice geometries hosting competing magnetic structures, and the variations introduced by distinct Wick decoupling schemes.

## Lagrange Multiplier Saddle-Point Formulations and Root-Finding Algorithms

### Microscopic Gauge Structure and the $SO(4)$ Constraint Algebra

The general extended Kitaev Hamiltonian on the tricoordinated honeycomb lattice is parametrized as:

$$H = \sum_{\langle ij \rangle_\gamma} \left[ K_\gamma S_i^\gamma S_j^\gamma + J \mathbf{S}_i \cdot \mathbf{S}_j + \Gamma (S_i^\alpha S_j^\beta + S_i^\beta S_j^\alpha) \right] - \sum_i \mathbf{h} \cdot \mathbf{S}_i$$

where $\gamma \in \{x, y, z\}$ designates the bond directional index, and $(\alpha, \ beta, \gamma)$ form an even permutation of Cartesian components. In the four-Majorana representation formulated by Kitaev and generalized by Seifert, Meng, and Vojta, each spin-$1/2$ degree of freedom is replaced by a matter Majorana fermion $c_i$ and three gauge Majorana fermions $\mathbf{b}_i = (b_i^x, b_i^y, b_i^z)^T$ satisfying the Clifford algebra $\{c_i, c_j\} = 2\delta_{ij}$, $\{b_i^\alpha, b_j^\beta\} = 2\delta^{\alpha\beta}\delta_{ij}$, and $\{c_i, b_j^\beta\} = 0$. Collecting these into a four-component Majorana multiplet $\boldsymbol{\eta}_i = (c_i, b_i^x, b_i^y, b_i^z)^T$, the local gauge redundancy forms an $SO(4) \cong [SU(2) \times SU(2)] / \mathbb{Z}_2$ group structure.

Physical spin states occupy the two-dimensional subspace defined by the local $\mathbb{Z}_2$ parity operator:

$$D_i = -i b_i^x b_i^y b_i^z c_i = +1$$

which corresponds to the eigenvalue constraint of the local fermion parity. In continuous gauge formulations, this physical sector corresponds to the singlet condition for three local gauge generators. Defining the three real skew-symmetric $4 \times 4$ generators $G^a$ ($a=1,2,3$) acting on $\boldsymbol{\eta}_i$:

$$G^1 = \begin{pmatrix} 0 & 1 & 0 & 0 \\ -1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & -1 & 0 \end{pmatrix}, \quad G^2 = \begin{pmatrix} 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & -1 \\ -1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \end{pmatrix}, \quad G^3 = \begin{pmatrix} 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \\ 0 & -1 & 0 & 0 \\ -1 & 0 & 0 & 0 \end{pmatrix}$$

the generators satisfy the $\mathfrak{so}(3)$ gauge Lie algebra $[G^a, G^b] = -2\epsilon^{abc} G^c$ and commute with physical spin operators. The local operator constraints enforce vanishing gauge charges:

$$\hat{g}_i^a \equiv -\frac{i}{4} \boldsymbol{\eta}_i^T G^a \boldsymbol{\eta}_i = 0, \quad a \in \{1, 2, 3\}$$

In the explicit mean-field implementations of Yilmaz, Kampf, and Yip and Ralko and Merino, this projection is formulated by requiring that the pairwise contractions between the three gauge Majoranas vanish on average:

$$\hat{\Gamma}_i^\gamma \equiv \frac{i}{2} \epsilon^{\alpha\beta\gamma} b_i^\alpha b_i^\beta = 0$$

These constraints eliminate unphysical on-site gauge-charge fluctuations while preserving the mean-field bilinear structure. To enforce these conditions at the saddle-point level, one introduces real local Lagrange multipliers $\boldsymbol{\lambda}_i = (\lambda_i^1, \lambda_i^2, \lambda_i^3)^T$ directly into the Hamiltonian:

$$H_\lambda = \sum_i \sum_{a=1}^3 \lambda_i^a \hat{g}_i^a = \frac{i}{2} \sum_i \sum_\gamma \lambda_i^\gamma \epsilon^{\alpha\beta\gamma} b_i^\alpha b_i^\beta$$

### Variational Free Energy and the Saddle-Point Min-Max Problem

Mean-field decoupling of all two-body spin interactions in the matter, gauge, and hybrid channels yields an effective quadratic Hamiltonian in terms of the Majorana operators:

$$H_{MF}(\boldsymbol{\Phi}, \boldsymbol{\lambda}) = \frac{i}{4} \sum_{i,j} \boldsymbol{\eta}_i^T \mathcal{A}_{ij}(\boldsymbol{\Phi}, \boldsymbol{\lambda}) \boldsymbol{\eta}_j + E_0(\boldsymbol{\Phi})$$

Here, $\boldsymbol{\Phi} = \{\chi_{ij}, u_{ij}^\gamma, \kappa_{ij}^{\alpha\beta}, \mathbf{m}_i\}$ represents the set of mean-field parameters, where $\chi_{ij} = \langle i c_i c_j \rangle$ denotes matter fermion hoppings, $u_{ij}^\gamma = \langle i b_i^\gamma b_j^\gamma \rangle$ are gauge link variables, $\kappa_{ij}^{\alpha\beta} = \langle i c_i b_j^\beta \rangle$ denote hybrid cross-correlations, and $\mathbf{m}_i = \langle \mathbf{S}_i \rangle$ are local magnetizations. The matrix $\mathcal{A}(\boldsymbol{\Phi}, \boldsymbol{\lambda})$ is real and skew-symmetric ($\mathcal{A}^T = -\mathcal{A}$), and $E_0(\boldsymbol{\Phi})$ denotes the classical energy shift resulting from Wick contractions.

At finite temperature $T = \beta^{-1}$, the grand canonical free energy is given by integrating out the Majorana fermions:

$$F(\boldsymbol{\Phi}, \boldsymbol{\lambda}) = E_0(\boldsymbol{\Phi}) - \frac{1}{2\beta} \sum_{k} \ln \left[ 2 \cosh\left( \frac{\beta E_k(\boldsymbol{\Phi}, \boldsymbol{\lambda})}{2} \right) \right]$$

where $\{\pm E_k\}$ ($E_k \ge 0$) are the eigenvalues of the Hermitian Bogoliubov–de Gennes (BdG) Hamiltonian $\mathcal{H}_{BdG} = i \mathcal{A}$. The physical ground state corresponds to a stationary saddle point defined by the min-max variational principle:

$$\min_{\boldsymbol{\Phi}} \max_{\boldsymbol{\lambda}} F(\boldsymbol{\Phi}, \boldsymbol{\lambda})$$

The stationarity conditions generate a set of coupled self-consistent equations:

$$\frac{\partial F}{\partial \Phi_m} = 0 \implies \Phi_m = \langle \hat{\mathcal{O}}_m \rangle_{\boldsymbol{\Phi}, \boldsymbol{\lambda}}$$

$$g_i^a(\boldsymbol{\Phi}, \boldsymbol{\lambda}) \equiv \frac{\partial F}{\partial \lambda_i^a} = -\frac{i}{4} \langle \boldsymbol{\eta}_i^T G^a \boldsymbol{\eta}_i \rangle_{\boldsymbol{\Phi}, \boldsymbol{\lambda}} = 0$$

### Numerical Update Equations and Acceleration Schemes

Solving the coupled system for $\boldsymbol{\Phi}$ and $\boldsymbol{\lambda}$ requires specialized root-finding algorithms because $\boldsymbol{\Phi}$ represents an energy minimum whereas $\boldsymbol{\lambda}$ represents an energy maximum.

In the damped gradient ascent scheme employed by Yilmaz et al., the Lagrange multipliers are updated dynamically between fixed-point iterations of the mean fields:

$$\Phi_m^{(n+1)} = (1 - \alpha_\Phi) \Phi_m^{(n)} + \alpha_\Phi \langle \hat{\mathcal{O}}_m \rangle^{(n)}$$

$$\lambda_i^{a, (n+1)} = \lambda_i^{a, (n)} + \alpha_\lambda g_i^a\left(\boldsymbol{\Phi}^{(n+1)}, \boldsymbol{\lambda}^{(n)}\right)$$

Here, $\alpha_\Phi \in [0.05, 0.2]$ and $\alpha_\lambda \in [0.01, 0.1]$ serve as empirical damping factors. While this method is straightforward to implement, it exhibits first-order asymptotic convergence and slows down significantly when the fermion gap closes.

Second-order convergence is achieved through the multidimensional Newton-Raphson scheme, where the constraint residuals $\mathbf{g}(\boldsymbol{\lambda}) = \mathbf{0}$ are resolved at each step using the exact Jacobian matrix:

$$\mathcal{J}_{(i,a), (j,b)} \equiv \frac{\partial g_i^a}{\partial \lambda_j^b} = \frac{\partial^2 F}{\partial \lambda_i^a \partial \lambda_j^b} = -\chi_{ia, jb}^{loc}(\omega = 0)$$

In terms of the positive BdG eigenvectors $\vert{}u_k\rangle$ and single-particle energies $E_k$, the static Majorana susceptibility is:

$$\mathcal{J}_{(i,a), (j,b)} = -\frac{1}{4} \sum_{k, l} \frac{\tanh\left(\frac{\beta E_k}{2}\right) + \tanh\left(\frac{\beta E_l}{2}\right)}{E_k + E_l} \langle u_k \vert{} G_i^a \vert{} u_l \rangle \langle u_l \vert{} G_j^b \vert{} u_k \rangle$$

The Newton update incorporates a Levenberg-Marquardt regularization parameter $\mu_{LM} \ge 0$ to ensure positive-definiteness away from convergence:

$$\boldsymbol{\lambda}^{(n+1)} = \boldsymbol{\lambda}^{(n)} - \left[ \mathcal{J}\left(\boldsymbol{\lambda}^{(n)}\right) + \mu_{LM} \mathbb{I} \right]^{-1} \mathbf{g}\left(\boldsymbol{\lambda}^{(n)}\right)$$

To accelerate the joint convergence of both mean fields and Lagrange multipliers without evaluating the full non-local susceptibility matrix at every iteration, the Direct Inversion in the Iterative Subspace (DIIS, or Anderson acceleration) is applied to the combined state vector $\mathbf{x} = (\boldsymbol{\Phi}, \boldsymbol{\lambda})^T$. Defining the iteration residual $\mathbf{R}(\mathbf{x}) \equiv \mathbf{x}_{calc} - \mathbf{x}$, the extrapolated update is formed from the history of the previous $M$ iterations as $\mathbf{x}^{(n+1)} = \sum_{m=1}^M c_m \mathbf{x}_{calc}^{(m)}$, where the coefficients $c_m$ minimize the residual norm subject to $\sum_{m=1}^M c_m = 1$. The coefficients are obtained by solving the linear Karush-Kuhn-Tucker system:

$$\begin{pmatrix} B_{11} & \dots & B_{1M} & -1 \\ \vdots & \ddots & \vdots & \vdots \\ B_{M1} & \dots & B_{MM} & -1 \\ -1 & \dots & -1 & 0 \end{pmatrix} \begin{pmatrix} c_1 \\ \vdots \\ c_M \\ \mu \end{pmatrix} = \begin{pmatrix} 0 \\ \vdots \\ 0 \\ -1 \end{pmatrix}, \quad B_{kl} = \langle \mathbf{R}^{(k)}, \mathbf{R}^{(l)} \rangle$$

Finally, when the saddle-point free energy exhibits negative curvature directions that destabilize standard root finders, an augmented Lagrangian penalty is introduced:

$$\mathcal{L}_\mu(\boldsymbol{\Phi}, \boldsymbol{\lambda}) = F(\boldsymbol{\Phi}, \boldsymbol{\lambda}) + \frac{\mu_P}{2} \sum_{i,a} \left[ g_i^a(\boldsymbol{\Phi}, \boldsymbol{\lambda}) \right]^2$$

This modifies the update equations to:

$$\boldsymbol{\Phi}^{(n+1)} = \arg\min_{\boldsymbol{\Phi}} \mathcal{L}_{\mu_P}\left(\boldsymbol{\Phi}, \boldsymbol{\lambda}^{(n)}\right)$$

$$\lambda_i^{a, (n+1)} = \lambda_i^{a, (n)} + \mu_P g_i^a\left(\boldsymbol{\Phi}^{(n+1)}, \boldsymbol{\lambda}^{(n)}\right)$$

where the penalty weight $\mu_P$ is adaptively increased if the residual norm reduction stagnates.

```
Algorithm 1: Stabilized Hybrid Newton-Raphson / Anderson MMFT Solver
Input: Initial parameters Phi^(0), lambda^(0); tolerances tol_Phi, tol_g; history depth M
Output: Self-consistent ground state Phi*, lambda*
Initialize iteration index n = 0, history buffer H = empty
while true do
    Construct real skew-symmetric matrix A(Phi^(n), lambda^(n))
    Diagonalize BdG Hamiltonian H_BdG = i * A to obtain eigenvalues E_k >= 0 and eigenvectors |u_k>
    Compute Fermi-Dirac occupation factors f_k = 1 / (exp(beta * E_k) + 1)
    Evaluate updated mean-field targets Phi_target = EvaluateExpectations(u, f)
    Evaluate constraint violation vector g_i^a = -i/4 * Tr[ G^a * <eta_i * eta_i^T> ]

    if max_norm(Phi_target - Phi^(n)) < tol_Phi and max_norm(g) < tol_g then
        Phi* = Phi^(n)
        lambda* = lambda^(n)
        return Phi*, lambda*
    end if

    if max_norm(g) > tol_g then
        Compute static Majorana susceptibility Jacobian J_(ia, jb) from eigenvectors |u_k>
        Solve linear system [ J + mu_LM * I ] * delta_lambda = -g
        lambda^(n+1) = lambda^(n) + delta_lambda
        Phi^(n+1) = (1 - alpha_Phi) * Phi^(n) + alpha_Phi * Phi_target
    else
        Define input vector x_in = [Phi^(n), lambda^(n)]
        Define output vector x_out = [Phi_target, lambda^(n) + alpha_lambda * g]
        Residual R = x_out - x_in
        Append (x_in, R) to history buffer H (retaining at most M entries)
        Build overlap matrix B_kl = dot(R^(k), R^(l))
        Solve Karush-Kuhn-Tucker system for mixing coefficients c_m
        Set [Phi^(n+1), lambda^(n+1)] = sum_m c_m * x_out^(m)
    end if

    n = n + 1
end while
```

### Convergence Instabilities and Empirical Adjustments

The numerical optimization of Majorana mean-field models across magnetic phase transitions is susceptible to four primary instability mechanisms:

|**Instability Mechanism**|**Physical and Mathematical Origin**|**Literature Benchmark**|**Empirical Resolution Strategy**|
|---|---|---|---|
|**Susceptibility Divergence**|Gap-closing at Dirac points ($E_k + E_l \to 0$), causing the Jacobian $\mathcal{J}$ to diverge.|Yilmaz et al.; Freitas et al.|Introduce Levenberg-Marquardt regularization $\mu_{LM} \approx 10^{-2}\vert{}K\vert{}$ and bound updates $\vert{}\delta \lambda\vert{} \le 0.05\vert{}K\vert{}$.|
|**Limit-Cycle Oscillations**|Bistability across first-order transitions between QSL and polarized states.|Yilmaz et al.|Apply thermal annealing from $T = 0.05\vert{}K\vert{} \to 0$ or increase linear damping $\alpha_\Phi \le 0.05$.|
|**Gauge-Ghost Trapping**|Convergence to unphysical saddle points corresponding to $D_i = -1$.|Ralko & Merino; Yilmaz et al.|Seed multipliers $\lambda_i^\gamma$ with a staggered configuration to break continuous unphysical symmetries.|
|**Subspace Collapse**|Eigenvector crossing causing discontinuous jumps in Anderson history vectors.|Yilmaz et al.|Flush DIIS history buffer upon residual norm spikes and revert to steepest descent line searches.|

## Competing Magnetic Orders and Multi-Sublattice Unit Cells

### Honeycomb Geometry, Bond Vectors, and Zigzag Ordering Domains

The pristine honeycomb lattice is described by a two-site Bravais basis ($A$ and $B$) with primitive translation vectors:

$$\mathbf{a}_1 = \sqrt{3} d_{NN} \, \hat{\mathbf{x}}, \quad \mathbf{a}_2 = d_{NN} \left( \frac{\sqrt{3}}{2} \hat{\mathbf{x}} + \frac{3}{2} \hat{\mathbf{y}} \right)$$

where $d_{NN}$ is the nearest-neighbor bond distance. Taking an $A$-sublattice site at $\mathbf{r}_A = (0, 0)$, the three adjacent $B$-sublattice sites are reached by directed bond vectors $\boldsymbol{\delta}_\gamma$:

$$\boldsymbol{\delta}_x = d_{NN} \left( -\frac{\sqrt{3}}{2}, -\frac{1}{2} \right), \quad \boldsymbol{\delta}_y = d_{NN} \left( \frac{\sqrt{3}}{2}, -\frac{1}{2} \right), \quad \boldsymbol{\delta}_z = d_{NN} (0, 1)$$

The bond vectors satisfy the lattice relations $\boldsymbol{\delta}_x - \boldsymbol{\delta}_z = -\mathbf{a}_2$ and $\boldsymbol{\delta}_y - \boldsymbol{\delta}_z = \mathbf{a}_1 - \mathbf{a}_2$. The reciprocal lattice vectors associated with this primitive frame are:

$$\mathbf{b}_1 = \frac{2\pi}{\sqrt{3} d_{NN}} \left( 1, -\frac{1}{\sqrt{3}} \right), \quad \mathbf{b}_2 = \frac{4\pi}{3 d_{NN}} (0, 1)$$

Collinear zigzag order spontaneously breaks the threefold rotational symmetry ($C_3$) of the honeycomb lattice, selecting one of three symmetry-equivalent domains characterized by an ordering wavevector $\mathbf{Q}_M$ at the Brillouin zone boundary:

$$\mathbf{Q}_{M_x} = \left( 0, \frac{2\pi}{3 d_{NN}} \right), \quad \mathbf{Q}_{M_y} = \left( \frac{\pi}{\sqrt{3} d_{NN}}, -\frac{\pi}{3 d_{NN}} \right), \quad \mathbf{Q}_{M_z} = \left( \frac{\pi}{\sqrt{3} d_{NN}}, \frac{\pi}{3 d_{NN}} \right)$$

In the $\text{ZZ}_z$ domain, spins are aligned ferromagnetically across $x$- and $y$-bonds and antiferromagnetically across $z$-bonds. In the complementary stripy phase ($\text{ST}_z$), spins are antiparallel across $x$- and $y$-bonds and parallel across $z$-bonds.

### The 4-Sublattice Rectangular Unit Cell

To capture uniform ferromagnetic (FM), Néel antiferromagnetic (AFM), and the $\text{ZZ}_z$ and $\text{ST}_z$ states without enlarging the unit cell beyond what is strictly necessary, a 4-sublattice orthorhombic unit cell is established. The supercell translation vectors are defined by:

$$\mathbf{A}_1 = \mathbf{a}_1 = \sqrt{3} d_{NN} \, \hat{\mathbf{x}}, \quad \mathbf{A}_2 = 2\mathbf{a}_2 - \mathbf{a}_1 = 3 d_{NN} \, \hat{\mathbf{y}}$$

The four atomic basis positions within the rectangular cell are given by:

$$\mathbf{r}_1 \, (A_1) = (0, 0)$$

$$\mathbf{r}_2 \, (B_1) = \boldsymbol{\delta}_z = (0, d_{NN})$$

$$\mathbf{r}_3 \, (A_2) = \boldsymbol{\delta}_z - \boldsymbol{\delta}_x = d_{NN} \left( \frac{\sqrt{3}}{2}, \frac{3}{2} \right)$$

$$\mathbf{r}_4 \, (B_2) = 2\boldsymbol{\delta}_z - \boldsymbol{\delta}_x = d_{NN} \left( \frac{\sqrt{3}}{2}, \frac{5}{2} \right)$$

Nearest-neighbor connectivity across sublattices proceeds as follows: $z$-bonds connect $(1 \leftrightarrow 2)$ and $(3 \leftrightarrow 4)$ entirely within the home cell; $x$-bonds connect $(3 \leftrightarrow 2)$ intra-cell and $(1 \leftrightarrow 4)$ across the translation $-\mathbf{A}_2$; and $y$-bonds connect $(3 \leftrightarrow 2)$ across the translation $+\mathbf{A}_1$ and $(1 \leftrightarrow 4)$ across the translation $\mathbf{A}_1 - \mathbf{A}_2$.

|**Magnetic Phase**|**Sublattice 1 (A1A_1A1​)**|**Sublattice 2 (B1B_1B1​)**|**Sublattice 3 (A2A_2A2​)**|**Sublattice 4 (B2B_2B2​)**|**Ordering Vector (Q\mathbf{Q}Q)**|**Literature References**|
|---|---|---|---|---|---|---|
|**Ferromagnetic (FM)**|$+\hat{\mathbf{n}}$|$+\hat{\mathbf{n}}$|$+\hat{\mathbf{n}}$|$+\hat{\mathbf{n}}$|$\mathbf{0}$|Rau et al.; Winter et al.|
|**Néel Antiferromagnet (AFM)**|$+\hat{\mathbf{n}}$|$-\hat{\mathbf{n}}$|$+\hat{\mathbf{n}}$|$-\hat{\mathbf{n}}$|$\mathbf{0}$ _(in 4-site cell)_|Rau et al.; Winter et al.|
|**Zigzag-$\mathbf{z}$ ($\text{ZZ}_z$)**|$+\hat{\mathbf{n}}$|$-\hat{\mathbf{n}}$|$-\hat{\mathbf{n}}$|$+\hat{\mathbf{n}}$|$\mathbf{Q}_{M_z} \equiv \frac{1}{2}\mathbf{A}_2^*$|Knolle et al.; Winter et al.|
|**Stripy-$\mathbf{z}$ ($\text{ST}_z$)**|$+\hat{\mathbf{n}}$|$+\hat{\mathbf{n}}$|$-\hat{\mathbf{n}}$|$-\hat{\mathbf{n}}$|$\mathbf{Q}_{M_z} \equiv \frac{1}{2}\mathbf{A}_2^*$|Rau et al.; Winter et al.|
|**Zigzag-$\mathbf{x}$ ($\text{ZZ}_x$)**|$+\hat{\mathbf{n}}$|$-\hat{\mathbf{n}}$|$+\hat{\mathbf{n}}$|$-\hat{\mathbf{n}}$ _(phase shifted)_|$\mathbf{Q}_{M_x}$ _(doubled along $\mathbf{A}_1$)_|Winter et al.|

In the presence of finite symmetric off-diagonal exchange $\Gamma > 0$, the classical spin direction $\hat{\mathbf{n}}$ cants away from the cubic axes toward the crystallographic axes, determined by the dominant eigenvector of the local exchange matrix:

$$\hat{\mathbf{n}} = \frac{1}{\sqrt{2 + r^2}} (1, 1, -r)^T, \quad r = \frac{-(K - \Gamma) + \sqrt{(K - \Gamma)^2 + 8\Gamma^2}}{2\Gamma}$$

For $K = \Gamma > 0$, this yields $r = 2$ and $\hat{\mathbf{n}} = \frac{1}{\sqrt{6}}(1, 1, -2)^T$, which aligns the magnetic moments directly along the in-plane zigzag direction.

### The 8-Sublattice Supercell for Non-Coplanar and Multi-$\mathbf{Q}$ Orders

While the 4-sublattice rectangular cell accommodates single-$\mathbf{Q}$ collinear orders, it is insufficient when external magnetic fields break the equivalence of the three zigzag domains or when competing exchanges stabilize non-coplanar structures. For magnetic fields tilted away from high-symmetry axes, the system favors multi-$\mathbf{Q}$ configurations, including the commensurate triple-$\mathbf{Q}$ "AF star" and "FM star" orders identified by Winter et al.. Describing these states, as well as preserving compatibility with the 4-sublattice Klein duality transformation $\mathcal{T}_4$, requires an 8-sublattice cell formed by doubling both primitive translations:

$$\mathbf{R}_{cell} = m_1 (2\mathbf{a}_1) + m_2 (2\mathbf{a}_2)$$

The cell comprises four $A$-type sites located at $\mathbf{r}_{A, \mu} \in \{\mathbf{0}, \mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_1 + \mathbf{a}_2\}$ and four $B$-type sites at $\mathbf{r}_{B, \mu} = \mathbf{r}_{A, \mu} + \boldsymbol{\delta}_z$. This doubled geometry nests all three zone-boundary ordering vectors $\mathbf{Q}_{M_x}, \mathbf{Q}_{M_y}, \mathbf{Q}_{M_z}$ simultaneously into the supercell $\Gamma$ point. Consequently, it can capture non-coplanar spin superpositions of the form $\mathbf{S}_i = \sum_{\gamma} \mathbf{S}_{M_\gamma} \cos(\mathbf{Q}_{M_\gamma} \cdot \mathbf{r}_i)$, as well as long-wavelength incommensurate modulations proximate to the ferromagnetic-to-zigzag phase boundaries.

### Momentum-Space BdG Bloch Matrices for the 4-Sublattice Cell

In the 4-sublattice cell, four Majorana species per site yield a 16-component Nambu-Majorana field:

$$\boldsymbol{\Psi}_{\mathbf{k}} = \left( \boldsymbol{\eta}_{1, \mathbf{k}}^T, \boldsymbol{\eta}_{2, \mathbf{k}}^T, \boldsymbol{\eta}_{3, \mathbf{k}}^T, \boldsymbol{\eta}_{4, \mathbf{k}}^T \right)^T$$

where $\boldsymbol{\eta}_{\mu, \mathbf{k}} = (c_{\mu, \mathbf{k}}, b_{\mu, \mathbf{k}}^x, b_{\mu, \mathbf{k}}^y, b_{\mu, \mathbf{k}}^z)^T$. The BdG Hamiltonian is represented in momentum space as:

$$H_{MF} = \frac{1}{2} \sum_{\mathbf{k} \in \text{BZ}} \boldsymbol{\Psi}_{\mathbf{k}}^\dagger \mathcal{H}_4(\mathbf{k}) \boldsymbol{\Psi}_{\mathbf{k}}$$

The $16 \times 16$ Hermitian Bloch matrix $\mathcal{H}_4(\mathbf{k})$ consists of sixteen $4 \times 4$ blocks:

$$\mathcal{H}_4(\mathbf{k}) = \begin{pmatrix} \mathcal{D}_1 & \mathcal{T}_{12}(\mathbf{k}) & 0 & \mathcal{T}_{14}(\mathbf{k}) \\ \mathcal{T}_{12}^\dagger(\mathbf{k}) & \mathcal{D}_2 & \mathcal{T}_{23}(\mathbf{k}) & 0 \\ 0 & \mathcal{T}_{23}^\dagger(\mathbf{k}) & \mathcal{D}_3 & \mathcal{T}_{34}(\mathbf{k}) \\ \mathcal{T}_{14}^\dagger(\mathbf{k}) & 0 & \mathcal{T}_{34}^\dagger(\mathbf{k}) & \mathcal{D}_4 \end{pmatrix}$$

The on-site blocks $\mathcal{D}_\mu$ are purely imaginary, skew-Hermitian matrices encoding local Zeeman interactions and the Lagrange multipliers enforcing the gauge constraints:

$$\mathcal{D}_\mu = \begin{pmatrix} 0 & -i h_x & -i h_y & -i h_z \\ i h_x & 0 & -i \lambda_\mu^z & i \lambda_\mu^y \\ i h_y & i \lambda_\mu^z & 0 & -i \lambda_\mu^x \\ i h_z & -i \lambda_\mu^y & i \lambda_\mu^x & 0 \end{pmatrix}$$

The inter-sublattice hopping blocks incorporate matter, gauge, and cross-hopping correlations:

$$\mathcal{T}_{12}(\mathbf{k}) = \begin{pmatrix} i \chi_{12}^z & 0 & 0 & i \kappa_{12}^{c b^z} \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ i \kappa_{21}^{b^z c} & 0 & 0 & i u_{12}^z \end{pmatrix}, \quad \mathcal{T}_{34}(\mathbf{k}) = \begin{pmatrix} i \chi_{34}^z & 0 & 0 & i \kappa_{34}^{c b^z} \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ i \kappa_{43}^{b^z c} & 0 & 0 & i u_{34}^z \end{pmatrix}$$

Inter-cell connectivity carries phase factors corresponding to the rectangular lattice vectors $\mathbf{A}_1$ and $\mathbf{A}_2$:

$$\mathcal{T}_{14}(\mathbf{k}) = \mathcal{T}_{14}^x \, e^{-i \mathbf{k} \cdot \mathbf{A}_2} + \mathcal{T}_{14}^y \, e^{i \mathbf{k} \cdot (\mathbf{A}_1 - \mathbf{A}_2)}$$

$$\mathcal{T}_{23}(\mathbf{k}) = \mathcal{T}_{23}^x + \mathcal{T}_{23}^y \, e^{-i \mathbf{k} \cdot \mathbf{A}_1}$$

where $\mathcal{T}_{\mu\nu}^\gamma$ contains the bond expectations along link $\gamma$. Diagonalization of $\mathcal{H}_4(\mathbf{k})$ produces 16 particle-hole symmetric bands $\pm E_n(\mathbf{k})$ ($n=1, \dots, 8$).

### Explicit Hopping Blocks and 8-Sublattice Bloch Matrix Formulation

In the presence of symmetric off-diagonal exchange $\Gamma$ and Heisenberg interaction $J$, each nearest-neighbor bond $\langle \mu\nu \rangle_\gamma$ contributes a generalized $4 \times 4$ hopping block $\mathcal{T}_{\mu\nu}^\gamma$:

$$\mathcal{T}_{\mu\nu}^\gamma = \begin{pmatrix} i \chi_{\mu\nu}^\gamma & i \kappa_{\mu\nu}^{c b^x} & i \kappa_{\mu\nu}^{c b^y} & i \kappa_{\mu\nu}^{c b^z} \\ i \kappa_{\nu\mu}^{b^x c} & i u_{\mu\nu}^{xx} & i u_{\mu\nu}^{xy} & i u_{\mu\nu}^{xz} \\ i \kappa_{\nu\mu}^{b^y c} & i u_{\mu\nu}^{yx} & i u_{\mu\nu}^{yy} & i u_{\mu\nu}^{yz} \\ i \kappa_{\nu\mu}^{b^z c} & i u_{\mu\nu}^{zx} & i u_{\mu\nu}^{zy} & i u_{\mu\nu}^{zz} \end{pmatrix}$$

where the bond parameters are defined as $\chi_{\mu\nu}^\gamma = \langle i c_\mu c_\nu \rangle$, $u_{\mu\nu}^{\alpha\beta} = \langle i b_\mu^\alpha b_\nu^\beta \rangle$, and $\kappa_{\mu\nu}^{c b^\beta} = \langle i c_\mu b_\nu^\beta \rangle$. Along an $x$-bond, the primary Kitaev channel engages the gauge component $u_{\mu\nu}^{xx}$, while the off-diagonal exchange $\Gamma$ populates the cross-gauge terms $u_{\mu\nu}^{yz}$ and $u_{\mu\nu}^{zy}$. Along a $y$-bond, the Kitaev channel populates $u_{\mu\nu}^{yy}$, with $\Gamma$ populating $u_{\mu\nu}^{xz}$ and $u_{\mu\nu}^{zx}$.

For the 8-sublattice supercell defined by doubled primitive vectors $\mathbf{A}_1' = 2\mathbf{a}_1$ and $\mathbf{A}_2' = 2\mathbf{a}_2$, the Nambu-Majorana multiplet expands to 32 components:

$$\boldsymbol{\Psi}_{\mathbf{k}}^{(8)} = \left( \boldsymbol{\eta}_{A_1, \mathbf{k}}^T, \boldsymbol{\eta}_{B_1, \mathbf{k}}^T, \dots, \boldsymbol{\eta}_{A_4, \mathbf{k}}^T, \boldsymbol{\eta}_{B_4, \mathbf{k}}^T \right)^T$$

The resulting $32 \times 32$ Bloch Hamiltonian matrix $\mathcal{H}_8(\mathbf{k})$ assumes a block bipartite structure:

$$\mathcal{H}_8(\mathbf{k}) = \begin{pmatrix} \mathbb{D}_A & \mathbb{T}_{AB}(\mathbf{k}) \\ \mathbb{T}_{AB}^\dagger(\mathbf{k}) & \mathbb{D}_B \end{pmatrix}$$

where $\mathbb{D}_A = \text{diag}(\mathcal{D}_{A_1}, \mathcal{D}_{A_2}, \mathcal{D}_{A_3}, \mathcal{D}_{A_4})$ and $\mathbb{D}_B = \text{diag}(\mathcal{D}_{B_1}, \mathcal{D}_{B_2}, \mathcal{D}_{B_3}, \mathcal{D}_{B_4})$ are block-diagonal on-site skew-Hermitian matrices containing local Zeeman terms and Lagrange multipliers. The $16 \times 16$ inter-sublattice block $\mathbb{T}_{AB}(\mathbf{k})$ contains the four $A$-to-$B$ sublattice interconnections across the three bond orientations:

$$\mathbb{T}_{AB}(\mathbf{k}) = \begin{pmatrix} \mathcal{T}^z & \mathcal{T}^x e^{-i \mathbf{k} \cdot \mathbf{a}_2} & \mathcal{T}^y e^{i \mathbf{k} \cdot (\mathbf{a}_1 - \mathbf{a}_2)} & 0 \\ \mathcal{T}^x & \mathcal{T}^z & 0 & \mathcal{T}^y e^{i \mathbf{k} \cdot (\mathbf{a}_1 - \mathbf{a}_2)} \\ \mathcal{T}^y & 0 & \mathcal{T}^z & \mathcal{T}^x e^{-i \mathbf{k} \cdot \mathbf{a}_2} \\ 0 & \mathcal{T}^y & \mathcal{T}^x & \mathcal{T}^z \end{pmatrix}$$

## Decoupling Channel Variations and Wick Contractions

### Algebraic Properties of the Three Majorana Representations

Fu, Knolle, and Perkins provided a unified analysis of the three operator mappings used to represent quantum spin-$1/2$ operators via Majorana fermions: the $SO(3)$ Tsvelik representation ($\mathbf{J}^\alpha$), the Kitaev representation ($\mathbf{N}^\alpha$), and the $SO(4)$ chiral representation ($\mathbf{M}^\alpha$).

|**Mathematical Property**|**SO(3)SO(3)SO(3) Representation (Jα\mathbf{J}^\alphaJα)**|**Kitaev Representation (Nα\mathbf{N}^\alphaNα)**|**SO(4)SO(4)SO(4) Chiral Representation (Mα\mathbf{M}^\alphaMα)**|
|---|---|---|---|
|**Operator Construction**|$J_i^\alpha = -\frac{i}{2}\epsilon^{\alpha\beta\gamma}\eta_i^\beta \eta_i^\gamma$|$N_i^\alpha = i b_i^\alpha c_i$|$M_i^\alpha = \frac{i}{2}(\eta_i^0 \eta_i^\alpha - \frac{1}{2}\epsilon^{\alpha\beta\gamma}\eta_i^\beta \eta_i^\gamma)$|
|**Majorana Multiplicity**|3 per site: $(\eta^x, \eta^y, \eta^z)$|4 per site: $(c, b^x, b^y, b^z)$|4 per site: $(\eta^0, \eta^x, \eta^y, \eta^z)$|
|**Local Hilbert Space**|Non-local pairing ($2^{3N/2}$)|Doubled ($2^{2N}$); local unphysical states|Doubled ($2^{2N}$); projects via left-chiral sector|
|**Projector Commutation**|No local projection required|$[N_i^\alpha, D_i] \ne 0$; mixes parity|$[M_i^\alpha, P_{i,L}] = 0$ commutes identically|
|**Zeeman Coupling Term**|Bilinear: $-\frac{i}{2}\sum_\alpha h_\alpha \epsilon^{\alpha\beta\gamma}\eta_i^\beta \eta_i^\gamma$|Bilinear: $-i \sum_\alpha h_\alpha b_i^\alpha c_i$|Mixed bilinear: $-\sum_\alpha h_\alpha M_i^\alpha$|
|**Integrability of Pure Kitaev**|Solvable via non-local Jordan-Wigner|Solvable via static gauge links $\hat{u}_{ij}^\gamma$|Solvable via chiral projection|

The $SO(3)$ representation has the advantage of introducing no unphysical local states, as three Majorana fermions per site generate an overall Hilbert space dimension of $2^{3N/2}$, exactly matching $2^{N/2}$ physical copies of the spin system. However, this non-local redundancy makes it difficult to decouple non-Kitaev terms without breaking translation invariance. The Kitaev representation simplifies the bond interactions into static $\mathbb{Z}_2$ gauge fields, but the spin operators $N_i^\alpha$ anticommute with the parity operator $D_i$, meaning that unconstrained mean-field decoupling mixes physical and unphysical sectors. The $SO(4)$ chiral representation resolves this issue by projecting the spin operators into the left-handed sector via $P_{i,L} = \frac{1+D_i}{2}$, which ensures that $[M_i^\alpha, P_{i,L}] = 0$ holds identically at the operator level.

### Wick Contraction Channels and Mean-Field Energetics

Decoupling the exchange term $H_{ij}^\gamma = -K_\gamma S_i^\gamma S_j^\gamma$ produces distinct physical channels depending on which pairings are retained. In the Kitaev representation, the four-Majorana term $H_{ij}^\gamma = \frac{K_\gamma}{4}(b_i^\gamma c_i)(b_j^\gamma c_j)$ yields three primary Wick pairings:

1. The **Kitaev pairing channel**, $\langle i b_i^\gamma b_j^\gamma \rangle (i c_i c_j) + \langle i c_i c_j \rangle (i b_i^\gamma b_j^\gamma)$, which preserves the separation between itinerant matter fermions and static $\mathbb{Z}_2$ gauge fields.

2. The **Fock cross-hopping channel**, $-\langle i c_i b_j^\gamma \rangle (i b_i^\gamma c_j) - \langle i b_i^\gamma c_j \rangle (i c_i b_j^\gamma)$, which couples matter and gauge sectors and generates hybridization bands.

3. The **Hartree magnetic channel**, $\langle i b_i^\gamma c_i \rangle (i b_j^\gamma c_j) + \langle i b_j^\gamma c_j \rangle (i b_i^\gamma c_i) \equiv -4 m_i^\gamma m_j^\gamma$, which drives conventional long-range magnetic ordering.


In the $SO(4)$ chiral representation, the composite nature of $M_i^\alpha$ expands each bond interaction into 16 distinct Wick contractions: four pure matter contractions $\langle \eta_i^0 \eta_j^0 \rangle (\eta_i^\gamma \eta_j^\gamma)$, four pure gauge contractions $\epsilon^{\gamma\alpha\beta}\epsilon^{\gamma\mu\nu} \langle \eta_i^\alpha \eta_j^\mu \rangle (\eta_i^\beta \eta_j^\nu)$, and eight mixed-chirality cross contractions $\epsilon^{\gamma\alpha\beta}\langle \eta_i^0 \eta_j^\alpha \rangle (\eta_i^\gamma \eta_j^\beta)$.

|**Decoupling Scheme**|**Ground-State Energy (E0/∣K∣E_0 / \vert{}K\vert{}E0​/∣K∣)**|**Critical Field (hc/∣K∣h_c / \vert{}K\vert{}hc​/∣K∣)**|**Physical Mechanism and Limitations**|
|---|---|---|---|
|**Exact Solution (Lieb's Theorem)**|$-0.39362$|$0.07 - 0.08$ _(DMRG benchmark)_|Flux-free ground state; topological gap opens as $\Delta \sim h_x h_y h_z / K^2$.|
|**Frozen-Link Channel ($u_{ij}^\gamma = 1$)**|$-0.39362$|$0.35 - 0.40$|Gauge fluxes are treated as static, overestimating QSL stability against high fields.|
|**Unconstrained Full Wick ($\boldsymbol{\lambda} = \mathbf{0}$)**|$-0.520$ to $-0.540$|$\approx 0.18$|Overbinds the ground state by mixing in unphysical negative-parity ($D_i = -1$) states.|
|**Hartree-Only Magnetic Channel**|$-0.250$ to $-0.300$|$0.12 - 0.15$|Underestimates quantum entanglement, driving an artificial sharp first-order transition.|
|**Constrained MMFT ($SO(4)$ Multipliers)**|$-0.394$ to $-0.410$|$0.06 - 0.08$|Enforcing $\langle \hat{\Gamma}_i^\gamma \rangle = 0$ restores the gap-closing transition at the $\Gamma$ point.|

The quantitative comparisons demonstrate that retaining hybrid cross-contractions while enforcing the gauge constraints via Lagrange multipliers is necessary to benchmark the stability of Kitaev spin liquids against competing magnetic orders. When the gauge fields are treated as frozen ($u_{ij}^\gamma = 1$), the absence of dynamic matter-gauge coupling artificially stabilizes the spin liquid, overestimating the critical breakdown field by a factor of nearly five.

Conversely, omitting the Lagrange multipliers allows unphysical states into the variational space, lowering the ground-state energy well below the exact Lieb bound. Full saddle-point optimization over $\boldsymbol{\lambda}_i$ restores the variational energy to within a few percent of the exact value and correctly reproduces the low-field gap-closing transition observed in numerical tensor network simulations.

## Benchmark Literature Survey and Reference Registry

- Lu Freitas, Tim Bauer, Reinhold Egger, and Rodrigo G. Pereira, _Electric polarization near vortices in the extended Kitaev model_, npj Quantum Materials 9, 33 (2024).

    - DOI: 10.1038/s41535-024-00643-5

    - arXiv: 2311.15750 [cond-mat.str-el]

- Lu Freitas, _Majorana mean-field theory with Z2 vortices: electrical tests in the extended Kitaev model_, Ph.D. thesis, Universidade Federal do Rio Grande do Norte (UFRN) / Heinrich-Heine-Universität Düsseldorf (2023).

- F. Yilmaz, A. P. Kampf, and S.-K. Yip, _Phase diagrams of the Kitaev-Heisenberg-$\Gamma$-$\Gamma'$ model in a magnetic field_, Phys. Rev. Research 4, 043024 (2022).

    - DOI: 10.1103/PhysRevResearch.4.043024

    - arXiv: 2202.06158 [cond-mat.str-el]

- F. Yilmaz, A. P. Kampf, and S.-K. Yip, _Interlayer exchange in bilayer Kitaev spin liquids_, Phys. Rev. B 109, 024439 (2024).

    - DOI: 10.1103/PhysRevB.109.024439

- A. Ralko and J. Merino, _Majorana Chiral Spin Liquid in a Model for Mott Insulating Cuprates_, Phys. Rev. Lett. 124, 217203 (2020).

    - DOI: 10.1103/PhysRevLett.124.217203

- Urban F. P. Seifert, Tobias Meng, and Matthias Vojta, _Fractionalized Fermi liquids and exotic superconductivity in the Kitaev-Kondo lattice_, Phys. Rev. B 97, 085118 (2018).

    - DOI: 10.1103/PhysRevB.97.085118

- Urban F. P. Seifert, J. Doneit, and Matthias Vojta, _Bilayer Kitaev models and fractionalized Fermi liquids_, Phys. Rev. B 98, 155101 (2018).

    - DOI: 10.1103/PhysRevB.98.155101

- Jianlong Fu, Johannes Knolle, and Natalia B. Perkins, _Majorana representations of spin and an alternative solution of the Kitaev honeycomb model_, Phys. Rev. B 97, 115142 (2018).

    - DOI: 10.1103/PhysRevB.97.115142

    - arXiv: 1801.09713 [cond-mat.str-el]

- Jeffrey G. Rau, Eric Kin-Ho Lee, and Hae-Young Kee, _Generic Spin Model for Bilayer and Honeycomb Materials with Strong Spin-Orbit Coupling_, Phys. Rev. Lett. 112, 077204 (2014).

    - DOI: 10.1103/PhysRevLett.112.077204

- Stephen M. Winter, Ying Li, Harald O. Jeschke, and Roser Valentí, _Challenges in design of Kitaev materials: Magnetic interactions in $\alpha$-RuCl$_3$_, Phys. Rev. B 93, 214431 (2016); and _Breakdown of Magnons in a Strongly Spin-Orbit Coupled Magnet_, Phys. Rev. Lett. 117, 277202 (2016).

    - DOI: 10.1103/PhysRevB.93.214431 and 10.1103/PhysRevLett.117.277202

    - arXiv: 1611.03436 [cond-mat.str-el]

- J. Knolle, G.-W. Chern, D. L. Kovrizhin, R. Moessner, and N. B. Perkins, _Magnon damping in the zigzag phase of the Kitaev-Heisenberg-$\Gamma$ model on a honeycomb lattice_, Phys. Rev. B 97, 134432 (2018).

    - DOI: 10.1103/PhysRevB.97.134432

    - arXiv: 1911.12829 [cond-mat.str-el]
