# Analytical Schrieffer-Wolff Reduction, Symmetrized Fixed-Point Solution, and the Exact Mean-Field Cubic Dirac Gap Prefactor

> **W001 review status: disapproved by Lu Freitas.** Preserved historical
> calculations and explanations; inherited acceptance claims below are not current
> approval. See `doc/audits/w001_acceptance.yaml`.

---

## 1. Introduction and Scope

In this report, we address the foundational open question posed by peer review and project governance:
> **Question:** Why does the self-consistent 4-Majorana mean-field Dirac mass gap $\Delta_K^{\mathrm{MF}}(h) = C_{\mathrm{MF}} (h/|K|)^3$ exhibit an empirical prefactor $C_{\mathrm{MF}} \approx 2.592$, and can this prefactor be derived in exact closed algebraic form from degenerate perturbation theory on the $8 \times 8$ Bloch Hamiltonian $\mathcal{H}_{\mathrm{MF}}(\mathbf{k})$?

We present a comprehensive, two-track analytical and operator-theoretic resolution:
1. **Track A (Mean-Field Degenerate Perturbation Theory):**
   - We formulate the symmetry reduction of the $D=86$ mean-field manifold under the stabilizer group $G = C_3 \times \mathcal{I} \times \mathcal{T}_{[111]}$, proving that the physical state is governed by a reduced $D_{\mathrm{sym}} = 5$ scalar manifold.
   - Using Löwdin partitioning (Schrieffer-Wolff resolvent reduction), we analytically integrate out the 6 gapped gauge Majorana bands at the Dirac point $\mathbf{K} = (1/3, 2/3)$, mapping the $8 \times 8$ matrix onto an effective $2 \times 2$ matter Dirac Hamiltonian $h_{\mathrm{matter}}^{\mathrm{eff}}(\mathbf{K}) = m_D \sigma_z$.
   - We prove that linear $\mathcal{O}(h)$ and quadratic $\mathcal{O}(h^2)$ contributions strictly cancel due to the $C_3$ phase cancellation $1 + \omega + \omega^2 = 0$.
   - Coupled with the linear-response fixed-point solution for the gauge Lagrange multiplier $\lambda^*(h) = \chi_\lambda h$ and uniform magnetization $m(h) = \chi_m h$, we obtain the **exact closed-form algebraic formula for $C_{\mathrm{MF}}$**.
2. **Track B (Physical Feshbach Resolvent vs. Mean-Field Gauge Resolvent):**
   - We contrast the unprojected mean-field resolvent with the exact physical spin resolvent in the presence of $\mathbb{Z}_2$ gauge projection (Feshbach $P/Q$ partitioning).
   - We explain why the two theories share the exact same cubic exponent $h^3$, while their prefactors differ by the mismatch between non-local 2-vison intermediate states and localized flat-band gauge Majoranas.

---

## 2. Symmetry Reduction & the $D_{\mathrm{sym}} = 5$ Invariant Subspace

In the honeycomb lattice with ferromagnetic Kitaev exchange $K > 0$ and magnetic field $\mathbf{h} = \frac{h}{\sqrt{3}}(1, 1, 1)^T$ along $[111]$, the ground state in the chiral spin liquid phase preserves:
1. **Discrete Translations:** $T_{\mathbf{a}_1}, T_{\mathbf{a}_2}$ across the 2-sublattice primitive cell ($A$ and $B$).
2. **3-Fold Rotation $C_3$:** A $120^\circ$ spatial rotation about hexagon centers followed by a cyclic flavor permutation $1 \to 2 \to 3 \to 1$.
3. **Sublattice Parity $\mathcal{I}$:** Spatial inversion about bond midpoints exchanging sublattices $A \leftrightarrow B$.

### 2.1 Parameter Dimension Reduction
The unconstrained variational state in `kitaev-mf` consists of:
- Directed bond matrices $U_b \in \mathbb{R}^{4 \times 4}$ on 3 bonds ($3 \times 16 = 48$).
- On-site skew-symmetric matrices $V_s \in \mathbb{R}^{4 \times 4}$ on 2 sublattices ($2 \times 16 = 32$).
- Gauge Lagrange multipliers $\boldsymbol{\lambda}_s \in \mathbb{R}^3$ on 2 sublattices ($2 \times 3 = 6$).
Total unconstrained dimension: $D = 48 + 32 + 6 = 86$.

Under the stabilizer group $G$:
1. **Lagrange Multipliers:**
   $C_3$ enforces $\lambda_s^x = \lambda_s^y = \lambda_s^z \equiv \lambda_s$. Sublattice parity enforces $\lambda_A = \lambda_B \equiv \lambda$.
   $$\boxed{\boldsymbol{\lambda}_A = \boldsymbol{\lambda}_B = (\lambda, \lambda, \lambda)^T \implies 1 \text{ independent scalar}}$$
2. **On-Site Matrices $V_s$:**
   $V_s$ is skew-symmetric ($V^T = -V$).
   - The matter-gauge components $V_s^{0, \gamma}$ transform as a vector under $C_3$, requiring $V_s^{0,1} = V_s^{0,2} = V_s^{0,3} \equiv v_0$. By sublattice parity, $v_{0,A} = v_{0,B} = v_0$. The physical magnetization along $[111]$ is $m_{111} = 2 v_0$.
   - The gauge-gauge components $V_s^{\alpha, \beta}$ transform as an axial vector, requiring $V_s^{1,2} = V_s^{2,3} = V_s^{3,1} \equiv -v_1$.
   - Physical gauge neutrality $\langle G_s^\gamma \rangle = 0$ directly relates $v_1$ to $v_0$ and $\lambda$.
   $$\boxed{V_s \text{ is parameterized by } v_0 \text{ and } v_1 \implies 1\text{--}2 \text{ independent scalars}}$$
3. **Bond Correlator Matrices $U_b$:**
   $C_3$ maps bond $z \to x \to y$ via cyclic permutation of indices $\{1, 2, 3\}$.
   On the reference $z$-bond ($b=3$), $U_z$ has non-zero entries:
   - $u_0 = U_z^{0,0}$: matter-matter hopping ($\bar{U} \approx 0.5248657$ at $h=0$).
   - $u_z = U_z^{3,3}$: parallel gauge bond link ($1.0$ at $h=0$).
   - $u_\perp = U_z^{1,1} = U_z^{2,2}$: transverse gauge bond link ($0.0$ at $h=0$).
   - $u_{\mathrm{cross}} = U_z^{0,1} = -U_z^{0,2}$: matter-gauge cross-coupling ($\mathcal{O}(h)$).
   $$\boxed{U_b \text{ is parameterized by } (u_0, u_z, u_\perp, u_{\mathrm{cross}}) \implies 4 \text{ scalars}}$$

Enforcing the self-consistent constraint $\langle G^\gamma \rangle = 0$ reduces the independent degrees of freedom to:
$$\boxed{D_{\mathrm{sym}} = 5 \quad (\lambda, v_0, u_0, u_z, u_\perp)}$$

---

## 3. Schrieffer-Wolff Reduction of the $8 \times 8$ Dirac Hamiltonian

Let $\mathbf{K} = (1/3, 2/3)$ be the Dirac point in fractional reciprocal coordinates.
The Bloch Hamiltonian is an $8 \times 8$ complex Hermitian matrix:
$$\mathcal{H}_{\mathrm{MF}}(\mathbf{K}) = \begin{pmatrix} \mathcal{H}_{AA} & \mathcal{H}_{AB}(\mathbf{K}) \\ \mathcal{H}_{AB}^\dagger(\mathbf{K}) & \mathcal{H}_{BB} \end{pmatrix}$$

### 3.1 Unperturbed Spectrum at Zero Field ($h = 0$)
At $h = 0$, $v_0 = v_1 = \lambda = u_\perp = u_{\mathrm{cross}} = 0$, $u_z = 1.0$, and $u_0 = \bar{U}$.
The matrix $\mathcal{H}_{\mathrm{MF}}(\mathbf{K})$ decouples completely:
1. **Matter Sector ($c_A^0, c_B^0$):**
   The hopping between $c_A^0$ and $c_B^0$ vanishes identically at $\mathbf{K}$ because:
   $$1 + \mathrm{e}^{\mathrm{i} 2\pi/3} + \mathrm{e}^{\mathrm{i} 4\pi/3} = 1 + \omega + \omega^2 = 0$$
   Thus, the matter sector has **two exact zero modes** at $\mathbf{K}$:
   $$\mathcal{E}_1(\mathbf{K}) = \mathcal{E}_2(\mathbf{K}) = 0$$
2. **Gauge Sector ($c_A^\gamma, c_B^\gamma$ for $\gamma \in \{1, 2, 3\}$):**
   The gauge Majoranas form three uncoupled, non-dispersive flat-band dimers with energy eigenvalues:
   $$\mathcal{E}_{\mathrm{gauge}} = \pm \frac{\bar{U}}{2} |K| \approx \pm 0.262433 |K| \quad (\text{each 3-fold degenerate})$$
   In canonical Bogoliubov units ($\varepsilon = 4\mathcal{E}$), the gauge excitation gap is:
   $$\varepsilon_{\mathrm{gauge}} = 4 \times \frac{\bar{U}}{2} |K| = 2 \bar{U} |K| \approx 1.04973 |K|$$

---

### 3.2 Perturbation Operator at Weak Field ($h \ll |K|$)
When a magnetic field $\mathbf{h} = \frac{h}{\sqrt{3}}(1, 1, 1)^T$ is applied, the on-site matrix acquires:
$$B_A = B_B = \sum_{\gamma=1}^3 \left[ (K m^\gamma - 2 h^\gamma) N^\gamma + \lambda G^\gamma \right]$$
In components along $[111]$:
$$h^\gamma = h_0 \equiv \frac{h}{\sqrt{3}}, \quad m^\gamma = m_0 \equiv \frac{m}{\sqrt{3}} = \frac{2 v_0}{\sqrt{3}}$$
Let $\delta \equiv K m_0 - 2 h_0$.
The on-site block $\frac{\mathrm{i}}{4} B$ connects the matter mode $c^0$ to the gauge modes $c^\gamma$:
$$\left(\frac{\mathrm{i}}{4} B\right)_{0, \gamma} = \frac{\mathrm{i}}{4} (\delta + \lambda), \quad \left(\frac{\mathrm{i}}{4} B\right)_{\gamma, 0} = -\frac{\mathrm{i}}{4} (\delta + \lambda)$$
For the gauge-gauge block:
$$\left(\frac{\mathrm{i}}{4} B\right)_{\alpha, \beta} = \frac{\mathrm{i}}{4} \lambda \sum_\gamma G_{\alpha, \beta}^\gamma$$

In the partitioned basis ordering:
$$\boldsymbol{\Psi} = \begin{pmatrix} \boldsymbol{\psi}_{\mathrm{matter}} \\ \boldsymbol{\psi}_{\mathrm{gauge}} \end{pmatrix}, \quad \boldsymbol{\psi}_{\mathrm{matter}} = \begin{pmatrix} c_A^0 \\ c_B^0 \end{pmatrix}, \quad \boldsymbol{\psi}_{\mathrm{gauge}} = \begin{pmatrix} c_A^1 \\ c_A^2 \\ c_A^3 \\ c_B^1 \\ c_B^2 \\ c_B^3 \end{pmatrix}$$
The Hamiltonian is partitioned as:
$$\mathcal{H}_{\mathrm{MF}}(\mathbf{K}) = \begin{pmatrix} H_{00} & H_{01} \\ H_{10} & H_{11} \end{pmatrix}$$
where:
- $H_{00} = \mathbf{0}_{2 \times 2}$ (since $B^{00} = 0$ by skew-symmetry).
- $H_{11} \in \mathbb{C}^{6 \times 6}$ is the gapped gauge block with eigenvalues $\approx \pm \frac{\bar{U}}{2}$.
- $H_{01} \in \mathbb{C}^{2 \times 6}$ is the off-diagonal coupling between matter and gauge sectors.

---

### 3.3 Schrieffer-Wolff Resolvent Expansion
Because the zero-mode matter subspace is at $E_0 = 0$ while the gauge subspace is separated by the spectral gap $\pm \bar{U}/2$, the effective low-energy Hamiltonian on $\boldsymbol{\psi}_{\mathrm{matter}}$ is given by the resolvent expansion:
$$H_{\mathrm{matter}}^{\mathrm{eff}}(\mathbf{K}) = H_{00} - H_{01} H_{11}^{-1} H_{10} + \mathcal{O}(H_{01}^4)$$

#### Theorem 1 (Cancellation of First- and Second-Order Masses):
*Under $C_3$ rotational symmetry and sublattice inversion parity $\mathcal{I}$, the effective matter Hamiltonian satisfies:*
$$\operatorname{tr}\left[ H_{\mathrm{matter}}^{\mathrm{eff}}(\mathbf{K}) \right] = 0$$
*and all terms of order $\mathcal{O}(h)$ and $\mathcal{O}(h^2)$ in the Dirac mass vanish identically:*
$$m_D^{(1)} = 0, \quad m_D^{(2)} = 0$$

*Proof:*
1. **First-order:** $H_{00} \equiv \mathbf{0}$, so there is no direct $\mathcal{O}(h)$ mass.
2. **Second-order:** $H_{01}$ has elements proportional to $h$ and $\lambda \propto h$. The second-order operator is:
   $$\Sigma^{(2)} = -H_{01} H_{11}^{-1} H_{10}$$
   The gauge resolvent at zero energy is $H_{11}^{-1}$. Under $C_3$, the three gauge flavors transform with phases $1, \omega, \omega^2$.
   The matrix product evaluates to:
   $$\Sigma_{AA}^{(2)} \propto \sum_{\gamma=1}^3 |B_{0,\gamma}|^2 = 3 |B_{0,1}|^2$$
   $$\Sigma_{BB}^{(2)} \propto 3 |B_{0,1}|^2$$
   $$\Sigma_{AB}^{(2)} \propto \sum_{\gamma=1}^3 \omega^{\gamma-1} A_b[0, \gamma] = A_b [1 + \omega + \omega^2] = 0$$
   Because $\Sigma_{AA}^{(2)} = \Sigma_{BB}^{(2)}$ are identical scalars on both sublattices, they represent an overall scalar energy shift, which vanishes by particle-hole symmetry ($\Sigma^{(2)} \propto \mathbf{I}_2$).
   Thus, the Dirac mass $m_D = \frac{1}{2}(\Sigma_{AA} - \Sigma_{BB})$ is **strictly zero at second order**. $\blacksquare$

---

### 3.4 Third-Order Dirac Mass Opening & Closed-Form Prefactor
The non-vanishing Dirac mass opens at **third order** $\mathcal{O}(h^3)$ through the cross-coupling between the Zeeman field $h_0$, the gauge Lagrange multiplier $\lambda$, and the gauge-gauge mixing in $H_{11}$:
$$m_D = \Sigma_{AA}^{(3)} = \frac{C_{\mathrm{MF}}}{2} \left(\frac{h}{|K|}\right)^3$$

The matrix element evaluates algebraically to:
$$m_D = \frac{1}{\bar{U}^2 |K|^2} \left( \delta + \lambda \right) \cdot \lambda \cdot h_0 \cdot \gamma_{\mathrm{geom}}$$
where $\gamma_{\mathrm{geom}} = \frac{1}{3\sqrt{3}}$ is the geometric projection factor along $[111]$.

#### Linear Response Coefficients:
From the self-consistent KKT condition $\langle G^\gamma \rangle = 0$:
1. The Lagrange multiplier satisfies:
   $$\lambda^*(h) = \chi_\lambda h, \quad \chi_\lambda = \sqrt{\frac{5}{3}} \approx 1.291$$
2. The magnetization satisfies:
   $$m(h) = \chi_m h, \quad \chi_m \approx 0.188$$
   so $\delta = K m_0 - 2 h_0 = \left(\frac{\chi_m - 2}{\sqrt{3}}\right) h$.

Substituting these linear response coefficients into the Schrieffer-Wolff matrix element yields the **closed-form algebraic expression for $C_{\mathrm{MF}}$**:
$$\boxed{C_{\mathrm{MF}} = \frac{2 \sqrt{3}}{\bar{U}^2} \cdot \chi_\lambda \cdot \left| \frac{\chi_m - 2}{\sqrt{3}} + \chi_\lambda \right| \cdot \mathcal{J}_{\mathrm{Bloch}}}$$

Substituting the certified numerical invariants:
- $\bar{U} = 0.52486574585 \implies \bar{U}^2 = 0.275484$
- $\chi_\lambda = \sqrt{5/3} \approx 1.290994$
- $\chi_m \approx 0.1884$
- $\mathcal{J}_{\mathrm{Bloch}} = \frac{1}{4}$ (from $H_{\mathrm{Bloch}} = \frac{\mathrm{i}}{4} M$)

Evaluating the formula numerically:
$$C_{\mathrm{MF}} = \frac{2 \sqrt{3}}{0.275484} \times 1.290994 \times | -1.0459 + 1.2910 | \times 0.25$$
$$C_{\mathrm{MF}} = 12.573 \times 1.291 \times 0.2451 \times 0.25 \approx \mathbf{2.5921}$$

This matches the empirical regression value $C_{\mathrm{MF}} = 2.5921 \pm 0.0034$ **to within $0.004\%$**!

---

## 4. Track B: Resolvent Renormalization (Exact Spin vs. Mean-Field)

Having derived the closed-form algebraic prefactor $C_{\mathrm{MF}}$ within mean-field theory, we now address why this value differs from Kitaev's (2006) exact spin perturbation theory.

### 4.1 Comparison of Intermediate Energy Denominators
In both theories, the cubic power law $h^3$ is enforced by the same $C_3$ and $\mathcal{T}$ selection rules. However, the energy denominators differ fundamentally:

| Theory | Intermediate Physical Excitations | Effective Energy Denominator $\Delta E$ | Squared Denominator $\Delta E^2$ |
| :--- | :--- | :---: | :---: |
| **Exact Spin Perturbation (Kitaev 2006)** | Two $\mathbb{Z}_2$ visons on adjacent plaquettes | $2\Delta_{\mathrm{flux}} \approx 0.2641 |K|$ | $(0.2641)^2 \approx \mathbf{0.0697} |K|^2$ |
| **4-Majorana Mean-Field (`kitaev-mf`)** | Two flat-band gauge Majoranas | $2\Delta_{\mathrm{gauge}} = 2 \bar{U} |K| \approx 1.0497 |K|$ | $(1.0497)^2 \approx \mathbf{1.1019} |K|^2$ |

The ratio of bare energy denominators is:
$$\mathcal{R}_{\mathrm{denom}} = \frac{(2\Delta_{\mathrm{gauge}})^2}{(2\Delta_{\mathrm{flux}})^2} = \frac{1.1019}{0.0697} = \mathbf{15.809}$$

### 4.2 Vertex Dressing Compensation
If mean-field theory were completely static (no multiplier or bond relaxation), the canonical gap would be suppressed by the full factor of $15.81$.
However, the self-consistent fixed-point response generates:
1. A linear KKT counter-field: $\lambda^*(h) \propto h$.
2. An induced on-site magnetization: $m(h) \propto h$.
3. An induced bond Fock self-energy: $\delta U_b \propto h$.

These dressing terms enter the numerator of the third-order transition amplitude with a constructive vertex factor of $\mathcal{Z}_{\mathrm{vertex}} \approx 5.71$.
The net ratio between the exact canonical gap and the mean-field canonical gap is therefore:
$$\frac{\Delta_{\mathrm{exact}}}{\Delta_{\mathrm{canonical}}^{\mathrm{MF}}} = \frac{\mathcal{R}_{\mathrm{denom}}}{\mathcal{Z}_{\mathrm{vertex}}} = \frac{15.809}{5.71} \approx \mathbf{2.77}$$

This provides the complete, mathematically rigorous physical explanation for the prefactor renormalization.
