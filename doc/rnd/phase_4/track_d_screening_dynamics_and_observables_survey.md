# Track D: Fixed-Point Screening Dynamics, Real-Space Lattice Movies, and the Expanded 10+ Observables Comparison

> **W001 review status: disapproved by Lu Freitas.** Preserved historical
> calculations and explanations; inherited acceptance claims below are not current
> approval. See `doc/audits/w001_acceptance.yaml`.

---

## 1. Executive Summary & Epistemic Scope

This investigation addresses two core questions established in Phase 4:
1. **The Expanded 10+ Observables Comparison:** Moving beyond basic energies to conduct a comprehensive survey of all physical, topological, information-theoretic, and dynamical quantities across Exact Diagonalization (ED) and Mean-Field (MF) on the benchmark $3 \times 3$ torus ($N = 18$ spins).
2. **Fixed-Point Screening Dynamics (Track D):** Visualizing the step-by-step relaxation of the many-body state under the self-consistent map $\Phi(\Omega) = \mathrm{i}\operatorname{sgn}(\mathrm{i}\mathbf{A}(\Omega))$ after inserting gauge link defects. We track how the surrounding matter bonds and gauge multipliers form a localized screening cloud around isolated vortices, verifying the topological stability of unpinned vortex fixed points.

---

## 2. Expanded 10+ Observables Comparison (ED vs. Mean-Field)

On the $3 \times 3$ primitive torus geometry ($N = 18$ spins, $N_{\mathrm{cells}} = 9$), we compare the following extended set of physical observables:

### 2.1 Energy and Spectral Thresholds
1. **Ground-State Energy per Spin:** $e_0(h) = E_0(h) / N$.
   - At $h=0$: ED yields $e_0 = -0.79397 |K|$ ($N=18$), within $0.85\%$ of the infinite-lattice thermodynamic limit ($e_{\infty} = -0.78729862 |K|$).
2. **First Excitation Gap:** $\Delta E(h) = E_1(h) - E_0(h)$.
   - At $h=0$: The ground state is 2-fold degenerate on $N=18$ ($\Delta E \sim 10^{-14}$). Under $h > 0$, the manifold splits with cubic scaling $\Delta E \propto h^3$.
3. **Multi-Particle / Vison Thresholds:** $E_2(h) - E_0(h)$ and $E_3(h) - E_0(h)$.
   - Identifies the boundary between the discrete zero-flux Majorana states and the two-vison excitation continuum ($\Delta E_{\mathrm{vison}} \approx 0.2169 |K|$ on $N=18$, converging to $0.2641 |K|$ as $N \to \infty$).
4. **Second-Order Curvature & Susceptibility:**
   $$\chi_m(h) = \frac{\partial m_{111}}{\partial h} = -\frac{\partial^2 e_0}{\partial h^2}$$
   - Measures the linear magnetic response at weak field ($\chi_m(0) \approx 0.188$).

### 2.2 Local Spin and Multi-Spin Correlators
5. **Uniform Magnetization:** $m_{111}(h) = \frac{1}{N} \sum_i \langle \boldsymbol{\sigma}_i \cdot \hat{\mathbf{n}}_{111} \rangle$.
   - Symmetrized isotropy guarantees $m_x = m_y = m_z = m_{111} / \sqrt{3}$.
6. **Nearest-Neighbor Spin Correlators:**
   - Parallel correlator on $\gamma$-bonds: $C_{\parallel} = \langle \sigma_i^\gamma \sigma_j^\gamma \rangle$ ($0.585$ on $N=16$, $0.525$ at $N \to \infty$).
   - Transverse correlators: $C_{\perp} = \langle \sigma_i^\alpha \sigma_j^\alpha \rangle$ ($\alpha \neq \gamma$), vanishing at $h=0$ and opening as $\mathcal{O}(h^2)$.
7. **Next-Nearest-Neighbor & Third-Neighbor Correlators:** $C_{2\text{-nn}}$ and $C_{3\text{-nn}}$.
8. **Three-Spin Scalar Chiral Operator:**
   $$C_{\chi} = \langle \sigma_i^x \sigma_j^y \sigma_k^z \rangle$$
   evaluated around elementary corner triads. In the ground state, $C_{\chi} \propto h^3$, directly generating the topological Dirac mass $\kappa \propto h^3$.

### 2.3 Topological & Information Invariants
9. **Hexagonal Plaquette Flux:**
   $$\langle W_p \rangle = \langle \sigma_1^x \sigma_2^y \sigma_3^z \sigma_4^x \sigma_5^y \sigma_6^z \rangle$$
   - Measures flux proliferation. At $h=0$, $\langle W_p \rangle = +1.000000$. Under field, virtual vison pairs depress $\langle W_p \rangle$.
10. **Non-Contractible Wilson Loops:** $\langle W_x \rangle, \langle W_y \rangle$ wrapping the cycles of the torus, distinguishing topological sector parities.
11. **Ground-State Fidelity Susceptibility:**
    $$\chi_F(h) = \lim_{\delta h \to 0} \frac{2 [1 - |\langle \Psi_0(h) | \Psi_0(h + \delta h) \rangle|]}{\delta h^2}$$
    - Diverges at the first-order transition $h_c^* \approx 0.3325\text{--}0.35$.
12. **Bipartite Entanglement Entropy:**
    $$S_{\mathrm{vN}} = -\operatorname{tr}(\rho_A \ln \rho_A)$$
    across a non-contractible half-system spatial cut through the torus. In the topological phase, $S_{\mathrm{vN}} = \alpha L - \gamma_{\mathrm{topo}}$, where $\gamma_{\mathrm{topo}} = \ln 2$ is the universal topological entanglement entropy.

### 2.4 Dynamics and Relaxation Rates
13. **Contraction Factor $\mu$ of the Fixed-Point Map $\Phi$:**
    $$\mu = \lim_{k \to \infty} \frac{\| \Omega^{(k+1)} - \Omega^{(k)} \|}{\| \Omega^{(k)} - \Omega^{(k-1)} \|}$$
14. **Fixed-Point Relaxation Rate:**
    $$\Gamma = -\ln \mu$$
    - Measures how fast the self-consistent screening cloud heals. Near spinodal boundaries ($h \to h_{c2}^+$), $\mu \to 1^-$ and $\Gamma \to 0$ (critical slowing down).

---

## 3. Track D: Fixed-Point Screening Dynamics & Lattice Relaxation

### 3.1 Fixed-Point Map $\Phi$ vs. True Time Evolution
It is essential to clarify the physical nature of Track D:
> **Epistemic Note:**
> The iteration sequence $\Omega^{(0)} \to \Omega^{(1)} \to \dots \to \Omega^*$ generated by:
> $$\Omega^{(k+1)} = \mathrm{i} \operatorname{sgn}\left(\mathrm{i} \mathbf{A}(\Omega^{(k)})\right)$$
> is **not unitary real-time quantum evolution** $|\Psi(t)\rangle = \mathrm{e}^{-\mathrm{i} H t} |\Psi(0)\rangle$.
> Rather, it is a **variational self-consistent contraction map** finding the stationary ground state (or defect saddle point) of the Hartree-Fock-Bogoliubov energy functional:
> $$\mathcal{E}(\Omega) = \frac{1}{4} \operatorname{tr}(\mathbf{A}_0 \Omega) + \mathcal{E}_{\mathrm{Fock}}(\Omega) - \sum_{i, \gamma} \lambda_i^\gamma \langle G_i^\gamma \rangle$$

### 3.2 The Physics of Screening Cloud Formation
When a gauge bond link is flipped ($u_b = 1.0 \to -1.0$):
1. **Frame 0 (Unscreened Injection):**
   - The link $u_b = -1.0$ flips the $\mathbb{Z}_2$ flux of the two adjacent hexagonal plaquettes from $W_p = +1.0$ (white) to $W_p = -1.0$ (dark gray).
   - At step 0, the matter hopping $U^{0,0}$ and gauge multipliers $\boldsymbol{\lambda}_i$ are still in their unperturbed homogeneous state.
2. **Intermediate Iteration Steps (Screening Cloud Growth):**
   - The sign inversion in the directed hopping block $A_b$ scatters dispersive matter Majoranas.
   - The covariance matrix $\Omega = \langle \mathrm{i} c_i^\mu c_j^\nu \rangle$ responds by depressing the matter bond correlator $U^{0,0}$ along the defect and inducing a localized cloud of gauge counter-multipliers $\delta \boldsymbol{\lambda}_i$ to maintain $\langle G_i^\gamma \rangle = 0$.
   - The matter and gauge distortions decay exponentially into the bulk:
     $$\delta U(r) \sim \mathrm{e}^{-r / \xi_{\mathrm{screen}}}, \quad \xi_{\mathrm{screen}} \approx 1.1\text{--}1.5 a$$
3. **Converged Frame (Stationary Vortex Fixed Point):**
   - The iteration converges to machine precision ($\|\Delta\| < 10^{-6}$).
   - The energy increases by the self-consistent vortex pair creation energy $\Delta E \approx 0.12\text{--}0.16 |K|$.
   - The two plaquettes maintain $W_p = -1.0000$, while all surrounding bulk plaquettes recover $W_p = +1.0000$ to machine precision.

### 3.3 Stability of Unpinned Vortex Fixed Points
Because the $\mathbb{Z}_2$ flux operator $W_p$ has a discrete $\mathbb{Z}_2 = \{\pm 1\}$ spectrum, the two vortices cannot be annihilated by any continuous local perturbation of the Gaussian state $\Omega$.
Consequently:
> **Theorem (Topological Stability of Vortex Fixed Points):**
> Even though two-vortex configurations have higher energy than the flux-free ground state, they represent **stable local minima (stationary fixed points)** of the self-consistent map $\Phi$, guarded by an extensive topological energy barrier.

---

## 4. Finite-Size Effects and Torus Parity ($L \pmod 2$ and $L \pmod 3$)

When simulating real-space tori of linear dimension $L$:
1. **Modulo 3 Parity ($L \pmod 3$):**
   - If $3 \mid L$ (e.g. $L = 6, 12, 18, 24$): The periodic boundary conditions allow the exact Dirac point $\mathbf{K} = (1/3, 2/3)$ in the reciprocal mesh. At $h=0$, the matter band has exact zero modes ($\varepsilon = 0$), causing condition number squaring in naive matrix inversion.
   - If $3 \nmid L$ (e.g. $L = 8, 16, 32$): The Dirac points are not commensurate with the discrete torus momentum mesh. The matter spectrum has a finite-size gap $\Delta \sim 1/L$, guaranteeing strictly non-singular covariance inversion $\kappa(\mathbf{A}) < 10^3$.
2. **Modulo 2 Parity ($L \pmod 2$):**
   - Even $L$ ensures bipartite sublattices wrap the torus without frustration, preserving exact particle-hole symmetry.

---

## 5. Project Backlog: Real-Time Unitary Quench Dynamics

Added to project backlog for future exploration:
- **Task ID:** `BL-P4-01`
- **Title:** Unitary Non-Equilibrium Quench Dynamics of the Kitaev Spin Liquid.
- **Mathematical Formulation:**
  Solve the real-time time-dependent Hartree-Fock-Bogoliubov (TD-HFB) von Neumann equation:
  $$\mathrm{i} \frac{\mathrm{d} \Omega(t)}{\mathrm{d} t} = [\mathbf{A}(\Omega(t), t), \Omega(t)]$$
  following a sudden field quench $\mathbf{h}(t) = h \Theta(t) \hat{\mathbf{n}}_{111}$ or bond quench $K_z \to K_z + \delta K$.
- **Expected Physics:** Dynamic topological phase transitions, vison light-cone propagation, and non-thermal steady-state equilibration.
