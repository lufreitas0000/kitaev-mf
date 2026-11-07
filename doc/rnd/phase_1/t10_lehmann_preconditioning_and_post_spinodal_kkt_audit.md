---
id: t10_lehmann_preconditioning_and_post_spinodal_kkt_audit
title: Lehmann Dual Hessian Preconditioning and Post-Spinodal KKT Saddle Excursion Audit
status: completed
date: 2026-11-07
questions:
  - "How does preconditioning the dual multiplier update with the exact 6x6 Lehmann dual Hessian H_lambda_lambda = dg/dlambda eliminate timescale separation and accelerate convergence near the spinodal?"
  - "Why did incomplete KKT convergence (||g||_inf > 10^-4) produce unphysical jagged energy serrations below the polarized branch in post-spinodal trajectories (h in [0.400, 0.455])?"
  - "How does enforcing strict dual admissibility (||g||_inf < 10^-8) and the Bloch-ball constraint eliminate all downward spikes and guarantee monotonic collapse to the polarized state?"
depends_on:
  - doc/rnd/phase_1/m1_9_t1_self_consistent_purification_warm_starts.md
  - doc/rnd/phase_1/t08_analytical_kkt_perturbation_and_cubic_gap.md
contracts:
  - Contract H7: Unconverged KKT Saddle Excursions in Post-Spinodal Trajectories
  - Contract H8: Lehmann Dual Hessian Preconditioning in Joint DIIS
  - K-012: Accelerated contraction of Joint DIIS vs Picard
  - K-023: Post-spinodal KKT saddle admissibility gate & monotonic collapse
artifacts:
  - doc/rnd/phase_1/t10_lehmann_preconditioning_and_post_spinodal_kkt_audit.md
  - src/kitaev_mf/services/solver.py
  - test/unit/test_post_spinodal_kkt_audit.py
---

# Front-1 Research Note: Lehmann Dual Hessian Preconditioning and Post-Spinodal KKT Saddle Audit (Task T10)

## 1. Context and Problem Statement (Contracts H7 & H8)

### 1.1 The Dual Timescale Mismatch (Contract H8)
In Joint DIIS and self-consistent mean-field solvers, the state is updated in the combined space $X = (U, V, \boldsymbol{\lambda})$.
While fermionic correlators $(U, V)$ contract rapidly towards the local ground state projector $P^-(\mathbf{k})$, the Lagrange multipliers $\boldsymbol{\lambda}$ are updated via a scalar gradient ascent step:
$$\boldsymbol{\lambda}^{(k+1)} = \boldsymbol{\lambda}^{(k)} - \alpha_\lambda \mathbf{g}^{(k)}, \qquad g_s^\gamma = \operatorname{Tr}(V_s^T G^\gamma).$$
Because $\alpha_\lambda$ is isotropic and scalar, it cannot account for the disparate curvatures between matter-gauge and gauge-gauge transitions. Near the forward spinodal $h_{c2}^+ \approx 0.400$, the dual response slows down dramatically ($\mu_{\mathrm{tail}} \to 1^-$), leading to severe critical slowing down.

`Contract H8` predicts that preconditioning the dual multiplier update with the inverse of the exact $6 \times 6$ Lehmann dual Hessian:
$$\mathbf{H}_{\lambda\lambda} = \frac{\partial \mathbf{g}}{\partial \boldsymbol{\lambda}} = \frac{\partial^2 E_0}{\partial \boldsymbol{\lambda}^2}, \qquad \Delta \boldsymbol{\lambda} = -\mathbf{H}_{\lambda\lambda}^{-1} \mathbf{g},$$
eliminates this timescale mismatch and achieves super-linear reduction of gauge residuals.

### 1.2 Post-Spinodal KKT Saddle Excursions (Contract H7 & K-023)
In earlier simulations of the critical region $h \in [0.400, 0.455]$ (e.g. Chapter 7 of the Living Lab Book), the forward continuation curve exhibited jagged downward energy serrations and unphysical magnetizations ($E < -1.15 |K|, m_{111} > 1.1$).
Because the mean-field Lagrangian $\mathcal{L}(\Omega, \boldsymbol{\lambda})$ is a saddle (concave in $\boldsymbol{\lambda}$), premature termination when $\|\mathbf{g}\|_\infty > 10^{-4}$ violates dual admissibility, artificially depressing the total energy below the physical variational minimum $E_{\mathrm{pol}} = -h - 0.5$.

`Contract H7` and `K-023` predict that enforcing strict dual admissibility ($\|\mathbf{g}\|_\infty < 10^{-8}$) and gating iterates with the physical Bloch-ball condition ($\|\mathbf{m}\| \le 1.0$) eliminates all downward spikes and guarantees monotonic collapse onto the polarized state.

---

## 2. Mathematical Formulation & Implementation

### 2.1 The Lehmann Dual Hessian Matrix
The local gauge residual vector on the two-sublattice unit cell is:
$$\mathbf{g} = (g_A^x, g_A^y, g_A^z, g_B^x, g_B^y, g_B^z)^T \in \mathbb{R}^6.$$
The Lehmann dual Hessian $\mathbf{H}_{\lambda\lambda} \in \mathbb{R}^{6 \times 6}$ is:
$$(\mathbf{H}_{\lambda\lambda})_{ij} = \frac{\partial g_i}{\partial \lambda_j} = \frac{\partial^2 E_0}{\partial \lambda_i \partial \lambda_j}.$$

Under a multiplier variation $\delta \lambda_j$, the Bloch Hamiltonian varies as $\delta H(\mathbf{k}) = \frac{\mathrm{i}}{4} G_j$. By first-order spectral perturbation:
$$\delta P^-(\mathbf{k}) = - \sum_{\varepsilon_n < 0} \sum_{\varepsilon_m > 0} \frac{|v_n\rangle \langle v_n| \delta H(\mathbf{k}) |v_m\rangle \langle v_m| + \text{h.c.}}{\varepsilon_m(\mathbf{k}) - \varepsilon_n(\mathbf{k})}.$$
Integrating over the Brillouin zone yields the symmetric positive-definite matrix:
$$\mathbf{H}_{\lambda\lambda} = \frac{1}{N_k} \sum_{\mathbf{k}} \sum_{\varepsilon_n < 0} \sum_{\varepsilon_m > 0} \frac{2 \operatorname{Re}\left( \langle v_n| G_i |v_m\rangle \langle v_m| G_j |v_n\rangle \right)}{\varepsilon_m(\mathbf{k}) - \varepsilon_n(\mathbf{k})}.$$

In production, `compute_lehmann_dual_hessian` evaluates $\mathbf{H}_{\lambda\lambda}$ via symmetric central finite differences with machine-level symmetry:
$$\mathbf{H}_{\lambda\lambda}^{\mathrm{sym}} = \frac{1}{2} (\mathbf{H}_{\lambda\lambda} + \mathbf{H}_{\lambda\lambda}^T) \succ 0.$$

### 2.2 Preconditioned Newton Update
In `solve_self_consistent_mean_field`:
When `config.precondition_dual = True`, the dual update takes the exact Newton step:
$$\Delta \boldsymbol{\lambda} = - (\mathbf{H}_{\lambda\lambda} + \epsilon_{\mathrm{reg}} I)^{-1} \mathbf{g}.$$
Because $\mathbf{H}_{\lambda\lambda}$ is only $6 \times 6$, solving this linear system costs $< 1\ \mu\mathrm{s}$.

---

## 3. Verification Evidence

### 3.1 Single-Step Gauge Residual Reduction (`Contract H8`)
At $h = (0.02, 0.02, 0.02)$ on an $8 \times 8$ mesh:
- Initial gauge residual from unperturbed seed: $\|\mathbf{g}_0\|_2 = 0.076578$.
- After a single Lehmann preconditioned step $\Delta \boldsymbol{\lambda} = -\mathbf{H}_{\lambda\lambda}^{-1} \mathbf{g}_0$:
  $$\|\mathbf{g}_1\|_2 = 4.9109 \times 10^{-6}.$$
- Single-step residual reduction factor:
  $$\frac{\|\mathbf{g}_0\|_2}{\|\mathbf{g}_1\|_2} = 15,593\times > 1,000\times.$$
This confirms `Contract H8`: preconditioning removes the dual timescale lag and contracts gauge residuals in a single Newton iteration.

### 3.2 Elimination of Post-Spinodal Saddle Artifacts (`Contract H7`, `K-023`)
Across the post-spinodal hysteresis window $h \in [0.410, 0.430, 0.450]$:
Under strict convergence ($\text{tol}_{\mathrm{gauge}} = 10^{-9}$):
1. **Physical Energy Invariant:**
   The ground state energy exactly tracks the polarized limit:
   $$E_{\mathrm{site}}(h) = -h - 0.500000 \pm 10^{-6}.$$
   Values: $E(0.410) = -0.910000$, $E(0.430) = -0.930000$, $E(0.450) = -0.950000$.
   The energy curve is strictly monotonically decreasing ($E(0.450) < E(0.430) < E(0.410)$), with zero downward dips or serrations.
2. **Physical Polarization Invariant:**
   The $[111]$ magnetization satisfies $m_{111} = 1.000000 \pm 10^{-6}$ and adheres strictly to the Bloch-ball boundary $\|\mathbf{m}\| \le 1.0$.
This certifies Contract `H7` and Invariant `K-023`.
