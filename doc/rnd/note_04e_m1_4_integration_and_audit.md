# Note 04e: Milestone M1.4 Integration, Numerical Anomalies, and Lessons Learned

**Status:** M1.4.4 R&D audit record
**Scope:** Document trial-and-error resolutions, algebraic sign traps, and physical domain boundaries discovered during the implementation of production observables and energetics services.

---

## 1. Incident 1: Element-wise Multiplication vs. Matrix Trace on Skew Matrices

### The Anomaly
In `src/kitaev_mf/services/hamiltonian.py`, the local Hartree exchange field was initially assembled using:
```python
tr_v_n[site, beta - 1] = torch.sum(v_torch[site].T * n_torch[beta])
```
While this passed all unit tests in the pure Kitaev limit ($V_s = \mathbf{0}$), testing under finite Heisenberg exchange $J = 0.04$ and field $\mathbf{h} = (0.02, 0.01, 0.03)$ revealed an energy discrepancy of $\Delta E \approx 5.08 \times 10^{-4}$ between the single-particle spectral energy $E_{\mathrm{quad}, 0} - C$ and the variational Wick functional $E_{\mathrm{total}}^{\mathrm{G}}$.

### Root Cause
For two matrices $A$ and $B$, the element-wise Hadamard sum evaluates to:
$$\sum_{\mu,\nu} A_{\mu\nu} B_{\mu\nu} = \operatorname{tr}(A^T B).$$
Setting $A = V_s^T = -V_s$ and $B = N^\beta$ gives:
$$\sum_{\mu,\nu} (V_s^T)_{\mu\nu} N_{\mu\nu}^\beta = \operatorname{tr}((V_s^T)^T N^\beta) = \operatorname{tr}(V_s N^\beta) = -\operatorname{tr}(V_s^T N^\beta).$$
Because both $V_s$ and $N^\beta$ are skew-symmetric ($(V_s)^T = -V_s$ and $(N^\beta)^T = -N^\beta$), taking the transpose of $V_s$ inside the element-wise sum inverted the overall sign of the Hartree field:
$$\text{torch.sum}(V_s^T * N^\beta) = -\operatorname{tr}(V_s^T N^\beta).$$
This inverted the on-site Hartree shift in $B_i$, producing a $2 \times E_{\mathrm{Hartree}}$ error between the tangent operator and the energy functional.

### The Resolution
Replaced the element-wise operation with the explicit matrix product trace:
```python
tr_v_n[site, beta - 1] = torch.trace(v_torch[site].T @ n_torch[beta])
```
This restored exact consistency $\vert{}(E_{\mathrm{quad}, 0} - C) - E_{\mathrm{total}}^{\mathrm{G}}\vert{} < 10^{-10}$ across all non-zero exchange couplings and magnetic fields.

---

## 2. Incident 2: Invariant vs. Dynamic $\mathbb{Z}_2$ Plaquette Fluxes (Wilson Loops)

### The Anomaly
An integration test asserted:
```python
assert np.allclose(flux_2, 1.0, atol=1e-10)
assert np.allclose(flux_4, 1.0, atol=1e-10)
```
for a state converged under finite perturbation $\Gamma = 0.04$ and $\mathbf{h} = (0.02, 0.02, 0.02)$, causing an immediate assertion failure ($W_p \approx 0.94$).

### Root Cause (Physical Misconception in Test Specification)
- In the pure Kitaev limit ($J = \Gamma = \Gamma' = D = \mathbf{h} = 0$), the local link operators $\hat{u}_{ij}^\gamma = i b_i^\gamma b_j^\gamma$ commute with the Hamiltonian ($[\hat{u}_{ij}^\gamma, H] = 0$). By Lieb's theorem, $u_b = +1.0$ is an exact conserved quantum number, so $W_p = \prod_{b \in \partial p} u_b \equiv +1.0$ to machine precision.
- Under finite $\Gamma$ and $\mathbf{h}$, the gauge Majoranas acquire dispersion and hybridize with matter Majoranas ($[\hat{u}_{ij}^\gamma, H] \neq 0$). The bond links fluctuate ($u_b < 1.0$), driving $W_p < 1.0$.
- Asserting $W_p = 1.0$ under finite perturbations is physically incorrect.

### The Resolution
1. Preserved exact $W_p = 1.0$ assertion strictly for the pure Kitaev state ($K = -1.0$).
2. For perturbed states, scale parity asserts that the **fluctuating flux and energy density of the 2-site cell match those of the 4-site supercell**:
   $$W_{p, \text{2-site}} = W_{p, \text{4-site}}, \qquad E_{\text{site}, \text{2-site}} = E_{\text{site}, \text{4-site}}.$$

---

## 3. Incident 3: Premature Solver Convergence and On-Site Magnetization Drift

### The Anomaly
The solver originally monitored only:
$$\text{diff\_u} = \max_b \vert{}U_b^{(n+1)} - U_b^{(n)}\vert{}, \qquad \text{max\_g} = \max_{s,\gamma} \vert{}g_s^\gamma\vert{}.$$
In magnetic phases (Zigzag, Stripy, field-polarized), the gauge residual $\text{max\_g}$ reached the tolerance threshold ($< 10^{-8}$) while on-site magnetization components $V_s^{0\gamma} = \langle \sigma_s^\gamma \rangle$ were still relaxing across iterations.

### Root Cause
$V_s \in \mathfrak{so}(4)$ contains 6 degrees of freedom: 3 gauge residuals $g_s^\gamma = \operatorname{tr}(G^\gamma V_s)$ and 3 physical dipole moments $\langle \sigma_s^\gamma \rangle = V_s^{0\gamma}$. Monitoring only $g_s$ left the magnetic relaxation unconstrained.

### The Resolution
Added `diff_v = max_s |V_s^{(n+1)} - V_s^{(n)}|` to `SolverMetrics` and required:
```python
if diff_u < tol and diff_v < tol and max_g < tol:
    converged = True
```
This guarantees that both matter hoppings, physical magnetic moments, and gauge constraints are simultaneously stationary at self-consistency.

---

## References

* **Freitas, L. R. D., Bauer, T., Egger, R., and Pereira, R. G.** (2024). *Electric polarization near vortices in the extended Kitaev model*. npj Quantum Materials **9**, 33. DOI: [10.1038/s41535-024-00643-5](https://doi.org/10.1038/s41535-024-00643-5).
* **Kitaev, A.** (2006). *Anyons in an exactly solved model and beyond*. Annals of Physics **321**, 2–111. DOI: [10.1016/j.aop.2005.10.005](https://doi.org/10.1016/j.aop.2005.10.005).
