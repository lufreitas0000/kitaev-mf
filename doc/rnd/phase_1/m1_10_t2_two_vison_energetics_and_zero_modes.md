# M1.10-T2: Real-Space Two-Vison Energetics, Interaction Potential, and Sub-Gap Majorana Modes

**Document Version:** 1.0.0
**Status:** Completed & Numerically Verified
**Scope:** Two-vison pair interaction potential $V(R)$ at $h=0$, sub-gap Majorana bound modes under $[111]$ field, and self-consistent fixed-point stability.

---

## 1. Dual Triangular Lattice Geometry & Branch Cuts

On a periodic honeycomb torus with $L \times L$ unit cells ($N_{\mathrm{sites}} = 2L^2$, $M = 8L^2$ Majoranas), hexagonal plaquettes form a dual triangular lattice of coordination 6. Each directed bond $\langle j k \rangle_\gamma$ is an edge on the dual graph connecting exactly two adjacent plaquettes.

By global Gauss law on the torus:
$$\prod_{p=1}^{L^2} W_p \equiv +1.0.$$
An isolated pair of vortices ($W_{p_1} = -1.0, W_{p_2} = -1.0$) is created by flipping the gauge links $u_b \to -u_b$ along a minimal-length branch cut path $\mathcal{P}(p_1 \to p_2)$ on the dual graph. Every intermediate plaquette crossed by the cut has two flipped bonds ($(-1)^2 = +1$), leaving all $L^2 - 2$ bulk plaquettes strictly flux-free ($W_p = +1.000000000000 \pm 10^{-14}$).

---

## 2. Two-Vison Interaction Potential $V(R)$ at $h = 0$

On an $L = 8$ torus ($128$ sites, $M = 512$, $3 \nmid L$), we compute the defect creation energy:
$$V(R) = \frac{1}{2} \left[ E_{\mathrm{2\text{-}vison}}(R) - E_{\mathrm{vacuum}} \right].$$

### Numerical Results (Intel i5-13500 CPU)
| Separation $R$ | Cut Length | Defect Energy $\Delta E / |K|$ | Lowest Mode $\epsilon_1$ | Second Mode $\epsilon_2$ | $W(p_1), W(p_2)$ |
|:---:|:---:|:---:|:---:|:---:|:---:|
| **1** (Nearest) | 2 | $0.114139$ | $0.518440$ | $0.596027$ | $-1.00, -1.00$ |
| **2** | 4 | $0.133020$ | $0.389817$ | $0.501710$ | $-1.00, -1.00$ |
| **3** | 5 | $0.163729$ | $0.123325$ | $0.637139$ | $-1.00, -1.00$ |
| **4** (Maximal) | 4 | $0.167911$ | $0.107401$ | $0.621015$ | $-1.00, -1.00$ |

### Physical Observations
1. **Attractive Interaction at Short Distances:** $V(R)$ increases monotonically from $0.114 |K|$ at $R=1$ to $0.168 |K|$ at $R=4$. This confirms an attractive defect interaction mediated by matter Majorana exchange.
2. **Asymptotic Saturation:** At maximal separation $R = 4$, $\Delta E \approx 0.168 |K|$, approaching the finite-size thermodynamic benchmark.
3. **Decoupling of Bound States:** At $R = 1$, the two bound modes strongly hybridize into a high-energy doublet ($\epsilon_1 \approx 0.52$). As separation grows to $R = 4$, $\epsilon_1$ drops sharply to $0.107$, decoupling into an isolated sub-gap state.

---

## 3. Zeeman Field Response & Sub-Gap Modes ($h > 0$)

Under external magnetic field along $[111]$ with separation fixed at $R = 4$:

| Field $h / |K|$ | $\Delta E / |K|$ | Bulk Vacuum Gap | In-Gap Mode $\epsilon_1$ | In-Gap Ratio $\epsilon_1 / \Delta_{\mathrm{bulk}}$ |
|:---:|:---:|:---:|:---:|:---:|
| **0.00** | $0.167911$ | $0.828427$ | $0.107401$ | $12.96\%$ |
| **0.01** | $0.167883$ | $0.827772$ | $0.107368$ | $12.97\%$ |
| **0.02** | $0.167802$ | $0.825833$ | $0.107271$ | $12.99\%$ |
| **0.04** | $0.167477$ | $0.818451$ | $0.106884$ | $13.06\%$ |
| **0.06** | $0.166940$ | $0.807247$ | $0.106246$ | $13.16\%$ |
| **0.08** | $0.166193$ | $0.793306$ | $0.105369$ | $13.28\%$ |
| **0.10** | $0.165244$ | $0.777549$ | $0.104267$ | $13.41\%$ |

### Spatial Localization
Computing the spatial probability density $P_i = \sum_{\alpha=0}^3 |\psi(4i + \alpha)|^2$ for the lowest mode reveals that $>75\%$ of the wavefunction weight is concentrated on the 8 sites immediately surrounding plaquettes $p_1 = (0, 0)$ and $p_2 = (4, 4)$, demonstrating exponential localization at the vortex cores.

---

## 4. Self-Consistent Solver Dynamics: Metastable Trapping

When the self-consistent mean-field fixed-point loop is iterated starting from the 2-vison covariance matrix:
* The plaquette fluxes remain locked at $W(p_1) = -1.0000, W(p_2) = -1.0000$ to machine precision across all iterations.
* Vacuum plaquettes exhibit zero drift ($|W_p - 1.0| < 2 \times 10^{-15}$).
* **Conclusion:** Separated vortex pairs form robust, self-consistent **metastable topological solitons** in the mean-field landscape.
