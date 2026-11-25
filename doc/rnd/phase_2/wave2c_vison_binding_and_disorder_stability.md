---
id: wave2c-vison-binding-and-disorder-stability
title: "Wave 2C (UOW-05 / Q-L2-1): Isolated Two-Vison Binding V(R, L), Exact Gutzwiller Parity Projection & Disorder/Vacancy KKT Stability"
questions:
  - Q-L2-1
contracts:
  - HC-02
  - HC-03
depends_on:
  - wave1-hc03-multi-l-visons
  - m1_9_t2_disordered_solver_scaling
---

# Wave 2C (`UOW-05` / `Q-L2-1`): Isolated Two-Vison Binding $V(R, L)$, Exact Gutzwiller Projection & Disorder/Vacancy Stability

**Provenance:**
- Services: `src/kitaev_mf/services/{vison_binding_wave2c,disorder_stability_wave2c,wave2c_benchmark}.py` + hardened `gutzwiller_{projection,stability}.py` (`UOW-04` `F1`–`F4` resolved)
- ACID Datasets (`# Git Dirty:` & `base_git_sha` tagged):
  - `data/processed/benchmarks/wave2/q_l2_1_two_vison_binding_projection.{dat,meta.json}` (95 rows: 14 exact Zero-Variance Gutzwiller `L=4..8` rows + 81 chiral `L=8..15` rows)
  - `data/processed/benchmarks/wave2/q_l2_1_disorder_vacancy_stability.{dat,meta.json}` (28 rows: 10 multi-seed bond-disordered two-vison rows + 15 inhomogeneous vacancy rows + 3 Willans–Chalker–Moessner vacancy flux-binding rows)
- **Human-Owned Hypothesis Locks Preserved:** `doc/hypotheses/HC-01.yaml`–`HC-03.yaml` (`TODO(owner)`) and `CONVENTIONS.md` remain completely untouched.

---

## 1. Physical Mental Image & Direct Answer to `Q-L2-1`

> **Physical Picture:** In the pure Kitaev spin liquid ($K_{\mathrm{Pauli}}=-1$), two $Z_2$ vortices (visons) connected by a dual branch cut are **topologically deconfined** quasiparticles whose pair cost saturates to $2\Delta_v \approx 0.3072\,|K|$ ($\Delta_v = 0.1536\,|K|$, Kitaev 2006; adjacent pair $V(1) = 0.263313\,|K|$, Panigrahi–Coleman–Tsvelik 2023) with attractive short-range binding $V_{\mathrm{bind}}(1) \approx -0.044\,|K|$. On finite periodic tori, however, **unprojected `scMFT`** at fixed periodic boundary conditions $bc=(+1,+1)$ suffers sharp Pedrocchi–Chesi–Loss (2011) parity distortions—even producing unphysical **negative** $V_{\mathrm{scMFT}}(R=2,3) < 0$ on resonant $3 \mid L$ tori ($L=6$) and odd-parity ghost vacua on $L=4, 7, 8$. Enforcing exact Gutzwiller $\mathcal{D}_j = +1$ projection via the 1D $yz$-loop holonomy law $S_x = -bc_y \prod_{b \in C_x} u_b^{y,z}$ plus canonical zero-mode parity pairing unlocks exact Zero-Variance many-body spin energies up to $L=8$ ($N=128$ spins, $2^{128}$ Hilbert space), restoring strict Lieb positivity everywhere. Under a chiral gap $\kappa_{\mathrm{static}} > 0$, intervalley Dirac interference $\cos(\frac{2\pi}{3} R + \phi)$ modulates the Majorana Zero Mode (MZM) hybridization envelope (Lahtinen 2011), decisively resolving Wave-1's monotone-fit ambiguity ($\Delta\mathrm{BIC} = +18.7 \text{ to } +64.2 \gg 10$). Finally, full inhomogeneous ($D = 80 L^2$) Kenney–Laub Padé `[3/2]` `scMFT` maintains 100% two-vison trapping and KKT gauge neutrality ($\max |g_s^\gamma| < 10^{-10}$) across quenched bond disorder ($\delta \le 0.20$) and sublattice-balanced vacancies ($n_{\mathrm{vac}} = 2$), while exposing a **sharp Lakatosian boundary**: an unpaired single vacancy ($|N_A - N_B| = 1$) on a gapped ($3 \nmid L$) or $C_3$-broken cluster at $h > 0$ hits a discontinuous dangling-Majorana level crossing at $\lambda_{j_\gamma}^\gamma \approx 2 h_\gamma$ ($\Delta \langle G_1^z \rangle > 1.0$, specifically `-0.468 -> +0.558`, jump `+1.026` in `probe_single_vacancy_kkt_level_crossing`), proving that pure $T=0$ Gaussian sign projection cannot satisfy KKT neutrality on an isolated spin-$1/2$ Curie vacancy moment without finite-$T$ / mixed-state occupation of the crossing doublet.

---

## 2. Exact Zero-Variance Gutzwiller $V_{\mathrm{proj}}(R, L)$ vs. Unprojected $V_{\mathrm{scMFT}}(R, L)$ ($h=0$)

Exploiting the `UOW-04 F1` exact 1D $yz$-loop holonomy theorem ($S_x = -bc_y \prod_{b \in C_x} u_b^{y,z}$), single-`eigh` column-swap Thouless evaluation (`16x` speedup), canonical $(v_0, \Xi v_0)$ zero-mode parity pairing (`allow_zero_mode=True` recovering exact $2^{18}$ Spin-Lanczos ED $E_0^{\mathrm{ED}}/N = -0.7621342368$ on $3\times 3$ `two_vison`), and Zero Variance at $h=0$ ($\mathrm{std}(E_{\mathrm{loc}}) \le 4.84 \times 10^{-13}$ across 4 independent $y,z$-bond-flipped bitstrings up to $N=128$):

| $L$ ($N$) | Class | $R$ | $bc_{2v}^{\mathrm{proj}}$ (`qp`) | $V_{\mathrm{proj}}(R, L)$ | $V_{\mathrm{scMFT}}(R, L)$ | Parity Shift $V_{\mathrm{proj}} - V_{\mathrm{scMFT}}$ | Probe $\mathrm{std}(E_{\mathrm{loc}})$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `4` (`32`) | `3!\|L` | `1, 2` | `(1,1)`, `(1,-1)` (`0`) | `+0.0520, +0.1791` | `+0.3402, +0.4673` | `-0.2882` (odd-$bc$ ghost vacuum) | `<= 1.37e-13` |
| `5` (`50`) | `3!\|L` | `1, 2` | `(1,1)`, `(1,1)` (`0`) | `+0.2991, +0.3219` | `+0.2784, +0.3012` | `+0.0207` | `<= 4.84e-13` |
| **`6` (`72`)** | **`3\|L`** | **`1, 2, 3`** | `(1,1)` (`0`), `(1,-1)` (`1`), `(-1,-1)` (`0`) | **`+0.5257, +0.5073, +0.2603`** | **`+0.0395, -0.0589, -0.1403`** | **`+0.486, +0.566, +0.401`** (restores Lieb $>0$!) | `<= 3.52e-14` |
| `7` (`98`) | `3!\|L` | `1, 2, 3` | `(1,1)`, `(1,1)`, `(1,-1)` (`0`) | `+0.1403, +0.2095, +0.2201` | `+0.2875, +0.3566, +0.3033` | `-0.1472, -0.1472, -0.0832` | `<= 1.45e-13` |
| **`8` (`128`)** | **`3!\|L`** | **`1..4`** | `(1,1)`, `(1,1)`, `(1,-1)`, `(1,-1)` (`0`) | **`+0.2511, +0.2947, +0.2727, +0.3148`** | `+0.2722, +0.3158, +0.2634, +0.3358` | `-0.0211, -0.0211, +0.0093, -0.0211` | `<= 6.74e-14` |

- **Bayesian Deconfinement vs. Linear Confinement (`L=8`, exact projected spin states):** Comparing $\mathcal{M}_1: V = V_\infty$ against $\mathcal{M}_0: V = 0.20\,N_{\mathrm{cut}} + V_0$ yields **$\Delta\mathrm{BIC}_{01} = +17.24 > 10$** ($B_{10} = 5{,}531 \gg 148$), vanishing unconstrained string tension $\sigma_{\mathrm{eff}} = 0.0169 \ll 0.20$, and $\Delta_v = V_\infty / 2 = 0.1417\,|K|$ ($7.8\%$ from Kitaev's $0.1536\,|K|$; $V_{\mathrm{proj}}(R=4, L=8) = 0.3148\,|K|$ within $2.5\%$ of $2\Delta_v = 0.3072\,|K|$). Pooling non-resonant $L \in \{5, 7, 8\}$ gives $\Delta\mathrm{BIC}_{01} = +21.69$ ($\sigma_{\mathrm{eff}} = 0.0195$).

---

## 3. Oscillatory Majorana Hybridization Envelope (Resolving Wave-1 `HC-03`)

Fitting the Lahtinen (2011, Sec. 3.1, Eq. 5) Dirac-intervalley oscillatory Majorana envelope $\mathcal{M}_{\mathrm{osc\text{-}exp}}: \omega_{\mathrm{MZM}}/\Delta_{\mathrm{can}} = A\,|\cos(\frac{2\pi}{3} R + \phi)|\,e^{-R/\xi_M}$ ($k_K = 2\pi/3$ fixed analytically, 1 phase parameter $\phi$) across $L \in \{8, 10, 12, 14, 15\}$ ($n=22$ points per $\kappa_{\mathrm{static}}$ for $R \ge 2$):

| $\kappa_{\mathrm{static}}$ | $\xi_{\mathrm{osc}}$ (lin / log) | $\Delta\mathrm{BIC}(\mathcal{M}_{\mathrm{osc\text{-}exp}} \text{ vs } \mathcal{M}_{\mathrm{mono\text{-}exp}})$ | $\Delta\mathrm{BIC}_{\log}(\mathcal{M}_{\mathrm{osc\text{-}exp}} \text{ vs } \mathcal{M}_{\mathrm{osc\text{-}pow}})$ | `leave_lmin` (osc vs mono) | `leave_lmax` (osc vs mono) | $V_{\mathrm{bind}}(R=1, L=15)$ |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `0.04` | `2.112` / `1.922` | **`+18.70`** | `+5.15` | `0.299` vs `0.512` | `0.421` vs `0.530` | `-0.1004 |K|` |
| `0.06` | `1.660` / `1.545` | **`+41.45`** | `+7.49` | `0.184` vs `0.475` | **`0.122`** vs `0.414` | `-0.1241 |K|` |
| `0.10` | `1.277` / `1.199` | **`+64.24`** | `+8.87` | `0.301` vs `0.612` | **`0.065`** vs `0.570` | `-0.1522 |K|` |

- Decisively resolves Wave-1's inconclusive monotone fit ($\Delta\mathrm{BIC} \in [2.39, 4.68]$): the oscillatory Majorana envelope beats monotone exponential by $\Delta\mathrm{BIC} \in [+18.7, +64.2] \gg 10$, cuts `leave_lmax` error by up to $8.7\times$ (`0.570 -> 0.065`), and restores strict $1/\kappa$ contraction of $\xi_M(\kappa)$.

---

## 4. Inhomogeneous Disorder, Dilute Vacancies, Padé `[3/2]` Advantage & Lakatosian Single-Vacancy Crossing

1. **Quenched Bond-Disordered Two-Vison Ensembles ($\delta \in \{0.00, 0.05, 0.10, 0.15, 0.20\}$, $\kappa_{\mathrm{static}}=0.06$, $L=6$ [8 seeds] & $L=8$ [6 seeds]):**
   - **100% convergence & 100% two-vison topological trapping** (`defs == 2`, $\min |W_p^{\mathrm{Pf}}| = 1.000000$, $\max |g_s^\gamma| < 10^{-11}$) across all 58 disordered realizations!
   - **Padé `[3/2]` vs. Newton–Schulz:** Padé `[3/2]` converges in `5.88–6.67` iterations vs. `15.88–20.00` for cubic Newton–Schulz (`2.7–3.0x` fewer iterations, `max |Omega_Pade - Omega_eigh| < 1e-11`).
   - Two-vison energy $V_{\mathrm{pair}}(\delta)$ remains robust (`0.4318 -> 0.4236 [0.3962, 0.4423]` BCa 95% CI on $L=8$ at $\delta=0.20$).
2. **Dilute Non-Magnetic Vacancies (`n_vac in {0, 1, 2}`, full inhomogeneous $D = 80 L^2$ KKT `scMFT`):**
   - **Sublattice-balanced clusters (`n_vac = 0, 2` on $L=4, \delta=0.05$ and $L=6, \delta=0.0$, $h \in \{0, 0.05, 0.10\}$):** 100% converged (`max_g <= 4.32e-11`, exact K-007 double-counting `< 1e-11`, site positivity diamond residual `<= 4.16e-17`). At $\kappa(A_{\mathrm{act}}) = 3.92 \times 10^5$ (`L=6, n_vac=2, h=0.10`), Padé `[3/2]` finishes in **11** iterations (`1.87e-13` accuracy) while cubic Newton–Schulz needs **38** iterations.
   - **Exact Odd-Active-Site Improper $O(4)$ Pfaffian Parity (`UOW-02 F2`):** On clean resonant `L=6, n_vac=1, h=0.10` ($N_{\mathrm{act}} = 71$ odd), inhomogeneous `scMFT` converges cleanly (`max_g = 7.59e-11`, `pade_it = 10` vs `ns_it = 33`) and verifies $\operatorname{Pf}(\mathbf{P}_0^{\oplus 71}\Omega_{\mathrm{act}}^{\mathrm{kmf}}\mathbf{P}_0^{\oplus 71}) = -\operatorname{Pf}(\Omega_{\mathrm{act}}^{\mathrm{kmf}}) = -1.00000$, while even-$N_{\mathrm{act}}$ clusters ($N_{\mathrm{act}} \in \{30, 32, 70, 72\}$) preserve sign (`+1`).
   - **Willans–Chalker–Moessner / Kao–Knolle Eq. (8) Vacancy-Pair Flux Binding ($h=0$):** Self-consistent A–B vacancy pairs lower total energy by binding fluxes to the merged 12-site vacancy super-plaquettes: $E_{\mathrm{bind}} = (E_{\mathrm{bound}} - E_{\mathrm{zero}})/2 = -0.0860\,|K|$ ($L=6$, resonant $3 \mid L$), **$-0.0227\,|K|$ ($L=8$)**, and **$-0.0173\,|K|$ ($L=10$)** (`max_g < 1e-11`; cf. thermodynamic limit $-0.0268\,|K|$, Kao–Knolle et al. PRX 11, 011034 (2021) Eq. 8).
   - **Lakatosian Counterexample — Pure-State $T=0$ KKT Breakdown on Unpaired Gapped/Disordered Vacancies:** On `L=4` ($3 \nmid L$, bulk gap $0.88$) or $C_3$-broken `delta = 0.05` single-vacancy clusters (`n_vac = 1`, $|N_A - N_B| = 1$) at $h > 0$, the dangling neighbor Majorana $b_{j_\gamma}^\gamma$ couples to $c^0$ via $\lambda_{j_\gamma}^\gamma - 2 h_\gamma + \delta_{\mathrm{Hartree}}$. Driving $\langle G_{j_\gamma}^\gamma \rangle \to 0$ pins $\lambda_{j_\gamma}^\gamma$ onto a zero-temperature level crossing (`sv_min ~ 6.3e-6`, $\kappa \sim 10^6\text{–}10^7$) with a discontinuous jump $\langle G_1^z \rangle: -0.468 \to +0.558$ ($\Delta \langle G_1^z \rangle = +1.026$ at $\lambda_1^z = 2 h_z \pm 0.01$ in `probe_single_vacancy_kkt_level_crossing`, reaching $+1.613$ at further DIIS cycling). Thus **no pure $T=0$ Gaussian projector ($\Omega^2 = -I$) can satisfy KKT neutrality** for an isolated Curie-like vacancy spin-$1/2$ moment in a gapped bulk; resolving it requires finite-$T$ Fermi-Dirac smearing (or mixed half-occupied zero-mode density matrices), linking `Q-L2-1` directly to finite-$T$ Wave-2 thermodynamics.

---

## 5. Primary Literature Verification

1. **Kitaev, Ann. Phys. 321, 2–111 (2006), App. A:** Isolated vison gap $\Delta_v \approx 0.1536\,|K|$ in Pauli normalization ($H = -|K|\sum \sigma_i^\alpha \sigma_j^\alpha$).
2. **Panigrahi, Coleman & Tsvelik, PRB 108, 045151 (2023), Eq. (28):** Closed-form adjacent vison-pair gap $V(1) = 0.263313(6)\,|K|$, attractive binding $V(1) - 2\Delta_v = -0.0439\,|K|$.
3. **Pedrocchi, Chesi & Loss, PRB 84, 165414 (2011), Sec. IV & Fig. 4:** Physical parity projection on $3 \mid L$ vs. $3 \nmid L$ tori, zero-mode parity survival, and finite-size two-vortex parity oscillations.
4. **Lahtinen, New J. Phys. 13, 075009 (2011), Sec. 3.1, Eq. (5):** Oscillatory MZM hybridization splitting $\mathcal{E}(R) = \Delta_f \cos(2\pi R/\lambda) e^{-R/\xi}$ with period-3 intervalley interference and $\xi \sim \kappa^{-1}$.
5. **Willans, Chalker & Moessner, PRL 104, 237203 (2010); Kao, Knolle, Halász, Moessner & Perkins, PRX 11, 011034 (2021), Sec. III & Eq. (8); Knolle, Moessner & Perkins, PRL 122, 047202 (2019):** Vacancy quasilocalized modes, merged 12-site vacancy flux-binding $E_{\mathrm{bind}} = -0.0268\,|K|$, dangling-Majorana Zeeman response, and bond-disorder stability.
