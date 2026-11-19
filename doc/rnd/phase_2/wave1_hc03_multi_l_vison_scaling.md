---
id: wave1-hc03-multi-l-visons
title: "Wave 1 HC-03: Stratified Multi-L Two-Vison Potential, Chiral NNN Bridge, and Intervalley MZM Oscillation"
questions:
  - L0-Q3
  - L2-Q4
contracts:
  - HC-03
depends_on:
  - m1_10_t2_two_vison_energetics_and_zero_modes
---

# Wave 1 `HC-03`: Multi-$L$ Stratified Two-Vison Scaling and Chiral `KWG` Cross-Bridge

**Provenance:** `scripts/wave1/run_hc03_multi_l_visons.py` @ commit `44692e54e8baf400e6b20ea3c2e0f064b607c9c0` (`data/processed/benchmarks/wave1/hc03_multi_l_vison_scaling.{dat,meta.json}`, SHA-256 `d588abb8cc45cea388c87214cce61ee1e243317a08be111a9b655ae76c0f37c7`). External chiral reference imported from `data/external/kitaev_walk_geometry/chiral_vison_tmatrix_scaling.dat`.

## 1. Cross-Repository Chiral Bridge & Topological Trapping Verification

- **Normalization Reconciliation:** In `kitaev-mf` canonical $4\mathcal{H}_{\mathrm{Bloch}}$ units, the gauge-covariant Haldane NNN element $c_{\mathrm{NNN}} = 2\,\kappa_{\mathrm{static}}$ yields bulk Dirac gap $\Delta_{\mathrm{bulk}}^{\mathrm{can}} = 12\sqrt{3}|\kappa_{\mathrm{static}}|$, which equals `kitaev-walk-geometry`'s $\Delta_K^{\mathrm{KWG}} = 6\sqrt{3}|\kappa_{\mathrm{KWG}}|$ at **$\kappa_{\mathrm{static}} = \kappa_{\mathrm{KWG}} / 2$** (disambiguating the `/ 6` arithmetic comment in `HC-03.yaml`). Direct comparison on adjacent bond-cut tori reproduces `chiral_vison_tmatrix_scaling.dat` $E_0/\Delta_K$ to **`< 4.4e-9`** across both resonant torus sizes tested:
  - $\kappa_{\mathrm{KWG}}=0.08$ ($\kappa_{\mathrm{static}}=0.04$): `0.28256246` ($L=12$), `0.31331511` ($L=18$).
  - $\kappa_{\mathrm{KWG}}=0.12$ ($\kappa_{\mathrm{static}}=0.06$): `0.24082086` ($L=12$), `0.24758574` ($L=18$).
- **Matrix Sign Purification vs `eigh` Cross-Check:** At $L=12, h=0.20, \kappa_{\mathrm{static}}=0$, `purify_newton_schulz` agrees with real symmetric `eigh(-A^2)` to `6.46e-12` across the full $(2L^2 \times 4)$-Majorana covariance matrix.
- **Stratified Trapping & Gauge Rigor:** All 94 chiral two-vison points ($\kappa_{\mathrm{static}} \in \{0.04, 0.06, 0.10\}$, $L \in \{8, 10, 12, 15, 16, 18, 20\}$ covering both $3 \mid L$ resonant and $3 \nmid L$ non-resonant tori) converge cleanly (`conv = 1`, `defs = 2`, $\min |W_p^{\mathrm{Pf}}| = 1.000000$, $\max |g_s^\gamma| \le 1.644 \times 10^{-11}$). Attractive two-vison binding saturates smoothly to `0.398 |K|` ($\kappa=0.04$), `0.451 |K|` ($\kappa=0.06$), and `0.537 |K|` ($\kappa=0.10$), while `[111]` $h=0.20$ at $\kappa_{\mathrm{static}}=0$ keeps strictly two visons ($\min |W_p^{\mathrm{Pf}}| \approx 0.897$) with strong $3 \mid L$ finite-size Dirac suppression ($V_{\mathrm{pair}}(L=12) = 0.1216\,|K|$ vs $V_{\mathrm{pair}}(L=16) = 0.3125\,|K|$).

## 2. Regime V Model Selection & Intervalley Oscillation Refutation / Refinement

Evaluating the pre-registered **monotone** decay model $\mathcal{M}_{\mathrm{exp}}: \omega/\Delta = A e^{-R/\xi_M}$ against power-law $\mathcal{M}_{\mathrm{pow}}: C R^{-\alpha}$ across raw mixed-phase separations $R \ge 2$ yields **inconclusive / non-decisive** evidence ($\Delta\mathrm{BIC} < 10$ everywhere) and large two-sided holdout errors:
- $\kappa_{\mathrm{static}} = 0.04$: $\Delta\mathrm{BIC} = +2.54$, `leave_lmin` $0.554$, `leave_lmax` $0.381$, $\xi_{\mathrm{fit}} = 4.02$.
- $\kappa_{\mathrm{static}} = 0.06$: $\Delta\mathrm{BIC} = +4.68$, `leave_lmin` $0.521$, `leave_lmax` $0.540$, $\xi_{\mathrm{fit}} = 3.52$.
- $\kappa_{\mathrm{static}} = 0.10$: $\Delta\mathrm{BIC} = +2.39$, `leave_lmin` $0.618$, `leave_lmax` $0.592$, $\xi_{\mathrm{fit}} = 4.78$.

**Scientific Diagnosis:** Inspection of raw rows proves why the raw monotone fit fails: hopping `p2` along the $\mathbf{a}_1$ zigzag axis modulates the overlap by **period-3 intervalley Dirac interference** $\cos(\mathbf{K}\cdot\mathbf{R}+\phi)$ (constructive jumps at $R = 3, 6$, destructive nodes at $R = 2, 4$, e.g. at $\kappa=0.10, L=12$: $\omega/\Delta = 0.3386, 0.0198, 0.0711, 0.0330, 0.00588, 0.00135$ for $R=1..6$). Filtering at **fixed phase** ($R=3 \to 6$ on $L=12$) recovers strict $1/\kappa$ contraction of the correlation length ($\xi_M^{\mathrm{fixed}} \approx 3.03, 1.36, 0.76$ for $\kappa=0.04, 0.06, 0.10$). Future PI threshold locking in `HC-03.yaml` should replace the monotone exponential null with the oscillatory Majorana envelope $|\cos(k_{\mathrm{eff}} R + \phi)| e^{-R/\xi_M}$ or stratify by $R \bmod 3$.
