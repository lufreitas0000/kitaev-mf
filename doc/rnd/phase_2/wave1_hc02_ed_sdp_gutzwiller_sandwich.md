---
id: wave1-hc02-ed-sdp-sandwich
title: "Wave 1 HC-02: Gutzwiller-Projected scMFT vs Lanczos ED and Clarabel SOS SDP Bounds"
questions:
  - L0-Q2
  - L2-Q3
contracts:
  - HC-02
depends_on:
  - t2_1_3_variational_overlaps_and_2rdm_bounds
---

# Wave 1 `HC-02`: Variational `scMFT` vs Exact ED and Certified `Clarabel` SDP Sandwich

**Provenance:** `scripts/wave1/run_hc02_variational_sandwich.py` @ commit `44692e54e8baf400e6b20ea3c2e0f064b607c9c0` (`data/processed/benchmarks/wave1/hc02_ed_sdp_scmft_sandwich.{dat,meta.json}`, SHA-256 `94e99ef5f0e9af1a9f66670f907ead83c2015455a84c726304e218c9ae52c150`). External ED & SDP datasets ingested via `data/external/kitaev_numerics/`.

## 1. Local Constraint Orientation & Exact Zero-Variance Principle

Mapping the $2N$-qubit Majorana-Jordan-Wigner Fock space onto the local parity sector $\mathcal{D}_s = c_{4s} c_{4s+1} c_{4s+2} c_{4s+3} = -Z_{2s} Z_{2s+1} = +1$ (spanned by $\{|01\rangle, |10\rangle\}$ per site with $\sigma^x = -Z$, $\sigma^y = +Y$, $\sigma^z = +X$, independently verified against all $2^{2N}$ Fock states to `3.84e-16` on $N=2$) eliminates the spurious positive energy (`+3.3376`) present in the unreviewed draft `t2_1_3`.

Furthermore, clean zero-flux periodic boundary conditions $(+1,+1)$ project identically to zero on finite honeycomb tori because the free-fermion ground state carries odd gauge parity $\pi_G = -1$. Selecting the even-parity antiperiodic/topological sectors `bc in {(-1,-1), (-1,+1), (+1,-1)}` restores non-vanishing variational states obeying the Zero-Variance Principle at $h=0$:
- **$N = 8$ ($2\times 2$):** $E_{\mathrm{proj}}^{\mathrm{exact}}(0) = E_{\mathrm{proj}}^{\mathrm{VMC}}(0) = -0.8090169944 = -(1+\sqrt{5})/4$ ($e_{\mathrm{loc}}\text{-std} = 2.22 \times 10^{-16}$).
- **$N = 12$ ($2\times 3$):** $E_{\mathrm{proj}}^{\mathrm{exact}}(0) = E_0^{\mathrm{ED}}(0) = -0.8166902250$ to all 10 reported digits ($e_{\mathrm{loc}}\text{-std} = 3.08 \times 10^{-16}$, correlation gap `0.00000000`).
- **$N = 18$ ($3\times 3$):** $E_{\mathrm{proj}}^{\mathrm{VMC}}(0) = -0.7939723679 = E_0^{\mathrm{ED}}(0)$ ($e_{\mathrm{loc}}\text{-std} = 3.09 \times 10^{-16}$).

## 2. Certified Sandwich Across $h \in [0, 0.50]$ (Regimes III & IV)

All **24 rows** across $N \in \{8, 12, 18\}$ and $h \in \{0.0, 0.02, 0.0485, 0.10, 0.25, 0.3325, 0.40, 0.50\}$ unconditionally satisfy $E_{\mathrm{SDP}}^{\mathrm{cert}}(h) \le E_0^{\mathrm{ED}}(h) \le E_{\mathrm{proj}}(h)$ (`sandwich_valid = 1`, dual slack $0.0$, 50 VMC chains with $\mathrm{BC}_a$ $B=2000$ 95% CIs and integrated autocorrelation $\tau_{\mathrm{int}}$):

- At $h = 0.25$ on $N = 12$: $E_{\mathrm{SDP}}^{\mathrm{cert}} = -1.031028$, $E_0^{\mathrm{ED}} = -0.934055$, $E_{\mathrm{proj}}^{\mathrm{exact}} = -0.896120$, $E_{\mathrm{proj}}^{\mathrm{VMC}} = -0.896210$ ($\mathrm{CI}_{95\%} = [-0.897951, -0.894243]$, $\tau_{\mathrm{int}} = 1.72$, $n_{\mathrm{eff}} = 1454.1$), gaining $-0.084349\,|K|$/site over unprojected `scMFT` ($-0.811771$).
- At $h = 0.3325$ on $N = 18$: $E_{\mathrm{proj}}^{\mathrm{VMC}} = -0.959371$ ($\mathrm{CI}_{95\%} = [-0.961796, -0.957123]$) vs $E_0^{\mathrm{ED}} = -0.999602$, reducing the unprojected `scMFT` discrepancy ($-0.832532$) from $0.167$ to $0.040\,|K|$/site.
- Above $h_{c2}^+ = 0.400$ (polarized product state), Gutzwiller projection yields zero lowering over `scMFT` ($E_{\mathrm{proj}}^{\mathrm{exact}} = -h - 0.5$, e.g. $-0.9000$ at $h=0.40$), isolating $\Delta E_{\mathrm{corr}} \approx 0.157\,|K|$/site of genuine transverse/flux quantum fluctuations captured by ED near the transition.
