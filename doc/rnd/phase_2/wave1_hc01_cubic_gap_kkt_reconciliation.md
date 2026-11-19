---
id: wave1-hc01-cubic-gap
title: "Wave 1 HC-01: Feshbach-Schur Cubic Dirac Gap and KKT Counter-Field Dressing"
questions:
  - L0-Q1
  - L2-Q2
contracts:
  - HC-01
depends_on:
  - m1_10_t1_k019a_second_order_erratum
  - m1_10_t1_wave2_cubic_gap_evidence_pilot
---

# Wave 1 `HC-01`: Feshbach–Schur Cubic Dirac Gap Reconciliation

**Provenance:** `scripts/wave1/run_hc01_cubic_gap.py` @ commit `44692e54e8baf400e6b20ea3c2e0f064b607c9c0` (`data/processed/benchmarks/wave1/hc01_cubic_gap_reconciliation.{dat,meta.json}`, SHA-256 `7d8a7557787d55909e3284b25cb7387a171fb5ea454512e23052cf82feced767`).

## 1. Epistemic Scope (Regime I: Deterministic Continuum $k$-Mesh Extrapolation)

Under [QSL Metrology Guidelines](../../../../contracts/qsl_metrology_guidelines.md) Regime I, continuum extrapolation uses quadratic Vandermonde/Richardson regression in $1/N_k$ across deterministic Monkhorst–Pack meshes $N_k \in \{24, 36, 48, 72\}$ in the pure ferromagnetic Kitaev `[111]` phase ($K=-1$, canonical $4\mathcal{H}_{\mathrm{Bloch}}$ units matching Kitaev's single-particle excitation convention).

## 2. Exact Feshbach–Schur Anatomy at the Dirac Point $\mathbf{K} = (1/3, 2/3)$

At $\mathbf{K}$, the unperturbed pure-Kitaev spectrum partitions into a 2D massless matter kernel $P$ ($\varepsilon_P = 0$) and a 6D gapped gauge/matter complement $Q$ bounded by $\Delta_{\mathrm{dimer}}^{\mathrm{can}} \approx 2.09943\,|K|$. Writing the self-consistent perturbation $\delta \mathbf{V}(h) = \mathbf{H}_{\mathrm{can}}(h) - \mathbf{H}_0$ in the $P \oplus Q$ eigenbasis yields three exact structural results across all $h \in [0.002, 0.030]$:

1. **Zero Diagonal Shift & Zero Bare-Zeeman Cubic Gap:**
   $\|\delta V_{PP}\|_\infty \le 7.73 \times 10^{-15}$ and $\Delta_{\mathrm{Schur}}^{\mathrm{bare\ Zeeman}} \le 1.60 \times 10^{-15}$ everywhere because $\mathbf{N}^\gamma$ vanishes identically on $Q$ ($Q \mathbf{V}_{\mathrm{Zeeman}} Q \equiv 0$).
2. **KKT Counter-Field Activation in $Q$-Space & $PQ$ Hartree Screening:**
   Setting $\boldsymbol{\lambda}^* = 0$ even while keeping full self-consistent Hartree bond/site shifts $(\delta \mathbf{U}, \delta \mathbf{V})$ produces $\Delta_{3\mathrm{rd}}^{\mathrm{no}\text{-}\lambda} \le 4.04 \times 10^{-15}$. Conversely, retaining the full $\delta V_{PQ}$ vertex with **strictly** $\delta V_{\lambda, QQ}$ inside $Q$ recovers $\Delta_{3\mathrm{rd}}^{\mathrm{dressed}} = +P\,\delta V_{PQ} D_0^{-1} \delta V_{QQ} D_0^{-1} \delta V_{QP} P$ to `< 1e-8` relative accuracy (`1.06453444e-05` for both columns at $h=0.01, N_k=48$). Meanwhile, freezing $(\mathbf{U}^{(0)}, \mathbf{V}^{(0)})$ while applying $\boldsymbol{\lambda}^*$ alone induces 95% destructive interference in $\delta V_{PQ}$, reducing the unscreened estimate to $5.008 \times 10^{-7}$ until self-consistent $(\delta \mathbf{U}, \delta \mathbf{V})$ screens $\lambda^* G_{PQ}$.
3. **Continuum Prefactor & Scaling Exponent:**
   - $c_{\mathrm{MF}}(N_k) = \Delta_{\mathrm{can}}(0.01) / h^3$: `10.654137` ($N_k=24$), `10.644855` ($36$), `10.642438` ($48$), `10.641100` ($72$).
   - Continuum limit ($1/N_k \to 0$): **$c_{\mathrm{MF}}^{\mathrm{cont}} = 10.641666$**, $\|\boldsymbol{\lambda}^*\|/h \to 2.237814 \approx \sqrt{5}$, and fitted small-$h$ power law across $h \in [0.002, 0.030]$ at $N_k = 48$ gives **$p = 2.999694$**. *(Note: $h=0.002, N_k=48$ reached iteration cap `80` with `conv=0` under strict $10^{-12}$ tolerance while converging smoothly in gap value $\Delta/h^3 = 10.643303$.)*
