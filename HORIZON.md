# Project Horizon — `kitaev-mf`

**Authority:** Intra-project ordered milestone horizon (`Phase -> Milestone -> Wave`) for `kitaev-mf`. Current active sprint gate lives in [`CURRENT_SPRINT.md`](CURRENT_SPRINT.md); incubating research questions live in [`BACKLOG.md`](BACKLOG.md); cross-repository QSL hierarchy (`L0–L3`) lives in [`../GUIDING_QUESTIONS.md`](../GUIDING_QUESTIONS.md).

---

## Horizon 1 (`RECON-20261116`, Phase 1 Foundation) — Reconciled on Linear `main` (**Completed**)

- Verified `A01`–`A20` Clifford/CAR algebra, `Gr(4, 8)` & `8N` real skew-symmetric Grassmannian manifold solvers, cubic Newton–Schulz matrix sign purification, Parlett–Reid 12-Majorana plaquette Pfaffians, KKT post-spinodal audit (`h_{c2}^+ = 0.400`), two-vison attractive branch-cut potential `K-020`, and Living Lab Book / Monograph (`K-001`–`K-023`).
- External Spin-ED and 2-RDM SDP responsibilities are delegated cleanly to [`kitaev-numerics`](../kitaev-numerics/) via provenance-stamped `data/external/` contracts (`ADR-003`).

---

## Horizon 2 (`M2.1` / `Wave 1` — `HC-01..03` Execution Complete, Awaiting PI Threshold Lock)

- **`Wave 1` Pre-Registered Checks (`doc/hypotheses/HC-01.yaml`–`HC-03.yaml`, `status: draft`):**
  1. `HC-01` (`K-024`, `IDEA-001` & `K-019-A`): Executed Feshbach–Schur continuum decomposition (`c_MF^cont = 10.641666`, `p = 2.999694`, `||lambda*||/h = 2.237814`), isolating KKT $\lambda^* G_{QQ}$ activation and $PQ$ Hartree screening (`<= 4.04e-15` without $\lambda^*$).
  2. `HC-02` (`K-025`, `IDEA-003` & `K-023`): Verified `24/24` $E_{\mathrm{SDP}}^{\mathrm{cert}} \le E_0^{\mathrm{ED}} \le E_{\mathrm{proj}}$ bounds across $N \in \{8, 12, 18\}$ with 50-chain $\mathrm{BC}_a$ 95% CIs and Zero-Variance $h=0$ exact projection (`-0.8166902250` on $N=12$).
  3. `HC-03` (`K-026`, `M2` Vison–Vison Potential & MZM Scaling): Verified `<4.4e-9` static chiral match to `kitaev-walk-geometry:HC-11` (`kappa_static = kappa_KWG/2`), converged 94 stratified chiral points ($L \in \{8..20\}$), and identified period-3 intervalley modulation explaining raw monotone $\Delta\mathrm{BIC} \in [2.39, 4.68]$ (`TODO(owner)` threshold lock open).

---

## Horizon 3 (`Wave 2–3` — Finite-Temperature, Dynamical Response & Extended $K\text{–}J\text{–}\Gamma\text{–}\Gamma'$)

- **`Wave 2` (`M2.1`–`M2.4` / Phase 2):** Static thermodynamics $C_v(T)/T$, $S(T)$ across $T_L \approx 0.012$ and $T_H \approx 0.375$, Gutzwiller projection sensitivity, and defect mobility.
- **`Wave 3` (`Phase 3` & `Phase 5` Extended Kitaev):** Port exact symbolic `SymPy` $(J, K, \Gamma, \Gamma')$ + fourth-order ring-exchange tensors from [`kitaev-derivation`](../kitaev-derivation/) into multi-orbital `scMFT` supercells and coordinate material-regime bilayers with [`ksl-sc-bilayer`](../ksl-sc-bilayer/).
