# Project Horizon — `kitaev-mf`

**Authority:** Intra-project ordered milestone horizon (`Phase -> Milestone -> Wave`) for `kitaev-mf`. Current active sprint gate lives in [`CURRENT_SPRINT.md`](CURRENT_SPRINT.md); incubating research questions live in [`BACKLOG.md`](BACKLOG.md); cross-repository QSL hierarchy (`L0–L3`) lives in [`../GUIDING_QUESTIONS.md`](../GUIDING_QUESTIONS.md).

---

## Horizon 1 (`RECON-20261116`, Phase 1 Foundation) — Reconciled on Linear `main` (**Completed**)

- Verified `A01`–`A20` Clifford/CAR algebra, `Gr(4, 8)` & `8N` real skew-symmetric Grassmannian manifold solvers, cubic Newton–Schulz matrix sign purification, Parlett–Reid 12-Majorana plaquette Pfaffians, KKT post-spinodal audit (`h_{c2}^+ = 0.400`), two-vison attractive branch-cut potential `K-020`, and Living Lab Book / Monograph (`K-001`–`K-023`).
- External Spin-ED and 2-RDM SDP responsibilities are delegated cleanly to [`kitaev-numerics`](../kitaev-numerics/) via provenance-stamped `data/external/` contracts (`ADR-003`).

---

## Horizon 2 (`Wave 1` — `IDEA-001`, `IDEA-003`, `K-019-A` & `M2` Preparation) — **Human-Gated at `HC-01..03.yaml`**

- **`Wave 1` Pre-Registered Checks (`doc/hypotheses/HC-01.yaml`–`HC-03.yaml`, `status: draft`):**
  1. `HC-01` (`IDEA-001` & `K-019-A`): Reconcile the empirical 4-flavor `scMFT` cubic Dirac gap prefactor ($c_{\mathrm{MF}} = 2.5921 \pm 0.0034$) against pure-Kitaev Feshbach–Schur full-$Q$ resolvent perturbation theory (`kitaev-walk-geometry`) and check KKT counter-field vertex dressing.
  2. `HC-02` (`IDEA-003` & `K-023`): Cross-check variational `scMFT` ground-state energy curve $E_{\mathrm{MF}}(h)$ and coexistence domain wall ($\sigma_{\mathrm{DW}} = +0.042|K|$) against unbiased `kitaev-numerics` ED ($N = 18, 24$) and `Clarabel` 2-RDM SDP lower bounds ($E_{\mathrm{SDP}}^{(2)} \le E_0 \le E_{\mathrm{proj}}$).
  3. `HC-03` (`M2` Vison–Vison Potential & MZM Scaling): Test multi-$L$ finite-size scaling of isolated vison-pair interaction $V(R)$ and hybridized MZM splitting $\omega_0(R)$ against `kitaev-walk-geometry` `HC-11` (`TODO(owner)` thresholds required prior to sweeps).

---

## Horizon 3 (`Wave 2–3` — Finite-Temperature, Dynamical Response & Extended $K\text{–}J\text{–}\Gamma\text{–}\Gamma'$)

- **`Wave 2` (`M2.1`–`M2.4` / Phase 2):** Static thermodynamics $C_v(T)/T$, $S(T)$ across $T_L \approx 0.012$ and $T_H \approx 0.375$, Gutzwiller projection sensitivity, and defect mobility.
- **`Wave 3` (`Phase 3` & `Phase 5` Extended Kitaev):** Port exact symbolic `SymPy` $(J, K, \Gamma, \Gamma')$ + fourth-order ring-exchange tensors from [`kitaev-derivation`](../kitaev-derivation/) into multi-orbital `scMFT` supercells and coordinate material-regime bilayers with [`ksl-sc-bilayer`](../ksl-sc-bilayer/).
