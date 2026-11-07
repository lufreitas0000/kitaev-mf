---
id: M1.9-T5
title: "Domain Wall Metastability and Spatial Coexistence inside the Spinodal Window"
status: Completed
date: "2026-11-06"
questions:
  - "Does a composite spatial domain wall between QSL and polarized phases possess positive surface tension sigma_DW > 0 inside the spinodal hysteresis window?"
  - "Are intermediate modulated or stripy states stabilized between the thermodynamic crossing h_c^* = 0.3325 and upper spinodal h_{c2}^+ = 0.400?"
depends_on:
  - M1.6-T5
  - M1.7-T1
contracts:
  - K-011
  - K-023
artifacts:
  - scripts/m1_9/audit_post_spinodal_and_master_table.py
  - data/processed/benchmarks/m1_9/post_spinodal_kkt_audit.dat
---

# Domain Wall Metastability and Spatial Coexistence inside the Spinodal Window

## 1. The Spinodal Coexistence Window
Between the thermodynamic crossing $h_c^* = 0.3325$ and the forward spinodal $h_{c2}^+ = 0.400$, both the topological QSL state and the field-polarized state are local minima of the Hartree-Fock energy functional.
To test whether a mixed spatial state is thermodynamically competitive, we assembled an $L_x \times L_y = 16 \times 8$ supercell ($N_{\mathrm{sub}} = 128$) and initialized a domain-wall seed:
- Left half ($x \in [0, 8)$): QSL ansatz ($W_p = 0.53$, $m=0.18$).
- Right half ($x \in [8, 16)$): Polarized ansatz ($W_p = 0.00$, $m=0.57$).

## 2. Relaxation and Domain Wall Profile
Under self-consistent DIIS iterations at $h = 0.335$:
1. The domain wall relaxes into a narrow boundary layer of width $\xi_{\mathrm{DW}} \approx 1.5$ lattice spacings.
2. The domain wall energy is strictly positive:
$$
\sigma_{\mathrm{DW}} = E_{\mathrm{DW}} - \frac{1}{2}(E_{\mathrm{QSL}} + E_{\mathrm{pol}}) = +0.042 |K| > 0.
$$
Because the surface tension $\sigma_{\mathrm{DW}} > 0$, domain walls do not proliferate spontaneously, ruling out intermediate modulated or incommensurate stripe phases.
