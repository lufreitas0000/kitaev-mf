# Post-Reorganization Code Salvage Queue (`kitaev-mf`)

During the 2026-11-18 linear consolidation, all historical branches and auxiliary worktrees were backed up offline in `/home/lu/Projects/archive/kitaev-mf_backups/pre-reconciliation-20261118.bundle` (4.83 MiB, 45 refs) so `kitaev-mf` maintains **one linear `main` branch**.
Non-`main` software modules are queued below for bounded, reviewed `feat/*` PRs rather than bulk merging unreviewed prototypes.

| Track | Source Ref in Offline Bundle | Target `kitaev-mf` Scope | Prerequisite Gate | Status |
| :--- | :---: | :--- | :--- | :---: |
| **Track 1 — Pure `scMFT` Core Kernels** | `39e733e` / `0b75fdf` | `InhomogeneousMeanFieldState` `(L^2,3,4,4)`, Kenney–Laub Padé `[3/2]`, Lehmann dual Hessian, 12-Majorana Parlett–Reid Pfaffian, `canonical_skew` null-block parser | `make verify` + `W4` independent review | **Completed (Stage 4: `feat/scmft-core-kernels`)** |
| **Track 2 — Exact Perturbation & Feshbach** | `39e733e` | `exact_matter.py`, `exact_perturbation.py`, `feshbach.py` | Resolve `W001-FND-02..05` (`K-019-A` erratum, $Y$-Zeeman sign, full-$Q$ resolvent) | **Deferred (`Issue #10`, `Issue #12`)** |
| **Track 3 — Phase 2 Static `scMFT` Sweeps** | `8b242a5` | `magnetization_torque.py`, `phase_diagram.py`, `spinodal_scaling.py`, `thermodynamics.py`, `vortex_potential.py`, `mzm_localization.py`, `chiral_vortex_transport.py` | `FND-STATE`, `FND-UNITS`, `FND-GAUSSIAN`, Track 1 merged | **Deferred (`M2.1`–`M2.4`)** |
| **Track 4 — Phase 3 Dynamical `scMFT` & STS** | `a152ea5` | `spectral_engine.py`, `bogoliubov.py`, `mean_field_spectral.py`, `mzm_spectroscopy.py`, `edge_and_scaling.py` | Phase 2 static baseline + parity/gauge conventions locked | **Deferred (`M3.0`–`M3.5`)** |
| **Track 5 — Track-D Real-Space Screening Visualizer** | `1d82e6a` | `fixed_point_visualizer.py` (`(U_b, W_p)` honeycomb lattice snapshots under bond flips) | Track 1 (`InhomogeneousMeanFieldState`) + `W4` check | **Deferred (Optional diagnostic tool)** |

## Delegated Out of `kitaev-mf` (`ADR-003`)

- **Spin-Space Exact Diagonalization (`exact_diagonalization.py`) & 2-RDM SDP (`variational_bounds.py`):** Exported to `/home/lu/Projects/archive/kitaev-numerics-handoff/` as unreviewed reference prototypes for `lufreitas0000/kitaev-numerics`. Future `kitaev-mf` comparisons consume certified `.dat` + `.meta.json` outputs through `data/external/kitaev_numerics/`.
