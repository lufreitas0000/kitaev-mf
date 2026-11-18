# Current sprint — RECON-20261116 completion & owner gate

Owning milestone: RECON-20261116
**Owning wave: RECON-20261116.** Active task status is owned by [PROJECT_STATE.yaml](PROJECT_STATE.yaml); this note explains current sprint scope.

## Reconciliation Status (Stages 0–5)

- **Git Topology & Family Boundaries ([ADR-003](doc/adr/ADR-003-linear-main-and-family-boundaries.md)):** Single `main` branch, single local worktree, offline pre-reconciliation bundle verified outside the repository, and clean 5-repo separation (`kitaev-mf` pure `scMFT` vs. `kitaev-numerics`, `kitaev-walk-geometry`, `ksl-sc-bilayer`, `Julia/KitaevDerivation.jl`) with `data/external/<source_repo>/` verified by `scripts/sync_external_data.py --verify`.
- **Merged Deliverables:** Wave 4 sampling/symmetry design (`PR #20`), R&D and W001 foundation absorption (`PR #22`), audited Track 1 pure-`scMFT` core kernels (`PR #23`, W4 reviewed), and PG-10 `native-v2` / LaTeX / CI-Pages overhaul (`governance/profile.toml`, `.github/workflows/docs.yml`, zero tracked PDFs).

## Active execution sprint on `main`: pending owner selection

Reconciliation changes no physical thresholds, hypotheses, or model scope. No Phase 2 physics campaign is running until Lu selects the next milestone and registers quantitative claims under [`doc/hypotheses/`](doc/hypotheses/README.md).

## Open Scientific Prerequisites (Human-Gated)

- Local Majorana representation and projector sign conventions: [K-019-A Issue #10](https://github.com/lufreitas0000/kitaev-mf/issues/10) (`m1_10_t1_k019a_second_order_erratum.md` resolves $PVP=0 \not\implies H_{\mathrm{eff}}^{(2)}=0$, leaving full-$Q$ resolvent and physical projection open).
- Exact finite-torus gauge/parity, flux-sector thresholds, and resolvent: [Issue #12](https://github.com/lufreitas0000/kitaev-mf/issues/12) (coordinate external resolvents with `kitaev-walk-geometry`).
- Analytical evaluation of the empirical `scMFT` cubic gap and linear counter-field regressions ($2.5921 \pm 0.0034$ and $\sqrt{5} h$): open [`IDEA-001`](BACKLOG.md).
- Independent finite-cluster spin-ED / 2-RDM SDP references (`H1`, `H3`): [validation Issue #14](https://github.com/lufreitas0000/kitaev-mf/issues/14) and [kitaev-numerics Issue #1](https://github.com/lufreitas0000/kitaev-mf/issues/1), imported strictly via `data/external/kitaev_numerics/`.
- Statistical design and geometry/KKT/vortex scoping: [Issue #7](https://github.com/lufreitas0000/kitaev-mf/issues/7), [Issue #8](https://github.com/lufreitas0000/kitaev-mf/issues/8), and [Issue #15](https://github.com/lufreitas0000/kitaev-mf/issues/15).

**Non-goals:** no unapproved Phase 2 solver campaign, no silent `HC-*.yaml` edits, and no promotion of empirical fits to analytical theorems.

## Next decision / stop conditions

Ask Lu to confirm whether the immediate post-reconciliation sprint targets (a) `IDEA-001` / `K-019-A` theoretical closure or (b) registering the first `doc/hypotheses/HC-01.yaml` contract for matched `scMFT` vs. `kitaev-numerics` finite-torus comparison. Longer-horizon items live in [BACKLOG.md](BACKLOG.md); obstructions and pivots live in [LOG_DISCOVERY_AND_PIVOTS.md](LOG_DISCOVERY_AND_PIVOTS.md).
