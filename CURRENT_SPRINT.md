# Current sprint — M2.1 Wave 1 (`HC-01..03`) Execution Complete

Owning milestone: M2.1
**Owning wave: Wave 1 (`HC-01`–`HC-03`).** Active task status is owned by [PROJECT_STATE.yaml](PROJECT_STATE.yaml); this note explains current sprint scope.

## Completed Wave 1 Deliverables (`44692e5`)

- **External Data Ingestion (`ADR-003`):** Verified 3 active upstream datasets under `data/external/` (`kitaev_numerics` Lanczos ED + `Clarabel` Pauli SOS SDP, and `kitaev_walk_geometry` chiral vison $T$-matrix) alongside 4 quarantined legacy datasets.
- **`HC-01` Feshbach–Schur Cubic Gap (`K-024`):** Proved continuum $c_{\mathrm{MF}}^{\mathrm{cont}} = 10.641666$, $\|\boldsymbol{\lambda}^*\|/h \to 2.237814$, slope $2.999694$, and established that 100% of the $Q$-space $O(h)$ gap-opening vertex comes from $\lambda^{*,\gamma} G^\gamma$ (`delta_bare_zeeman <= 1.60e-15`, `delta_no_lambda <= 4.04e-15`).
- **`HC-02` Certified ED / SDP / Gutzwiller Sandwich (`K-025`):** Reconciled $\mathcal{D}_s = -Z_{2s}Z_{2s+1}=+1$ (`<3.84e-16` on $N=2$ Fock space), verified even-parity Zero-Variance Principle on $N=12$ (`-0.8166902250`, $e_{\mathrm{loc}}\text{-std}=3.08 \times 10^{-16}$) and $N=18$ (`-0.7939723679`), and confirmed `24/24` sandwich inequalities with 50-chain $\mathrm{BC}_a$ 95% CIs across $h \in [0, 0.50]$.
- **`HC-03` Multi-$L$ Vison Scaling & KWG Bridge (`K-026`):** Reconciled $\kappa_{\mathrm{static}} = \kappa_{\mathrm{KWG}}/2$ (`<4.4e-9` match to `kitaev-walk-geometry`), converged 94 chiral two-vison points across $L \in \{8, 10, 12, 15, 16, 18, 20\}$, and uncovered period-3 intervalley $\cos(\mathbf{K}\cdot\mathbf{R})$ modulation (`Delta_BIC` in `[2.39, 4.68]` for raw monotone fit vs crisp $\xi_M \approx \{3.03, 1.36, 0.76\}$ on fixed-phase $R=3 \to 6$).

## Open Scientific Prerequisites (Human-Gated)

- Present Wave 1 numerical results for PI `TODO(owner)` threshold locking on `doc/hypotheses/HC-01.yaml`–`HC-03.yaml` (including replacing the monotone null in `HC-03.yaml` with the oscillatory intervalley envelope or $R \bmod 3$ stratification).
- Analytical closed-form BZ evaluation of `IDEA-001` (`c_MF = 10.641666`, `sqrt(5)`) remains open.

## Next decision / stop conditions

Await PI lock on `HC-01..03.yaml` thresholds and authorization for Wave 2 finite-temperature thermodynamics ($C_v(T)/T$, $S(T)$ across $T_L, T_H$). Longer-horizon items live in [BACKLOG.md](BACKLOG.md) and [HORIZON.md](HORIZON.md); obstructions live in [LOG_DISCOVERY_AND_PIVOTS.md](LOG_DISCOVERY_AND_PIVOTS.md).
