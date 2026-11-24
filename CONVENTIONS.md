# CONVENTIONS.md — 4-Majorana Gaussian `scMFT` (`kitaev-mf`)

**Status:** `kitaev-mf` internal conventions (`Phase 1–2`, locked) + opt-in QSL `conventions-v2.0` cross-repository bridge (`lufreitas0000/kitaev-mf#25`, synchronized with `research/qsl/contracts/qsl_conventions.json`).

> **Rule for All Agents:** Read this file and [`doc/shared/NOTATION.md`](doc/shared/NOTATION.md) before modifying algebra, solvers, or cross-repository adapters. Never alter internal Majorana normalizations, signs, bond orientations, or `TODO(owner)` thresholds in `doc/hypotheses/HC-*.yaml` to make a test pass. Use [`src/kitaev_mf/domain/qsl_bridge.py`](src/kitaev_mf/domain/qsl_bridge.py) at every cross-QSL boundary.

---

## 1. Locked Internal `kitaev-mf` Conventions (Zero Breaking Changes)

| Quantity | Internal `kitaev-mf` Symbol / Code | Exact Definition |
| :--- | :--- | :--- |
| **Pauli & spin convention** | $\sigma_i^\alpha$ (`S_i^\alpha = \sigma_i^\alpha / 2`) | Eigenvalues $\pm 1$; pure FM Kitaev uses $K_{\mathrm{Pauli}} = -1.0$ in $H = +K_{\mathrm{Pauli}}\sum_{\langle ij\rangle_\alpha} \sigma_i^\alpha \sigma_j^\alpha - \frac{h}{\sqrt{3}}\sum_i(\sigma_i^x+\sigma_i^y+\sigma_i^z)$. |
| **Local Majorana basis** | $\mathbf{c}_i = (c_i^0, c_i^1, c_i^2, c_i^3)^T$ | Related to Kitaev's $(c_i, b_i^x, b_i^y, b_i^z)$ by the odd-parity map $(c_i^0, c_i^1, c_i^2, c_i^3) = (-c_i, b_i^x, b_i^y, b_i^z)$ with $\{\,c_i^\mu, c_j^\nu\,\} = 2\delta_{ij}\delta^{\mu\nu}\mathbb{1}$. |
| **Local parity constraint** | $\mathcal{D}_i$ | $\mathcal{D}_i = +c_i^0 c_i^1 c_i^2 c_i^3 = +b_i^x b_i^y b_i^z c_i = -(-1)^{n_{1,i}+n_{2,i}} = +1$ on the physical spin subspace. |
| **Legacy $\mathrm{SO}(4)$ matrices** | `G_MATRICES[a]`, `M_MATRICES[a]`, `N_MATRICES[a]` (`a in {1,2,3}`) | `constants.py`: $[M^a, M^b] = +2\epsilon_{abc}M^c$, $[G^a, G^b] = -2\epsilon_{abc}G^c$, $[M^a, G^b] = 0$, $N^a = \frac{1}{2}(M^a + G^a)$. |
| **Quadratic Hamiltonian & covariance** | $H = \frac{\mathrm{i}}{4}\mathbf{c}^T \mathbf{A} \mathbf{c}$, $\Omega_{ab} = \langle \mathrm{i} c_a c_b \rangle_{a\neq b}$ | Real skew-symmetric $\mathbf{A}^T = -\mathbf{A}$, $\Omega = \mathrm{i}\operatorname{sgn}(\mathrm{i}\mathbf{A})$, $\langle H \rangle = -\frac{1}{4}\operatorname{Tr}(\mathbf{A}\Omega)$. Bond matrix $U_b^{\mu\nu} = \langle \mathrm{i} c_i^\mu c_j^\nu \rangle$ ($i \in A \to j \in B$), on-site $V_s^{\mu\nu} = \langle \mathrm{i} c_s^\mu c_s^\nu \rangle$. |
| **Lagrange multiplier term** | `state.lambdas` ($\boldsymbol{\lambda}_s^{\mathrm{kmf}}$) | Enters on-site skew block as $B_s \supset +\sum_{a=1}^3 \lambda_s^{\mathrm{kmf},a} G^a$, i.e. $H_\lambda = +2\sum_{s,a} \lambda_s^{\mathrm{kmf},a} T_s^a$ where $T_s^a = \frac{\mathrm{i}}{8}\mathbf{c}_s^T G^a \mathbf{c}_s$. |
| **Energy normalization** | `EnergyBreakdown.per_site` | Intensive bulk energy strictly **per honeycomb site** $e = E / N_{\mathrm{sites}}$ ($N_{\mathrm{sites}} = 2 L^2$); localized vison energies extensive $V(R) = E_{2v}(R) - E_0$. |
| **Directed bonds & plaquettes** | $b = (i \in A \to j \in B, \gamma_b)$, $W_p$ | Primitive $A \to B$ offsets: $z\,(0,0)$, $x\,(0,-1)$, $y\,(-1,0)$. Plaquette Wilson loop $W_p = \prod_{b \in \partial p} u_b$ (and 12-Majorana Pfaffian $\mathcal{W}_p = \operatorname{Pf}(\Omega|_{\partial p})$). |

---

## 2. Opt-In Cross-QSL `conventions-v2.0` Bridge (`src/kitaev_mf/domain/qsl_bridge.py`)

To interoperate losslessly with `ksl-sc-bilayer` (`conventions-v2.0`), `kitaev-numerics`, and `kitaev-walk-geometry` without disturbing verified `kitaev-mf` Wave-1 physics (`HC-01..03`), [`qsl_bridge.py`](src/kitaev_mf/domain/qsl_bridge.py) implements the exact boundary adapter (`QSLConventionBridgePort` / `DefaultQSLConventionBridge`):

### 2.1 Universal $(\mathsf{L}^\alpha, \mathsf{R}^\alpha)$ Generators & $P_0$ Parity Involution

On the canonical `conventions-v2.0` quartet $\boldsymbol{\gamma}_j = (+c_j, b_j^x, b_j^y, b_j^z)^T = P_0 \mathbf{c}_j$ with $P_0 = \operatorname{diag}(-1, 1, 1, 1)$ (`P0_PARITY_BRIDGE`, $P_0^T = P_0 = P_0^{-1}$, $\det P_0 = -1$):

- **Universal Cartesian-keyed generators (`SO4_L`, `SO4_R`, strictly $\alpha \in \{\text{'x'},\text{'y'},\text{'z'}\}$ — numeric keys forbidden):**
  $$(\mathsf{L}^\alpha)_{0,\alpha} = -\frac{1}{2},\quad (\mathsf{L}^\alpha)_{\alpha,0} = +\frac{1}{2},\quad (\mathsf{L}^\alpha)_{\beta,\gamma} = -\frac{1}{2}\epsilon_{\alpha\beta\gamma},$$
  $$(\mathsf{R}^\alpha)_{0,\alpha} = +\frac{1}{2},\quad (\mathsf{R}^\alpha)_{\alpha,0} = -\frac{1}{2},\quad (\mathsf{R}^\alpha)_{\beta,\gamma} = -\frac{1}{2}\epsilon_{\alpha\beta\gamma},$$
  satisfying $[\mathsf{L}^\alpha, \mathsf{L}^\beta] = \epsilon_{\alpha\beta\gamma}\mathsf{L}^\gamma$, $[\mathsf{R}^\alpha, \mathsf{R}^\beta] = \epsilon_{\alpha\beta\gamma}\mathsf{R}^\gamma$, $[\mathsf{L}^\alpha, \mathsf{R}^\beta] = 0$.
- **Exact Left $\leftrightarrow$ Right Swap & Legacy Matrix Dictionary:**
  $$P_0 \mathsf{L}^\alpha P_0 = \mathsf{R}^\alpha, \qquad P_0 \mathsf{R}^\alpha P_0 = \mathsf{L}^\alpha,$$
  $$+\frac{1}{2} M^a = \mathsf{R}^\alpha \quad (\text{physical spin } S_j^\alpha \text{ in } (-c, \mathbf{b})), \qquad -\frac{1}{2} G^a = \mathsf{L}^\alpha \quad (\text{gauge } G_j^{\mathrm{v2},\alpha} = -T_j^a \text{ in } (-c, \mathbf{b})),$$
  $$N^a = \mathsf{R}^\alpha - \mathsf{L}^\alpha, \qquad P_0 N^a P_0 = \mathsf{L}^\alpha - \mathsf{R}^\alpha = -N^a.$$
- **Operator & Expectation Bridge (Exact on 16-State Two-Site Fock Space):**
  - Physical spin: $S_j^\alpha = \frac{\mathrm{i}}{4}\boldsymbol{\gamma}_j^T \mathsf{L}^\alpha \boldsymbol{\gamma}_j = \frac{\mathrm{i}}{4}\mathbf{c}_j^T \mathsf{R}^\alpha \mathbf{c}_j = P_+ \frac{\sigma_j^\alpha}{2} P_+$.
  - Standard-sign parton gauge: $G_j^{\mathrm{v2},\alpha} = \frac{\mathrm{i}}{4}\boldsymbol{\gamma}_j^T \mathsf{R}^\alpha \boldsymbol{\gamma}_j = \frac{\mathrm{i}}{4}\mathbf{c}_j^T \mathsf{L}^\alpha \mathbf{c}_j = -T_j^a = -P_- \frac{\sigma_j^\alpha}{2} P_-$, with $G_j^{\mathrm{v2},\alpha} P_+ = 0$.
  - Unprojected Pauli bilinear: $\sigma_j^\alpha = \mathrm{i} b_j^\alpha c_j = \mathrm{i} c_j^0 c_j^a = 2(S_j^\alpha - G_j^{\mathrm{v2},\alpha})$.
  - Lagrange multipliers: $\boldsymbol{\lambda}_s^{\mathrm{v2}} = -2\,\boldsymbol{\lambda}_s^{\mathrm{kmf}}$ (`lambdas_kmf_to_v2` / `lambdas_v2_to_kmf`), preserving $H_\lambda = +\boldsymbol{\lambda}_s^{\mathrm{v2}} \cdot \mathbf{G}_s^{\mathrm{v2}} = +2\boldsymbol{\lambda}_s^{\mathrm{kmf}} \cdot \mathbf{T}_s$ identically.
  - Covariance & bond blocks: $V_s^{\mathrm{v2}} = P_0 V_s^{\mathrm{kmf}} P_0$, $U_b^{\mathrm{v2}} = P_0 U_b^{\mathrm{kmf}} P_0$, $\Omega^{\mathrm{v2}} = \mathbf{P}_0^{\oplus N}\Omega^{\mathrm{kmf}}\mathbf{P}_0^{\oplus N}$ (diagonal $(0,0)$ and $(\alpha,\beta)$ blocks invariant; off-diagonal $(0,\alpha)$ and $(\alpha,0)$ matter–gauge blocks negate).

### 2.2 Exclusive Symbols: $W_p$ and Pure-Kitaev $(u_{ij}^\alpha, v_{ij}^\alpha)$

- **Pure-Kitaev directed bond scalars ($i \in A \to j \in B$):**
  $$u_{ij}^\alpha \equiv \langle \mathrm{i} b_i^\alpha b_j^\alpha \rangle = U_b^{\gamma_b,\gamma_b}, \qquad v_{ij}^\alpha \equiv \langle \mathrm{i} c_i c_j \rangle = U_b^{0,0}$$
  (at pure FM Kitaev $u = +1$, $v \approx -0.5248657458506$ identically in both `kitaev-mf` and `ksl-sc-bilayer`).
- **Plaquette Wilson loop ($W_p$, reserved exclusively for Kitaev plaquettes):**
  On the physical 6-site hexagon subspace ($\mathcal{D}_s = +1$), the three $A \to B$ and three $B \to A$ perimeter steps yield $(-\mathrm{i})^3(+\mathrm{i})^3 = +1$, giving the exact projected operator identity (verified on the $2^{12} = 4096$-state sparse Fock space):
  $$P_{\mathrm{phys}}\, W_p^{\mathrm{spin}}\, P_{\mathrm{phys}} = P_{\mathrm{phys}}\left(+\prod_{b \in \partial p} \mathrm{i} b_{A_b}^{\alpha_b} b_{B_b}^{\alpha_b}\right) P_{\mathrm{phys}}, \qquad W_p^{\mathrm{MF}} = \prod_{b \in \partial p} u_b.$$

### 2.3 Cross-QSL Coupling & Chiral $\kappa$ Port Conversions

| Source $\leftrightarrow$ Target | Exact Port Equation | Function in `qsl_bridge.py` |
| :--- | :--- | :--- |
| `kitaev-mf` / `kitaev-numerics` $\leftrightarrow$ `ksl-sc-bilayer` | $K_{\mathrm{Pauli}} = -K_{\mathrm{bilayer}}$ (`-1.0 <-> +1.0` FM) | `kitaev_coupling_pauli_to_bilayer` / `_bilayer_to_pauli` |
| `kitaev-walk-geometry` $\leftrightarrow$ `kitaev-mf` | $\kappa_{\mathrm{static}} = \frac{1}{2}\kappa_{\mathrm{KWG}}$, $\Delta_{\mathrm{can}}^{\mathrm{MF}} = 12\sqrt{3}|\kappa_{\mathrm{static}}| = 6\sqrt{3}|\kappa_{\mathrm{KWG}}| = \Delta_K^{\mathrm{KWG}}$ | `chiral_kappa_kwg_to_kmf_static` / `_kmf_static_to_kwg` |
| `kitaev-mf` $\leftrightarrow$ `ksl-sc-bilayer` | $\kappa_{\mathrm{KSLSC}} = \kappa_{\mathrm{static}}$ (full BdG gap $2\Delta_{\mathrm{KSL}} = \Delta_{\mathrm{can}}^{\mathrm{MF}} = 12\sqrt{3}|\kappa_{\mathrm{static}}|$) | `chiral_kappa_kmf_static_to_kslsc` / `chiral_full_gap_canonical` |

---

## 3. Open Owner Decisions (`TODO(owner)`)

- [ ] `TODO(owner)` (`lufreitas0000/kitaev-mf#25`): Decide whether a future major refactor of `kitaev-mf` should switch internal `(c^0, c^1, c^2, c^3)` to `(+c, b^x, b^y, b^z)` or permanently retain `(-c, b^x, b^y, b^z)` with `qsl_bridge.py` (`P_0`) as the canonical boundary adapter.
- [ ] `TODO(owner)`: Lock pre-registered numerical `support_if` / `refute_if` thresholds in `doc/hypotheses/HC-01.yaml`–`HC-03.yaml`.
