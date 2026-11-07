---
id: T14
title: "Feshbach Projectors, Resolvent Domain, and Eliminated-Block Operator Bounds"
status: completed
date: 2026-11-07
questions:
  - "How are the exact physical Feshbach projectors P_ff and Q_flux defined on the physical spin Hilbert space and Majorana sector space?"
  - "What are the exact physical flux thresholds E_0, E_Q^min, and Delta_flux on finite geometries?"
  - "What is the analytical proof that external low-energy windows E <= E_0 + Lambda lie in the resolvent set rho(Q H_K Q) with bounded operator norm?"
  - "Why is the global operator norm condition ||Q V Q|| < Delta_flux fallacious in the thermodynamic limit, and how does linked-cluster locality restore perturbative control?"
  - "How do the resolvent bounds scale across finite honeycomb tori L in {4, 6, 8, 10, 12}?"
depends_on:
  - "doc/rnd/exact-perturbation-theory/t12_matter_parity_selection_rule.md"
  - "doc/rnd/exact-perturbation-theory/t13_exact_matter_solver_and_vison_extrapolation.md"
  - "doc/rnd/exact-perturbation-theory/notes/03_low_energy_projector_and_feshbach.md"
contracts:
  - "B-030"
  - "B-031"
  - "B-032"
  - "B-033"
artifacts:
  - "src/kitaev_mf/services/feshbach.py"
  - "test/unit/test_feshbach_bounds.py"
---

# T14: Feshbach Projectors, Resolvent Domain, and Eliminated-Block Operator Bounds

## 1. Physical Projectors and Flux Thresholds (B-030, B-031)

### 1.1 Projector Decomposition
Let $\mathcal{H}_{\text{phys}}$ be the physical spin Hilbert space of dimension $2^{2N}$. The mutually orthogonal, complete projectors onto the flux-free and vortex sectors are defined by
$$
P_{\text{ff}} = \prod_p \frac{1 + W_p}{2}, \qquad Q_{\text{flux}} = \mathbf{I} - P_{\text{ff}},
$$
satisfying
$$
P_{\text{ff}}^2 = P_{\text{ff}}, \quad Q_{\text{flux}}^2 = Q_{\text{flux}}, \quad P_{\text{ff}} Q_{\text{flux}} = 0, \quad P_{\text{ff}} + Q_{\text{flux}} = \mathbf{I}.
$$
Because the Kitaev Hamiltonian strictly commutes with all hexagonal plaquette operators, $[H_K, W_p] = 0$, the unperturbed Hamiltonian is strictly block-diagonal:
$$
P_{\text{ff}} H_K Q_{\text{flux}} = 0, \qquad Q_{\text{flux}} H_K P_{\text{ff}} = 0.
$$

### 1.2 Spectral Edges and Physical Flux Gap
Define the physical flux-free ground-state energy:
$$
E_0 = \inf \sigma(P_{\text{ff}} H_K P_{\text{ff}}),
$$
and the bottom of the eliminated vortex block:
$$
E_Q^{\min} = \inf \sigma(Q_{\text{flux}} H_K Q_{\text{flux}}).
$$
The exact physical flux excitation gap is
$$
\Delta_{\text{flux}} = E_Q^{\min} - E_0 > 0.
$$
On finite honeycomb tori with physical parity projection across all topological boundary twist sectors $\boldsymbol{\tau} \in \{\pm 1\}^2$, $\Delta_{\text{flux}}(L) \in [0.179, 0.315] \, K$ for $L \in \{4, 6, 8, 10, 12\}$, extrapolating to $\Delta_{\text{flux}}^\infty \approx 0.300 K$ in the infinite-separation limit.

---

## 2. Resolvent Domain and Operator Norm Bounds (B-032)

### 2.1 Resolvent Set Invariant
Consider an external energy window $E \le E_0 + \Lambda$ with cutoff $\Lambda < \Delta_{\text{flux}}$.
The spectral distance from $E$ to the spectrum of the eliminated self-adjoint block $Q_{\text{flux}} H_K Q_{\text{flux}}$ is:
$$
\operatorname{dist}(E, \sigma(Q_{\text{flux}} H_K Q_{\text{flux}})) = \inf_{\lambda \in \sigma(Q H_K Q)} |\lambda - E| = E_Q^{\min} - E \ge \Delta_{\text{flux}} - \Lambda = \delta > 0.
$$
Consequently:
1. $E \in \rho(Q_{\text{flux}} H_K Q_{\text{flux}})$ (the resolvent set) unconditionally for all $E \le E_0 + \Lambda$.
2. The eliminated-block resolvent $R_Q(E) = (E - Q_{\text{flux}} H_K Q_{\text{flux}})^{-1}$ is a bounded, self-adjoint operator on the $Q$-subspace.
3. Its exact operator norm equals the spectral radius:
$$
\|R_Q(E)\| = \frac{1}{\operatorname{dist}(E, \sigma(Q_{\text{flux}} H_K Q_{\text{flux}}))} = \frac{1}{E_Q^{\min} - E} \le \frac{1}{\delta}.
$$
At the ground-state energy $E = E_0$ ($\Lambda = 0$), the bound is exact and tight:
$$
\|R_Q(E_0)\| = \frac{1}{\Delta_{\text{flux}}}.
$$

---

## 3. Demarcation of the Extensive Operator-Norm Fallacy (B-033)

### 3.1 The Fallacy of Global $\|Q V Q\|$ Smallness
In naive text-book presentations of Schrieffer–Wolff perturbation theory, one asserts that the perturbed resolvent
$$
R_Q(E, V) = (E - Q(H_K + V)Q)^{-1} = R_Q(E) \sum_{n=0}^\infty \left[ Q V Q R_Q(E) \right]^n
$$
converges if the operator norm satisfies the smallness condition:
$$
\|Q V Q R_Q(E)\| \le \|V\| \|R_Q(E)\| < 1 \iff \|V\| < \Delta_{\text{flux}}.
$$
However, for a uniform Zeeman magnetic field $\mathbf{h}$, the perturbation is extensive:
$$
V = -\sum_{j=1}^{2N} \mathbf{h} \cdot \boldsymbol{\sigma}_j \implies \|V\| = 2N \|\mathbf{h}\|.
$$
In the thermodynamic limit $N \to \infty$, $\|V\| \to \infty$ for any arbitrarily small $h > 0$. Therefore:
$$
\lim_{N\to\infty} \frac{\|V\|}{\Delta_{\text{flux}}} = \infty.
$$
The global operator-norm smallness condition is **mathematically impossible** to satisfy for an extensive many-body system.

### 3.2 Resolution via Linked-Cluster Locality
The perturbative series for low-energy effective operators $H_{\text{eff}}$ is valid not because of extensive operator-norm bounds, but because of **quantum locality and linked-cluster expansion**:
1. External flux-free states $|\Psi_0\rangle$ have zero vortices.
2. A single spin flip $\sigma_j^\alpha$ is strictly local, creating exactly two adjacent visons.
3. In third order, $P_{\text{ff}} V R_Q(E_0) V R_Q(E_0) V P_{\text{ff}}$, the intermediate states can only return to the zero-flux sector if the three flipped bonds form a closed loop (elementary triangle of the honeycomb lattice).
4. Each localized intermediate state has a finite local defect energy denominator:
   $$E_{\text{intermediate}} - E_0 \ge \Delta_{\text{flux}} > 0.$$
5. Consequently, the expansion parameter is the **intensive, local ratio**:
   $$\lambda_{\text{local}} = \frac{h}{\Delta_{\text{flux}}} \ll 1,$$
   which remains strictly finite and small ($h/\Delta_{\text{flux}} \approx 0.16 \ll 1$ at $h=0.05 K$) as $N \to \infty$.
