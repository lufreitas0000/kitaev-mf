---
id: M1.10-T1
title: "Exact Pure-Kitaev Spin Perturbation Theory, Normalization Map, and Mean-Field Reconciliation"
status: Completed
date: "2026-11-06"
questions:
  - "What is the exact operator-level mapping of the third-order Zeeman perturbation onto the zero-flux Majorana subspace?"
  - "How does the 1/4 Bloch generator normalization and the self-consistent KKT linear response reconcile the mean-field Dirac mass gap with third-order spin perturbation theory?"
depends_on:
  - M1.1-T1
  - M1.4-T1
  - M1.6-T7
contracts:
  - K-010
  - K-019
  - Contract H6
artifacts:
  - doc/labbook/chapters/ch_m1_6_solvers_and_basins.tex
  - doc/manuscript/sections/sec03_field_response_and_multistability.tex
  - test/unit/test_post_spinodal_kkt_audit.py
---

# M1.10-T1: Exact Pure-Kitaev Spin Perturbation Theory, Normalization Map, and Mean-Field Reconciliation

## 1. Hamiltonian & Normalization Conventions

We consider the ferromagnetic isotropic Kitaev honeycomb model:
$$H = H_K + V,$$
with:
$$H_K = -K \sum_{\langle i j \rangle_\gamma} \sigma_i^\gamma \sigma_j^\gamma, \qquad K > 0,$$
and Zeeman field perturbation:
$$V = -\sum_i \sum_{\gamma=x,y,z} h_\gamma \sigma_i^\gamma.$$

### Pauli vs. Spin-$1/2$ Mappings
If another literature convention uses canonical spin-$1/2$ operators $S_i^\gamma = \frac{1}{2} \sigma_i^\gamma$:
$$H_{\mathrm{spin}} = -J \sum_{\langle i j \rangle_\gamma} S_i^\gamma S_j^\gamma - \sum_i \mathbf{h}_S \cdot \mathbf{S}_i,$$
then:
$$J S_i^\gamma S_j^\gamma = \frac{J}{4} \sigma_i^\gamma \sigma_j^\gamma \implies K = \frac{J}{4},$$
$$\mathbf{h}_S \cdot \mathbf{S}_i = \frac{1}{2} \mathbf{h}_S \cdot \boldsymbol{\sigma}_i \implies \mathbf{h} = \frac{1}{2} \mathbf{h}_S.$$

---

## 2. Third-Order Perturbation Theory in Spin Space

Let $P$ project onto the unperturbed zero-flux physical ground-state subspace $P_{\mathrm{ff}}$, with ground-state energy $E_0$. Let $Q = 1 - P$.

Since a single spin operator $\sigma_i^\gamma$ flips the two plaquette fluxes sharing the $\gamma$-bond at site $i$, the expectation value of any single spin in the flux-free sector vanishes identically:
$$P V P \equiv 0.$$

Therefore, in Rayleigh-Schrödinger / Feshbach perturbation theory, all linear and quadratic terms vanish ($H_{\mathrm{eff}}^{(1)} = 0, H_{\mathrm{eff}}^{(2)} = 0$). The lowest non-vanishing effective Hamiltonian in the zero-flux sector appears at third order:
$$H_{\mathrm{eff}}^{(3)} = P V \frac{Q}{E_0 - H_K} V \frac{Q}{E_0 - H_K} V P.$$

### Non-Vanishing Intermediate Operator Paths
To return to the zero-flux sector ($P \sigma_j^\alpha \sigma_k^\gamma \sigma_l^\beta P \ne 0$), the three spin operations must form an elementary closed flux-annihilation loop.

On a triad of three sites $\langle\langle j, k, l \rangle\rangle$ where site $k$ is connected to $j$ via an $\alpha$-bond and to $l$ via a $\beta$-bond, the three operators $\sigma_j^\alpha$, $\sigma_k^\gamma$, and $\sigma_l^\beta$ (with $\alpha \ne \beta \ne \gamma$ being the three bond types meeting at $k$) flip each of the three plaquettes meeting at vertex $k$ exactly twice ($(-1)^2 = +1$).

In terms of the Kitaev four-Majorana representation ($\sigma_i^\mu = \mathrm{i} b_i^\mu c_i$ with local gauge constraints $\mathcal{D}_i = -c_i^0 c_i^x c_i^y c_i^z = b_i^x b_i^y b_i^z c_i = +1$):
$$b_k^\gamma = \epsilon^{\alpha\beta\gamma} b_k^\alpha b_k^\beta c_k.$$
Substituting into the three-spin operator:
$$\sigma_j^\alpha \sigma_k^\gamma \sigma_l^\beta = (\mathrm{i} b_j^\alpha c_j) (\mathrm{i} b_k^\gamma c_k) (\mathrm{i} b_l^\beta c_l) = \epsilon^{\alpha\beta\gamma} (\mathrm{i} b_j^\alpha b_k^\alpha)(\mathrm{i} b_k^\beta b_l^\beta) (\mathrm{i} c_j c_l).$$

In the zero-flux sector where the static gauge link variables are $u_{jk} = \mathrm{i} b_j^\alpha b_k^\alpha = +1$ and $u_{kl} = \mathrm{i} b_k^\beta b_l^\beta = +1$:
$$\mathbf{P (\sigma_j^\alpha \sigma_k^\gamma \sigma_l^\beta) P = \epsilon^{\alpha\beta\gamma} (\mathrm{i} c_j c_l)}.$$

The three-spin interaction projects **exactly to a next-nearest-neighbor matter Majorana hopping** $\mathrm{i} c_j c_l$ between sites on the same sublattice.

---

## 3. The Chiral Majorana Dirac Mass Gap

Summing over the time orderings and field components yields the effective Hamiltonian:
$$H_{\mathrm{eff}}^{(3)} = -\mathrm{i} \kappa \sum_{\langle\langle j, l \rangle\rangle} \nu_{jl} c_j c_l,$$
where $\nu_{jl} = \pm 1$ is the chiral orientation along the next-nearest-neighbor bond.

In Kitaev's static single-mode approximation (replacing the full two-vison + matter Feshbach resolvent $Q(E_0 - H_K)^{-1}Q$ by a single static energy denominator $\Delta E$), for magnetic field along $[111]$ ($\mathbf{h} = \frac{h}{\sqrt{3}}(1, 1, 1)^T$), the effective chiral coupling takes the form:
$$\kappa_{\mathrm{static}} = \frac{h_x h_y h_z}{\Delta E^2} = \frac{h^3}{3\sqrt{3} \Delta E^2}.$$

At the Dirac point $\mathbf{K} = (1/3, 2/3)$, the unperturbed nearest-neighbor dispersion vanishes:
$$f(\mathbf{K}) = 1 + \mathrm{e}^{-\mathrm{i} 2\pi/3} + \mathrm{e}^{\mathrm{i} 2\pi/3} = 0.$$
The next-nearest-neighbor chiral structure factor evaluates to:
$$g(\mathbf{K}) = \sum_{m=1}^6 \nu_m \sin(\mathbf{K} \cdot \mathbf{R}_m) = 3\sqrt{3}.$$

The canonical mass matrix acting on the two-sublattice spinor $\psi = (c_A, c_B)^T$ is:
$$\mathcal{M}(\mathbf{K}) = \begin{pmatrix} 2 \kappa g(\mathbf{K}) & 0 \\ 0 & -2 \kappa g(\mathbf{K}) \end{pmatrix} = \begin{pmatrix} 6\sqrt{3} \kappa & 0 \\ 0 & -6\sqrt{3} \kappa \end{pmatrix}.$$

The physical single-particle excitation gap is:
$$\mathbf{\Delta_{\mathrm{canonical}} = 2 \times (6\sqrt{3} \kappa) = 12\sqrt{3} \kappa}.$$

---

## 4. Reconciling Spin Perturbation Theory with 4-Flavor Mean-Field Theory

In the production 4-Majorana HFB solver, numerical evaluation across $h \in [0.005, 0.150]$ yields $\Delta_K^{\mathrm{MF}}(h) = (2.5921 \pm 0.0034) (h/|K|)^3$. Comparing this result with third-order spin perturbation theory involves three exact structural elements:

### Element A: Bloch Matrix Normalization ($1/4$)
In `src/kitaev_mf/services/hamiltonian.py`, the single-particle Bloch Hamiltonian is defined as:
$$H_{\mathrm{Bloch}}(\mathbf{k}) = \frac{\mathrm{i}}{4} A(\mathbf{k}).$$
Its eigenvalues $w$ are related to the canonical excitation frequencies $\varepsilon_{\mathrm{phys}}$ (eigenvalues of $\mathrm{i} A$) by Proposition 3.1 (`prop:m1_2_spectral_symmetries`):
$$w = \frac{1}{4} \varepsilon_{\mathrm{phys}} \implies \Delta_{\mathrm{canonical}}^{\mathrm{MF}}(h) = 4 \times \Delta_K^{\mathrm{MF}}(h) = 4 \times 2.5921 \left(\frac{h}{|K|}\right)^3 = \mathbf{10.3684} \, \frac{h^3}{|K|^2}.$$

### Element B: Intermediate Two-Vison Continuum vs. Flat Gauge-Band Energy Scales
1. **Spin Perturbation Theory (Two-Vison + Matter Continuum):** Acting with $\sigma_i^\gamma = \mathrm{i} c_i^0 c_i^\gamma$ on the physical ground state flips a $\mathbb{Z}_2$ bond link and excites a non-local vison pair together with dispersive matter Majoranas. Using the isolated two-vison threshold $\Delta_{\mathrm{flux}} \approx 0.264 |K|$ in the static single-mode approximation yields the reference scale:
   $$\Delta_{\mathrm{static}} = \frac{4}{\Delta_{\mathrm{flux}}^2} h^3 \approx 57.40 \, \frac{h^3}{|K|^2},$$
   whereas the exact continuum-integrated Feshbach value $\kappa_{\mathrm{exact}}$ is evaluated explicitly over the two-vison + matter spectrum in Track 4 (`T12`--`T15`).
2. **Unprojected 4-Flavor HFB Mean-Field Theory (Localized Gauge Bands):** In Hartree-Fock-Bogoliubov factorization, background $\mathbb{Z}_2$ fluxes remain uniform, and the intermediate single-particle excitations are localized flat-band gauge Majoranas $c_i^\alpha$ with two-particle excitation energy:
   $$2\Delta_{\mathrm{gauge}} = 2 |K| \bar{U} \approx 1.0497315 |K|.$$
   The bare squared-denominator ratio between these two intermediate scales is:
   $$\left(\frac{2\Delta_{\mathrm{gauge}}}{\Delta_{\mathrm{flux}}}\right)^2 = \left(\frac{1.0497315}{0.2640}\right)^2 = 15.8106.$$

### Element C: Self-Consistent KKT Multiplier and Fock Vertex Dressing (`Contract H6`)
In 4-flavor HFB theory, the bare Zeeman operator $-\sum_\gamma h_\gamma N^\gamma$ only couples the dispersive matter flavor $c^0$ to the gauge flavors $c^\gamma$ ($N^\gamma_{0\gamma} = -1$) and contains no direct $c^\alpha \leftrightarrow c^\beta$ matrix elements. Therefore, one cannot simply substitute $\Delta E \to 2\Delta_{\mathrm{gauge}}$ into the spin-perturbation formula $\frac{4}{\Delta E^2} h^3$. Instead, the $\mathcal{O}(h^3)$ chiral Dirac mass in HFB theory is generated by the **self-consistent linear response** at $\mathcal{O}(h)$:
- The external field induces a linear KKT Lagrange multiplier counter-field $\lambda_s^\gamma = \mathcal{O}(h)$ (`K-009`), which couples to the dual gauge generators $G^\gamma$ ($G^\gamma_{\alpha\beta} = -\epsilon_{\gamma\alpha\beta}$) that directly mix distinct gauge flavors $c^\alpha \leftrightarrow c^\beta$ at the intermediate vertex.
- Simultaneously, the field induces off-diagonal Fock bond expectations $\delta U_b^{0\gamma}, \delta U_b^{\alpha\beta} = \mathcal{O}(h)$.

Expanding the KKT stationarity condition $[\mathcal{H}(\mathbf{k}), \Pi^-(\mathbf{k})] = 0$ analytically to $\mathcal{O}(h^3)$ to derive the closed-form prefactor $c_{\mathrm{MF}} = 2.59214\dots$ is the formal target of `Contract H6` (Task `T08`).
