# Roadmap — Exact perturbation theory of the pure Kitaev model

## Objective

Develop a mathematically controlled perturbation theory for
$$H=H_K+V,\quad H_K=-K\sum_{\langle ij\rangle_\alpha}\sigma_i^\alpha\sigma_j^\alpha,\quad V=-\sum_{i,\alpha}h_\alpha\sigma_i^\alpha,$$
with $K>0$ and each bond counted once. Analytical derivations are primary; Python is an independent validation tool. The initial target is a controlled low-energy effective description through third order, with physical projection, intermediate flux sectors, and resolvent domains explicit. Do not assume thermodynamic convergence from finite-volume perturbation theory.

## Repository policy

- Working derivations remain Markdown under this directory; defer final LaTeX.
- Exploratory source and tests are added only after mathematical contracts are agreed.
- Generated arrays, logs, plots, caches, and build products go under untracked \`build/rnd/exact-perturbation-theory/\`.
- Keep this track isolated from production mean-field source.
- Reuse production interfaces only after an explicit compatibility review.

## Phase I — Conventions and unperturbed baseline

| ID | Deliverable | Exit criterion |
|---|---|---|
| R0.1 | Hamiltonian normalization map | Derive relation between Pauli-spin K, spin-1/2 J, and main-branch Bloch convention from operator definitions and bond counting |
| R0.2 | Flux-free spectrum | Bloch and finite-torus singular-value spectra agree under the derived map |
| R0.3 | Vacuum energy and Fock width | Independent converged quadrature and finite-size sums validate analytic values |
| R0.4 | Geometry/gauge contract | Bond orientation, plaquette ordering, link flips, Wilson loops, and twists are fixed |

Existing M0.1–M0.3 notes use $H_K=-K\sum\sigma_i^\alpha\sigma_j^\alpha$ and give $\epsilon_{\max}=6K$, $e_{\rm vac}\simeq-1.5746K$ per two-site cell, and width density $-2e_{\rm vac}$. These are analytical notes, not independently validated project outputs.

**Normalization audit:** the coupling map and real-space bond coefficient are derived in `notes/02_normalization_and_bloch_conventions.md`. The production Bloch matrix is scaled by one quarter relative to the canonical matter excitation generator because its builder returns $iA(\mathbf k)/4$. The production solver uses the negative-energy spectral projector, so this positive scale does not change the ideal occupied subspace or its fixed-point map; the numerical zero-mode threshold and public matrix contract remain to be audited under B-001. Do not interpret the factor-four eigenvalue difference alone as a physical disagreement.

## Phase II — Physical projection and parity

| ID | Deliverable | Exit criterion |
|---|---|---|
| R1.1 | Product of local constraints | Derive finite-torus projection from $D_i=+1$ with fixed site/bond ordering |
| R1.2 | Gauge-dependent matter parity | Required parity for each gauge representative, twist, and Wilson sector |
| R1.3 | Gaussian vacuum parity | Pfaffian/covariance prescription with fixed signs and zero-mode treatment |
| R1.4 | Spin-space ED | Small-cluster physical energies and flux labels agree with projected Majorana results |

Physical flux thresholds cannot be claimed from Gaussian vacuum energies alone. Complete this phase before interpreting vortex energies physically.

## Phase III — Flux-sector energetics

| ID | Deliverable | Exit criterion |
|---|---|---|
| R2.1 | Adjacent vortex pair | One-link flip, two adjacent plaquette fluxes, Gaussian and projected energies reported separately |
| R2.2 | Separated pair | Controlled endpoint separations; string-path gauge invariance |
| R2.3 | Finite-size study | Multiple shapes, sizes, modulo-three classes, twists, and Wilson sectors |
| R2.4 | Thermodynamic extrapolation | Controlled estimates and uncertainty |

Kitaev (2006), Appendix A, values are approximate thermodynamic benchmarks, not exact finite-size targets: $E_v\simeq0.1536K$, NN pair $\simeq0.2672K$, NNN pair $\simeq0.2372K$, and far-pair limit $\simeq0.3072K$.

## Phase IV — Projectors and resolvent control

| ID | Deliverable | Exit criterion |
|---|---|---|
| R3.1 | Physical projectors | Define $P_{\rm ff}$ and $Q_{\rm flux}=1-P_{\rm ff}$, including Wilson-loop policy |
| R3.2 | Physical thresholds | Compute $E_0$, $E_Q^{\min}$, and $\Delta_{\rm flux}$ on identical finite systems |
| R3.3 | Feshbach domain | State exact resolvent-set condition and prove only justified finite-size margins |
| R3.4 | Perturbed eliminated block | Analyze $Q_{\rm flux}VQ_{\rm flux}$ without assuming an extensive global norm is small |
| R3.5 | Matter-window limitation | Show what $P_\Lambda$ does and does not isolate; no unjustified SW reduction across gapless matter states |

## Phase V — Perturbation theory

| ID | Deliverable | Exit criterion |
|---|---|---|
| R4.1 | Spin-flip selection rules | Derive flux changes caused by each $\sigma_i^\alpha$ |
| R4.2 | Perturbative conventions | Fix projection, energy reference, degeneracy, and resolvent conventions |
| R4.3 | Lowest nonzero orders | Derive allowed terms and prove vanishing terms by selection rules |
| R4.4 | Third-order interaction | Derive all allowed paths and anisotropic field dependence; do not import coefficients |
| R4.5 | Error and validity | State finite-volume error estimates and limits of thermodynamic conclusions |

## Phase VI — Independent numerical validation

Implement incrementally after Phase I contracts are fixed: periodic geometry; link configurations and flux observables; bipartite T and antisymmetric A; singular values and Gaussian energy; gauge-equivalence and Bloch comparisons; parity/Pfaffian; small-cluster spin ED; flux-pair finite-size sweeps. Start with a transparent NumPy/SciPy kernel independent of the mean-field backend. Consider an adapter to the production LinalgBackendPort only after convention reconciliation.

## Completion criteria

Final exposition begins only after conventions are reconciled, parity is derived and checked, physical flux energies are reproducible, the Feshbach domain is justified, perturbative coefficients include all allowed intermediate sectors, and numerical tests/provenance are complete.