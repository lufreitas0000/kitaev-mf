# ADR-001: Algebraic Gauge Projection, Multi-Sublattice Geometry, and PyTorch Array Engine

## Context
Majorana mean-field theory (MMFT) on the extended Kitaev model introduces local auxiliary gauge redundancies. In the thesis of Lucas Freitas (2024), enforcing the gauge constraint $\langle c_i^T G^\gamma c_i \rangle = 0$ required an empirical linear fit $\delta\bm{\lambda} = \alpha(\Gamma)\mathbf{h}$. Furthermore, mapping competing magnetic orders (zigzag and stripy) requires expanding beyond the standard 2-sublattice primitive cell.

## Decisions
1. **Algebraic Rigor:** Adopt the exact Clifford / CAR $C^*$-algebra identity $\mathbf{T}_i^2 = \frac{3}{8}(\mathbb{1} - D_i)$, proving that the physical Hilbert space coincides with the gauge singlet sector.
2. **Coupled Root-Finding Solver:** Abandon decoupled Picard updates for $\bm{\lambda}$. Formulate the problem as a unified min-max saddle-point system $\min_\Omega \max_{\bm{\lambda}} \mathcal{L}(\Omega, \bm{\lambda})$ solved via Joint Anderson Acceleration (DIIS) and Augmented Lagrangian (ALM) penalty methods.
3. **Decoupling Channel Parameter:** Parameterize the spin representation as $\mathrm{S}_i^\alpha(\eta) = \frac{i}{4} c_i^T (\mathrm{M}^\alpha - \eta \mathrm{G}^\alpha) c_i$, retaining $\eta \in [0, 1]$ as an explicit simulation parameter.
4. **Unit Cell Hierarchy:** Implement both the primitive 2-sublattice cell and the 4-sublattice rectangular cell ($\mathbf{A}_1 = \sqrt{3}d_{NN}\hat{\mathbf{x}}$, $\mathbf{A}_2 = 3d_{NN}\hat{\mathbf{y}}$).
5. **Numerical Stack:** Standardize on PyTorch (`torch.linalg`) running double precision (`complex128`), allowing seamless execution between local hardware (RTX 3050 6GB / 64GB DDR5) and Colab Pro GPUs (L4 22.5GB).

## Consequences
- Gauge constraint enforcement becomes automated and robust without empirical parameters.
- Momentum-space calculations can explore both fractionalized spin liquids and collinear zigzag orders.
- The codebase remains clean, strictly typed, and hardware-independent.
