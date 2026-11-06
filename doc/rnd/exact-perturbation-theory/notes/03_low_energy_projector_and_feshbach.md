# 03 — Low-energy external states and Feshbach reduction

## Projectors and thresholds

Let $P_{\mathrm{ff}}$ project onto all physical states with zero
plaquette flux. Unless a separate superselection restriction is explicitly
imposed, this includes all global Wilson-loop sectors. Define
$$
Q_{\mathrm{flux}}=1-P_{\mathrm{ff}}.
$$
Thus $Q_{\mathrm{flux}}$ contains physical states with nonzero plaquette
flux. If one instead restricts to a fixed global Wilson-loop sector,
the Hilbert space and the complementary projector must be redefined
consistently; other zero-plaquette-flux Wilson sectors must not silently
be counted as vortex sectors.

Define the physical flux-free threshold
$$
E_0=\inf\sigma(P_{\mathrm{ff}}H_KP_{\mathrm{ff}})
$$
and the bottom of the eliminated block
$$
E_Q^{\min}=\inf\sigma(Q_{\mathrm{flux}}H_KQ_{\mathrm{flux}}).
$$
The flux threshold is
$$
\Delta_{\mathrm{flux}}=E_Q^{\min}-E_0.
$$
These definitions refer to the same finite geometry, boundary conditions,
and physical projection. Gaussian matter-vacuum energies alone are not
sufficient if the vacuum has the wrong physical parity.

A low-energy external-state window is
$$
P_\Lambda=P_{\mathrm{ff}}\,
\mathbf 1_{[E_0,E_0+\Lambda]}(H_K).
$$
It selects external states inside the flux-free sector. It does not isolate
a spectral block from its full orthogonal complement: that complement
contains flux-free matter states immediately above the cutoff. In a
gapless thermodynamic matter sector, no uniform gap separates these
states. Therefore $P_\Lambda$ alone does not justify a global
Schrieffer–Wolff block diagonalization.

## Feshbach operator

For $H=H_K+V$, the exact Feshbach operator associated with
$P_{\mathrm{ff}}+Q_{\mathrm{flux}}=1$ is
$$
H_{\mathrm F}(E)=P_{\mathrm{ff}}HP_{\mathrm{ff}}
+P_{\mathrm{ff}}HQ_{\mathrm{flux}}
(E-Q_{\mathrm{flux}}HQ_{\mathrm{flux}})^{-1}
Q_{\mathrm{flux}}HP_{\mathrm{ff}}.
$$
It is defined only when
$$
E\in\rho(Q_{\mathrm{flux}}HQ_{\mathrm{flux}}),
$$
the resolvent set of the eliminated block. For a self-adjoint eliminated
block, a sufficient real-energy condition is
$$
E<\inf\sigma(Q_{\mathrm{flux}}HQ_{\mathrm{flux}}).
$$

For the unperturbed eliminated block, the lower edge is
$E_0+\Delta_{\mathrm{flux}}$. Hence external energies in
$$
E\in[E_0,E_0+\Lambda],\qquad \Lambda<\Delta_{\mathrm{flux}},
$$
lie below its spectrum, with margin
$$
\delta=\Delta_{\mathrm{flux}}-\Lambda>0.
$$
For real $E\leq E_0+\Lambda$,
$$
\left\|(E-Q_{\mathrm{flux}}H_KQ_{\mathrm{flux}})^{-1}\right\|
\leq\delta^{-1}.
$$
This bound assumes that $E_Q^{\min}$ is the true bottom of the chosen
physical eliminated block.

For nonzero $V$, the exact resolvent contains
$Q_{\mathrm{flux}}VQ_{\mathrm{flux}}$. The unperturbed margin does not
by itself prove that the perturbed block remains invertible. A sufficient
bounded-operator condition is
$$
\|Q_{\mathrm{flux}}VQ_{\mathrm{flux}}\|<\delta,
$$
which permits a Neumann-series estimate. For an extensive field
perturbation this global norm is generally extensive, so this condition
may be too strong in the thermodynamic limit. A sector-resolved,
local, or linked-cluster argument is then needed.

A single $\sigma_i^\alpha$ insertion changes the two plaquette fluxes
adjacent to the $\alpha$-bond. Higher perturbative orders can access
other vortex configurations; their thresholds must be included when
estimating intermediate-state denominators.

## Interpretation

The interval above controls the *unperturbed flux-sector resolvent* at
the chosen external energies. It does not establish a globally isolated
low-energy band of the full many-body Hamiltonian or prove convergence
of a thermodynamic Schrieffer–Wolff series. The Feshbach operator is
energy dependent; projected eigenvalues satisfy its corresponding
energy-dependent eigenvalue equation.
