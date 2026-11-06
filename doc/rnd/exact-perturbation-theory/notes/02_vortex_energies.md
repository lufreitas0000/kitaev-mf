# 02 — Vortex-pair energies and finite-size targets

For the isotropic convention $J_x=J_y=J_z=K$, Kitaev's Appendix A
reports an isolated-vortex energy of approximately
$$
E_v\simeq0.1536K.
$$
This is a thermodynamic estimate extrapolated from finite systems. A torus
requires an even number of vortices.

The appendix also reports attractive short-distance interactions. The
quoted binding-energy reductions are approximately $0.04K$ for
nearest-neighbour vortices and $0.07K$ for next-nearest-neighbour
vortices. These are reductions relative to $2E_v$, not total pair
energies. The corresponding rough targets are
$$
E_{\mathrm{pair}}^{\mathrm{NN}}\simeq
2(0.1536)-0.04=0.2672K,
$$
$$
E_{\mathrm{pair}}^{\mathrm{NNN}}\simeq
2(0.1536)-0.07=0.2372K,
$$
and
$$
E_{\mathrm{pair}}^{\mathrm{far}}\longrightarrow
2E_v\simeq0.3072K.
$$

These are approximate reference targets, not exact constants. Kitaev
notes appreciable finite-size effects in the pair-interaction estimates
and does not apply a finite-size correction. The numerical project must
compute finite-torus values directly, state the physical parity convention,
and perform size and boundary-twist checks before comparing extrapolations
with these targets.

## Required comparisons

1. Flux-free sector against a one-link-flipped gauge configuration.
2. Gauge strings whose endpoints realize several vortex separations.
3. Several torus shapes and linear sizes, including different values modulo
   three because the isotropic dispersion has Dirac points.
4. Different global Wilson-loop sectors or boundary twists.
5. Gaussian-vacuum energies and physical-sector energies with the correct
   matter parity, reported separately.
