# Finite-torus conventions

## Scope and lattice

Use a periodic honeycomb torus with $N_1\times N_2$ primitive cells, where $N_1,N_2\geq 2$. Cell coordinates are $r=(m,n)\in\mathbb Z_{N_1}\times\mathbb Z_{N_2}$, and all coordinate arithmetic is modulo $(N_1,N_2)$. Each cell contains $A_r$ and $B_r$.

The three oriented A-to-B links are

- $z_r: A_r\to B_r$,
- $x_r: A_r\to B_{r+n_1}$,
- $y_r: A_r\to B_{r+n_2}$,

with $n_1=(1,0)$ and $n_2=(0,1)$. Store one real link variable per bond,
$u^\alpha_r\in\{\pm1\}$, always in this canonical A-to-B convention. In particular, traversing a link in the reverse direction does not replace its stored variable by $-u^\alpha_r$; the stored variable labels the unoriented bond with a fixed sublattice orientation.

The vortex-free reference representative is $u^x_r=u^y_r=u^z_r=+1$.

## Plaquette boundary and flux

A six-edge plaquette based at $A_r$ has the cycle

$$
A_r\xrightarrow{x}B_{r+n_1}
\xrightarrow{y}A_{r+n_1-n_2}
\xrightarrow{z}B_{r+n_1-n_2}
\xrightarrow{x}A_{r-n_2}
\xrightarrow{y}B_r
\xrightarrow{z}A_r.
$$

The arrows label the bond type; reverse traversals still use the canonical A-to-B link variable. The corresponding Majorana gauge flux is defined by the product of the six stored link variables,

$$
W_r
=
u^x_r\,
u^y_{r+n_1-n_2}\,
u^z_{r+n_1-n_2}\,
u^x_{r-n_2}\,
u^y_{r-n_2}\,
u^z_r.
$$

This convention gives $W_r=+1$ in the all-positive reference gauge. It is gauge invariant under
$u^\alpha_r\mapsto \eta_A(r)u^\alpha_r\eta_B(r+\delta_\alpha)$, where
$\eta_A,\eta_B\in\{\pm1\}$ and $\delta_z=0$, $\delta_x=n_1$, $\delta_y=n_2$.

Every link belongs to exactly two plaquettes. Consequently,

$$
\prod_{r}W_r=+1,
$$

so the number of plaquettes with $W_r=-1$ is even. Flipping one link changes exactly the two plaquette fluxes incident on that link. These are exact finite-graph identities, independent of the matter spectrum.

The symbol $W_r$ here denotes the Majorana link-flux product. Its precise operator relation, including any convention-dependent overall sign, to the spin-space Kitaev plaquette operator must be derived explicitly before comparing with spin-space ED; do not silently identify the two definitions.

## Noncontractible Wilson loops

Choose the two primitive noncontractible cycles represented by alternating $z$ and reverse-$x$ links, and by alternating $z$ and reverse-$y$ links. At fixed transverse coordinate, define

$$
\mathcal W_1(n)
=
\prod_{m=0}^{N_1-1}
u^z_{(m,n)}u^x_{(m-1,n)},
\qquad
\mathcal W_2(m)
=
\prod_{n=0}^{N_2-1}
u^z_{(m,n)}u^y_{(m,n-1)}.
$$

All indices are periodic. Each path closes after one winding and each vertex on the path occurs twice in the product, so these quantities are gauge invariant. Loops in the same homology class can be deformed into one another by multiplying the plaquette fluxes in the strip between them. In the flux-free sector, the loop value is therefore independent of the chosen transverse representative.

The plaquette fluxes subject to $\prod_rW_r=1$, together with two independent noncontractible Wilson-loop values, label gauge-equivalence classes of link configurations on the torus.

## Boundary twists

A convenient representative of the two global twists is obtained by changing signs on the two periodic seams:

$$
u^x_{(N_1-1,n)}\mapsto \tau_1u^x_{(N_1-1,n)},
\qquad
u^y_{(m,N_2-1)}\mapsto \tau_2u^y_{(m,N_2-1)},
\qquad
\tau_1,\tau_2\in\{\pm1\}.
$$

These seam changes leave every local plaquette flux unchanged and flip the corresponding primitive Wilson loop. They are link-variable twists; do not additionally impose matter-Majorana boundary phases unless an explicit mapping between the two conventions has been established. In a fixed gauge, the twists are represented by global holonomies, not by local vortices.

## Exact combinatorial checks required for B-004

The geometry implementation must test, for several tori including the smallest allowed sizes:

1. There are $2N_1N_2$ sites, $3N_1N_2$ links, and $N_1N_2$ plaquettes.
2. Every plaquette cycle has six distinct vertices and six distinct links.
3. Every link is incident on exactly two distinct plaquettes.
4. Every link appears exactly twice in the product of all plaquette fluxes.
5. The product of all plaquette fluxes is $+1$ for arbitrary link signs.
6. Flipping any one link changes precisely its two incident plaquette fluxes.
7. Arbitrary local gauge transformations preserve all plaquette fluxes and Wilson loops.
8. The two seam twists preserve all plaquette fluxes and flip the corresponding Wilson-loop values.
9. In a flux-free configuration, Wilson loops are invariant under changing the transverse representative within each homology class.

These are exact sign/counting tests; floating-point tolerances are not relevant. The later matter-matrix tests should separately use a scale-aware numerical tolerance.

## Remaining conventions to settle before parity and spin-space ED

- Fix a deterministic ordering of sites, bonds, plaquettes, and Majoranas. All Pfaffian and physical-parity signs depend on it.
- Derive the spin plaquette operator in the same site/component ordering and verify its relation to $W_r$ from $D_i=+1$.
- State explicitly whether a calculation fixes $(\mathcal W_1,\mathcal W_2)$ or sums over all four Wilson-loop sectors.
- Treat the physical projection and matter-fermion boundary-condition mapping as a separate derivation; the link-sector labels alone do not determine the allowed matter parity.
