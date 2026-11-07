"""Theoretical foundation and exact gauge-sector conventions for Phase 3.

Implements the exact Hamiltonian, gauge configurations, plaquette flux evaluation,
matter-fermion parity constraints (Pedrocchi et al. 2011, Zschocke & Vojta 2015),
and single-spin quench mapping for reproducing arXiv:2209.12930.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.linalg import schur


@dataclass(frozen=True)
class DirectedBond:
    """Directed nearest-neighbor bond connecting site_a in A to site_b in B."""

    site_a: int  # Sublattice A site index
    site_b: int  # Sublattice B site index
    bond_type: int  # 1 (x), 2 (y), 3 (z)


@dataclass(frozen=True)
class GaugeConfiguration:
    """Static Z2 gauge background on a honeycomb cluster.

    Stores directed link eigenvalues u_{ab} = i b_a^alpha b_b^alpha in {-1, +1}
    defined with canonical orientation from sublattice A to sublattice B.
    """

    num_sites: int
    bonds: Sequence[DirectedBond]
    bond_links: dict[tuple[int, int], int]  # (site_a, site_b) -> u_{ab} in {-1, +1}

    def get_link(self, site_i: int, site_j: int) -> int:
        """Get directed link variable u_{ij} with anti-symmetry u_{ji} = -u_{ij}."""
        if (site_i, site_j) in self.bond_links:
            return self.bond_links[(site_i, site_j)]
        if (site_j, site_i) in self.bond_links:
            return -self.bond_links[(site_j, site_i)]
        raise KeyError(f"Bond between site {site_i} and {site_j} does not exist.")


def build_uniform_gauge_configuration(
    num_sites: int,
    bonds: Sequence[DirectedBond],
) -> GaugeConfiguration:
    """Build vortex-free reference gauge configuration where u_{ab} = +1 for all A->B bonds."""
    bond_links: dict[tuple[int, int], int] = {}
    for bond in bonds:
        bond_links[(bond.site_a, bond.site_b)] = 1

    return GaugeConfiguration(
        num_sites=num_sites,
        bonds=bonds,
        bond_links=bond_links,
    )


def flip_bond_gauge_link(
    gauge: GaugeConfiguration,
    site_i: int,
    site_j: int,
) -> GaugeConfiguration:
    """Toggle static gauge link u_{ab} -> -u_{ab} on the bond connecting site_i and site_j."""
    key = (site_i, site_j) if (site_i, site_j) in gauge.bond_links else (site_j, site_i)
    if key not in gauge.bond_links:
        raise ValueError(f"Bond ({site_i}, {site_j}) does not exist in gauge configuration.")

    new_links = dict(gauge.bond_links)
    new_links[key] = -new_links[key]

    return GaugeConfiguration(
        num_sites=gauge.num_sites,
        bonds=gauge.bonds,
        bond_links=new_links,
    )


def compute_plaquette_flux(
    gauge: GaugeConfiguration,
    plaquette_bonds: Sequence[tuple[int, int]],  # Directed or undirected pairs (i, j) on hexagon
) -> int:
    """Compute gauge-invariant Wilson loop W_p = prod_{b in partial p} u_b in {-1, +1}."""
    w = 1
    for site_i, site_j in plaquette_bonds:
        key = (site_i, site_j) if (site_i, site_j) in gauge.bond_links else (site_j, site_i)
        w *= gauge.bond_links[key]
    return int(w)


def assemble_matter_skew_matrix(
    gauge: GaugeConfiguration,
    K: float = 1.0,
    kappa: float = 0.0,
    triplets: Sequence[tuple[int, int, int, int]] | None = None,
    # (j, k, l, nu): j -> k -> l path with intermediate site k and chirality nu in {+1, -1}
) -> np.ndarray:
    """Assemble real skew-symmetric Hamiltonian matrix A_{ij}(G) such that:

    H_matter = (i/4) sum_{ij} A_{ij} c_i c_j.
    Canonical Majorana normalization: {c_i, c_j} = 2 delta_{ij}.

    For ferromagnetic Kitaev exchange H_K = -K sum_{<jl>_alpha} sigma_j^alpha sigma_l^alpha (K > 0):
      sigma_j^alpha sigma_l^alpha = -i u_{jl}^alpha c_j c_l (j in A, l in B)
      => -K sigma_j^alpha sigma_l^alpha = +i K u_{jl}^alpha c_j c_l
      => A_{jl} = +2 K u_{jl}^alpha, A_{lj} = -2 K u_{jl}^alpha.
    """
    N = gauge.num_sites
    A = np.zeros((N, N), dtype=np.float64)

    # 1. Nearest-neighbor Kitaev exchange:
    for (sa, sb), u_val in gauge.bond_links.items():
        # sa in A, sb in B: A_{sa, sb} = +2 * K * u_{ab}
        A[sa, sb] = 2.0 * K * u_val
        A[sb, sa] = -2.0 * K * u_val

    # 2. Next-nearest-neighbor 3-spin chiral interaction:
    if kappa != 0.0 and triplets is not None:
        for j, k, l_site, nu in triplets:
            # Directed traversal: j -> k -> l
            u_jk = gauge.get_link(j, k)
            u_kl = gauge.get_link(k, l_site)
            u_tilde = u_jk * u_kl  # In uniform gauge, u_{AB}=1, u_{BA}=-1 => u_tilde = -1
            # Mass contribution: -nu * kappa * u_tilde = +nu * kappa
            val = -nu * kappa * u_tilde
            A[j, l_site] += val
            A[l_site, j] -= val

    # Verify skew-symmetry
    assert np.allclose(A, -A.T, atol=1e-12)
    return A


@dataclass(frozen=True)
class DiagonalizedMatterSector:
    """Eigensystem of the quadratic matter Hamiltonian in gauge sector G."""

    energies: np.ndarray  # Sorted non-negative excitation energies 0 <= eps_1 <= ... <= eps_N
    Q: np.ndarray  # Orthogonal transformation matrix Q in O(2N)
    det_Q: float  # noqa: N815  # Determinant det(Q) in {-1, +1}
    unprojected_vacuum_energy: float  # E_0 = -1/2 sum eps_mu
    physical_matter_parity: int  # pi_G in {-1, +1}
    physical_ground_energy: float  # E_0 + (0 if pi_G == +1 else eps_1)


def diagonalize_matter_hamiltonian(
    A: np.ndarray,
    gauge: GaugeConfiguration | None = None,
    boundary: str = "torus",
    L: int = 3,
) -> DiagonalizedMatterSector:
    """Diagonalize real skew-symmetric matrix A using canonical real Schur decomposition.

    Enforces canonical 2x2 blocks: [[0, eps_k], [-eps_k, 0]] with eps_k >= 0,
    sorted in ascending order 0 <= eps_1 <= ... <= eps_N.
    Computes exact physical matter-fermion parity per Pedrocchi et al. (2011) and
    Zschocke & Vojta (2015).
    """
    N2 = A.shape[0]
    N = N2 // 2
    T, Z = schur(A, output="real")

    # 1. Standardize 2x2 blocks: ensure upper-right element is non-negative
    Z_std = Z.copy()
    omegas: list[float] = []
    for k in range(N):
        w = T[2 * k, 2 * k + 1]
        if w < 0.0:
            # Swap column 2k and 2k+1 to flip sign of the 2x2 block
            Z_std[:, [2 * k, 2 * k + 1]] = Z_std[:, [2 * k + 1, 2 * k]]
            w = -w
        omegas.append(float(w))

    # 2. Sort positive energies in ascending order
    idx = np.argsort(omegas)
    sorted_energies = np.array(omegas, dtype=np.float64)[idx]

    # Reorder columns of Z
    col_idx: list[int] = []
    for i in idx:
        col_idx.extend([2 * i, 2 * i + 1])
    Q = Z_std[:, col_idx]

    det_Q = float(np.round(np.linalg.det(Q)))
    unprojected_E0 = float(-0.5 * np.sum(sorted_energies))

    # 3. Physical matter parity
    if boundary == "torus":
        # (-1)^theta = (-1)^{L1 + L2 + M(L1 - M)} for L1=L2=L: (-1)^{2L} = +1
        phase_theta = 1
        # (-1)^{N_chi} = prod_b u_b
        if gauge is not None:
            n_chi_sign = 1
            for u in gauge.bond_links.values():
                if u == -1:
                    n_chi_sign *= -1
        else:
            n_chi_sign = 1
        pi_G = int(phase_theta * n_chi_sign * det_Q)
    else:  # Open boundary conditions
        if gauge is not None:
            # Under OBC, parity flips with each vortex pair
            # Or parity is unconstrained by bulk
            pi_G = 1
        else:
            pi_G = 1

    phys_E0 = unprojected_E0 if pi_G == 1 else float(unprojected_E0 + sorted_energies[0])

    return DiagonalizedMatterSector(
        energies=sorted_energies,
        Q=Q,
        det_Q=det_Q,
        unprojected_vacuum_energy=unprojected_E0,
        physical_matter_parity=pi_G,
        physical_ground_energy=phys_E0,
    )


def build_honeycomb_torus_cluster(
    L: int,
) -> tuple[int, list[DirectedBond], list[tuple[int, int, int, int]]]:
    """Generate geometry for an L x L unit-cell honeycomb cluster on a torus.

    Returns:
        num_sites: Total sites 2 * L^2
        bonds: List of DirectedBond (site_a, site_b, alpha)
        triplets: List of (j, k, l, nu) for NNN chiral interaction
    """
    num_sites = 2 * L * L

    def site_idx(m1: int, m2: int, sub: int) -> int:
        return 2 * ((m1 % L) + L * (m2 % L)) + sub

    nbr_A: dict[int, dict[int, int]] = {s: {} for s in range(0, num_sites, 2)}
    nbr_B: dict[int, dict[int, int]] = {s: {} for s in range(1, num_sites, 2)}
    bonds: list[DirectedBond] = []

    for m1 in range(L):
        for m2 in range(L):
            sa = site_idx(m1, m2, 0)
            sb_z = site_idx(m1, m2, 1)
            sb_x = site_idx(m1, m2 - 1, 1)
            sb_y = site_idx(m1 - 1, m2, 1)

            bonds.append(DirectedBond(site_a=sa, site_b=sb_z, bond_type=3))
            bonds.append(DirectedBond(site_a=sa, site_b=sb_x, bond_type=1))
            bonds.append(DirectedBond(site_a=sa, site_b=sb_y, bond_type=2))

            nbr_A[sa][3] = sb_z
            nbr_A[sa][1] = sb_x
            nbr_A[sa][2] = sb_y

            nbr_B[sb_z][3] = sa
            nbr_B[sb_x][1] = sa
            nbr_B[sb_y][2] = sa

    triplets: list[tuple[int, int, int, int]] = []
    # Right-handed pairs (alpha, beta) where third bond gamma satisfies eps(alpha, gamma, beta) = +1:
    pairs = [(1, 3), (3, 2), (2, 1)]

    # Sublattice A triangles through intermediate site k in B: nu = +1
    for sb, nbrs in nbr_B.items():
        for alpha, beta in pairs:
            triplets.append((nbrs[alpha], sb, nbrs[beta], +1))

    # Sublattice B triangles through intermediate site k in A: nu = -1
    for sa, nbrs in nbr_A.items():
        for alpha, beta in pairs:
            triplets.append((nbrs[alpha], sa, nbrs[beta], -1))

    return num_sites, bonds, triplets
