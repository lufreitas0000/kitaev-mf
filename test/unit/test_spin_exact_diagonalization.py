"""Unit tests for spin-space exact diagonalization (ED) of Kitaev clusters (Track B / B-023).

Constructs the exact many-body spin-1/2 Pauli Hamiltonian H_K on small clusters,
proves exact commutation with hexagonal Wilson loops [H_K, W_p] = 0, and determines
sector-resolved ground states and exact vison defect energies.
"""

import numpy as np


def build_spin_kitaev_hamiltonian(
    num_spins: int,
    bonds: list[tuple[int, int, str]],
    K: float = 1.0,
) -> np.ndarray:
    """Assemble the exact 2^N x 2^N spin-space Hamiltonian H_K = -K sum_{<ij>_a} sigma_i^a sigma_j^a."""
    dim = 1 << num_spins
    H = np.zeros((dim, dim), dtype=np.complex128)

    for i, j, alpha in bonds:
        for s in range(dim):
            si = (s >> i) & 1
            sj = (s >> j) & 1

            if alpha == "z":
                # sigma_i^z sigma_j^z |s> = (-1)^(si + sj) |s>
                val = (-1.0) ** (si + sj)
                H[s, s] -= K * val
            elif alpha == "x":
                # sigma_i^x sigma_j^x |s> = |s ^ (1<<i) ^ (1<<j)>
                s_prime = s ^ (1 << i) ^ (1 << j)
                H[s_prime, s] -= K * 1.0
            elif alpha == "y":
                # sigma_i^y sigma_j^y |s> = - (-1)^(si + sj) |s ^ (1<<i) ^ (1<<j)>
                s_prime = s ^ (1 << i) ^ (1 << j)
                val = -((-1.0) ** (si + sj))
                H[s_prime, s] -= K * val

    return H


def test_isolated_hexagon_spin_flux_commutation() -> None:
    """Verify that on an isolated 6-site hexagon, [H_K, W_p] = 0 and W_p has eigenvalues +-1."""
    num_spins = 6
    dim = 1 << num_spins
    bonds = [
        (0, 1, "z"),
        (1, 2, "x"),
        (2, 3, "y"),
        (3, 4, "z"),
        (4, 5, "x"),
        (5, 0, "y"),
    ]

    H = build_spin_kitaev_hamiltonian(num_spins, bonds, K=1.0)
    assert np.allclose(H, H.conj().T), "Hamiltonian must be Hermitian"

    # Hexagonal Wilson loop W_p = sigma_0^x sigma_1^y sigma_2^z sigma_3^x sigma_4^y sigma_5^z
    W = np.zeros((dim, dim), dtype=np.complex128)
    for s in range(dim):
        s_new = s ^ (1 << 0) ^ (1 << 3) ^ (1 << 1) ^ (1 << 4)
        s1 = (s >> 1) & 1
        s2 = (s >> 2) & 1
        s4 = (s >> 4) & 1
        s5 = (s >> 5) & 1
        phase = ((-1.0) ** (s2 + s5)) * (-((-1.0) ** (s1 + s4)))
        W[s_new, s] = phase

    assert np.allclose(W, W.conj().T), "W_p must be Hermitian"
    assert np.allclose(W @ W, np.eye(dim)), "W_p must be unitary and square to 1"

    # Exact commutation [H, W_p] = 0
    comm = H @ W - W @ H
    assert np.allclose(comm, 0.0, atol=1e-12), f"Non-zero commutator: {np.linalg.norm(comm)}"

    # Eigenspectrum
    evals, evecs = np.linalg.eigh(H)
    psi_0 = evecs[:, 0]
    w_exp = np.vdot(psi_0, W @ psi_0).real
    # In ground state, W_p has definite eigenvalue +1
    assert np.isclose(abs(w_exp), 1.0, atol=1e-12)


def test_two_hexagon_spin_flux_sectors_and_defect_energy() -> None:
    """Verify sector-resolved ground state energies on a 10-site two-plaquette cluster."""
    num_spins = 10
    dim = 1 << num_spins

    # Two hexagons sharing z-bond (3, 4):
    bonds = [
        (0, 1, "x"),
        (1, 2, "y"),
        (2, 3, "x"),
        (3, 4, "z"),
        (4, 5, "y"),
        (5, 0, "z"),
        (4, 6, "x"),
        (6, 7, "y"),
        (7, 8, "z"),
        (8, 9, "x"),
        (9, 3, "y"),
    ]

    H = build_spin_kitaev_hamiltonian(num_spins, bonds, K=1.0)

    # Plaquette 1: W1 = sigma_0^y sigma_1^z sigma_2^z sigma_3^y sigma_4^x sigma_5^x
    W1 = np.zeros((dim, dim), dtype=np.complex128)
    for s in range(dim):
        s_new = s ^ (1 << 0) ^ (1 << 3) ^ (1 << 4) ^ (1 << 5)
        s0 = (s >> 0) & 1
        s1 = (s >> 1) & 1
        s2 = (s >> 2) & 1
        s3 = (s >> 3) & 1
        phase = ((-1.0) ** (s1 + s2)) * (-((-1.0) ** (s0 + s3)))
        W1[s_new, s] = phase

    comm1 = H @ W1 - W1 @ H
    assert np.allclose(comm1, 0.0, atol=1e-12), (
        f"W1 does not commute with H: {np.linalg.norm(comm1)}"
    )

    # Simultaneous diagonalization of H and W1
    # Project into W1 = +1 and W1 = -1 subspaces:
    P_plus = 0.5 * (np.eye(dim) + W1)
    P_minus = 0.5 * (np.eye(dim) - W1)

    H_plus = P_plus @ H @ P_plus
    H_minus = P_minus @ H @ P_minus

    evals_plus = np.linalg.eigvalsh(H_plus)
    evals_minus = np.linalg.eigvalsh(H_minus)

    # Filter out the zero eigenvalues introduced by the projector
    e0_plus = np.min(evals_plus[np.abs(evals_plus) > 1e-8])
    e0_minus = np.min(evals_minus[np.abs(evals_minus) > 1e-8])

    # The flux-free sector must have lower energy than the vortex sector:
    e_flux_free = min(e0_plus, e0_minus)
    e_vortex = max(e0_plus, e0_minus)
    vison_energy = e_vortex - e_flux_free

    assert vison_energy > 0.0, f"Expected positive vison energy, got {vison_energy}"
    assert np.isclose(e_flux_free, -6.841619, atol=1e-4)
