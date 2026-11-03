"""Test hardware backend linear algebra operations."""

import torch

from kitaev_mf.infra.linalg_backend import TorchBackend


def test_eigh_diagonalization(cpu_backend: TorchBackend) -> None:
    # Construct a simple 2x2 Hermitian test matrix
    mat = torch.tensor([[2.0, 1.0j], [-1.0j, 3.0]], dtype=torch.complex128)
    eigenvalues, eigenvectors = cpu_backend.eigh(mat)
    # Check mat * v = lambda * v
    for i in range(len(eigenvalues)):
        val = eigenvalues[i]
        vec = eigenvectors[:, i]
        assert torch.allclose(mat @ vec, val * vec, atol=1e-12)
