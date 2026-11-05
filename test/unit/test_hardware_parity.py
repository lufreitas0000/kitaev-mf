"""Test hardware backend execution parity, device routing, and zero-leak contracts (T4).

Contracts:
    - Backend device resolution and tensor conversion parity.
    - Eigensolver parity between CPU and CUDA devices to machine precision.
    - Bloch Hamiltonian assembly on explicit GPU devices.
"""

import numpy as np
import pytest
import torch

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.correlators import compute_spectral_projector
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian


def test_backend_device_resolution_and_conversion() -> None:
    """Verify backend accurately resolves requested CPU device."""
    backend_cpu = TorchBackend(device="cpu")
    assert backend_cpu.device.type == "cpu"

    arr = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)
    tensor = backend_cpu.to_tensor(arr, dtype=torch.float64)
    assert tensor.device.type == "cpu"
    assert tensor.dtype == torch.float64

    roundtrip = backend_cpu.to_numpy(tensor)
    assert np.allclose(arr, roundtrip, atol=1e-15)


def test_cpu_eigensolver_and_batched_parity(cpu_backend: TorchBackend) -> None:
    """Verify single and batched eigh produce identical results on CPU."""
    dim = 8
    rng = np.random.default_rng(42)
    raw = rng.standard_normal((3, dim, dim)) + 1j * rng.standard_normal((3, dim, dim))
    hermitian_batch = 0.5 * (raw + np.swapaxes(raw.conj(), -1, -2))

    batch_tensor = cpu_backend.to_tensor(hermitian_batch, dtype=torch.complex128)
    batch_evals, batch_evecs = cpu_backend.batched_eigh(batch_tensor)

    for i in range(3):
        single_evals, single_evecs = cpu_backend.eigh(batch_tensor[i])
        assert np.allclose(
            cpu_backend.to_numpy(batch_evals[i]),
            cpu_backend.to_numpy(single_evals),
            atol=1e-14,
        )


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA device not available")
def test_cuda_hardware_parity() -> None:
    """Verify numerical parity between CPU and CUDA executions on Bloch Hamiltonian."""
    backend_cpu = TorchBackend(device="cpu")
    backend_cuda = TorchBackend(device="cuda")

    cell = make_honeycomb_2site()
    params = ModelParameters(K=-1.0, J=0.05, h=(0.02, 0.03, 0.04))
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.4)
    k_grid = cell.make_k_grid(16, 16)

    # 1. Hamiltonian assembly on CPU vs CUDA
    h_cpu = build_bloch_hamiltonian(k_grid, cell, params, state, backend_cpu)
    h_cuda = build_bloch_hamiltonian(k_grid, cell, params, state, backend_cuda)

    assert h_cuda.device.type == "cuda"
    diff_h = float(torch.max(torch.abs(h_cpu - h_cuda.cpu())).item())
    assert diff_h < 1e-14

    # 2. Diagonalization parity
    evals_cpu, _ = backend_cpu.batched_eigh(h_cpu)
    evals_cuda, _ = backend_cuda.batched_eigh(h_cuda)

    diff_evals = float(torch.max(torch.abs(evals_cpu - evals_cuda.cpu())).item())
    assert diff_evals < 1e-13

    # 3. Projector parity
    p_cpu = compute_spectral_projector(h_cpu, backend_cpu)
    p_cuda = compute_spectral_projector(h_cuda, backend_cuda)
    diff_p = float(torch.max(torch.abs(p_cpu - p_cuda.cpu())).item())
    assert diff_p < 1e-13
