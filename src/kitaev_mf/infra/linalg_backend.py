"""Hardware-agnostic linear algebra adapter implementing PyTorch/CPU/GPU dispatch."""

import numpy as np
import torch

from kitaev_mf.domain.ports import LinalgBackendPort


class TorchBackend(LinalgBackendPort):
    """PyTorch implementation of LinalgBackendPort supporting CUDA and CPU seamlessly."""

    def __init__(self, device: str | None = None) -> None:
        if device is None:
            self._device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            self._device = torch.device(device)

    @property
    def device(self) -> torch.device:
        return self._device

    def to_tensor(self, array: np.ndarray, dtype: torch.dtype = torch.complex128) -> torch.Tensor:
        return torch.tensor(array, dtype=dtype, device=self._device)

    def to_numpy(self, tensor: torch.Tensor) -> np.ndarray:
        return tensor.detach().cpu().numpy()

    def eigh(self, matrix: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        with torch.no_grad():
            eigenvalues, eigenvectors = torch.linalg.eigh(matrix)
        return eigenvalues, eigenvectors

    def batched_eigh(self, matrices: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        with torch.no_grad():
            eigenvalues, eigenvectors = torch.linalg.eigh(matrices)
        return eigenvalues, eigenvectors
