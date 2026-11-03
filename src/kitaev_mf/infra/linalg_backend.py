"""Hardware-agnostic linear algebra service implementing PyTorch/CPU/GPU dispatch."""
from typing import Protocol, Tuple
import torch
import numpy as np

class LinalgBackend(Protocol):
    """Protocol defining the linear algebra interface for eigensolvers and root finders."""
    def eigh(self, matrix: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]: ...
    def batched_eigh(self, matrices: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]: ...

class TorchBackend:
    """PyTorch implementation of LinalgBackend supporting CUDA and CPU seamlessly."""
    def __init__(self, device: str | None = None) -> None:
        if device is None:
            self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        else:
            self.device = torch.device(device)

    def to_tensor(self, array: np.ndarray, dtype: torch.dtype = torch.complex128) -> torch.Tensor:
        return torch.tensor(array, dtype=dtype, device=self.device)

    def to_numpy(self, tensor: torch.Tensor) -> np.ndarray:
        return tensor.detach().cpu().numpy()

    def eigh(self, matrix: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """Eigenvalues and eigenvectors for a single Hermitian matrix."""
        with torch.no_grad():
            eigenvalues, eigenvectors = torch.linalg.eigh(matrix)
        return eigenvalues, eigenvectors

    def batched_eigh(self, matrices: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """Eigenvalues and eigenvectors for a batch of Hermitian matrices (e.g., k-space)."""
        with torch.no_grad():
            eigenvalues, eigenvectors = torch.linalg.eigh(matrices)
        return eigenvalues, eigenvectors
