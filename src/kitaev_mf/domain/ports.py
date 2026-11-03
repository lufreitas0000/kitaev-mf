"""Hexagonal Architecture Port Protocols for hardware linear algebra and physics builders."""

from typing import Protocol, runtime_checkable

import numpy as np
import torch


@runtime_checkable
class LinalgBackendPort(Protocol):
    """Protocol defining hardware acceleration operations for eigensolvers and array handling."""

    @property
    def device(self) -> torch.device:
        """Target execution device (CPU or CUDA)."""
        ...

    def to_tensor(self, array: np.ndarray, dtype: torch.dtype = torch.complex128) -> torch.Tensor:
        """Convert a NumPy array to a PyTorch tensor on the target hardware device."""
        ...

    def to_numpy(self, tensor: torch.Tensor) -> np.ndarray:
        """Convert a PyTorch tensor back to a NumPy array on the host CPU."""
        ...

    def eigh(self, matrix: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """Compute eigenvalues and eigenvectors for a single Hermitian matrix."""
        ...

    def batched_eigh(self, matrices: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """Compute eigenvalues and eigenvectors for a batch of Hermitian matrices."""
        ...


@runtime_checkable
class HamiltonianBuilderPort(Protocol):
    """Protocol for momentum-space Bloch-Majorana matrix assembly."""

    def __call__(
        self,
        k_pts: np.ndarray | torch.Tensor,
        backend: LinalgBackendPort | None = None,
    ) -> torch.Tensor:
        """Build the Bloch-Majorana Hamiltonian at the specified reciprocal points."""
        ...
