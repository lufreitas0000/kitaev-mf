"""Hexagonal Architecture Port Protocols for hardware linear algebra and physics builders."""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np
import torch
from numpy.typing import NDArray

if TYPE_CHECKING:
    from kitaev_mf.domain.lattice import UnitCell
    from kitaev_mf.domain.qsl_bridge import (
        MajoranaBasis,
        PureKitaevBondScalars,
        QSLMeanFieldPortState,
        SiteSU2Decomposition,
    )
    from kitaev_mf.domain.state import MeanFieldState


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


@runtime_checkable
class QSLConventionBridgePort(Protocol):
    """Protocol for lossless cross-QSL ``conventions-v2.0`` state and observable bridging."""

    @property
    def conventions_version(self) -> str:
        """Locked cross-QSL contract version identifier."""
        ...

    def to_v2_state(self, state: MeanFieldState) -> QSLMeanFieldPortState:
        """Convert ``MeanFieldState`` from ``kitaev-mf`` ``(-c, b)`` to ``(+c, b)`` basis."""
        ...

    def from_v2_state(self, port_state: QSLMeanFieldPortState) -> MeanFieldState:
        """Convert ``QSLMeanFieldPortState`` from ``(+c, b)`` back to ``kitaev-mf`` ``(-c, b)`` basis."""
        ...

    def extract_bond_scalars(
        self,
        state: MeanFieldState | QSLMeanFieldPortState,
        cell: UnitCell,
    ) -> PureKitaevBondScalars:
        """Extract directed ``A -> B`` pure-Kitaev scalars ``(u_ij^alpha, v_ij^alpha)``."""
        ...

    def compute_plaquette_wilson_loops(
        self,
        u_gauge: NDArray[np.float64],
        cell: UnitCell,
    ) -> NDArray[np.float64]:
        """Compute Kitaev plaquette Wilson loops ``W_p`` from directed gauge links ``u_ij^alpha``."""
        ...

    def extract_site_su2(
        self,
        v_site: NDArray[np.float64],
        *,
        basis: MajoranaBasis = "kitaev_mf",
    ) -> SiteSU2Decomposition:
        """Extract basis-invariant ``(<S_j>, <G_j^{v2}>, <sigma_j>_unproj)`` expectations."""
        ...
