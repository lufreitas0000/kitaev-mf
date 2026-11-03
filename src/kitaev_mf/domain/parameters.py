"""Immutable parameter dataclasses for the extended Kitaev model."""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class ModelParameters:
    """Coupling constants for the extended J-K-Gamma Kitaev model."""
    J: float = 0.0
    K: float = -1.0
    Gamma: float = 0.0
    h: tuple[float, float, float] = (0.0, 0.0, 0.0)  # Magnetic field (hx, hy, hz)
    D: float = 0.0                                  # Dzyaloshinskii-Moriya coupling
    V0: float = 0.0                                 # Local tip electrostatic potential
    eta: float = 1.0                                # Decoupling channel: 1.0 (Kitaev N), 0.0 (Chiral M)

    @property
    def h_vec(self) -> np.ndarray:
        return np.array(self.h, dtype=np.float64)
