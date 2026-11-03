"""Pure physics services and numerical solvers."""

from kitaev_mf.services.correlators import (
    compute_correlators_from_bloch,
    compute_gauge_residuals,
    compute_spectral_projector,
)
from kitaev_mf.services.hamiltonian import (
    build_bloch_hamiltonian,
    build_bloch_hamiltonian_2site,
    build_bloch_hamiltonian_4site,
)
from kitaev_mf.services.solver import (
    JointAndersonSolver,
    solve_self_consistent_mean_field,
)

__all__ = [
    "build_bloch_hamiltonian",
    "build_bloch_hamiltonian_2site",
    "build_bloch_hamiltonian_4site",
    "compute_spectral_projector",
    "compute_correlators_from_bloch",
    "compute_gauge_residuals",
    "JointAndersonSolver",
    "solve_self_consistent_mean_field",
]
