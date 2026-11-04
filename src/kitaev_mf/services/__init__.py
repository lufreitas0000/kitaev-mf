"""Pure physics services and numerical solvers."""

from kitaev_mf.services.correlators import (
    compute_correlators_from_bloch,
    compute_gauge_residuals,
    compute_spectral_projector,
)
from kitaev_mf.services.energetics import (
    check_energy_consistency,
    compute_gaussian_energy,
    compute_spectral_ground_state_energy,
)
from kitaev_mf.services.hamiltonian import (
    build_bloch_hamiltonian,
    build_bloch_hamiltonian_2site,
    build_bloch_hamiltonian_4site,
)
from kitaev_mf.services.observables import (
    compute_bond_spin_correlators,
    compute_plaquette_fluxes,
    compute_sublattice_magnetizations,
)
from kitaev_mf.services.quadrature import (
    midpoint_convergence,
    midpoint_form_factor_average,
)
from kitaev_mf.services.solver import (
    JointAndersonSolver,
    solve_self_consistent_mean_field,
)

__all__ = [
    "build_bloch_hamiltonian",
    "build_bloch_hamiltonian_2site",
    "build_bloch_hamiltonian_4site",
    "check_energy_consistency",
    "compute_bond_spin_correlators",
    "compute_correlators_from_bloch",
    "compute_gauge_residuals",
    "compute_gaussian_energy",
    "compute_plaquette_fluxes",
    "compute_spectral_ground_state_energy",
    "compute_spectral_projector",
    "compute_sublattice_magnetizations",
    "JointAndersonSolver",
    "midpoint_convergence",
    "midpoint_form_factor_average",
    "solve_self_consistent_mean_field",
]
