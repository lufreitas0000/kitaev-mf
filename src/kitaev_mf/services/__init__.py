"""Pure physics services and numerical solvers."""

from kitaev_mf.services.canonical_skew import ScalableSchurResult, canonical_skew_schur
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
from kitaev_mf.services.purification import (
    PurificationMetrics,
    purify_newton_schulz,
    purify_pade_sign,
)
from kitaev_mf.services.quadrature import (
    midpoint_convergence,
    midpoint_form_factor_average,
)
from kitaev_mf.services.real_space import RealSpaceTorusCluster, compute_skew_pfaffian
from kitaev_mf.services.response import (
    HellmannFeynmanResult,
    KleinDualResult,
    compute_riemannian_gradient_norm,
    evaluate_finite_difference_magnetization,
    evaluate_klein_dual_pair,
    evaluate_thermodynamic_susceptibility,
)
from kitaev_mf.services.solver import (
    JointAndersonSolver,
    compute_lehmann_dual_hessian,
    solve_self_consistent_mean_field,
)

__all__ = [
    "ScalableSchurResult",
    "build_bloch_hamiltonian",
    "build_bloch_hamiltonian_2site",
    "build_bloch_hamiltonian_4site",
    "canonical_skew_schur",
    "check_energy_consistency",
    "compute_bond_spin_correlators",
    "compute_correlators_from_bloch",
    "compute_gauge_residuals",
    "compute_gaussian_energy",
    "compute_lehmann_dual_hessian",
    "compute_plaquette_fluxes",
    "compute_riemannian_gradient_norm",
    "compute_skew_pfaffian",
    "compute_spectral_ground_state_energy",
    "compute_spectral_projector",
    "compute_sublattice_magnetizations",
    "evaluate_finite_difference_magnetization",
    "evaluate_klein_dual_pair",
    "evaluate_thermodynamic_susceptibility",
    "HellmannFeynmanResult",
    "JointAndersonSolver",
    "KleinDualResult",
    "midpoint_convergence",
    "PurificationMetrics",
    "purify_newton_schulz",
    "purify_pade_sign",
    "RealSpaceTorusCluster",
    "midpoint_form_factor_average",
    "solve_self_consistent_mean_field",
]
