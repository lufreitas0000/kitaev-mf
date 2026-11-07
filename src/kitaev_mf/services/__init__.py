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
from kitaev_mf.services.exact_matter import (
    build_antisymmetric_majorana_matrix_A,
    build_matter_matrix_T,
    build_reference_links,
    compute_matter_state_parity,
    compute_matter_vacuum_parity,
    compute_physical_ground_state_energy,
    compute_thermodynamic_vacuum_quadrature,
    compute_torus_gauge_sign,
    compute_torus_plaquette_fluxes,
    compute_torus_wilson_loops,
    create_vison_pair_links,
    evaluate_vison_pair_energy,
    extrapolate_vison_gap_infinite_separation,
    find_dual_path,
    get_plaquette_neighbors,
    is_physical_state,
)
from kitaev_mf.services.exact_perturbation import (
    ExactPerturbationResult,
    audit_perturbation_validity_bounds,
    build_spin_zeeman_perturbation,
    compute_exact_third_order_coupling,
    compute_spin_resolvent_third_order,
)
from kitaev_mf.services.feshbach import (
    FeshbachThresholds,
    audit_extensive_norm_fallacy,
    compute_spin_feshbach_thresholds,
    compute_torus_feshbach_thresholds,
    evaluate_eliminated_block_resolvent_norm,
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
)
from kitaev_mf.services.quadrature import (
    midpoint_convergence,
    midpoint_form_factor_average,
)
from kitaev_mf.services.real_space import RealSpaceTorusCluster
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
    solve_self_consistent_mean_field,
)

__all__ = [
    "ExactPerturbationResult",
    "FeshbachThresholds",
    "HellmannFeynmanResult",
    "JointAndersonSolver",
    "KleinDualResult",
    "PurificationMetrics",
    "RealSpaceTorusCluster",
    "audit_extensive_norm_fallacy",
    "audit_perturbation_validity_bounds",
    "build_antisymmetric_majorana_matrix_A",
    "build_bloch_hamiltonian",
    "build_bloch_hamiltonian_2site",
    "build_bloch_hamiltonian_4site",
    "build_matter_matrix_T",
    "build_reference_links",
    "build_spin_zeeman_perturbation",
    "check_energy_consistency",
    "compute_bond_spin_correlators",
    "compute_correlators_from_bloch",
    "compute_exact_third_order_coupling",
    "compute_gauge_residuals",
    "compute_gaussian_energy",
    "compute_matter_state_parity",
    "compute_matter_vacuum_parity",
    "compute_physical_ground_state_energy",
    "compute_plaquette_fluxes",
    "compute_riemannian_gradient_norm",
    "compute_spectral_ground_state_energy",
    "compute_spectral_projector",
    "compute_spin_feshbach_thresholds",
    "compute_spin_resolvent_third_order",
    "compute_sublattice_magnetizations",
    "compute_thermodynamic_vacuum_quadrature",
    "compute_torus_feshbach_thresholds",
    "compute_torus_gauge_sign",
    "compute_torus_plaquette_fluxes",
    "compute_torus_wilson_loops",
    "create_vison_pair_links",
    "evaluate_eliminated_block_resolvent_norm",
    "evaluate_finite_difference_magnetization",
    "evaluate_klein_dual_pair",
    "evaluate_thermodynamic_susceptibility",
    "evaluate_vison_pair_energy",
    "extrapolate_vison_gap_infinite_separation",
    "find_dual_path",
    "get_plaquette_neighbors",
    "is_physical_state",
    "midpoint_convergence",
    "midpoint_form_factor_average",
    "purify_newton_schulz",
    "solve_self_consistent_mean_field",
    "JointAndersonSolver",
    "KleinDualResult",
    "midpoint_convergence",
    "PurificationMetrics",
    "purify_newton_schulz",
    "RealSpaceTorusCluster",
    "midpoint_form_factor_average",
    "solve_self_consistent_mean_field",
]
