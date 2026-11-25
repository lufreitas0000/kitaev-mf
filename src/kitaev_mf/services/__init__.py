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
from kitaev_mf.services.gutzwiller_benchmark import export_wave2b_gutzwiller_stability_benchmark
from kitaev_mf.services.gutzwiller_stability import (
    AutocorrelationDiagnostics,
    FluxSectorStabilityResult,
    NDecouplingBreakingResult,
    SU2GaugeStabilityResult,
    ThoulessSU2BlockCheck,
    build_flux_sector_bonds,
    build_site_su2_gauge_orthogonal,
    compute_chain_diagnostics,
    evaluate_flux_sector_gutzwiller,
    evaluate_n_decoupling_su2_breaking,
    evaluate_su2_gauge_stability,
    exact_energy_variance,
    make_staggered_su2_patterns,
    verify_thouless_su2_block_diagonality,
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
    "AutocorrelationDiagnostics",
    "FluxSectorStabilityResult",
    "HellmannFeynmanResult",
    "JointAndersonSolver",
    "KleinDualResult",
    "NDecouplingBreakingResult",
    "PurificationMetrics",
    "RealSpaceTorusCluster",
    "SU2GaugeStabilityResult",
    "ScalableSchurResult",
    "ThoulessSU2BlockCheck",
    "build_bloch_hamiltonian",
    "build_bloch_hamiltonian_2site",
    "build_bloch_hamiltonian_4site",
    "build_flux_sector_bonds",
    "build_site_su2_gauge_orthogonal",
    "canonical_skew_schur",
    "check_energy_consistency",
    "compute_bond_spin_correlators",
    "compute_chain_diagnostics",
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
    "evaluate_flux_sector_gutzwiller",
    "evaluate_klein_dual_pair",
    "evaluate_n_decoupling_su2_breaking",
    "evaluate_su2_gauge_stability",
    "evaluate_thermodynamic_susceptibility",
    "exact_energy_variance",
    "export_wave2b_gutzwiller_stability_benchmark",
    "make_staggered_su2_patterns",
    "midpoint_convergence",
    "midpoint_form_factor_average",
    "purify_newton_schulz",
    "purify_pade_sign",
    "solve_self_consistent_mean_field",
    "verify_thouless_su2_block_diagonality",
]
