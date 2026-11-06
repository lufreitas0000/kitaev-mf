"""Domain entities defining solver configuration contracts and execution metrics."""

from dataclasses import dataclass, field
from typing import Literal


@dataclass(frozen=True)
class SolverConfig:
    """Immutable configuration parameters for the self-consistent mean-field solver.

    Attributes:
        k_grid_shape: Number of Monkhorst-Pack subdivisions (N_k1, N_k2) in fractional coordinates.
        max_iterations: Maximum number of outer self-consistent iterations.
        tol_correlator: Convergence threshold on the max-norm of correlator updates ||Delta U||_inf.
        tol_gauge: Convergence threshold on the max-norm of gauge residuals ||g||_inf.
        mixing_alpha: Base linear mixing parameter alpha in (0, 1] for Picard state updates.
        lambda_step: Step length alpha_lambda for updating Lagrange multipliers from gauge residuals.
        diis_history_depth: Maximum number of historical iterations M retained for DIIS acceleration.
        diis_reg: Diagonal Tikhonov regularization epsilon parameter for the DIIS overlap matrix.
        method: Fixed-point acceleration algorithm ('fp', 'picard', 'diis', 'broyden', 'jfnk').
        broyden_alpha: Initial scaling / step parameter for Broyden inverse Jacobian update.
        jfnk_gmres_tol: Relative convergence tolerance for inner GMRES iterations in JFNK.
        jfnk_gmres_maxiter: Maximum number of inner Krylov subspace iterations per Newton step.
        jfnk_eps: Finite-difference directional perturbation scale for matrix-free Jacobian products.
        skew_mode: On-site skew-symmetry handling mode ('raw', 'symmetrize_output', 'retract_step', 'intrinsic').
    """

    k_grid_shape: tuple[int, int] = (24, 24)
    max_iterations: int = 200
    tol_correlator: float = 1e-10
    tol_gauge: float = 1e-10
    mixing_alpha: float = 0.5
    lambda_step: float = 0.2
    diis_history_depth: int = 6
    diis_reg: float = 1e-12
    method: Literal["diis", "fp", "picard", "broyden", "jfnk"] = "diis"
    broyden_alpha: float = 0.5
    jfnk_gmres_tol: float = 1e-4
    jfnk_gmres_maxiter: int = 15
    jfnk_eps: float = 1e-7
    skew_mode: Literal["raw", "symmetrize_output", "retract_step", "intrinsic"] = "retract_step"

    def __post_init__(self) -> None:
        if len(self.k_grid_shape) != 2 or any(n <= 0 for n in self.k_grid_shape):
            raise ValueError(
                f"k_grid_shape must be a pair of positive integers, got {self.k_grid_shape}"
            )
        if self.max_iterations <= 0:
            raise ValueError(f"max_iterations must be positive, got {self.max_iterations}")
        if self.tol_correlator <= 0.0 or self.tol_gauge <= 0.0:
            raise ValueError("Convergence tolerances must be strictly positive.")
        if not (0.0 < self.mixing_alpha <= 1.0):
            raise ValueError(f"mixing_alpha must be in (0, 1], got {self.mixing_alpha}")
        if self.lambda_step <= 0.0:
            raise ValueError(f"lambda_step must be strictly positive, got {self.lambda_step}")
        if self.diis_history_depth < 1:
            raise ValueError(
                f"diis_history_depth must be at least 1, got {self.diis_history_depth}"
            )
        if self.diis_reg < 0.0:
            raise ValueError(f"diis_reg must be non-negative, got {self.diis_reg}")
        if not (0.0 < self.broyden_alpha <= 1.0):
            raise ValueError(f"broyden_alpha must be in (0, 1], got {self.broyden_alpha}")
        if self.jfnk_gmres_tol <= 0.0:
            raise ValueError(f"jfnk_gmres_tol must be positive, got {self.jfnk_gmres_tol}")
        if self.jfnk_gmres_maxiter < 1:
            raise ValueError(f"jfnk_gmres_maxiter must be >= 1, got {self.jfnk_gmres_maxiter}")
        if self.jfnk_eps <= 0.0:
            raise ValueError(f"jfnk_eps must be positive, got {self.jfnk_eps}")


@dataclass(frozen=True)
class SolverMetrics:
    """Immutable telemetry recording the trajectory and final state of the solver loop.

    Attributes:
        converged: Whether correlator, on-site, and gauge residuals dropped below tolerance.
        iterations: Total number of outer self-consistent iterations executed.
        final_correlator_residual: Final max-norm difference ||Delta U||_inf across bonds.
        final_onsite_residual: Final max-norm difference ||Delta V||_inf across sublattices.
        final_gauge_residual: Final max-norm gauge violation ||g||_inf across sublattices.
        residual_history: Chronological sequence of residual norms max(||Delta U||, ||Delta V||, ||g||).
        contraction_factors: Sequence of step-to-step contraction ratios ||R_{k+1}|| / ||R_k||.
        comm_residuals: Sequence of BZ-averaged commutator stationarity residuals eps_comm.
        map_evaluations: Total number of self-consistency map Phi evaluations incurred.
    """

    converged: bool
    iterations: int
    final_correlator_residual: float
    final_onsite_residual: float
    final_gauge_residual: float
    residual_history: tuple[float, ...] = field(default_factory=tuple)
    contraction_factors: tuple[float, ...] = field(default_factory=tuple)
    comm_residuals: tuple[float, ...] = field(default_factory=tuple)
    map_evaluations: int = 0

    def __post_init__(self) -> None:
        if self.iterations < 0:
            raise ValueError(f"iterations cannot be negative, got {self.iterations}")
        if (
            self.final_correlator_residual < 0.0
            or self.final_onsite_residual < 0.0
            or self.final_gauge_residual < 0.0
        ):
            raise ValueError("Residual norms cannot be negative.")
        if self.map_evaluations < 0:
            raise ValueError(f"map_evaluations cannot be negative, got {self.map_evaluations}")
