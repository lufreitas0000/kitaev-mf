"""Domain entities defining solver configuration contracts and execution metrics."""

from dataclasses import dataclass, field


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
    """

    k_grid_shape: tuple[int, int] = (24, 24)
    max_iterations: int = 200
    tol_correlator: float = 1e-10
    tol_gauge: float = 1e-10
    mixing_alpha: float = 0.5
    lambda_step: float = 0.2
    diis_history_depth: int = 6
    diis_reg: float = 1e-12

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


@dataclass(frozen=True)
class SolverMetrics:
    """Immutable telemetry recording the trajectory and final state of the solver loop.

    Attributes:
        converged: Whether correlator, on-site, and gauge residuals dropped below tolerance.
        iterations: Total number of iterations executed.
        final_correlator_residual: Final max-norm difference ||Delta U||_inf across bonds.
        final_onsite_residual: Final max-norm difference ||Delta V||_inf across sublattices.
        final_gauge_residual: Final max-norm gauge violation ||g||_inf across sublattices.
        residual_history: Chronological sequence of residual norms max(||Delta U||, ||Delta V||, ||g||).
    """

    converged: bool
    iterations: int
    final_correlator_residual: float
    final_onsite_residual: float
    final_gauge_residual: float
    residual_history: tuple[float, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if self.iterations < 0:
            raise ValueError(f"iterations cannot be negative, got {self.iterations}")
        if (
            self.final_correlator_residual < 0.0
            or self.final_onsite_residual < 0.0
            or self.final_gauge_residual < 0.0
        ):
            raise ValueError("Residual norms cannot be negative.")
