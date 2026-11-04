"""Test domain entities and validation rules for SolverConfig and SolverMetrics."""

from dataclasses import FrozenInstanceError

import pytest

from kitaev_mf.domain.solver import SolverConfig, SolverMetrics


def test_solver_config_defaults_and_immutability() -> None:
    config = SolverConfig()
    assert config.k_grid_shape == (24, 24)
    assert config.max_iterations == 200
    assert config.tol_correlator == 1e-10
    assert config.tol_gauge == 1e-10
    assert config.mixing_alpha == 0.5
    assert config.lambda_step == 0.2
    assert config.diis_history_depth == 6
    assert config.diis_reg == 1e-12

    # Immutability check: modifying a frozen dataclass raises FrozenInstanceError
    with pytest.raises(FrozenInstanceError):
        config.max_iterations = 100  # type: ignore[misc]


def test_solver_config_custom_instances_for_scaling() -> None:
    """Verify that custom parameterizations for finite-size scaling instantiate cleanly."""
    custom_config = SolverConfig(
        k_grid_shape=(96, 96),
        max_iterations=500,
        tol_correlator=1e-12,
        tol_gauge=1e-12,
        mixing_alpha=0.25,
        lambda_step=0.05,
        diis_history_depth=10,
        diis_reg=1e-10,
    )
    assert custom_config.k_grid_shape == (96, 96)
    assert custom_config.max_iterations == 500
    assert custom_config.tol_correlator == 1e-12
    assert custom_config.diis_history_depth == 10


def test_solver_config_validations() -> None:
    with pytest.raises(ValueError, match="k_grid_shape"):
        SolverConfig(k_grid_shape=(0, 24))

    with pytest.raises(ValueError, match="max_iterations"):
        SolverConfig(max_iterations=0)

    with pytest.raises(ValueError, match="tolerances"):
        SolverConfig(tol_correlator=-1e-5)

    with pytest.raises(ValueError, match="mixing_alpha"):
        SolverConfig(mixing_alpha=0.0)

    with pytest.raises(ValueError, match="mixing_alpha"):
        SolverConfig(mixing_alpha=1.5)

    with pytest.raises(ValueError, match="lambda_step"):
        SolverConfig(lambda_step=-0.1)

    with pytest.raises(ValueError, match="diis_history_depth"):
        SolverConfig(diis_history_depth=0)

    with pytest.raises(ValueError, match="diis_reg"):
        SolverConfig(diis_reg=-1e-10)


def test_solver_metrics_defaults_and_validations() -> None:
    metrics = SolverMetrics(
        converged=True,
        iterations=15,
        final_correlator_residual=1.2e-11,
        final_onsite_residual=3.4e-12,
        final_gauge_residual=4.5e-12,
        residual_history=(0.1, 0.01, 1e-3, 1.2e-11),
    )
    assert metrics.converged is True
    assert metrics.iterations == 15
    assert metrics.final_onsite_residual == 3.4e-12
    assert len(metrics.residual_history) == 4

    with pytest.raises(FrozenInstanceError):
        metrics.iterations = 20  # type: ignore[misc]

    with pytest.raises(ValueError, match="iterations"):
        SolverMetrics(
            converged=False,
            iterations=-1,
            final_correlator_residual=0.1,
            final_onsite_residual=0.1,
            final_gauge_residual=0.1,
        )

    with pytest.raises(ValueError, match="Residual norms"):
        SolverMetrics(
            converged=False,
            iterations=1,
            final_correlator_residual=-0.1,
            final_onsite_residual=0.1,
            final_gauge_residual=0.1,
        )

    with pytest.raises(ValueError, match="Residual norms"):
        SolverMetrics(
            converged=False,
            iterations=1,
            final_correlator_residual=0.1,
            final_onsite_residual=-0.1,
            final_gauge_residual=0.1,
        )

    with pytest.raises(ValueError, match="Residual norms"):
        SolverMetrics(
            converged=False,
            iterations=1,
            final_correlator_residual=0.1,
            final_onsite_residual=0.1,
            final_gauge_residual=-0.1,
        )
