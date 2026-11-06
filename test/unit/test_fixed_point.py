"""Unit tests verifying mathematical contracts of fixed-point accelerators."""

from __future__ import annotations

import numpy as np

from kitaev_mf.services.fixed_point import (
    BroydenSolver,
    DampedPicard,
    JFNKSolver,
    JointAndersonSolver,
    PlainFixedPoint,
)


def test_contractive_linear_system_all_solvers() -> None:
    """All accelerators must converge on a strictly contractive linear system."""
    dim = 6
    rng = np.random.default_rng(42)

    # Random contractive matrix with spectral radius rho ~ 0.5
    mat = rng.standard_normal((dim, dim))
    mat = 0.5 * mat / float(np.linalg.norm(mat, 2))
    b = rng.standard_normal(dim)

    # Fixed point x* = (I - mat)^{-1} b
    x_exact = np.linalg.solve(np.eye(dim) - mat, b)

    def phi_map(x_vec: np.ndarray) -> np.ndarray:
        return np.asarray(mat @ x_vec + b, dtype=np.float64)

    solvers = [
        PlainFixedPoint(),
        DampedPicard(alpha=0.6),
        JointAndersonSolver(history_depth=4),
        BroydenSolver(alpha=0.6),
        JFNKSolver(gmres_tol=1e-5),
    ]

    for solver in solvers:
        solver.reset()
        x_curr = np.zeros(dim)
        for _ in range(50):
            f_val = phi_map(x_curr)
            if np.linalg.norm(f_val - x_curr) < 1e-10:
                break
            x_curr, _ = solver.step(x_curr, f_val, residual_fn=lambda x: phi_map(x) - x)

        np.testing.assert_allclose(x_curr, x_exact, atol=1e-8)


def test_non_contractive_saddle_system() -> None:
    """Accelerators (DIIS, Broyden, JFNK) must converge on coupled saddle systems where FP diverges."""
    dim = 8
    rng = np.random.default_rng(101)

    k_primal = rng.standard_normal((4, 4))
    k_primal = k_primal.T @ k_primal + np.eye(4)
    k_dual = rng.standard_normal((4, 4))
    k_dual = k_dual.T @ k_dual + np.eye(4)
    c_block = 0.4 * rng.standard_normal((4, 4))

    mat_saddle = np.block([[k_primal, c_block], [-c_block.T, k_dual]])
    step_size = 0.7 / float(np.linalg.norm(mat_saddle, 2))
    x_target = rng.standard_normal(dim)
    b_vec = -mat_saddle @ x_target

    def phi_map(x_vec: np.ndarray) -> np.ndarray:
        return np.asarray(x_vec - step_size * (mat_saddle @ x_vec + b_vec), dtype=np.float64)

    solvers = [
        JointAndersonSolver(history_depth=5, reg=1e-12),
        BroydenSolver(alpha=0.5),
        JFNKSolver(gmres_tol=1e-5, gmres_maxiter=10),
    ]

    for solver in solvers:
        solver.reset()
        x_curr = np.zeros(dim)
        converged = False
        for _ in range(80):
            f_val = phi_map(x_curr)
            res_norm = float(np.linalg.norm(f_val - x_curr))
            if res_norm < 1e-9:
                converged = True
                break
            x_curr, _ = solver.step(x_curr, f_val, residual_fn=lambda x: phi_map(x) - x)

        assert converged is True
        np.testing.assert_allclose(x_curr, x_target, atol=1e-7)


def test_jfnk_map_evaluation_accounting() -> None:
    """Verify that JFNKSolver accurately reports extra map evaluations from GMRES."""
    dim = 4

    def phi_map(x_vec: np.ndarray) -> np.ndarray:
        return np.asarray(0.2 * x_vec + 1.0, dtype=np.float64)

    solver = JFNKSolver(gmres_tol=1e-4)
    x0 = np.zeros(dim)
    f0 = phi_map(x0)

    call_count = 0

    def residual_fn(x_in: np.ndarray) -> np.ndarray:
        nonlocal call_count
        call_count += 1
        res = phi_map(x_in) - x_in
        return np.asarray(res, dtype=np.float64)

    x1, extra_evals = solver.step(x0, f0, residual_fn=residual_fn)
    assert extra_evals > 0
    assert call_count == extra_evals

    # Without residual_fn, JFNKSolver falls back to direct fixed-point output
    x_fallback, evals_fallback = solver.step(x0, f0, residual_fn=None)
    np.testing.assert_allclose(x_fallback, f0)
    assert evals_fallback == 0


def test_accelerator_input_validations() -> None:
    """Verify that invalid parameters in fixed-point accelerators raise ValueError."""
    import pytest

    with pytest.raises(ValueError, match="alpha"):
        DampedPicard(alpha=0.0)

    with pytest.raises(ValueError, match="history_depth"):
        JointAndersonSolver(history_depth=0)

    with pytest.raises(ValueError, match="reg"):
        JointAndersonSolver(reg=-1.0)

    with pytest.raises(ValueError, match="alpha"):
        BroydenSolver(alpha=0.0)
