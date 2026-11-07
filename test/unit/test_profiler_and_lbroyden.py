"""Unit tests for ResourceProfiler, LimitedMemoryBroydenSolver, and symmetric seed perturbations."""

from __future__ import annotations

import numpy as np
import pytest

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.infra.profiler import ResourceProfiler
from kitaev_mf.services.fixed_point import LimitedMemoryBroydenSolver
from kitaev_mf.services.seeds import add_symmetric_perturbation, make_qsl_seed


def test_resource_profiler_flop_counts() -> None:
    """Verify analytical FLOP accounting for eigensolvers and GEMM operations."""
    profiler = ResourceProfiler()

    # Real GEMM: 2 * M * N * K
    profiler.count_gemm(10, 20, 30, is_complex=False)
    assert profiler.total_flops == pytest.approx(2.0 * 10 * 20 * 30)

    # Complex GEMM: 8 * M * N * K
    profiler.reset()
    profiler.count_gemm(10, 20, 30, is_complex=True)
    assert profiler.total_flops == pytest.approx(8.0 * 10 * 20 * 30)

    # Real eigensolver: (4/3) M^3
    profiler.reset()
    profiler.count_eigh(12, is_complex=False)
    assert profiler.total_flops == pytest.approx((4.0 / 3.0) * (12**3))

    # Complex eigensolver: (16/3) M^3
    profiler.reset()
    profiler.count_eigh(12, is_complex=True)
    assert profiler.total_flops == pytest.approx((16.0 / 3.0) * (12**3))

    # Purification step: Newton-Schulz 4 M^3
    profiler.reset()
    profiler.count_purification_step(16, step_type="newton_schulz")
    assert profiler.total_flops == pytest.approx(4.0 * (16**3))


def test_limited_memory_broyden_linear_system() -> None:
    """Verify that L-Broyden accurately solves a linear root-finding problem G(x) = A x - b."""
    rng = np.random.default_rng(42)
    dim = 15
    # Construct well-conditioned positive definite matrix
    q, _ = np.linalg.qr(rng.standard_normal((dim, dim)))
    diag = np.linspace(1.0, 5.0, dim)
    a_mat = q @ np.diag(diag) @ q.T
    b_vec = rng.standard_normal(dim)
    x_true = np.linalg.solve(a_mat, b_vec)

    # Target residual G(x) = A x - b = 0 => Phi(x) = x - (A x - b)
    def phi_fn(x: np.ndarray) -> np.ndarray:
        return np.asarray(x - (a_mat @ x - b_vec), dtype=np.float64)

    solver = LimitedMemoryBroydenSolver(history_size=10, alpha=0.5)
    x = np.zeros(dim)
    for _ in range(60):
        phi_x = phi_fn(x)
        g = phi_x - x
        if np.linalg.norm(g) < 1e-11:
            break
        x, _ = solver.step(x, phi_x)

    assert np.allclose(x, x_true, atol=1e-8)

    # Check that default alpha is 0.5
    default_solver = LimitedMemoryBroydenSolver()
    assert default_solver.alpha == 0.5
    assert default_solver.history_size == 10
    assert default_solver.tol_singular == 1e-12

    # Robustness under 30 independent perturbed trials (Contract H9)
    success_count = 0
    for trial in range(30):
        t_rng = np.random.default_rng(100 + trial)
        x_init = x_true + t_rng.uniform(-0.35, 0.35, size=dim)
        s_trial = LimitedMemoryBroydenSolver(history_size=10, alpha=0.5)
        curr_x = x_init.copy()
        for _ in range(50):
            p_val = phi_fn(curr_x)
            res_norm = float(np.linalg.norm(p_val - curr_x))
            if res_norm < 1e-9:
                success_count += 1
                break
            curr_x, _ = s_trial.step(curr_x, p_val)

    # 100% convergence across all 30 perturbed trials
    assert success_count == 30


def test_add_symmetric_perturbation_invariance() -> None:
    """Verify that add_symmetric_perturbation preserves C3 bond permutation and sublattice equivalence."""
    cell = make_honeycomb_2site()
    qsl = make_qsl_seed(cell)
    rng = np.random.default_rng(123)

    pert_state = add_symmetric_perturbation(qsl, delta=0.05, rng=rng)

    # Check sublattice equivalence on on-site Majoranas
    assert np.allclose(pert_state.V[0], pert_state.V[1], atol=1e-15)
    # Check on-site skew-symmetry
    assert np.allclose(pert_state.V[0], -pert_state.V[0].T, atol=1e-15)
    # Check lambda equivalence and C3 isotropy
    assert np.allclose(pert_state.lambdas[0], pert_state.lambdas[1], atol=1e-15)
    assert np.isclose(pert_state.lambdas[0, 0], pert_state.lambdas[0, 1], atol=1e-15)
    assert np.isclose(pert_state.lambdas[0, 1], pert_state.lambdas[0, 2], atol=1e-15)
