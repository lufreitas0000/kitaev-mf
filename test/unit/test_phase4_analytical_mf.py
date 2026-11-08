"""Unit tests for Phase 4 symmetry reduction, Schrieffer-Wolff block-diagonalization, and symmetrized solver."""

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.phase_4.symmetrized_solver import solve_symmetrized_fixed_point
from kitaev_mf.services.phase_4.symmetry_reduction import (
    SymmetricParameters,
    block_diagonalize_dirac_point,
    build_dirac_bloch_hamiltonian,
    extract_symmetric_parameters,
    pack_symmetric_state,
)


def test_symmetric_parameter_packing_and_projection() -> None:
    """Verify round-trip vector packing and C3 parameter extraction."""
    sym = SymmetricParameters(
        lambda_val=0.0123,
        v0=0.0045,
        v1=0.0015,
        u0=-0.5248,
        uz=0.9998,
        uperp=-0.0001,
        ucross=-0.0002,
    )
    vec = sym.to_vector()
    sym_rec = SymmetricParameters.from_vector(vec)
    assert np.allclose(sym.to_vector(), sym_rec.to_vector())

    cell = make_honeycomb_2site()
    state = pack_symmetric_state(sym, cell)
    sym_extracted = extract_symmetric_parameters(state)

    assert np.isclose(sym.lambda_val, sym_extracted.lambda_val, atol=1e-12)
    assert np.isclose(sym.v0, sym_extracted.v0, atol=1e-12)
    assert np.isclose(sym.v1, sym_extracted.v1, atol=1e-12)
    assert np.isclose(sym.u0, sym_extracted.u0, atol=1e-12)
    assert np.isclose(sym.uz, sym_extracted.uz, atol=1e-12)
    assert np.isclose(sym.uperp, sym_extracted.uperp, atol=1e-12)
    assert np.isclose(sym.ucross, sym_extracted.ucross, atol=1e-12)


def test_pure_kitaev_hamiltonian_and_eigenvalues() -> None:
    """Verify that H_MF(K) at h=0 possesses exactly 2 zero modes and 6-fold degenerate flat gauge bands."""
    sym = SymmetricParameters.pure_kitaev()
    params = ModelParameters(K=-1.0)
    h8 = build_dirac_bloch_hamiltonian(params, sym)

    evals = np.linalg.eigvalsh(h8)
    # 2 zero modes at K
    assert np.isclose(evals[3], 0.0, atol=1e-12)
    assert np.isclose(evals[4], 0.0, atol=1e-12)

    # 6 gauge eigenvalues at +/- |u_matter| / 2
    u_gauge_expected = abs(sym.u0) / 2.0
    for idx in (0, 1, 2):
        assert np.isclose(evals[idx], -u_gauge_expected, atol=1e-5)
    for idx in (5, 6, 7):
        assert np.isclose(evals[idx], u_gauge_expected, atol=1e-5)


def test_schrieffer_wolff_matter_reduction() -> None:
    """Verify that Schrieffer-Wolff reduction accurately reproduces the 8x8 Dirac gap."""
    sym = SymmetricParameters(
        lambda_val=0.018288,
        v0=0.008716,
        v1=0.002905,
        u0=-0.524736,
        uz=0.999790,
        uperp=-7.38e-07,
        ucross=-0.000103,
    )
    inv_sqrt3 = 1.0 / np.sqrt(3.0)
    h_val = 0.02
    h_vec = (h_val * inv_sqrt3, h_val * inv_sqrt3, h_val * inv_sqrt3)
    params = ModelParameters(K=-1.0, h=h_vec)

    h8 = build_dirac_bloch_hamiltonian(params, sym)
    evals8 = np.linalg.eigvalsh(h8)
    gap8 = float(evals8[4] - evals8[3])

    W, h_mat, h_gauge = block_diagonalize_dirac_point(params, sym)
    evals_mat = np.linalg.eigvalsh(h_mat)
    gap_mat = float(evals_mat[1] - evals_mat[0])

    # Relative difference between full 8x8 gap and SW 2x2 reduced gap < 1%
    rel_err = abs(gap8 - gap_mat) / gap8
    assert rel_err < 0.01


def test_symmetrized_solver_convergence_and_cubic_gap() -> None:
    """Verify that the DIIS-accelerated symmetrized solver converges and recovers the certified gap."""
    h_val = 0.02
    inv_sqrt3 = 1.0 / np.sqrt(3.0)
    h_vec = (h_val * inv_sqrt3, h_val * inv_sqrt3, h_val * inv_sqrt3)
    params = ModelParameters(K=-1.0, h=h_vec)

    res = solve_symmetrized_fixed_point(params, nk=18, max_iter=25, tol=1e-7)
    assert res.metrics.converged is True
    # Shared SolverMetrics container contract; this does not accept the physics claim.
    assert isinstance(res.metrics.residual_history, tuple)
    assert isinstance(res.metrics.contraction_factors, tuple)
    assert isinstance(res.metrics.comm_residuals, tuple)
    # The certified benchmark value at h=0.02 is 2.130461e-05
    assert np.isclose(res.dirac_gap_mf, 2.130461e-05, rtol=0.01)
    # The ratio gap / h^3 should match ~ 2.66
    c_ratio = res.dirac_gap_mf / (h_val**3)
    assert 2.5 <= c_ratio <= 2.8
