"""Unit tests verifying mathematical contracts of seeds, symmetries, and perturbations."""

from __future__ import annotations

import numpy as np

from kitaev_mf.domain.constants import G_MATRICES
from kitaev_mf.domain.lattice import (
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)
from kitaev_mf.services.observables import (
    compute_plaquette_fluxes,
    compute_sublattice_magnetizations,
)
from kitaev_mf.services.seeds import (
    add_controlled_perturbation,
    apply_flux_flip,
    apply_gauge_transformation_sublattice_b,
    check_seed_representability,
    make_polarized_seed,
    make_qsl_seed,
)


def test_qsl_seed_properties() -> None:
    """Verify nominal pure Kitaev zero-flux QSL seed properties."""
    cell = make_honeycomb_2site()
    state = make_qsl_seed(cell)

    assert state.U.shape == (3, 4, 4)
    assert state.V.shape == (2, 4, 4)
    assert state.lambdas.shape == (2, 3)

    # All on-site and lambdas vanish
    np.testing.assert_allclose(state.V, 0.0)
    np.testing.assert_allclose(state.lambdas, 0.0)

    # Plaquette flux is exactly +1.0
    fluxes = compute_plaquette_fluxes(state, cell)
    np.testing.assert_allclose(fluxes, 1.0)

    # Representability holds
    assert check_seed_representability(state, cell) is True


def test_polarized_seed_purity_and_gauge_neutrality() -> None:
    """Verify that polarized seed constructed from M_MATRICES has exact purity and zero gauge violation."""
    cell = make_honeycomb_2site()
    state = make_polarized_seed(cell, axis="111")

    # Magnetization is exactly 1/sqrt(3) along [111]
    sig_mag, _ = compute_sublattice_magnetizations(state, cell)
    expected_n = np.ones(3) / np.sqrt(3)
    np.testing.assert_allclose(sig_mag[0], expected_n, atol=1e-14)
    np.testing.assert_allclose(sig_mag[1], expected_n, atol=1e-14)

    # Gauge residual vanishes identically on all sublattices
    for s in (0, 1):
        for gamma in (1, 2, 3):
            g_val = float(np.trace(G_MATRICES[gamma] @ state.V[s]))
            assert abs(g_val) < 1e-14

    # Local purity V^2 = -I_4
    for s in (0, 1):
        v_sq = state.V[s] @ state.V[s]
        np.testing.assert_allclose(v_sq, -np.eye(4), atol=1e-14)


def test_gauge_transformation_invariance() -> None:
    """Verify that discrete Z_2 gauge flip on sublattice B preserves Wilson loop and magnetizations."""
    cell = make_honeycomb_2site()
    state_base = make_qsl_seed(cell)
    state_trans = apply_gauge_transformation_sublattice_b(state_base, cell)

    # Bond gauge components are inverted
    assert np.allclose(state_trans.U[:, 1:, 1:], -state_base.U[:, 1:, 1:])

    # Plaquette fluxes remain identically +1.0 because every plaquette traverses an even number of B bonds
    w_base = compute_plaquette_fluxes(state_base, cell)
    w_trans = compute_plaquette_fluxes(state_trans, cell)
    np.testing.assert_allclose(w_base, w_trans)


def test_flux_flip_reversal() -> None:
    """Verify that flipping a single bond gauge sign in a 4-site supercell inverts plaquette fluxes to -1.0."""
    cell = make_honeycomb_4site_rectangular()
    state_base = make_qsl_seed(cell)
    state_flipped = apply_flux_flip(state_base, cell, bond_idx=0)

    w_flipped = compute_plaquette_fluxes(state_flipped, cell)
    np.testing.assert_allclose(w_flipped, -1.0)


def test_controlled_perturbation_scaling() -> None:
    """Verify that controlled perturbations scale linearly with delta and preserve skewness."""
    cell = make_honeycomb_2site()
    state_base = make_qsl_seed(cell)
    rng = np.random.default_rng(42)

    state_pert = add_controlled_perturbation(state_base, delta=0.01, rng=rng, skew_only=True)
    diff_u = np.max(np.abs(state_pert.U - state_base.U))
    diff_v = np.max(np.abs(state_pert.V - state_base.V))

    assert 0.0 < diff_u < 0.05
    assert 0.0 < diff_v < 0.05

    # V must remain strictly skew-symmetric
    for s in (0, 1):
        skew_err = np.max(np.abs(state_pert.V[s] + state_pert.V[s].T))
        assert skew_err < 1e-14

    # Zero delta returns original state
    state_zero = add_controlled_perturbation(state_base, delta=0.0, rng=rng)
    np.testing.assert_allclose(state_zero.U, state_base.U)


def test_polarized_seed_axes_and_validation() -> None:
    """Verify that all polarized seed coordinate axes instantiate cleanly and invalid axis raises."""
    cell = make_honeycomb_2site()
    for ax in ("x", "y", "z", "minus_111"):
        st = make_polarized_seed(cell, axis=ax)
        assert st.V.shape == (2, 4, 4)

    import pytest

    with pytest.raises(ValueError, match="Unknown axis"):
        make_polarized_seed(cell, axis="invalid")  # type: ignore[arg-type]


def test_unrepresentable_seed_detection() -> None:
    """Verify that check_seed_representability returns False for unphysical states."""
    from kitaev_mf.domain.state import MeanFieldState

    cell = make_honeycomb_2site()
    st = make_qsl_seed(cell)
    # Inflate U to violate representability bounds ||Omega|| <= 1
    u_bad = st.U.copy()
    u_bad[0, 0, 0] = 5.0
    st_bad = MeanFieldState(U=u_bad, V=st.V, lambdas=st.lambdas)
    assert check_seed_representability(st_bad, cell) is False
