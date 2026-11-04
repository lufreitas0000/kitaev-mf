"""TDD contract test suite for production observables and Wilson loops (M1.4.4)."""

import numpy as np

from kitaev_mf.domain.constants import N_MATRICES
from kitaev_mf.domain.lattice import (
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.observables import (
    compute_bond_spin_correlators,
    compute_plaquette_fluxes,
    compute_sublattice_magnetizations,
)


def test_sublattice_magnetization_pure_kitaev_and_field() -> None:
    """Verify that pure Kitaev state has vanishing magnetization, while Zeeman fields polarize."""
    cell = make_honeycomb_2site()
    state_qsl = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5248657)

    sigma_mag, spin_mag = compute_sublattice_magnetizations(state_qsl, cell)
    assert sigma_mag.shape == (2, 3)
    assert spin_mag.shape == (2, 3)
    assert np.allclose(sigma_mag, 0.0, atol=1e-14)
    assert np.allclose(spin_mag, 0.0, atol=1e-14)

    # State with finite on-site magnetization V_s^{0, 1} = 0.4
    v_pol = np.zeros((2, 4, 4), dtype=np.float64)
    v_pol[0, 0, 1] = 0.4
    v_pol[0, 1, 0] = -0.4
    state_pol = MeanFieldState(U=state_qsl.U, V=v_pol, lambdas=state_qsl.lambdas)

    sigma_pol, spin_pol = compute_sublattice_magnetizations(state_pol, cell)
    assert np.isclose(sigma_pol[0, 0], 0.4, atol=1e-14)
    assert np.isclose(spin_pol[0, 0], 0.2, atol=1e-14)
    assert np.isclose(sigma_pol[1, 0], 0.0, atol=1e-14)


def test_bond_spin_correlators_pure_kitaev() -> None:
    """Verify pure Kitaev bond correlators reproduce <sigma_i^gamma sigma_j^gamma> = U_b^{00} * U_b^{gg}."""
    cell = make_honeycomb_2site()
    u_matter_exact = -0.5248657458506
    state_qsl = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=u_matter_exact)

    corrs = compute_bond_spin_correlators(state_qsl, cell)
    assert corrs.shape == (3, 3, 3)

    for b_idx, bond in enumerate(cell.bonds):
        gamma_idx = bond.bond_type - 1
        # Parallel spin correlation on bond gamma
        assert np.isclose(corrs[b_idx, gamma_idx, gamma_idx], -u_matter_exact, atol=1e-12)
        # All off-diagonal spin components vanish
        for a in range(3):
            for b in range(3):
                if not (a == gamma_idx and b == gamma_idx):
                    assert np.isclose(corrs[b_idx, a, b], 0.0, atol=1e-14)


def test_plaquette_fluxes_zero_flux_ground_state() -> None:
    """Verify Wilson loop expectation values W_p = +1.0 across 2-site and 4-site honeycomb cells."""
    cell_2 = make_honeycomb_2site()
    state_2 = MeanFieldState.from_pure_kitaev(cell_2, u_gauge=1.0, u_matter=-0.52)
    fluxes_2 = compute_plaquette_fluxes(state_2, cell_2)

    assert fluxes_2.shape == (1,)
    assert np.allclose(fluxes_2, 1.0, atol=1e-14)

    cell_4 = make_honeycomb_4site_rectangular()
    state_4 = MeanFieldState.from_pure_kitaev(cell_4, u_gauge=1.0, u_matter=-0.52)
    fluxes_4 = compute_plaquette_fluxes(state_4, cell_4)

    assert fluxes_4.shape == (2,)
    assert np.allclose(fluxes_4, 1.0, atol=1e-14)


def test_production_observables_match_fock_oracle() -> None:
    """Verify production observable service agrees with independent Wick contractions."""
    cell = make_honeycomb_2site()
    rng = np.random.default_rng(2026)
    u_rand = rng.standard_normal((3, 4, 4))
    v_raw = rng.standard_normal((2, 4, 4))
    v_skew = 0.5 * (v_raw - np.swapaxes(v_raw, -1, -2))
    state = MeanFieldState(U=u_rand, V=v_skew, lambdas=np.zeros((2, 3)))

    sigma_mag, _ = compute_sublattice_magnetizations(state, cell)
    for s in range(2):
        for a in (1, 2, 3):
            expected = float(0.5 * np.trace(state.V[s].T @ N_MATRICES[a]))
            assert np.isclose(sigma_mag[s, a - 1], expected, atol=1e-14)

    corrs = compute_bond_spin_correlators(state, cell)
    for b_idx in range(3):
        for a in (1, 2, 3):
            for b in (1, 2, 3):
                m_i = float(np.trace(state.V[0].T @ N_MATRICES[a]))
                m_j = float(np.trace(state.V[1].T @ N_MATRICES[b]))
                contract = float(
                    np.trace(N_MATRICES[a].T @ state.U[b_idx] @ N_MATRICES[b] @ state.U[b_idx].T)
                )
                expected_corr = 0.25 * m_i * m_j - 0.5 * contract
                assert np.isclose(corrs[b_idx, a - 1, b - 1], expected_corr, atol=1e-14)
