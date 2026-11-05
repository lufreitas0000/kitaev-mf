"""Unit tests for Tier-2 solvable limits and dimensional reductions.

Contracts:
    TST-RED-01: Decoupled 1D zigzag-chain limit
    TST-RED-02: Isolated-dimer flat-band saturation limit
    TST-ASY-01: High-field classical energy and polarization convergence
"""

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.infra.linalg_backend import TorchBackend
from kitaev_mf.services.correlators import compute_correlators_from_bloch
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.hamiltonian import build_bloch_hamiltonian
from kitaev_mf.services.observables import compute_sublattice_magnetizations


def test_tst_red_01_zigzag_chain_reduction(cpu_backend: TorchBackend) -> None:
    """TST-RED-01: Decoupled 1D zigzag chain limit (K_x = K_y = -1, K_z = 0).

    The 2D honeycomb lattice decomposes into decoupled 1D zigzag chains.
    With K_z = 0, there are 2 active bonds (x and y). The sum of matter hoppings
    is -(4 / pi), so each active bond receives:
        U_x^{00} = U_y^{00} = -2 / pi approx -0.636619772368.
    Decoupled z-bond correlators vanish identically: U_z == 0.
    """
    cell = make_honeycomb_2site()
    # K_x = -1.0, K_y = -1.0, K_z = 0.0
    params = ModelParameters(K=(-1.0, -1.0, 0.0))
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.4)

    # 48x48 Riemann mesh balances quadrature accuracy and execution budget (< 20ms)
    k_grid = cell.make_k_grid(48, 48)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)
    u_out, v_out = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)

    expected_u_1d = -2.0 / np.pi

    # In make_honeycomb_2site(): bond 0 is z-bond, bond 1 is x-bond, bond 2 is y-bond
    u_z = u_out[0]
    u_x = u_out[1]
    u_y = u_out[2]

    # Active x- and y-bonds: U_b^{00} == -2 / pi
    assert np.isclose(u_x[0, 0], expected_u_1d, atol=2e-3)
    assert np.isclose(u_y[0, 0], expected_u_1d, atol=2e-3)

    # Symmetry between x and y chains
    assert np.isclose(u_x[0, 0], u_y[0, 0], atol=1e-12)

    # Decoupled z-bond has vanishing hopping
    assert np.allclose(u_z, 0.0, atol=1e-14)

    # Gauge residuals and on-site correlators remain neutral
    assert np.allclose(v_out, 0.0, atol=1e-14)


def test_tst_red_02_isolated_dimer_flat_bands(cpu_backend: TorchBackend) -> None:
    """TST-RED-02: Isolated dimer flat-band limit (K_x = K_y = 0, K_z = -1).

    The lattice decomposes into independent two-site dimers on z-bonds.
    1. Flavors c^0 (matter) and c^3 (z-gauge) form flat dimer bands at +/- 0.5.
       Flavors c^1 and c^2 (x- and y-gauge) are completely decoupled at E == 0.
       Therefore, the 4 non-negative bands are [0.0, 0.0, 0.5, 0.5].
    2. Matter hopping saturates to unity on the dimer: U_z^{00} = -1.0, U_x = U_y == 0.
    3. Energy per cell saturates to E_cell = -|K_z| = -1.0 (E_site = -0.5).
    """
    cell = make_honeycomb_2site()
    # K_x = 0.0, K_y = 0.0, K_z = -1.0
    params = ModelParameters(K=(0.0, 0.0, -1.0))
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-1.0)

    k_grid = cell.make_k_grid(12, 12)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)
    eigenvalues, _ = cpu_backend.batched_eigh(h_batch)
    e_np = eigenvalues.detach().cpu().numpy()

    # The 4 positive/zero bands: 2 decoupled gauge zero modes at 0.0, 2 dimer bands at 0.5
    pos_energies = e_np[:, 4:]
    assert np.allclose(pos_energies[:, :2], 0.0, atol=1e-14)
    assert np.allclose(pos_energies[:, 2:], 0.5, atol=1e-14)

    u_out, v_out = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)
    u_z = u_out[0]
    u_x = u_out[1]
    u_y = u_out[2]

    # Dimer z-bond saturates
    assert np.isclose(u_z[0, 0], -1.0, atol=1e-14)
    assert np.isclose(u_z[3, 3], 1.0, atol=1e-14)

    # Decoupled x- and y-bonds vanish
    assert np.allclose(u_x, 0.0, atol=1e-14)
    assert np.allclose(u_y, 0.0, atol=1e-14)

    # Energetics check
    state_saturated = MeanFieldState(U=u_out, V=v_out, lambdas=state.lambdas)
    energy = compute_gaussian_energy(state_saturated, cell, params)
    assert np.isclose(energy.per_cell, -1.0, atol=1e-14)
    assert np.isclose(energy.per_site, -0.5, atol=1e-14)


def test_tst_asy_01_high_field_classical_energy(cpu_backend: TorchBackend) -> None:
    """TST-ASY-01: High-field classical energy and polarization convergence.

    Under dominant Zeeman field along [1, 1, 1], spins align with the field:
    ||<sigma_s>||_2 -> 1.0, and energy density asymptotically approaches -h.
    """
    cell = make_honeycomb_2site()
    h_mag = 50.0
    h_unit = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
    h_vec = (float(h_mag * h_unit[0]), float(h_mag * h_unit[1]), float(h_mag * h_unit[2]))

    params = ModelParameters(K=-1.0, J=0.05, h=h_vec)
    state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.3)

    k_grid = cell.make_k_grid(8, 8)
    h_batch = build_bloch_hamiltonian(k_grid, cell, params, state, cpu_backend)
    u_out, v_out = compute_correlators_from_bloch(h_batch, cell, k_grid, cpu_backend)

    state_polarized = MeanFieldState(U=u_out, V=v_out, lambdas=state.lambdas)
    sigma_mag, _ = compute_sublattice_magnetizations(state_polarized, cell)

    for s in range(cell.num_sublattices):
        # Direction aligns with h
        proj = float(np.dot(sigma_mag[s], h_unit))
        norm = float(np.linalg.norm(sigma_mag[s]))
        # Polarization magnitude approaches saturation (within O(K/h) ~ 2%)
        assert norm > 0.98
        assert np.isclose(proj, norm, atol=1e-3)

    energy = compute_gaussian_energy(state_polarized, cell, params)
    # Energy per site E_site / h -> -1.0
    relative_energy = energy.per_site / h_mag
    assert np.isclose(relative_energy, -1.0, atol=0.03)
