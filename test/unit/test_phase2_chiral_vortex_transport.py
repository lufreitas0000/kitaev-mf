"""Finite chiral-matter transport contracts, separate from finite-field mean field."""

import numpy as np
import pytest

from kitaev_mf.services.phase_2.chiral_vortex_transport import (
    assemble_chiral_bdg_generator,
    bond_endpoints,
    build_dual_graph,
    build_homology_braiding_loop,
    compute_vortex_braiding_holonomy,
    generate_vortex_pair_cut,
)


def test_generator_is_locally_gauge_covariant() -> None:
    """Independent random site signs must conjugate NN and NNN terms alike."""
    length = 5
    rng = np.random.default_rng(20261108)
    gauge = rng.choice([-1.0, 1.0], 2 * length * length)
    _, incidence = build_dual_graph(length)
    signs = {bond: float(rng.choice([-1, 1])) for bond in incidence}
    transformed = {
        bond: sign
        * gauge[bond_endpoints(length, *bond)[0]]
        * gauge[bond_endpoints(length, *bond)[1]]
        for bond, sign in signs.items()
    }
    original = assemble_chiral_bdg_generator(length, bond_signs=signs)
    changed = assemble_chiral_bdg_generator(length, bond_signs=transformed)
    # Only sign flips and additions occur, so these identities are exact here.
    assert np.array_equal(original, -original.T)
    assert np.array_equal(changed, gauge[:, None] * original * gauge[None, :])


def test_cut_has_exactly_two_flux_endpoints() -> None:
    length = 5
    first, second = (0, 0), (2, 3)
    cut = set(generate_vortex_pair_cut(length, first, second))
    _, incidence = build_dual_graph(length)
    flux = {p: 1 for pair in incidence.values() for p in pair}
    for bond in cut:
        for plaquette in incidence[bond]:
            flux[plaquette] *= -1
    assert {p for p, value in flux.items() if value < 0} == {first, second}


@pytest.mark.parametrize("fixed, expected", [((4, 4), -1), ((0, 0), 1)])
def test_gauge_closed_transport_and_refinement(fixed: tuple[int, int], expected: int) -> None:
    coarse = compute_vortex_braiding_holonomy(p_fixed=fixed, n_interp_steps=6)
    fine = compute_vortex_braiding_holonomy(p_fixed=fixed, n_interp_steps=10)
    for result in (coarse, fine):
        assert result.enclosing == (expected == -1)
        assert result.gauge_closure_error == 0
        assert result.min_subspace_gap > 0.5
        assert result.min_overlap_singular_value > 0.7
        # Finite model-specific numerical criterion, not a universal braid theorem.
        assert np.allclose(result.holonomy_matrix, expected * np.eye(2), atol=1e-2, rtol=0)
        assert np.allclose(
            result.holonomy_matrix.conj().T @ result.holonomy_matrix,
            np.eye(2),
            atol=1e-12,
            rtol=0,
        )
    assert np.allclose(coarse.holonomy_matrix, fine.holonomy_matrix, atol=1e-3, rtol=0)


def test_transport_rejects_invalid_or_unresolved_paths() -> None:
    with pytest.raises(ValueError, match="interpolation"):
        compute_vortex_braiding_holonomy(n_interp_steps=0)
    with pytest.raises(ValueError, match="nonzero kappa"):
        compute_vortex_braiding_holonomy(kappa=0)
    with pytest.raises(ValueError, match="both"):
        compute_vortex_braiding_holonomy(custom_loop_bonds=[])
    with pytest.raises(ValueError, match="refine"):
        compute_vortex_braiding_holonomy(n_interp_steps=2, min_overlap=0.999999999)
    bonds, sites, start = build_homology_braiding_loop(10, (4, 4))
    with pytest.raises(ValueError, match="Stationary vortex"):
        compute_vortex_braiding_holonomy(p_fixed=start)
    with pytest.raises(ValueError, match="ordered closed loop"):
        compute_vortex_braiding_holonomy(custom_loop_bonds=bonds[:-1], custom_enclosed_sites=sites)
    with pytest.raises(ValueError, match="gauge boundary"):
        compute_vortex_braiding_holonomy(custom_loop_bonds=bonds, custom_enclosed_sites=frozenset())
