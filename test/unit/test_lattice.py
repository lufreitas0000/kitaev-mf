"""Test honeycomb lattice graph Lambda = (V, E) and unit cell geometries."""

import numpy as np

from kitaev_mf.domain.lattice import (
    DEFAULT_D_NN,
    make_honeycomb_2site,
    make_honeycomb_4site_rectangular,
)


def test_2site_lattice_geometry() -> None:
    cell = make_honeycomb_2site()
    assert cell.num_sublattices == 2
    assert cell.vertices == (0, 1)
    assert len(cell.edges) == 3

    # Verify unit lattice spacing |a1| = |a2| = 1.0
    assert np.isclose(np.linalg.norm(cell.a1), 1.0, atol=1e-15)
    assert np.isclose(np.linalg.norm(cell.a2), 1.0, atol=1e-15)

    # Verify all 3 nearest-neighbor bonds have length d_nn = 1/sqrt(3)
    for bond in cell.edges:
        vec = cell.bond_vector(bond)
        assert np.isclose(np.linalg.norm(vec), DEFAULT_D_NN, atol=1e-15)

    b1, b2 = cell.reciprocal_vectors
    assert np.isclose(np.dot(cell.a1, b1), 2.0 * np.pi, atol=1e-15)
    assert np.isclose(np.dot(cell.a1, b2), 0.0, atol=1e-15)
    assert np.isclose(np.dot(cell.a2, b1), 0.0, atol=1e-15)
    assert np.isclose(np.dot(cell.a2, b2), 2.0 * np.pi, atol=1e-15)


def test_4site_rectangular_lattice_geometry_and_scale_parity() -> None:
    cell_2 = make_honeycomb_2site()
    cell_4 = make_honeycomb_4site_rectangular()

    assert cell_4.num_sublattices == 4
    assert cell_4.vertices == (0, 1, 2, 3)
    assert len(cell_4.edges) == 6

    # Check bond type counts: 2 x-bonds, 2 y-bonds, 2 z-bonds
    types = [b.bond_type for b in cell_4.edges]
    assert types.count(1) == 2
    assert types.count(2) == 2
    assert types.count(3) == 2

    # Check supercell relation: A1 = a1 - a2 = (1, 0), A2 = a1 + a2 = (0, sqrt(3))
    assert np.allclose(cell_4.a1, cell_2.a1 - cell_2.a2, atol=1e-15)
    assert np.allclose(cell_4.a2, cell_2.a1 + cell_2.a2, atol=1e-15)

    # Every bond of type gamma in the 4-site cell must match the 2-site bond vector
    vec_by_type_2 = {b.bond_type: cell_2.bond_vector(b) for b in cell_2.edges}
    for bond in cell_4.edges:
        vec_4 = cell_4.bond_vector(bond)
        assert np.isclose(np.linalg.norm(vec_4), DEFAULT_D_NN, atol=1e-15)
        assert np.allclose(vec_4, vec_by_type_2[bond.bond_type], atol=1e-15)

    b1, b2 = cell_4.reciprocal_vectors
    assert np.isclose(np.dot(cell_4.a1, b1), 2.0 * np.pi, atol=1e-15)
    assert np.isclose(np.dot(cell_4.a2, b2), 2.0 * np.pi, atol=1e-15)
    assert np.isclose(np.dot(cell_4.a1, b2), 0.0, atol=1e-15)
    assert np.isclose(np.dot(cell_4.a2, b1), 0.0, atol=1e-15)
