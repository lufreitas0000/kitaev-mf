"""Test honeycomb lattice and unit cell geometries."""
import numpy as np
from kitaev_mf.domain.lattice import make_honeycomb_2site, make_honeycomb_4site_rectangular


def test_2site_lattice_geometry():
    cell = make_honeycomb_2site()
    assert cell.num_sublattices == 2
    assert len(cell.bonds) == 3
    
    b1, b2 = cell.reciprocal_vectors
    assert np.isclose(np.dot(cell.a1, b1), 2.0 * np.pi)
    assert np.isclose(np.dot(cell.a1, b2), 0.0)
    assert np.isclose(np.dot(cell.a2, b1), 0.0)
    assert np.isclose(np.dot(cell.a2, b2), 2.0 * np.pi)


def test_4site_rectangular_lattice_geometry():
    cell = make_honeycomb_4site_rectangular()
    assert cell.num_sublattices == 4
    # 4 sites with coordination z=3 implies (4 * 3) / 2 = 6 bonds
    assert len(cell.bonds) == 6
    
    # Check bond type counts: exactly 2 x-bonds, 2 y-bonds, and 2 z-bonds
    types = [b.bond_type for b in cell.bonds]
    assert types.count(1) == 2
    assert types.count(2) == 2
    assert types.count(3) == 2
    
    b1, b2 = cell.reciprocal_vectors
    assert np.isclose(np.dot(cell.a1, b1), 2.0 * np.pi)
    assert np.isclose(np.dot(cell.a2, b2), 2.0 * np.pi)
    assert np.isclose(np.dot(cell.a1, b2), 0.0)
    assert np.isclose(np.dot(cell.a2, b1), 0.0)
