"""Unit tests for Phase 4 Track C: Exact Diagonalization in the spin formalism."""

from kitaev_mf.services.phase_4.exact_diagonalization import (
    ExactDiagonalizationSpinSolver,
    HoneycombCluster,
)


def test_cluster_topology_and_bonds() -> None:
    """Verify that HoneycombCluster correctly generates periodic torus bonds and plaquettes."""
    c8 = HoneycombCluster.make_torus(2, 2)
    assert c8.num_spins == 8
    assert len(c8.bonds) == 12  # 3 bonds per cell * 4 cells = 12

    # Check that each site has coordination number 3 (1 x-bond, 1 y-bond, 1 z-bond)
    for site in range(8):
        bonds_on_site = [b for b in c8.bonds if b[0] == site or b[1] == site]
        assert len(bonds_on_site) == 3

    c16 = HoneycombCluster.make_torus(2, 4)
    assert c16.num_spins == 16
    assert len(c16.bonds) == 24
    assert len(c16.plaquettes) == 8


def test_sparse_hamiltonian_hermiticity() -> None:
    """Verify that the assembled sparse spin Hamiltonian is exactly Hermitian."""
    c8 = HoneycombCluster.make_torus(2, 2)
    solver = ExactDiagonalizationSpinSolver(c8, K=1.0, h=0.05)
    h_mat = solver.build_sparse_hamiltonian()

    diff = h_mat - h_mat.conj().T
    assert abs(diff).max() < 1e-14


def test_pure_kitaev_ground_state_energy_bounds() -> None:
    """Verify that pure Kitaev ED energy per spin satisfies known variational bounds."""
    c8 = HoneycombCluster.make_torus(2, 2)
    solver = ExactDiagonalizationSpinSolver(c8, K=1.0, h=0.0)
    obs = solver.compute_observables(num_states=2)

    # For K=1.0, ground state energy per site is negative, around -0.8 to -0.9
    assert -1.0 <= obs.energy_per_spin <= -0.7
    # Magnetization is exactly zero at h=0
    assert abs(obs.magnetization) < 1e-12


def test_weak_field_magnetization_and_energy_shift() -> None:
    """Verify that applying a field h > 0 induces non-zero magnetization and lowers ground-state energy."""
    c8 = HoneycombCluster.make_torus(2, 2)
    solver0 = ExactDiagonalizationSpinSolver(c8, K=1.0, h=0.0)
    obs0 = solver0.compute_observables(num_states=2)

    solver_field = ExactDiagonalizationSpinSolver(c8, K=1.0, h=0.05)
    obs_field = solver_field.compute_observables(num_states=2)

    # Energy with field is lower: E(h) < E(0)
    assert obs_field.energy_ground < obs0.energy_ground
    # Magnetization is positive along [111]
    assert obs_field.magnetization > 0.0
