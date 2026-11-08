# ruff: noqa: E402
# Repository-path bootstrap must precede local imports in this standalone entrypoint.
"""Unit test suite for Phase 2 Milestone M2.1 Task T2.1.3: Variational Projection Bounds.

Verifies:
- Exact spin-space diagonalization on N=6 hexagon (E=-4.0) and N=8 torus (E=-4*sqrt(3)).
- Rayleigh-Ritz variational inequality: E_proj(h) >= E_exact(h).
- 2-RDM semidefinite programming lower bound: E_SDP(h) <= E_exact(h).
- Strict variational bracket: E_SDP(h) <= E_exact(h) <= E_proj(h) for all h in [0.0, 0.40].
- Gutzwiller norm overlap exactness N_G(0) = 2^-N at zero field.
- Bayesian model selection on projection decay.
- ACID provenance header and numerical integrity of exported benchmark dataset.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from kitaev_mf.services.phase_2.variational_bounds import (
    BayesianDecayFitter,
    ExactDiagonalizationKitaev,
    GutzwillerProjector,
    TwoRDM_SDP_Solver,
)

DATASET_PATH = (
    REPO_ROOT / "data" / "processed" / "benchmarks" / "m2_1" / "variational_bounds_benchmark.dat"
)


def test_t2_1_3_ed_exactness_and_commutation() -> None:
    """Verify analytical exactness of ED on 6-site hexagon and 8-site honeycomb torus."""
    # 1. 6-site isolated hexagon
    hex_bonds = [
        (0, 1, "z"),
        (1, 2, "x"),
        (2, 3, "y"),
        (3, 4, "z"),
        (4, 5, "x"),
        (5, 0, "y"),
    ]
    ed6 = ExactDiagonalizationKitaev(num_spins=6, bonds=hex_bonds, K=-1.0, h=(0.0, 0.0, 0.0))
    e6, _ = ed6.solve_ground_state()
    assert np.isclose(e6, -4.000000, atol=1e-10), (
        f"Expected E_exact = -4.000000 on hexagon, got {e6:.6f}"
    )

    # Verify field lowering of ground state energy
    ed6_field = ExactDiagonalizationKitaev(
        num_spins=6,
        bonds=hex_bonds,
        K=-1.0,
        h=(0.1 / np.sqrt(3.0), 0.1 / np.sqrt(3.0), 0.1 / np.sqrt(3.0)),
    )
    e6_f, _ = ed6_field.solve_ground_state()
    assert e6_f < e6, "Ground state energy must decrease monotonically with field"

    # 2. 8-site periodic honeycomb torus
    torus_bonds = [
        (0, 1, "z"),
        (0, 5, "x"),
        (0, 3, "y"),
        (4, 5, "z"),
        (4, 1, "x"),
        (4, 7, "y"),
        (2, 3, "z"),
        (2, 7, "x"),
        (2, 1, "y"),
        (6, 7, "z"),
        (6, 3, "x"),
        (6, 5, "y"),
    ]
    ed8 = ExactDiagonalizationKitaev(num_spins=8, bonds=torus_bonds, K=-1.0, h=(0.0, 0.0, 0.0))
    e8, _ = ed8.solve_ground_state()
    expected_e8 = -4.0 * np.sqrt(3.0)  # -6.92820323
    assert np.isclose(e8, expected_e8, atol=1e-6), (
        f"Expected analytical E_exact = -4*sqrt(3) = {expected_e8:.6f}, got {e8:.6f}"
    )


def test_t2_1_3_sdp_semidefinite_lower_bound() -> None:
    """Verify that convex 2-RDM SDP satisfies E_SDP <= E_exact rigorously."""
    hex_bonds = [
        (0, 1, "z"),
        (1, 2, "x"),
        (2, 3, "y"),
        (3, 4, "z"),
        (4, 5, "x"),
        (5, 0, "y"),
    ]
    for h_mag in (0.0, 0.10, 0.20):
        h_comp = float(h_mag / np.sqrt(3.0))
        h_vec = (h_comp, h_comp, h_comp)
        ed = ExactDiagonalizationKitaev(num_spins=6, bonds=hex_bonds, K=-1.0, h=h_vec)
        e_exact, _ = ed.solve_ground_state()

        sdp = TwoRDM_SDP_Solver(num_spins=6, bonds=hex_bonds, K=-1.0, h=h_vec)
        e_sdp, met = sdp.solve(max_iter=200, tol=1e-5)

        assert met["psd_min_eigenvalue"] >= -1e-10, (
            "Moment matrix must remain positive semidefinite"
        )
        assert e_sdp <= e_exact + 1e-10, (
            f"SDP lower bound violated at h={h_mag}: E_SDP={e_sdp:.6f} > E_exact={e_exact:.6f}"
        )


def test_t2_1_3_gutzwiller_rayleigh_ritz_bound() -> None:
    """Verify that Gutzwiller projected state satisfies E_proj >= E_exact unconditionally."""
    hex_bonds = [
        (0, 1, "z"),
        (1, 2, "x"),
        (2, 3, "y"),
        (3, 4, "z"),
        (4, 5, "x"),
        (5, 0, "y"),
    ]
    ed = ExactDiagonalizationKitaev(num_spins=6, bonds=hex_bonds, K=-1.0, h=(0.0, 0.0, 0.0))
    e_exact, _ = ed.solve_ground_state()

    # Any physical state evaluated through Gutzwiller projector must satisfy Rayleigh-Ritz
    proj6 = GutzwillerProjector(num_spins=6)
    # Test random skew covariance matrix
    np.random.seed(42)
    R = np.random.randn(24, 24)
    A = R - R.T
    w, V = np.linalg.eigh(1j * A)
    omega = np.real(-1j * (V * np.sign(w)) @ V.conj().T)

    e_proj = proj6.evaluate_projected_energy(omega, ed)
    assert e_proj >= e_exact - 1e-10, (
        f"Rayleigh-Ritz bound violated: E_proj={e_proj:.6f} < E_exact={e_exact:.6f}"
    )


def test_t2_1_3_zero_field_gutzwiller_norm_exactness() -> None:
    """Verify exact theoretical zero-field Gutzwiller overlap N_G(0) on Kitaev hexagon and torus."""
    # 1. N=6 Hexagon: Boundary spectral zero modes eliminate the total parity monomial,
    # giving strictly S = empty (1.0), so N_G = 1/64 = 0.015625 exactly.
    proj6 = GutzwillerProjector(num_spins=6)
    bonds_num = [(0, 1, 3), (1, 2, 1), (2, 3, 2), (3, 4, 3), (4, 5, 1), (5, 0, 2)]
    A = np.zeros((24, 24), dtype=np.float64)
    for i, j, gamma in bonds_num:
        A[4 * i, 4 * j] += 1.0
        A[4 * j, 4 * i] -= 1.0
        A[4 * i + gamma, 4 * j + gamma] -= 1.0
        A[4 * j + gamma, 4 * i + gamma] += 1.0

    w, V = np.linalg.eigh(1j * A + 1e-7 * np.eye(24))
    omega0_6 = np.real(-1j * (V * np.sign(w)) @ V.conj().T)

    expected_ng6 = 1.0 / 64.0  # 2^-6
    ng0_6 = proj6.compute_norm_overlap(omega0_6)
    assert np.isclose(ng0_6, expected_ng6, atol=1e-8), (
        f"Expected N_G(0) = 2^-6 = {expected_ng6:.6f}, got {ng0_6:.6f}"
    )

    # 2. N=8 Honeycomb Torus: Fully gapped spectrum with PBC where both S = empty (1.0)
    # and S = all (1.0) contribute equally, giving N_G = 2 / 256 = 1 / 128 = 0.0078125 exactly.
    from kitaev_mf.domain.lattice import make_honeycomb_2site
    from kitaev_mf.domain.parameters import ModelParameters
    from kitaev_mf.services.real_space import RealSpaceTorusCluster
    from kitaev_mf.services.seeds import make_qsl_seed

    cluster8 = RealSpaceTorusCluster(L=2)
    cell8 = make_honeycomb_2site()
    seed8 = make_qsl_seed(cell8)
    params = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))
    omega8, _, met8 = cluster8.solve_real_space_purification(params, seed8, max_iter=20)
    assert met8["converged"], "Torus purification failed at h=0"

    proj8 = GutzwillerProjector(num_spins=8)
    ng0_8 = proj8.compute_norm_overlap(omega8)
    expected_ng8 = 1.0 / 128.0  # 2 / 256
    assert np.isclose(ng0_8, expected_ng8, atol=1e-8), (
        f"Expected N_G(0) = 1/128 = {expected_ng8:.6f}, got {ng0_8:.6f}"
    )


def test_t2_1_3_dataset_provenance_and_integrity() -> None:
    """Verify ACID provenance header and strict variational bracket in exported benchmark table."""
    assert DATASET_PATH.exists(), f"Benchmark dataset not found: {DATASET_PATH}"
    assert DATASET_PATH.stat().st_size > 1000, "Dataset file must be non-empty"

    content = DATASET_PATH.read_text(encoding="utf-8")
    assert "# ACID Provenance Header" in content
    assert "Git Commit:" in content
    assert "Variational Bracket Guarantee" in content
    assert "Timestamp:" in content

    data = np.loadtxt(DATASET_PATH, comments="#")
    assert len(data) == 17, f"Expected 17 field points in sweep, got {len(data)}"

    _h = data[:, 0]
    e_sdp_6 = data[:, 1]
    e_exact_6 = data[:, 2]
    e_proj_6 = data[:, 3]
    bracket_ok_6 = data[:, 8]

    # Every single point must satisfy E_SDP <= E_exact <= E_proj
    assert np.all(bracket_ok_6 == 1), (
        "Variational bracket must be strictly valid at all field points"
    )
    assert np.all(e_sdp_6 <= e_exact_6 + 1e-10), "SDP lower bound violated in benchmark"
    assert np.all(e_proj_6 >= e_exact_6 - 1e-10), "Rayleigh-Ritz upper bound violated in benchmark"

    # Energy monotonicity with field:
    assert np.all(np.diff(e_exact_6) < 0), (
        "Exact ground state energy must strictly decrease with field"
    )


def test_t2_1_3_bayesian_model_selection_decay() -> None:
    """Verify Bayesian model selection comparing quadratic decay M1 vs linear null M0."""
    data = np.loadtxt(DATASET_PATH, comments="#")
    h = data[:, 0]
    ng8 = data[:, 7]

    # Fit on un-clamped physical QSL regime (h <= 0.35, N_G > 1e-10)
    mask = (h <= 0.35 + 1e-5) & (ng8 > 1e-10)
    bayes = BayesianDecayFitter.fit_and_compare(h[mask], ng8[mask], num_sites=8)

    assert bayes["num_points"] >= 10
    assert bayes["kappa_fit"] > 0, "Decay rate kappa must be positive"
    assert bayes["rss1_log"] < bayes["rss0_log"], (
        f"Candidate quadratic model M1 must have lower RSS than linear M0 (RSS1={bayes['rss1_log']:.4f}, RSS0={bayes['rss0_log']:.4f})"
    )
    assert bayes["delta_bic_01"] > 10.0, (
        f"Expected decisive Delta_BIC_01 > 10.0, got {bayes['delta_bic_01']:.2f}"
    )
    assert bayes["bayes_factor_10"] > 148.0, (
        f"Expected Bayes factor > 148.0, got {bayes['bayes_factor_10']:.2e}"
    )
    assert bayes["candidate_favored"], "Candidate quadratic decay model must be decisively favored"
