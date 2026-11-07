"""Unit tests for Task T2.1.2: Plaquette Flux Representation and Bayesian Model Selection.

Verifies:
1. Zero-field exact equivalence Pf(Omega|_{partial p}) == W_p^fact == 1.000000.
2. Parlett-Reid skew elimination congruence Pf(M)^2 == det(M) to machine precision (< 1e-14).
3. Decisive Bayesian model selection (Delta_BIC_01 > 10, B_10 > 148) supporting quadratic
   trans-plaquette polarization M1 over linear relaxation M0.
4. Quantum fluctuation softening: Pf(Omega|_{partial p}) > W_p^fact for all 0 < h < h_c^*.
5. ACID dataset schema, metadata headers, and 100% solver convergence.
"""

from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

# Ensure repository root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np  # noqa: E402

from kitaev_mf.domain import ModelParameters, make_honeycomb_2site  # noqa: E402
from kitaev_mf.services.real_space import RealSpaceTorusCluster, compute_skew_pfaffian  # noqa: E402
from kitaev_mf.services.seeds import make_qsl_seed  # noqa: E402
from scripts.phase_2.m2_1.benchmark_plaquette_flux_bayesian import compute_bic_and_aicc  # noqa: E402

DATASET_PATH = Path("data/processed/benchmarks/m2_1/plaquette_flux_bayesian_benchmark.dat")


def test_t2_1_2_zero_field_exact_equivalence() -> None:
    """Verify that at h=0, factored Wilson loop and 12-Majorana Pfaffian agree to machine precision."""
    cluster = RealSpaceTorusCluster(L=4)
    cell2 = make_honeycomb_2site()
    seed = make_qsl_seed(cell2)
    params = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))

    omega, state, met = cluster.solve_real_space_purification(
        params=params,
        initial_state=seed,
        max_iter=30,
        tol=1e-8,
    )
    assert met["converged"]

    # Factored 6-link Wilson loop
    u_z = state.U[0, 3, 3]
    u_x = state.U[1, 1, 1]
    u_y = state.U[2, 2, 2]
    w_fact = float((u_x * u_y * u_z) ** 2)

    # Exact 12-Majorana Pfaffian
    p0 = (0, 0)
    w_pfaff = cluster.compute_plaquette_pfaffian_flux(omega, p0)

    # Machine-precision identity
    assert np.isclose(w_fact, 1.0, atol=1e-12), f"Expected W_fact=1.0, got {w_fact}"
    assert np.isclose(w_pfaff, 1.0, atol=1e-12), f"Expected W_pfaff=1.0, got {w_pfaff}"
    assert abs(w_fact - w_pfaff) < 1e-12, (
        f"Zero-field discrepancy must vanish: {abs(w_fact - w_pfaff)}"
    )

    # Check Parlett-Reid congruence: Pf(M)^2 == det(M)
    indices = cluster.compute_plaquette_gauge_majoranas(p0)
    sub_omega = omega[np.ix_(indices, indices)]
    pf_val = compute_skew_pfaffian(sub_omega)
    det_val = float(np.linalg.det(sub_omega))
    assert np.isclose(pf_val**2, det_val, atol=1e-13), (
        f"Pf(M)^2 = {pf_val**2:.15e} must equal det(M) = {det_val:.15e}"
    )


def test_t2_1_2_bayesian_model_selection_quadratic_scaling() -> None:
    """Verify decisive Bayesian model selection (Delta_BIC > 10) supporting quadratic model M_1."""
    assert DATASET_PATH.exists(), f"Benchmark dataset not found: {DATASET_PATH}"

    # Load dataset
    data = np.loadtxt(DATASET_PATH, comments="#")
    # Columns: h w_fact w_pfaff diff pf_gt_fact kappa_A idempotency_res lambda_norm conv iters time_sec
    h_vals = data[:, 0]
    diffs = data[:, 3]

    # Filter QSL phase (h <= 0.35)
    mask = h_vals <= 0.35
    h_qsl = h_vals[mask]
    diff_qsl = diffs[mask]
    n_samples = len(h_qsl)
    assert n_samples >= 10, f"Need at least 10 points for robust model selection, got {n_samples}"

    # Candidate Model M_1: delta W_p = c_W * h^2
    h2 = h_qsl**2
    c_W = float(np.sum(diff_qsl * h2) / np.sum(h2**2))
    fit_M1 = c_W * h2
    rss_1 = float(np.sum((diff_qsl - fit_M1) ** 2))
    bic_1, aicc_1 = compute_bic_and_aicc(rss_1, n_samples, n_params=1)

    # Null Model M_0: delta W_p = c_1 * h
    c_1 = float(np.sum(diff_qsl * h_qsl) / np.sum(h_qsl**2))
    fit_M0 = c_1 * h_qsl
    rss_0 = float(np.sum((diff_qsl - fit_M0) ** 2))
    bic_0, aicc_0 = compute_bic_and_aicc(rss_0, n_samples, n_params=1)

    delta_bic = bic_0 - bic_1
    delta_aicc = aicc_0 - aicc_1
    b10 = math.exp(min(0.5 * delta_bic, 500.0))

    # Pre-registered acceptance criteria:
    # 1. Decisive support: Delta_BIC > 10.0 and Bayes factor B_10 > 148
    assert delta_bic > 10.0, f"Expected decisive Delta_BIC > 10.0, got {delta_bic:.2f}"
    assert delta_aicc > 10.0, f"Expected decisive Delta_AICc > 10.0, got {delta_aicc:.2f}"
    assert b10 > 148.0, f"Expected Bayes factor B_10 > 148, got {b10:.2e}"

    # 2. Physics parameter check: c_W ~ 1.23 in ferromagnetic Kitaev
    assert 1.00 < c_W < 1.50, f"Expected c_W in [1.00, 1.50], got {c_W:.4f}"


def test_t2_1_2_quantum_fluctuation_softening() -> None:
    """Verify that quantum gauge fluctuations soften flux degradation: Pf > W_fact for 0 < h <= 0.35."""
    assert DATASET_PATH.exists()
    data = np.loadtxt(DATASET_PATH, comments="#")
    h_vals = data[:, 0]
    w_fact = data[:, 1]
    w_pfaff = data[:, 2]
    pf_gt_fact = data[:, 4]

    # In QSL phase for h > 0, Pfaffian flux MUST be strictly greater than factored link product
    mask = (h_vals > 0.0) & (h_vals <= 0.35)
    assert np.all(pf_gt_fact[mask] == 1), (
        "Quantum fluctuations must soften flux decay (Pf > W_fact) for all h in (0, 0.35]"
    )
    # Check that at h=0.20, flux retention is higher by at least 0.04
    idx_020 = np.argmin(np.abs(h_vals - 0.20))
    flux_gain = w_pfaff[idx_020] - w_fact[idx_020]
    assert flux_gain >= 0.04, f"Expected flux retention gain >= 0.04 at h=0.20, got {flux_gain:.4f}"


def test_t2_1_2_acid_provenance_and_dataset_integrity() -> None:
    """Verify ACID provenance header and numerical integrity of benchmark table."""
    assert DATASET_PATH.exists()
    assert DATASET_PATH.stat().st_size > 1000, "Dataset file must be non-empty (>1 kB)"

    content = DATASET_PATH.read_text(encoding="utf-8")
    assert "# ACID Provenance Header" in content
    assert "Git Commit:" in content
    assert "Geometry: Honeycomb Torus L=6" in content
    assert "Bayesian Selection" in content
    assert "Delta_BIC_01" in content
    assert "Timestamp:" in content

    data = np.loadtxt(DATASET_PATH, comments="#")
    assert len(data) == 19, f"Expected 19 field points, got {len(data)}"

    # All solver runs must be 100% converged
    conv_col = data[:, 8]
    assert np.all(conv_col == 1), "All solver runs in the benchmark must have converged (conv == 1)"

    # Idempotency residual:
    # At h=0 on 3 | L (L=6) torus, the 4 gapless Dirac zero modes at K and K' give
    # ||Omega^2 + I||_F / sqrt(2N) = sqrt(4) / 12 = 1/6 = 0.16667 exactly.
    # For all h > 0, the topological chiral gap eliminates all zero modes, yielding machine precision.
    idemp_col = data[:, 6]
    assert np.isclose(idemp_col[0], 1.0 / 6.0, atol=1e-4), (
        f"Expected exact 4-Dirac-mode zero-mode signature 1/6 at h=0, got {idemp_col[0]:.6f}"
    )
    assert np.all(idemp_col[1:] < 1e-12), (
        f"Max gapped idempotency residual exceeded: {np.max(idemp_col[1:]):.4e}"
    )


def test_t2_1_2_cli_entrypoint(tmp_path: Path) -> None:
    """Verify that the benchmark script CLI executes properly with custom flags."""
    out_file = tmp_path / "cli_test_benchmark.dat"
    cmd = [
        sys.executable,
        "scripts/phase_2/m2_1/benchmark_plaquette_flux_bayesian.py",
        "--L",
        "4",
        "--num-h",
        "4",
        "--h-max",
        "0.10",
        "--out",
        str(out_file),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    assert out_file.exists()
    assert "Bayesian Model Selection Results" in res.stdout
    assert "Delta_BIC_01" in res.stdout
