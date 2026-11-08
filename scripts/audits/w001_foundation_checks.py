"""Collect finite W001 foundation diagnostics without repairing legacy services.

Exit success means observations were collected, not that the audited services
satisfy the reference contracts. The JSON records mismatches explicitly.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

from kitaev_mf.services.exact_perturbation import build_spin_zeeman_perturbation
from kitaev_mf.services.feshbach import compute_spin_feshbach_thresholds

ROOT = Path(__file__).resolve().parents[2]


def sector_minimum(H: np.ndarray, P: np.ndarray) -> float:
    """Restrict to Ran(P), retaining genuine zero eigenvalues."""
    vals, vecs = np.linalg.eigh(P)
    basis = vecs[:, vals > 0.5]
    return float(np.linalg.eigvalsh(basis.conj().T @ H @ basis)[0])


def observe_legacy_thresholds(H: np.ndarray, W: np.ndarray) -> dict[str, Any]:
    try:
        result = compute_spin_feshbach_thresholds(H, W)
    except (ValueError, ZeroDivisionError) as exc:
        return {"outcome": "exception", "exception": type(exc).__name__, "message": str(exc)}
    return {
        "outcome": "returned",
        "E0": result.E0,
        "E_Q_min": result.E_Q_min,
        "delta_flux": result.delta_flux,
    }


def site_operator(op: np.ndarray, site: int, size: int) -> np.ndarray:
    result: np.ndarray = np.array([[1.0]], dtype=np.complex128)
    for current in reversed(range(size)):
        result = np.kron(result, op if current == site else np.eye(2))
    return result


def collect() -> dict[str, Any]:
    X = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    Y = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    Z = np.diag([1, -1]).astype(np.complex128)
    assert np.array_equal(X @ Y - Y @ X, 2j * Z)
    H = np.diag([-3.0, -1.0])
    W = np.diag([1.0, -1.0])
    P = (np.eye(2) + W) / 2
    shifts = []
    for shift in (0.0, 3.0):
        shifted = H + shift * np.eye(2)
        ref0 = sector_minimum(shifted, P)
        refq = sector_minimum(shifted, np.eye(2) - P)
        shifts.append(
            {
                "shift": shift,
                "reference_E0": ref0,
                "reference_E_Q_min": refq,
                "reference_gap": refq - ref0,
                "legacy": observe_legacy_thresholds(shifted, W),
            }
        )
    H4 = np.diag([-4.0, -3.0, -2.0, -1.0])
    W1 = np.diag([1.0, 1.0, -1.0, -1.0])
    W2 = np.diag([1.0, -1.0, 1.0, -1.0])
    joint = (np.eye(4) + W1) @ (np.eye(4) + W2) / 4
    joint0 = sector_minimum(H4, joint)
    jointq = sector_minimum(H4, np.eye(4) - joint)
    observed_y = build_spin_zeeman_perturbation(1, (0.0, 1.0, 0.0))
    pauli = {"x": X, "y": Y, "z": Z}
    size = 6
    bonds = [(0, 1, "z"), (1, 2, "x"), (2, 3, "y"), (3, 4, "z"), (4, 5, "x"), (5, 0, "y")]
    H6 = np.zeros((2**size, 2**size), dtype=np.complex128)
    for i, j, kind in bonds:
        H6 -= site_operator(pauli[kind], i, size) @ site_operator(pauli[kind], j, size)
    V6 = np.zeros_like(H6)
    for site in range(size):
        for op in pauli.values():
            V6 -= 0.02 * site_operator(op, site, size)
    energies, vectors = np.linalg.eigh(H6)
    E0 = float(energies[0])
    ground = np.abs(energies - E0) < 1e-10
    frame = vectors[:, ground]
    P6 = frame @ frame.conj().T
    excited = vectors[:, ~ground]
    denominators = E0 - energies[~ground]
    R6 = (excited * (1.0 / denominators)) @ excited.conj().T
    second = frame.conj().T @ V6 @ R6 @ V6 @ frame
    scalar = np.trace(second) / second.shape[0]
    source_hashes = {}
    for name in (
        "src/kitaev_mf/services/feshbach.py",
        "src/kitaev_mf/services/exact_perturbation.py",
        "scripts/audits/w001_foundation_checks.py",
    ):
        source_hashes[name] = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
    return {
        "scope": "finite diagnostic observations; no production repair or bulk certification",
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "git_dirty": bool(
            subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip()
        ),
        "source_sha256": source_hashes,
        "numpy_version": np.__version__,
        "machine_epsilon": float(np.finfo(float).eps),
        "energy_shift": shifts,
        "joint_projector": {
            "reference_gap": jointq - joint0,
            "single_W_legacy": observe_legacy_thresholds(H4, W1),
        },
        "sector_order": {
            "fixed_P_E0": -1.0,
            "fixed_Q_E_min": -3.0,
            "legacy": observe_legacy_thresholds(H, -W),
        },
        "pauli_Y": {
            "matches_standard_Zeeman_minus_Y": bool(np.array_equal(observed_y, -Y)),
            "equals_plus_Y": bool(np.array_equal(observed_y, Y)),
            "difference_frobenius_norm": float(np.linalg.norm(observed_y + Y)),
        },
        "hexagon_second_order": {
            "geometry": "isolated 6-spin ring, six typed bonds; not a periodic bulk torus",
            "field_components": [0.02, 0.02, 0.02],
            "E0": E0,
            "ground_rank": int(np.count_nonzero(ground)),
            "gap_above_ground": float(energies[~ground][0] - E0),
            "PVP_norm": float(np.linalg.norm(P6 @ V6 @ P6)),
            "second_order_norm": float(np.linalg.norm(second)),
            "second_order_scalar_residual": float(
                np.linalg.norm(second - scalar * np.eye(second.shape[0]))
            ),
            "second_order_eigenvalues": np.linalg.eigvalsh(second).tolist(),
            "projector_residual": float(np.linalg.norm(P6 @ P6 - P6)),
        },
    }


def main() -> None:
    output = ROOT / "build/rnd/exact-perturbation-theory/w001_t4_observations.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    observations = collect()
    output.write_text(json.dumps(observations, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(observations, indent=2))


if __name__ == "__main__":
    main()
