"""Opt-in reproduction of the retired, failing four-flavor pi-holonomy guarantee.

Run explicitly with pytest; intentionally excluded from normal test_*.py discovery.
PI authorized retirement on 2026-11-08; assertions are unchanged, not xfailed.
"""

import numpy as np

from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.phase_2.vortex_braiding import VortexBraidingCalculator


def topological_params() -> ModelParameters:
    component = 0.05 / np.sqrt(3.0)
    return ModelParameters(K=-1.0, h=(component, component, component))


def test_t2_2_3_enclosing_loop_pi_holonomy() -> None:
    """Verify that winding around a vortex imparts topological pi statistical holonomy."""
    calc = VortexBraidingCalculator(L=6)
    res = calc.evaluate_braiding_holonomy(topological_params(), N_sub=4)

    # Both independent phases must be physically meaningful:
    dist_non_to_zero = min(res.gamma_non_enclosing, 2.0 * np.pi - res.gamma_non_enclosing)
    assert dist_non_to_zero < 0.05, f"gamma_non should be near 0, got {res.gamma_non_enclosing:.4f}"

    dist_enc_to_pi = min(abs(res.gamma_enclosing - np.pi), abs(res.gamma_enclosing - 3.0 * np.pi))
    assert dist_enc_to_pi < 0.50, f"gamma_enc should be near pi, got {res.gamma_enclosing:.4f}"

    assert res.relative_error_from_pi <= 0.18, (
        f"L=6: Holonomy error {res.relative_error_from_pi * 100:.1f}% exceeded 18% threshold "
        f"(Delta_gamma = {res.delta_gamma:.4f} rad)"
    )
