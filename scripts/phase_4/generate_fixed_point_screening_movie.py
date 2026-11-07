"""Execution script generating animated screening movies and frames for Phase 4 Track D.

Simulates the relaxation under the self-consistent map Phi for:
1. Single bond flip at h = 0.0
2. Single bond flip at h = 0.02
3. L/2 gauge string (maximally separated vortices) at h = 0.0
4. L/2 gauge string (maximally separated vortices) at h = 0.02
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure local package takes precedence
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.state import InhomogeneousMeanFieldState, MeanFieldState
from kitaev_mf.services.phase_4.fixed_point_visualizer import FixedPointScreeningVisualizer
from kitaev_mf.services.real_space import RealSpaceTorusCluster


def run_all_movie_scenarios(L: int = 8, max_steps: int = 12) -> None:
    """Execute all 4 movie scenarios on an L x L real-space cluster."""
    cell = make_honeycomb_2site()
    cluster = RealSpaceTorusCluster(L=L, d_nn=1.0 / np.sqrt(3.0))
    print(
        f"\n[Movie Generator] Initializing L={L} Torus ({2 * L * L} sites, M={8 * L * L} Majoranas)..."
    )

    # Pure Kitaev converged reference state
    hom_state = MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5248657)
    base_inhom = InhomogeneousMeanFieldState.from_homogeneous(hom_state, cluster.num_cells)

    base_output_dir = (
        Path(__file__).resolve().parents[2] / "data" / "processed" / "movies" / "phase_4"
    )

    # =========================================================================
    # Scenario 1: Single Bond Flip at h = 0.0
    # =========================================================================
    print("\n--- Running Scenario 1: Single Bond Flip (h = 0.0) ---")
    s1_dir = base_output_dir / "scenario_1_single_flip_h0"
    params_h0 = ModelParameters(K=-1.0, h=(0.0, 0.0, 0.0))
    vis_s1 = FixedPointScreeningVisualizer(cluster, params_h0, s1_dir)

    # Flip 1 z-bond near the center of the lattice: cell (L//2, L//2), bond 0 (z-bond)
    flipped_s1 = {(L // 2, L // 2, 0)}
    vis_s1.run_screening_movie(
        initial_state=base_inhom,
        flipped_bonds=flipped_s1,
        max_steps=max_steps,
        gif_filename="screening_single_flip_h0.gif",
        title_prefix="Single Bond Flip (h=0)",
    )

    # =========================================================================
    # Scenario 2: Single Bond Flip at h = 0.02
    # =========================================================================
    print("\n--- Running Scenario 2: Single Bond Flip (h = 0.02) ---")
    s2_dir = base_output_dir / "scenario_2_single_flip_h002"
    inv_sqrt3 = 1.0 / np.sqrt(3.0)
    h_vec = (0.02 * inv_sqrt3, 0.02 * inv_sqrt3, 0.02 * inv_sqrt3)
    params_h002 = ModelParameters(K=-1.0, h=h_vec)
    vis_s2 = FixedPointScreeningVisualizer(cluster, params_h002, s2_dir)

    vis_s2.run_screening_movie(
        initial_state=base_inhom,
        flipped_bonds=flipped_s1,
        max_steps=max_steps,
        gif_filename="screening_single_flip_h002.gif",
        title_prefix="Single Bond Flip (h=0.02)",
    )

    # =========================================================================
    # Scenario 3: L/2 Gauge String (Maximally Separated Vortices) at h = 0.0
    # =========================================================================
    print(f"\n--- Running Scenario 3: L/2 Gauge String (d = {L // 2}, h = 0.0) ---")
    s3_dir = base_output_dir / "scenario_3_string_h0"
    vis_s3 = FixedPointScreeningVisualizer(cluster, params_h0, s3_dir)

    # Cut along line m1 = L//2 for m2 in range(L//2)
    # Using z-bonds (b_idx=0) along the line creates two vortices separated by distance L//2
    flipped_s3 = {(L // 2, m2, 0) for m2 in range(L // 2)}
    vis_s3.run_screening_movie(
        initial_state=base_inhom,
        flipped_bonds=flipped_s3,
        max_steps=max_steps,
        gif_filename="screening_string_h0.gif",
        title_prefix=f"Gauge String d={L // 2} (h=0)",
    )

    # =========================================================================
    # Scenario 4: L/2 Gauge String at h = 0.02
    # =========================================================================
    print(f"\n--- Running Scenario 4: L/2 Gauge String (d = {L // 2}, h = 0.02) ---")
    s4_dir = base_output_dir / "scenario_4_string_h002"
    vis_s4 = FixedPointScreeningVisualizer(cluster, params_h002, s4_dir)

    vis_s4.run_screening_movie(
        initial_state=base_inhom,
        flipped_bonds=flipped_s3,
        max_steps=max_steps,
        gif_filename="screening_string_h002.gif",
        title_prefix=f"Gauge String d={L // 2} (h=0.02)",
    )

    print(f"\n[Movie Generator] All 4 scenarios completed! Artifacts saved to: {base_output_dir}")


if __name__ == "__main__":
    run_all_movie_scenarios(L=8, max_steps=10)
