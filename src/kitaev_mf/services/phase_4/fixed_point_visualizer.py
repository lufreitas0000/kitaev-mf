"""Visualization engine for self-consistent fixed-point screening dynamics on real-space honeycomb tori.

Renders publication-quality frames and animated GIFs illustrating the relaxation of
the state under the self-consistent map Phi after:
1. Flipping a single gauge bond (creating a pair of nearest-neighbor vortices).
2. Flipping an L/2 gauge string (creating a pair of maximally separated vortices).

Features:
- Hexagonal plaquette polygons shaded by Pfaffian flux W_p in [-1, +1] (+1: white, -1: dark gray).
- Bond edges with thickness and color mapped to U_b bond correlators.
- Highlighted defect / branch-cut edges.
- Information header with Frame index, Energy E(k), and max gauge residual.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Sequence
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from PIL import Image

from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.domain.state import InhomogeneousMeanFieldState


class FixedPointScreeningVisualizer:
    """Renders step-by-step lattice frames of the self-consistent fixed-point relaxation."""

    def __init__(
        self,
        cluster: RealSpaceTorusCluster,
        params: ModelParameters,
        output_dir: Path | str,
    ) -> None:
        """Initialize visualizer.

        Args:
            cluster: Real-space torus cluster of linear size L.
            params: Model parameters (exchange K, magnetic field h).
            output_dir: Directory where PNG frames and GIFs will be saved.
        """
        self.cluster = cluster
        self.L = cluster.L
        self.params = params
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Lattice geometry constants
        self.d_nn = 1.0 / np.sqrt(3.0)
        self.a1 = np.array([0.5, np.sqrt(3.0) / 2.0], dtype=np.float64)
        self.a2 = np.array([-0.5, np.sqrt(3.0) / 2.0], dtype=np.float64)

        # Precompute 2D Cartesian coordinates of all sites
        self.site_coords = np.zeros((cluster.num_sites, 2), dtype=np.float64)
        for m2 in range(self.L):
            for m1 in range(self.L):
                r_cell = m1 * self.a1 + m2 * self.a2
                s_A = cluster.site_index(m1, m2, 0)
                s_B = cluster.site_index(m1, m2, 1)
                self.site_coords[s_A] = r_cell
                self.site_coords[s_B] = r_cell + np.array([0.0, self.d_nn], dtype=np.float64)

    def render_frame(
        self,
        frame_idx: int,
        inhom_state: InhomogeneousMeanFieldState,
        plaquette_fluxes: dict[tuple[int, int], float],
        flipped_bonds: set[tuple[int, int, int]],
        energy: float,
        residual: float,
        title_prefix: str = "Fixed-Point Screening",
    ) -> Path:
        """Render a single high-resolution PNG frame of the lattice state."""
        fig, ax = plt.subplots(figsize=(10, 9), dpi=150)
        ax.set_aspect("equal")
        ax.set_facecolor("#1A1A1A")  # Dark background to highlight contrast

        L = self.L
        patches = []
        flux_colors = []

        # 1. Draw hexagonal plaquettes shaded by W_p
        for m2 in range(L):
            for m1 in range(L):
                p = (m1, m2)
                w_val = float(plaquette_fluxes.get(p, 1.0))
                # Normalized color: w_val = +1 -> 1.0 (white), w_val = -1 -> 0.2 (dark gray)
                # Linear interpolation in [0.2, 1.0]
                gray_val = float(0.6 + 0.4 * w_val)
                gray_val = min(max(gray_val, 0.15), 1.0)

                # Plaquette vertices
                # Surrounding hexagon around cell (m1, m2):
                # 6 vertices in clockwise/counter-clockwise order
                v_sites = [
                    self.cluster.site_index(m1, m2, 0),
                    self.cluster.site_index(m1, m2, 1),
                    self.cluster.site_index(m1, (m2 + 1) % L, 0),
                    self.cluster.site_index((m1 - 1) % L, (m2 + 1) % L, 1),
                    self.cluster.site_index((m1 - 1) % L, (m2 + 1) % L, 0),
                    self.cluster.site_index((m1 - 1) % L, m2, 1),
                ]
                coords = self.site_coords[v_sites]

                # Check for PBC wrapping in polygon; if wrapped, use unwrap around centroid
                center = coords[0]
                unwrapped_coords = []
                for pt in coords:
                    diff = pt - center
                    # If wrapped across torus, shift
                    # Simple centroid polygon
                    unwrapped_coords.append(pt)

                poly = Polygon(coords, closed=True)
                patches.append(poly)
                flux_colors.append((gray_val, gray_val, gray_val))

        # Add plaquette collection
        p_coll = PatchCollection(patches, facecolors=flux_colors, edgecolors="#444444", linewidths=0.5)
        ax.add_collection(p_coll)

        # 2. Draw bonds with color and thickness denoting U_b
        disp = [(0, 0, 3), (0, -1, 1), (-1, 0, 2)]
        for m2 in range(L):
            for m1 in range(L):
                c = m1 + L * m2
                s_A = self.cluster.site_index(m1, m2, 0)
                pt_A = self.site_coords[s_A]

                for b_idx, (d1, d2, gamma) in enumerate(disp):
                    # Check if bond wraps boundary
                    is_boundary = (m1 + d1 < 0) or (m1 + d1 >= L) or (m2 + d2 < 0) or (m2 + d2 >= L)
                    s_B = self.cluster.site_index(m1 + d1, m2 + d2, 1)
                    pt_B = self.site_coords[s_B]

                    if is_boundary:
                        continue  # Skip boundary wraps for clean visualization

                    # Diagonal matter hopping and gauge link
                    u_matter = abs(float(inhom_state.U[c, b_idx, 0, 0]))
                    u_gauge = float(inhom_state.U[c, b_idx, gamma, gamma])

                    # Check if bond is in flipped gauge cut
                    bond_key = (m1, m2, b_idx)
                    is_flipped = bond_key in flipped_bonds

                    if is_flipped:
                        # Highlighted defect bond in bright cyan / orange
                        ax.plot(
                            [pt_A[0], pt_B[0]],
                            [pt_A[1], pt_B[1]],
                            color="#00FFFF",  # Cyan highlight for flipped link
                            linewidth=3.0,
                            zorder=4,
                        )
                    else:
                        # Color mapped to matter hopping strength
                        lw = max(0.8, 4.0 * u_matter)
                        # Regular bond in cool blue / white
                        ax.plot(
                            [pt_A[0], pt_B[0]],
                            [pt_A[1], pt_B[1]],
                            color="#4A90E2",
                            linewidth=lw,
                            alpha=0.85,
                            zorder=2,
                        )

        # 3. Information and metrics overlay
        h_norm = float(np.linalg.norm(self.params.h_vec))
        min_w = float(min(plaquette_fluxes.values()))
        info_text = (
            f"{title_prefix}\n"
            f"Frame: {frame_idx:02d} | Field h = {h_norm:.4f} | Lattice: {L}x{L} (N = {2*L*L})\n"
            f"Energy E = {energy:.6f} | Step Residual ||Δ|| = {residual:.2e} | Min(W_p) = {min_w:.4f}"
        )
        ax.text(
            0.02,
            0.98,
            info_text,
            transform=ax.transAxes,
            fontsize=11,
            verticalalignment="top",
            fontfamily="monospace",
            color="#FFFFFF",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#000000", alpha=0.8, edgecolor="#666666"),
        )

        # Fit display boundaries
        all_x = self.site_coords[:, 0]
        all_y = self.site_coords[:, 1]
        pad = 0.8
        ax.set_xlim(all_x.min() - pad, all_x.max() + pad)
        ax.set_ylim(all_y.min() - pad, all_y.max() + pad)
        ax.axis("off")

        frame_file = self.output_dir / f"frame_{frame_idx:03d}.png"
        fig.savefig(frame_file, bbox_inches="tight", facecolor=fig.get_facecolor(), pad_inches=0.1)
        plt.close(fig)
        return frame_file

    def run_screening_movie(
        self,
        initial_state: InhomogeneousMeanFieldState,
        flipped_bonds: set[tuple[int, int, int]],
        max_steps: int = 15,
        tol: float = 1e-6,
        gif_filename: str = "screening_evolution.gif",
        title_prefix: str = "Screening Relaxation",
    ) -> Path:
        """Run step-by-step fixed-point iteration, rendering frames and compiling an animated GIF."""
        curr_state = initial_state
        frame_files: list[Path] = []

        print(f"[Visualizer] Starting movie capture in {self.output_dir} (max {max_steps} frames)...")

        # Frame 0: Initial state right after defect insertion
        a0 = self.cluster.assemble_skew_hamiltonian(self.params, curr_state, flipped_bonds=flipped_bonds)
        omega0, e0 = self.cluster.solve_covariance_real_eigh(a0)
        fluxes0 = self.cluster.compute_all_plaquette_pfaffian_fluxes(omega0)

        f0_path = self.render_frame(
            frame_idx=0,
            inhom_state=curr_state,
            plaquette_fluxes=fluxes0,
            flipped_bonds=flipped_bonds,
            energy=e0,
            residual=1.0,
            title_prefix=f"{title_prefix} (Defect Injected)",
        )
        frame_files.append(f0_path)

        prev_omega = omega0.copy()

        for step in range(1, max_steps + 1):
            # 1. Update inhomogeneous state from covariance
            curr_state = self.cluster.extract_inhomogeneous_state(prev_omega, curr_state.lambdas)

            # 2. Re-assemble Hamiltonian
            a_k = self.cluster.assemble_skew_hamiltonian(self.params, curr_state, flipped_bonds=flipped_bonds)

            # 3. Diagonalize covariance
            omega_k, e_k = self.cluster.solve_covariance_real_eigh(a_k)

            # 4. Compute residual
            step_res = float(np.max(np.abs(omega_k - prev_omega)))
            prev_omega = omega_k.copy()

            # 5. Compute plaquette fluxes
            fluxes_k = self.cluster.compute_all_plaquette_pfaffian_fluxes(omega_k)

            f_path = self.render_frame(
                frame_idx=step,
                inhom_state=curr_state,
                plaquette_fluxes=fluxes_k,
                flipped_bonds=flipped_bonds,
                energy=e_k,
                residual=step_res,
                title_prefix=f"{title_prefix} (Relaxation)",
            )
            frame_files.append(f_path)
            print(f"  Frame {step:02d} rendered: E = {e_k:.6f}, ||Δ|| = {step_res:.2e}")

            if step_res < tol:
                print(f"[Visualizer] Converged to stationary fixed point in {step} steps!")
                break

        # Compile animated GIF using PIL
        gif_path = self.output_dir / gif_filename
        pil_images = [Image.open(f) for f in frame_files]
        if pil_images:
            pil_images[0].save(
                gif_path,
                save_all=True,
                append_images=pil_images[1:],
                duration=400,  # 400 ms per frame
                loop=0,
            )
            print(f"[Visualizer] Animated GIF compiled successfully: {gif_path}")

        return gif_path
