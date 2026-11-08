"""Service for Task T2.3.1: Global Field-Directional Magnetization Surface & Torque Discontinuities.

Computes anisotropic magnetization surfaces m(theta, phi, h), magnetic torque
vectors tau = m x h, and performs Bayesian model comparison against isotropic response.
"""

from dataclasses import dataclass

import numpy as np

from kitaev_mf.domain.lattice import UnitCell, make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.observables import compute_sublattice_magnetizations
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field


@dataclass(frozen=True)
class TorquePoint:
    """Numerical observables at a specific field magnitude and angle."""

    h_mag: float
    theta_deg: float
    theta_rad: float
    h_dir: tuple[float, float, float]
    m_vector: tuple[float, float, float]
    m_magnitude: float
    tau_vector: tuple[float, float, float]
    tau_projection: float
    tau_norm: float
    energy_per_site: float


@dataclass(frozen=True)
class TorqueBayesianResult:
    """Bayesian model selection results comparing anisotropic torque vs isotropic null."""

    n_samples: int
    rss_null: float
    rss_candidate: float
    bic_null: float
    bic_candidate: float
    delta_bic: float
    bayes_factor: float
    favors_anisotropic_torque: bool
    amplitude_c: float


class MagnetizationTorqueCalculator:
    """Calculates field-directional magnetization surfaces and magnetic torques."""

    def __init__(self, cell: UnitCell | None = None, k_grid: tuple[int, int] = (16, 16)) -> None:
        self.cell = cell or make_honeycomb_2site()
        self.k_grid = k_grid
        self.solver_config = SolverConfig(
            k_grid_shape=self.k_grid,
            max_iterations=60,
            method="diis",
        )

    @staticmethod
    def field_direction_in_110_plane(theta_rad: float) -> np.ndarray:
        """Field unit vector in the plane connecting [001] (theta=0) and [110] (theta=pi/2).

        At theta = arccos(1/sqrt(3)) (~54.74 deg), h is aligned with [111].
        """
        return np.array(
            [
                np.sin(theta_rad) / np.sqrt(2.0),
                np.sin(theta_rad) / np.sqrt(2.0),
                np.cos(theta_rad),
            ],
            dtype=np.float64,
        )

    def evaluate_point(
        self,
        h_mag: float,
        theta_deg: float,
        warm_seed: MeanFieldState | None = None,
    ) -> tuple[TorquePoint, MeanFieldState]:
        """Compute magnetization and torque at specified field and tilt angle."""
        th_rad = float(np.radians(theta_deg))
        h_dir = self.field_direction_in_110_plane(th_rad)
        h_vec = h_mag * h_dir

        params = ModelParameters(K=-1.0, h=tuple(h_vec))
        seed = warm_seed or make_qsl_seed(self.cell)

        state, _ = solve_self_consistent_mean_field(self.cell, params, seed, self.solver_config)

        sigma_mag, _ = compute_sublattice_magnetizations(state, self.cell)
        m_mean = np.mean(sigma_mag, axis=0)
        m_mag = float(np.linalg.norm(m_mean))

        # Torque tau = m x h
        tau_vec = np.cross(m_mean, h_vec)
        tau_norm = float(np.linalg.norm(tau_vec))
        # Projection along (1, -1, 0)/sqrt(2)
        tau_proj = float((tau_vec[0] - tau_vec[1]) / np.sqrt(2.0))

        eb = compute_gaussian_energy(state, self.cell, params)
        e_site = float(eb.E_total / self.cell.num_sublattices)

        pt = TorquePoint(
            h_mag=h_mag,
            theta_deg=theta_deg,
            theta_rad=th_rad,
            h_dir=(float(h_dir[0]), float(h_dir[1]), float(h_dir[2])),
            m_vector=(float(m_mean[0]), float(m_mean[1]), float(m_mean[2])),
            m_magnitude=m_mag,
            tau_vector=(float(tau_vec[0]), float(tau_vec[1]), float(tau_vec[2])),
            tau_projection=tau_proj,
            tau_norm=tau_norm,
            energy_per_site=e_site,
        )
        return pt, state

    def sweep_angles(
        self,
        h_mag: float,
        thetas_deg: np.ndarray | None = None,
    ) -> list[TorquePoint]:
        """Sweep field orientation in the (1, 1, 0) plane at fixed field magnitude."""
        if thetas_deg is None:
            thetas_deg = np.linspace(0.0, 90.0, 13)

        points: list[TorquePoint] = []
        curr_state = None
        for th in thetas_deg:
            pt, curr_state = self.evaluate_point(h_mag, float(th), warm_seed=curr_state)
            points.append(pt)
        return points


class TorqueBayesianFitter:
    """Performs Bayesian model comparison between isotropic response and anisotropic Kitaev torque."""

    @staticmethod
    def fit_and_compare(
        thetas_rad: np.ndarray,
        tau_projections: np.ndarray,
        h_magnitudes: np.ndarray | None = None,
    ) -> TorqueBayesianResult:
        """Compare isotropic null model tau=0 against anisotropic torque candidate.

        When h_magnitudes are provided, fits field-dependent amplitude C(h) across field slices
        accounting for the non-linear growth of torque near the spinodal breakdown.
        """
        n = len(tau_projections)
        assert n > 2, "Need at least 3 points for model comparison"

        # Null model M0: tau = 0
        rss_null = float(np.sum(tau_projections**2))

        # Candidate model M1: tau = C(h) * basis(theta)
        # Basis vanishes at theta=0 ([001]), theta=arccos(1/sqrt(3)) ([111]), theta=pi/2 ([110])
        basis = np.sin(2.0 * thetas_rad) * (3.0 * np.cos(thetas_rad) ** 2 - 1.0)

        if h_magnitudes is not None and len(np.unique(h_magnitudes)) > 1:
            unique_h = np.unique(h_magnitudes)
            k_params = len(unique_h)
            fit = np.zeros_like(tau_projections)
            amps = []
            for h_val in unique_h:
                mask = h_magnitudes == h_val
                b_slice = basis[mask]
                t_slice = tau_projections[mask]
                norm_sq = float(np.dot(b_slice, b_slice))
                c_val = float(np.dot(t_slice, b_slice) / norm_sq) if norm_sq > 1e-12 else 0.0
                fit[mask] = c_val * b_slice
                amps.append(c_val)
            c_rep = float(np.mean(amps))
        else:
            k_params = 1
            norm_sq = float(np.dot(basis, basis))
            c_rep = float(np.dot(tau_projections, basis) / norm_sq) if norm_sq > 1e-12 else 0.0
            fit = c_rep * basis

        rss_cand = float(np.sum((tau_projections - fit) ** 2))

        # BIC calculation
        bic_0 = float(n * np.log(max(rss_null / n, 1e-30)))
        bic_1 = float(n * np.log(max(rss_cand / n, 1e-30)) + float(k_params) * np.log(n))

        delta_bic = bic_0 - bic_1
        bf = float(np.exp(min(0.5 * delta_bic, 50.0)))

        return TorqueBayesianResult(
            n_samples=n,
            rss_null=rss_null,
            rss_candidate=rss_cand,
            bic_null=bic_0,
            bic_candidate=bic_1,
            delta_bic=delta_bic,
            bayes_factor=bf,
            favors_anisotropic_torque=(delta_bic > 10.0 and bf > 148.0),
            amplitude_c=c_rep,
        )
