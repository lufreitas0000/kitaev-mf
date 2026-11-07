"""Service for Task T2.3.3: Global Phase Diagram Mapping & First-Order Boundary Anisotropy.

Maps the thermodynamic first-order phase boundary h_c*(theta) between the non-Abelian QSL
and the field-polarized phase across field angles in the (1, 1, 0) plane, and evaluates
Bayesian model comparison against isotropic phase boundary hypotheses.
"""

from dataclasses import dataclass

import numpy as np

from kitaev_mf.domain.lattice import UnitCell, make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.seeds import make_polarized_seed, make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field


@dataclass(frozen=True)
class PhaseBoundaryPoint:
    """Critical thermodynamic transition field at a specific field tilt angle."""

    theta_deg: float
    theta_rad: float
    h_c_star: float
    energy_crossing: float
    qsl_energy_330: float
    pol_energy_330: float
    qsl_energy_340: float
    pol_energy_340: float


@dataclass(frozen=True)
class PhaseBoundaryBayesianResult:
    """Bayesian model selection results comparing anisotropic phase boundary against isotropic null."""

    n_samples: int
    h_base: float
    delta_h: float
    theta_max_deg: float
    rss_null: float
    rss_candidate: float
    bic_null: float
    bic_candidate: float
    delta_bic: float
    bayes_factor: float
    favors_anisotropic_boundary: bool


class PhaseBoundaryCalculator:
    """Calculates thermodynamic first-order critical fields h_c*(theta)."""

    def __init__(
        self,
        cell: UnitCell | None = None,
        k_grid: tuple[int, int] = (16, 16),
    ) -> None:
        self.cell = cell or make_honeycomb_2site()
        self.k_grid = k_grid
        self.solver_config = SolverConfig(
            k_grid_shape=self.k_grid,
            max_iterations=100,
            method="diis",
        )

    @staticmethod
    def field_direction_in_110_plane(theta_rad: float) -> np.ndarray:
        """Unit vector in the (1, 1, 0) plane connecting [001] (theta=0) and [110] (theta=pi/2)."""
        return np.array([
            np.sin(theta_rad) / np.sqrt(2.0),
            np.sin(theta_rad) / np.sqrt(2.0),
            np.cos(theta_rad),
        ], dtype=np.float64)

    def find_critical_field(self, theta_deg: float) -> PhaseBoundaryPoint:
        """Find thermodynamic crossing h_c^*(theta) where E_qsl(h) = E_pol(h)."""
        th_rad = float(np.radians(theta_deg))
        h_dir = self.field_direction_in_110_plane(th_rad)

        # Evaluate at brackets h=0.330 and h=0.340
        p1 = ModelParameters(K=-1.0, h=(0.330 * h_dir[0], 0.330 * h_dir[1], 0.330 * h_dir[2]))
        s_q1, m_q1 = solve_self_consistent_mean_field(self.cell, p1, make_qsl_seed(self.cell), self.solver_config)
        assert m_q1.converged, f"QSL solver failed to converge at h=0.330, theta={theta_deg}"
        e_q1 = float(compute_gaussian_energy(s_q1, self.cell, p1).E_total / self.cell.num_sublattices)
        s_p1, m_p1 = solve_self_consistent_mean_field(self.cell, p1, make_polarized_seed(self.cell), self.solver_config)
        assert m_p1.converged, f"Polarized solver failed to converge at h=0.330, theta={theta_deg}"
        e_p1 = float(compute_gaussian_energy(s_p1, self.cell, p1).E_total / self.cell.num_sublattices)
        d1 = e_q1 - e_p1

        p2 = ModelParameters(K=-1.0, h=(0.340 * h_dir[0], 0.340 * h_dir[1], 0.340 * h_dir[2]))
        s_q2, m_q2 = solve_self_consistent_mean_field(self.cell, p2, make_qsl_seed(self.cell), self.solver_config)
        assert m_q2.converged, f"QSL solver failed to converge at h=0.340, theta={theta_deg}"
        e_q2 = float(compute_gaussian_energy(s_q2, self.cell, p2).E_total / self.cell.num_sublattices)
        s_p2, m_p2 = solve_self_consistent_mean_field(self.cell, p2, make_polarized_seed(self.cell), self.solver_config)
        assert m_p2.converged, f"Polarized solver failed to converge at h=0.340, theta={theta_deg}"
        e_p2 = float(compute_gaussian_energy(s_p2, self.cell, p2).E_total / self.cell.num_sublattices)
        d2 = e_q2 - e_p2

        # Verify bracketing
        assert d1 < 0.0 and d2 > 0.0, f"Root not bracketed at theta={theta_deg}: d1={d1}, d2={d2}"

        # Linear root interpolation: d(h) = d1 + (h - 0.330)*(d2 - d1)/0.010 = 0
        h_c = float(0.330 + 0.010 * (-d1) / (d2 - d1))
        # Energy at crossing:
        e_cross = float(e_p1 + (h_c - 0.330) * (e_p2 - e_p1) / 0.010)

        return PhaseBoundaryPoint(
            theta_deg=theta_deg,
            theta_rad=th_rad,
            h_c_star=h_c,
            energy_crossing=e_cross,
            qsl_energy_330=e_q1,
            pol_energy_330=e_p1,
            qsl_energy_340=e_q2,
            pol_energy_340=e_p2,
        )

    def sweep_boundary(self, thetas_deg: np.ndarray | None = None) -> list[PhaseBoundaryPoint]:
        """Sweep the thermodynamic phase boundary across angles."""
        if thetas_deg is None:
            thetas_deg = np.linspace(0.0, 90.0, 13)

        points: list[PhaseBoundaryPoint] = []
        for th in thetas_deg:
            pt = self.find_critical_field(float(th))
            points.append(pt)
        return points


class PhaseBoundaryBayesianFitter:
    """Performs Bayesian model comparison between isotropic and anisotropic phase boundaries."""

    @staticmethod
    def fit_and_compare(
        thetas_deg: np.ndarray,
        h_c_values: np.ndarray,
    ) -> PhaseBoundaryBayesianResult:
        """Compare isotropic null model h_c = h_0 against candidate h_c = h_base + Delta_h * basis(theta)."""
        n = len(h_c_values)
        assert n > 2, "Need at least 3 points for model comparison"

        # Model 0: Isotropic boundary (k=1)
        h_0 = float(np.mean(h_c_values))
        rss_null = float(np.sum((h_c_values - h_0)**2))

        # Model 1: Anisotropic boundary (k=2)
        # Basis peaks at theta = arccos(1/sqrt(3)) ~ 54.74 deg ([111] direction):
        # sin^4(th)*cos^2(th) = (2/3)^2 * (1/3) = 4/27.
        th_rad = np.radians(thetas_deg)
        basis = (np.sin(th_rad)**4) * (np.cos(th_rad)**2)
        max_basis = 4.0 / 27.0
        basis_norm = basis / max_basis

        x_mat = np.column_stack([np.ones(n, dtype=np.float64), basis_norm])
        coeffs, _, _, _ = np.linalg.lstsq(x_mat, h_c_values, rcond=None)
        h_base = float(coeffs[0])
        delta_h = float(coeffs[1])

        fit_cand = x_mat @ coeffs
        rss_cand = float(np.sum((h_c_values - fit_cand)**2))

        bic_null = float(n * np.log(max(rss_null / n, 1e-30)) + 1.0 * np.log(n))
        bic_cand = float(n * np.log(max(rss_cand / n, 1e-30)) + 2.0 * np.log(n))

        delta_bic = bic_null - bic_cand
        bf = float(np.exp(min(0.5 * delta_bic, 50.0)))

        # Find angle of maximum critical field
        max_idx = int(np.argmax(h_c_values))
        th_max = float(thetas_deg[max_idx])

        return PhaseBoundaryBayesianResult(
            n_samples=n,
            h_base=h_base,
            delta_h=delta_h,
            theta_max_deg=th_max,
            rss_null=rss_null,
            rss_candidate=rss_cand,
            bic_null=bic_null,
            bic_candidate=bic_cand,
            delta_bic=delta_bic,
            bayes_factor=bf,
            favors_anisotropic_boundary=(delta_bic > 10.0 and bf > 148.0),
        )
