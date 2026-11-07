"""Service for Task T2.3.2: Critical Exponents & Scaling at the Polarized Breakdown Transition.

Computes forward continuation magnetization curves m_111(h), differential susceptibility
chi_111(h) = dm/dh, and performs Bayesian model comparison for the mean-field spinodal exponent gamma_sp = 0.50.
"""

from dataclasses import dataclass

import numpy as np

from kitaev_mf.domain.lattice import UnitCell, make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.observables import compute_sublattice_magnetizations
from kitaev_mf.services.seeds import make_qsl_seed
from kitaev_mf.services.solver import solve_self_consistent_mean_field


@dataclass(frozen=True)
class SpinodalPoint:
    """Numerical observables at a specific field magnitude along [111]."""

    h_mag: float
    m_111: float
    converged: bool
    iterations: int


@dataclass(frozen=True)
class SusceptibilityPoint:
    """Differential susceptibility evaluated between adjacent field points."""

    h_mid: float
    distance_to_spinodal: float
    chi_111: float


@dataclass(frozen=True)
class SpinodalScalingResult:
    """Bayesian model selection results comparing spinodal divergence against regular crossover."""

    n_samples: int
    fitted_gamma: float
    gamma_relative_error: float
    rss_null: float
    rss_candidate: float
    bic_null: float
    bic_candidate: float
    delta_bic: float
    bayes_factor: float
    favors_spinodal_divergence: bool
    amplitude_a: float


class SpinodalScalingCalculator:
    """Performs high-resolution forward continuation approaching the spinodal breakdown."""

    def __init__(
        self,
        cell: UnitCell | None = None,
        k_grid: tuple[int, int] = (24, 24),
        h_spinodal: float = 0.400,
    ) -> None:
        self.cell = cell or make_honeycomb_2site()
        self.k_grid = k_grid
        self.h_spinodal = h_spinodal
        self.solver_config = SolverConfig(
            k_grid_shape=self.k_grid,
            max_iterations=60,
            tol_correlator=1e-8,
            tol_gauge=1e-8,
            method="diis",
        )

    def run_forward_continuation(
        self,
        h_values: np.ndarray | None = None,
    ) -> list[SpinodalPoint]:
        """Trace forward adiabatic continuation from weak field up to the spinodal limit."""
        if h_values is None:
            h_values = np.array([
                0.350, 0.360, 0.370, 0.375, 0.380,
                0.385, 0.390, 0.392, 0.394, 0.396, 0.398,
            ])

        curr_state: MeanFieldState = make_qsl_seed(self.cell)
        points: list[SpinodalPoint] = []
        h_dir = (1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0), 1.0 / np.sqrt(3.0))

        for h_mag in h_values:
            params = ModelParameters(K=-1.0, h=(h_mag * h_dir[0], h_mag * h_dir[1], h_mag * h_dir[2]))
            state, metrics = solve_self_consistent_mean_field(
                self.cell, params, curr_state, self.solver_config
            )
            curr_state = state

            sigma_mag, _ = compute_sublattice_magnetizations(state, self.cell)
            m_111 = float(np.mean(sigma_mag))

            points.append(
                SpinodalPoint(
                    h_mag=float(h_mag),
                    m_111=m_111,
                    converged=metrics.converged,
                    iterations=metrics.iterations,
                )
            )

        return points

    def compute_susceptibility_series(
        self,
        points: list[SpinodalPoint],
    ) -> list[SusceptibilityPoint]:
        """Compute differential susceptibility chi_111 = dm/dh from continuation curve."""
        h_arr = np.array([pt.h_mag for pt in points])
        m_arr = np.array([pt.m_111 for pt in points])

        dh = np.diff(h_arr)
        dm = np.diff(m_arr)
        chi = dm / dh
        h_mid = 0.5 * (h_arr[:-1] + h_arr[1:])

        series: list[SusceptibilityPoint] = []
        for hm, c in zip(h_mid, chi, strict=True):
            t = float(self.h_spinodal - hm)
            series.append(
                SusceptibilityPoint(
                    h_mid=float(hm),
                    distance_to_spinodal=t,
                    chi_111=float(c),
                )
            )
        return series


class SpinodalBayesianFitter:
    """Performs Bayesian model comparison for the spinodal susceptibility exponent."""

    @staticmethod
    def fit_and_compare(
        susceptibility_series: list[SusceptibilityPoint],
    ) -> SpinodalScalingResult:
        """Compare constant crossover chi=c0 against mean-field spinodal chi=A*(h_sp-h)^(-0.5)."""
        n = len(susceptibility_series)
        assert n > 2, "Need at least 3 points for model comparison"

        t_arr = np.array([p.distance_to_spinodal for p in susceptibility_series])
        chi_arr = np.array([p.chi_111 for p in susceptibility_series])

        # Power-law linear regression in log-log space to extract unconstrained exponent:
        # ln(chi) = -gamma * ln(t) + ln(A)
        log_t = np.log(t_arr)
        log_chi = np.log(chi_arr)
        slope, _ = np.polyfit(log_t, log_chi, 1)
        fitted_gamma = float(-slope)
        rel_error = float(abs(fitted_gamma - 0.50) / 0.50)

        # Null Model M0: Regular constant crossover chi(h) = c0 (k=1)
        c0 = float(np.mean(chi_arr))
        rss_null = float(np.sum((chi_arr - c0)**2))

        # Candidate Model M1: Mean-field spinodal scaling with exact gamma=0.50 (k=1)
        # chi(h) = A * t^(-0.50)
        basis = t_arr**(-0.50)
        norm_sq = float(np.dot(basis, basis))
        amp_a = float(np.dot(chi_arr, basis) / norm_sq) if norm_sq > 1e-12 else 0.0
        fit_1 = amp_a * basis
        rss_cand = float(np.sum((chi_arr - fit_1)**2))

        bic_null = float(n * np.log(max(rss_null / n, 1e-30)) + 1.0 * np.log(n))
        bic_cand = float(n * np.log(max(rss_cand / n, 1e-30)) + 1.0 * np.log(n))

        delta_bic = bic_null - bic_cand
        bf = float(np.exp(min(0.5 * delta_bic, 50.0)))

        return SpinodalScalingResult(
            n_samples=n,
            fitted_gamma=fitted_gamma,
            gamma_relative_error=rel_error,
            rss_null=rss_null,
            rss_candidate=rss_cand,
            bic_null=bic_null,
            bic_candidate=bic_cand,
            delta_bic=delta_bic,
            bayes_factor=bf,
            favors_spinodal_divergence=(delta_bic > 10.0 and bf > 148.0 and rel_error <= 0.15),
            amplitude_a=amp_a,
        )
