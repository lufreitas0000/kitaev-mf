"""Thermodynamics service for Phase 2: Specific heat, fractionalization entropy stages, and Bayesian model comparison."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


@dataclass(frozen=True)
class ThermodynamicPoint:
    """Thermodynamic observables at a given temperature."""

    temperature: float
    c_v_vison: float
    c_v_maj: float
    c_v_total: float
    s_vison: float
    s_maj: float
    s_total: float


@dataclass(frozen=True)
class SpecificHeatResult:
    """Summary of thermodynamic curve across temperature."""

    points: list[ThermodynamicPoint]
    t_low: float
    c_v_low: float
    t_high: float
    c_v_high: float
    t_mid: float
    c_v_mid: float
    s_mid: float
    t_max: float
    s_max: float
    peak_ratio: float


@dataclass(frozen=True)
class ThermodynamicBayesianResult:
    """Bayesian model comparison between two-peak fractionalized and single-peak models."""

    n_samples: int
    rss_null: float
    rss_cand: float
    bic_null: float
    bic_cand: float
    delta_bic: float
    bayes_factor: float
    favors_two_peak: bool


class SpecificHeatCalculator:
    """Calculates finite-temperature specific heat and fractionalization entropy."""

    def __init__(
        self,
        delta_v: float = 0.137,
        k_grid: tuple[int, int] = (32, 32),
    ) -> None:
        self.delta_v = float(delta_v)
        self.k_grid = k_grid
        self._eps_k = self._compute_majorana_dispersion()

    def _compute_majorana_dispersion(self) -> np.ndarray:
        """Compute the itinerant matter Majorana dispersion across the 2D Brillouin zone.

        For pure Kitaev exchange K=-1.0, the single-particle dispersion is
        epsilon(k) = |1 + exp(i 2pi k1) + exp(i 2pi k2)| / 2.
        """
        kx = np.linspace(0.0, 1.0, self.k_grid[0], endpoint=False)
        ky = np.linspace(0.0, 1.0, self.k_grid[1], endpoint=False)
        kx_grid, ky_grid = np.meshgrid(kx, ky)
        k_pts = np.column_stack([kx_grid.ravel(), ky_grid.ravel()])

        phase1 = 2.0 * np.pi * k_pts[:, 0]
        phase2 = 2.0 * np.pi * k_pts[:, 1]
        f_k = 1.0 + np.exp(1j * phase1) + np.exp(1j * phase2)
        # Scaled single-particle matter dispersion (bandwidth 3.0 |K|):
        return np.abs(f_k)

    @staticmethod
    def _schottky_cv(x: np.ndarray | float) -> np.ndarray | float:
        """Compute x^2 * exp(-x) / (1 + exp(-x))^2 with numerical stability."""
        safe_x = np.asarray(x)
        # For large x, exp(-x) vanishes
        mask = safe_x < 50.0
        res = np.zeros_like(safe_x, dtype=np.float64)
        x_valid = safe_x[mask]
        exp_neg = np.exp(-x_valid)
        res[mask] = (x_valid**2) * exp_neg / ((1.0 + exp_neg) ** 2)
        return float(res) if np.isscalar(x) else res

    @staticmethod
    def _fermionic_entropy(x: np.ndarray | float) -> np.ndarray | float:
        """Compute ln(1 + exp(-x)) + x * exp(-x) / (1 + exp(-x)) with numerical stability."""
        safe_x = np.asarray(x)
        mask = safe_x < 50.0
        res = np.zeros_like(safe_x, dtype=np.float64)
        x_valid = safe_x[mask]
        exp_neg = np.exp(-x_valid)
        res[mask] = np.log(1.0 + exp_neg) + (x_valid * exp_neg) / (1.0 + exp_neg)
        return float(res) if np.isscalar(x) else res

    def compute_thermodynamic_point(self, temp: float) -> ThermodynamicPoint:
        """Evaluate vison and Majorana thermodynamic observables at temperature T."""
        t_val = float(temp)
        assert t_val > 0.0, "Temperature must be positive"

        # 1. Vison contribution: 1/2 of a 2-level system per spin (N_p = N_spins / 2)
        xv = self.delta_v / t_val
        cv_v = float(0.5 * self._schottky_cv(xv))
        s_v = float(0.5 * self._fermionic_entropy(xv))

        # 2. Itinerant Majorana contribution: 1 matter band per 2 spins
        xm = self._eps_k / t_val
        cv_m_k = np.asarray(self._schottky_cv(xm))
        cv_m = float(0.5 * np.mean(cv_m_k))
        s_m_k = np.asarray(self._fermionic_entropy(xm))
        s_m = float(0.5 * np.mean(s_m_k))

        cv_tot = cv_v + cv_m
        s_tot = s_v + s_m

        return ThermodynamicPoint(
            temperature=t_val,
            c_v_vison=cv_v,
            c_v_maj=cv_m,
            c_v_total=cv_tot,
            s_vison=s_v,
            s_maj=s_m,
            s_total=s_tot,
        )

    def compute_thermodynamic_curve(self, temps: np.ndarray | None = None) -> SpecificHeatResult:
        """Compute thermodynamic observables across temperatures and extract peak structure."""
        if temps is None:
            temps = np.logspace(-2.5, 0.8, 40)

        points = [self.compute_thermodynamic_point(float(t)) for t in temps]
        cv_total = np.array([p.c_v_total for p in points])

        # Identify peaks
        peaks, _ = find_peaks(cv_total)
        assert len(peaks) >= 2, f"Expected at least two peaks in C_v(T), found {len(peaks)}"

        idx_low = peaks[0]
        idx_high = peaks[-1]

        t_low = float(points[idx_low].temperature)
        cv_low = float(points[idx_low].c_v_total)
        t_high = float(points[idx_high].temperature)
        cv_high = float(points[idx_high].c_v_total)
        peak_ratio = float(t_high / t_low)

        # Identify intermediate dip between peaks
        mid_slice = cv_total[idx_low : idx_high + 1]
        idx_mid = idx_low + int(np.argmin(mid_slice))

        t_mid = float(points[idx_mid].temperature)
        cv_mid = float(points[idx_mid].c_v_total)
        s_mid = float(points[idx_mid].s_total)

        t_max = float(points[-1].temperature)
        s_max = float(points[-1].s_total)

        return SpecificHeatResult(
            points=points,
            t_low=t_low,
            c_v_low=cv_low,
            t_high=t_high,
            c_v_high=cv_high,
            t_mid=t_mid,
            c_v_mid=cv_mid,
            s_mid=s_mid,
            t_max=t_max,
            s_max=s_max,
            peak_ratio=peak_ratio,
        )


class ThermodynamicBayesianFitter:
    """Performs Bayesian model selection between single-peak and two-peak models."""

    @staticmethod
    def _single_schottky(temp: np.ndarray, a: float, delta: float) -> np.ndarray:
        x = delta / temp
        safe_x = np.asarray(x)
        mask = safe_x < 50.0
        res = np.zeros_like(safe_x, dtype=np.float64)
        x_valid = safe_x[mask]
        exp_neg = np.exp(-x_valid)
        res[mask] = a * (x_valid**2) * exp_neg / ((1.0 + exp_neg) ** 2)
        return res

    @classmethod
    def _double_schottky(
        cls,
        temp: np.ndarray,
        a1: float,
        d1: float,
        a2: float,
        d2: float,
    ) -> np.ndarray:
        return np.asarray(
            cls._single_schottky(temp, a1, d1) + cls._single_schottky(temp, a2, d2),
            dtype=np.float64,
        )

    @classmethod
    def fit_and_compare(
        cls,
        temps: np.ndarray,
        cv_values: np.ndarray,
    ) -> ThermodynamicBayesianResult:
        """Compare single-peak null model (k=2) against two-peak candidate (k=4)."""
        n = len(cv_values)
        assert n > 4, "Need at least 5 points for model comparison"

        # Model 0: Single Schottky (k=2)
        p0_null = [0.2, 0.5]
        popt_null, _ = curve_fit(
            cls._single_schottky,
            temps,
            cv_values,
            p0=p0_null,
            bounds=([0.0, 0.001], [1.0, 5.0]),
        )
        fit_null = cls._single_schottky(temps, *popt_null)
        rss_null = float(np.sum((cv_values - fit_null) ** 2))
        bic_null = float(n * np.log(rss_null / n) + 2 * np.log(n))

        # Model 1: Double Schottky (k=4)
        p0_cand = [0.2, 0.1, 0.2, 1.0]
        popt_cand, _ = curve_fit(
            cls._double_schottky,
            temps,
            cv_values,
            p0=p0_cand,
            bounds=([0.0, 0.01, 0.0, 0.2], [1.0, 0.5, 1.0, 5.0]),
        )
        fit_cand = cls._double_schottky(temps, *popt_cand)
        rss_cand = float(np.sum((cv_values - fit_cand) ** 2))
        bic_cand = float(n * np.log(rss_cand / n) + 4 * np.log(n))

        delta_bic = float(bic_null - bic_cand)
        # Avoid overflow in exp(0.5 * delta_bic)
        bayes_factor = float(np.exp(min(0.5 * delta_bic, 700.0)))
        favors_two_peak = delta_bic > 10.0 and bayes_factor > 148.0

        return ThermodynamicBayesianResult(
            n_samples=n,
            rss_null=rss_null,
            rss_cand=rss_cand,
            bic_null=bic_null,
            bic_cand=bic_cand,
            delta_bic=delta_bic,
            bayes_factor=bayes_factor,
            favors_two_peak=favors_two_peak,
        )
