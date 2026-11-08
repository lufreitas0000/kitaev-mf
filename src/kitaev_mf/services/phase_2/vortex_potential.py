"""Pure physics service for vortex separation potential V(R), vison deconfinement and screening."""

from dataclasses import dataclass

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


@dataclass(frozen=True)
class VortexPairPoint:
    """Numerical result for a single vortex-antivortex separation configuration."""

    p1: tuple[int, int]
    p2: tuple[int, int]
    R: float
    cut_length: int
    E_flux_free: float
    E_2vortex: float
    V_R: float
    w_p1: float
    w_p2: float


@dataclass(frozen=True)
class DeconfinementBayesianResult:
    """Statistical outcome of Bayesian model selection between confinement and deconfinement."""

    n_samples: int
    rss_null_confinement: float
    rss_candidate_deconfinement: float
    bic_null_confinement: float
    bic_candidate_deconfinement: float
    delta_bic: float
    bayes_factor: float
    unconstrained_sigma: float
    v_inf: float
    vison_gap: float
    asymptotic_flatness: float
    vison_gap_relative_error: float
    favors_deconfinement: bool


class VortexSeparationPotentialCalculator:
    """Calculates the inter-vortex interaction potential V(R) on a honeycomb torus cluster."""

    def __init__(self, L: int = 7, d_nn: float = 1.0) -> None:
        """Initialize cluster geometry and precompute flux-free reference state.

        Args:
            L: Linear size of the torus (L x L unit cells).
            d_nn: Nearest-neighbor bond distance.
        """
        self.L = L
        self.d_nn = d_nn
        self.geom = RealSpaceTorusCluster(L=L, d_nn=d_nn)
        self.cell = make_honeycomb_2site(d_nn)
        self.qsl_seed = make_qsl_seed(self.cell)

    def compute_plaquette_distance(
        self,
        p1: tuple[int, int],
        p2: tuple[int, int],
    ) -> float:
        """Compute the gauge-invariant Euclidean separation on the dual triangular lattice.

        Args:
            p1: Coordinate (m1, m2) of first plaquette.
            p2: Coordinate (m1, m2) of second plaquette.

        Returns:
            Euclidean distance R taking periodic minimum image into account.
        """
        L = self.L
        dm1 = (p2[0] - p1[0]) % L
        dm2 = (p2[1] - p1[1]) % L

        min_r2 = float("inf")
        for s1 in (-L, 0, L):
            for s2 in (-L, 0, L):
                d1 = dm1 + s1
                d2 = dm2 + s2
                dx = d1 * np.sqrt(3.0) + d2 * (np.sqrt(3.0) / 2.0)
                dy = d2 * 1.5
                r2 = dx**2 + dy**2
                if r2 < min_r2:
                    min_r2 = r2

        return float(np.sqrt(min_r2))

    def evaluate_flux_free_energy(self, params: ModelParameters) -> float:
        """Evaluate reference ground state energy in the flux-free sector W_p = +1."""
        A0 = self.geom.assemble_skew_hamiltonian(params, self.qsl_seed)
        _, E0 = self.geom.solve_covariance_real_eigh(A0)
        return float(E0)

    def evaluate_vortex_pair(
        self,
        p1: tuple[int, int],
        p2: tuple[int, int],
        params: ModelParameters,
        E_flux_free: float | None = None,
    ) -> VortexPairPoint:
        """Compute ground state energy and interaction potential for a vortex pair at p1 and p2.

        Args:
            p1: Coordinate of first vortex.
            p2: Coordinate of second vortex.
            params: Exchange couplings and external field.
            E_flux_free: Optional precomputed reference energy.

        Returns:
            VortexPairPoint containing energy, potential, and flux verification.
        """
        if E_flux_free is None:
            E_flux_free = self.evaluate_flux_free_energy(params)

        cut = self.geom.generate_vortex_pair_cut(p1, p2)
        Av = self.geom.assemble_skew_hamiltonian(params, self.qsl_seed, flipped_bonds=cut)
        omegav, Ev = self.geom.solve_covariance_real_eigh(Av)

        fluxes = self.geom.compute_plaquette_fluxes(omegav)
        w_p1 = fluxes.get(p1, 1.0)
        w_p2 = fluxes.get(p2, 1.0)

        R = self.compute_plaquette_distance(p1, p2)
        V_R = float(Ev - E_flux_free)

        return VortexPairPoint(
            p1=p1,
            p2=p2,
            R=R,
            cut_length=len(cut),
            E_flux_free=float(E_flux_free),
            E_2vortex=float(Ev),
            V_R=V_R,
            w_p1=float(w_p1),
            w_p2=float(w_p2),
        )

    def scan_separation_profile(
        self,
        params: ModelParameters,
        max_separation: int | None = None,
    ) -> list[VortexPairPoint]:
        """Scan vortex separation profile V(R) across distinct non-equivalent lattice separations.

        Args:
            params: Exchange couplings and field.
            max_separation: Optional maximum coordinate displacement along axes.

        Returns:
            Sorted list of VortexPairPoint ordered by increasing Euclidean separation R.
        """
        p1 = (0, 0)
        E0 = self.evaluate_flux_free_energy(params)
        max_m = max_separation if max_separation is not None else self.L // 2

        seen_distances: dict[float, VortexPairPoint] = {}

        for m2 in range(max_m + 1):
            for m1 in range(max_m + 1):
                if m1 == 0 and m2 == 0:
                    continue
                p2 = (m1, m2)
                point = self.evaluate_vortex_pair(p1, p2, params, E_flux_free=E0)
                key = round(point.R, 3)
                if key not in seen_distances or point.cut_length < seen_distances[key].cut_length:
                    seen_distances[key] = point

        results = list(seen_distances.values())
        results.sort(key=lambda pt: pt.R)
        return results


class DeconfinementBayesianFitter:
    """Performs rigorous Bayesian model comparison between confining string tension and deconfined visons."""

    KITAEV_EXACT_VISON_GAP = 0.1536
    CONFINEMENT_MIN_STRING_TENSION = 0.20

    def fit_and_compare(
        self,
        r_vals: np.ndarray,
        v_vals: np.ndarray,
        cut_lengths: np.ndarray | None = None,
        sigma_conf: float = CONFINEMENT_MIN_STRING_TENSION,
    ) -> DeconfinementBayesianResult:
        """Compare Null Model M0 (Confining String) vs Candidate Model M1 (Topological Deconfinement).

        Args:
            r_vals: Array of Euclidean vortex separations.
            v_vals: Array of measured interaction potentials V(R).
            cut_lengths: Optional branch cut lengths N_cut. If None, uses r_vals / sqrt(3).
            sigma_conf: Minimal string tension under confinement hypothesis M0.

        Returns:
            DeconfinementBayesianResult containing RSS, BIC, Bayes factor, and vison gap metrics.
        """
        n = len(v_vals)
        if n < 3:
            raise ValueError(f"Need at least 3 data points for model comparison, got {n}")

        cuts = cut_lengths if cut_lengths is not None else np.asarray(r_vals / np.sqrt(3.0))

        # 1. Null Model M0 (Confining String): V(N_cut) = sigma_conf * N_cut + V0 (k = 1 free parameter V0)
        v0_conf = float(np.mean(v_vals - sigma_conf * cuts))
        pred_null = sigma_conf * cuts + v0_conf
        rss_null = float(np.sum((v_vals - pred_null) ** 2))

        # 2. Candidate Model M1 (Topological Deconfinement): V(R) = V_inf (k = 1 free parameter V_inf)
        v_inf = float(np.mean(v_vals))
        rss_cand = float(np.sum((v_vals - v_inf) ** 2))

        # 3. Unconstrained linear slope for verification
        poly = np.polyfit(cuts, v_vals, 1)
        unconstrained_sigma = float(poly[0])

        # 4. Bayesian Information Criterion: BIC = n * ln(RSS / n) + k * ln(n)
        # Both models have k = 1 free parameter (V0 and V_inf), enabling direct likelihood ratio comparison
        k = 1
        bic_null = n * np.log(max(rss_null / n, 1e-12)) + k * np.log(n)
        bic_cand = n * np.log(max(rss_cand / n, 1e-12)) + k * np.log(n)
        delta_bic = float(bic_null - bic_cand)

        # Bayes Factor: B_10 = exp(Delta_BIC / 2)
        bayes_factor = float(np.exp(np.clip(delta_bic / 2.0, -50.0, 50.0)))

        # 5. Asymptotic plateau flatness: |V(R_max) - V(R_max-1)| / V(R_max)
        v_max = float(v_vals[-1])
        v_prev = float(v_vals[-2]) if n >= 2 else v_max
        asymptotic_flatness = float(abs(v_max - v_prev) / (abs(v_max) + 1e-12))

        # 6. Vison gap evaluation
        vison_gap = float(v_inf / 2.0)
        rel_error = float(
            abs(vison_gap - self.KITAEV_EXACT_VISON_GAP) / self.KITAEV_EXACT_VISON_GAP
        )

        favors_deconfinement = bool(delta_bic > 10.0 and bayes_factor > 148.0 and rel_error < 0.15)

        return DeconfinementBayesianResult(
            n_samples=n,
            rss_null_confinement=rss_null,
            rss_candidate_deconfinement=rss_cand,
            bic_null_confinement=float(bic_null),
            bic_candidate_deconfinement=float(bic_cand),
            delta_bic=delta_bic,
            bayes_factor=bayes_factor,
            unconstrained_sigma=unconstrained_sigma,
            v_inf=v_inf,
            vison_gap=vison_gap,
            asymptotic_flatness=asymptotic_flatness,
            vison_gap_relative_error=rel_error,
            favors_deconfinement=favors_deconfinement,
        )
