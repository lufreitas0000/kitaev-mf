"""Pure physics service for Majorana Zero Modes (MZMs) bound to isolated vortices in the Kitaev QSL."""

from collections import deque
from dataclasses import dataclass

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


@dataclass(frozen=True)
class MZMModeMetrics:
    """Spectral and spatial metrics of a vortex-bound Majorana zero mode."""

    L: int
    h: float
    p1: tuple[int, int]
    p2: tuple[int, int]
    R: float
    eps_mzm: float
    eps_bulk: float
    in_gap_ratio: float
    w_core: float
    w_uniform: float
    enrichment_ratio: float
    hop_profile: dict[int, float]


@dataclass(frozen=True)
class MZMBayesianSplittingResult:
    """Bayesian model selection between constant splitting and topological exponential decoupling."""

    n_samples: int
    rss_null_constant: float
    rss_candidate_exponential: float
    bic_null_constant: float
    bic_candidate_exponential: float
    delta_bic: float
    bayes_factor: float
    decay_length_xi: float
    favors_exponential_decay: bool


class MZMLocalizationCalculator:
    """Analyzes the spectral isolation and wavefunction spatial localization of vortex-bound MZMs."""

    def __init__(self, L: int = 8, d_nn: float = 1.0) -> None:
        """Initialize honeycomb torus cluster."""
        self.L = L
        self.d_nn = d_nn
        self.geom = RealSpaceTorusCluster(L=L, d_nn=d_nn)
        self.cell = make_honeycomb_2site(d_nn)
        self.qsl_seed = make_qsl_seed(self.cell)

    def compute_bulk_gap(self, params: ModelParameters) -> float:
        """Compute lowest positive excitation energy in the flux-free sector."""
        A0 = self.geom.assemble_skew_hamiltonian(params, self.qsl_seed)
        S0 = -A0 @ A0
        w0 = np.sqrt(np.maximum(np.linalg.eigvalsh(0.5 * (S0 + S0.T)), 0.0))
        eps_sorted = np.sort(w0)[::2]
        return float(eps_sorted[0])

    def build_lattice_adjacency(self) -> dict[int, list[int]]:
        """Construct the site-to-site adjacency graph for graph hop distances."""
        adj: dict[int, list[int]] = {s: [] for s in range(self.geom.num_sites)}
        disp = [(0, 0, 3), (0, -1, 1), (-1, 0, 2)]
        L = self.L
        for m2 in range(L):
            for m1 in range(L):
                sA = self.geom.site_index(m1, m2, 0)
                for d1, d2, _ in disp:
                    sB = self.geom.site_index(m1 + d1, m2 + d2, 1)
                    adj[sA].append(sB)
                    adj[sB].append(sA)
        return adj

    def compute_hop_distance_profile(
        self,
        psi_site: np.ndarray,
        p1: tuple[int, int],
        p2: tuple[int, int],
    ) -> dict[int, float]:
        """Compute the average wavefunction density as a function of graph hop distance from vortex cores."""
        adj = self.build_lattice_adjacency()
        vortex_boundary = self.geom.compute_plaquette_boundary_sites(
            p1
        ) | self.geom.compute_plaquette_boundary_sites(p2)

        hop_dist: dict[int, int] = dict.fromkeys(range(self.geom.num_sites), -1)
        q: deque[int] = deque()
        for s in vortex_boundary:
            hop_dist[s] = 0
            q.append(s)

        while q:
            curr = q.popleft()
            curr_d = hop_dist[curr]
            for nbr in adj[curr]:
                if hop_dist[nbr] == -1:
                    hop_dist[nbr] = curr_d + 1
                    q.append(nbr)

        max_d = max(hop_dist.values())
        profile: dict[int, float] = {
            d: float(np.mean([psi_site[s] for s, dist in hop_dist.items() if dist == d]))
            for d in range(max_d + 1)
        }
        return profile

    def evaluate_mzm_pair(
        self,
        params: ModelParameters,
        p1: tuple[int, int] | None = None,
        p2: tuple[int, int] | None = None,
    ) -> MZMModeMetrics:
        """Compute MZM excitation energy, in-gap ratio, core weight, and spatial profile."""
        L = self.L
        pos1 = p1 if p1 is not None else (0, 0)
        pos2 = p2 if p2 is not None else (L // 2, L // 2)

        # 1. Bulk gap
        eps_bulk = self.compute_bulk_gap(params)

        # 2. Assemble vortex pair Hamiltonian
        cut = self.geom.generate_vortex_pair_cut(pos1, pos2)
        Av = self.geom.assemble_skew_hamiltonian(params, self.qsl_seed, flipped_bonds=cut)
        Sv = -Av @ Av

        # 3. Diagonalize Sv = Av Av^T
        w2, v = np.linalg.eigh(0.5 * (Sv + Sv.T))
        w2 = np.maximum(w2, 0.0)
        eps_mzm = float(np.sqrt(w2[0]))

        # Lowest mode eigenvector
        psi0 = v[:, 0]
        total_sites = self.geom.num_sites
        psi_site = np.sum(np.abs(psi0.reshape((total_sites, 4))) ** 2, axis=1)
        total_norm = float(np.sum(psi_site))
        if total_norm > 1e-14:
            psi_site /= total_norm

        # 4. Core localization weight
        w_core, w_unif, ratio = self.geom.compute_vortex_mode_core_weight(psi0, pos1, pos2)

        # 5. Hop distance decay profile
        hop_profile = self.compute_hop_distance_profile(psi_site, pos1, pos2)

        # Euclidean separation on dual triangular lattice with periodic image minimization
        dm1 = (pos2[0] - pos1[0]) % L
        dm2 = (pos2[1] - pos1[1]) % L
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
        R = float(np.sqrt(min_r2))

        h_mag = float(np.linalg.norm(params.h))
        in_gap_ratio = float(eps_mzm / (eps_bulk + 1e-12))

        return MZMModeMetrics(
            L=L,
            h=h_mag,
            p1=pos1,
            p2=pos2,
            R=R,
            eps_mzm=eps_mzm,
            eps_bulk=eps_bulk,
            in_gap_ratio=in_gap_ratio,
            w_core=w_core,
            w_uniform=w_unif,
            enrichment_ratio=ratio,
            hop_profile=hop_profile,
        )


class MZMSplittingBayesianFitter:
    """Performs Bayesian model selection between constant splitting and topological exponential decoupling."""

    def fit_and_compare(
        self,
        l_vals: np.ndarray,
        eps_vals: np.ndarray,
    ) -> MZMBayesianSplittingResult:
        """Compare Null Model M0 (constant splitting) vs Candidate Model M1 (exponential decay).

        Args:
            l_vals: Array of linear cluster dimensions L (e.g. [4, 5, 7, 8]).
            eps_vals: Array of lowest MZM hybridization energies eps_0(L).

        Returns:
            MZMBayesianSplittingResult containing RSS, BIC, Bayes factor, and decay length xi.
        """
        n = len(eps_vals)
        if n < 3:
            raise ValueError(f"Need at least 3 points for model comparison, got {n}")

        # 1. Null Model M0 (Constant splitting): eps(L) = eps_bar (k = 1 free parameter)
        eps_mean = float(np.mean(eps_vals))
        rss_null = float(np.sum((eps_vals - eps_mean) ** 2))
        k_null = 1
        bic_null = n * np.log(max(rss_null / n, 1e-12)) + k_null * np.log(n)

        # 2. Candidate Model M1 (Exponential decoupling): ln(eps) = -L / xi + ln A (k = 2 free parameters)
        poly = np.polyfit(l_vals, np.log(np.maximum(eps_vals, 1e-12)), 1)
        slope, intercept = poly
        pred_cand = np.exp(slope * l_vals + intercept)
        rss_cand = float(np.sum((eps_vals - pred_cand) ** 2))
        k_cand = 2
        bic_cand = n * np.log(max(rss_cand / n, 1e-12)) + k_cand * np.log(n)

        delta_bic = float(bic_null - bic_cand)
        bayes_factor = float(np.exp(np.clip(delta_bic / 2.0, -50.0, 50.0)))
        decay_length = float(-1.0 / slope) if slope < 0.0 else float("inf")

        favors_exponential = bool(delta_bic > 10.0 and bayes_factor > 148.0 and slope < 0.0)

        return MZMBayesianSplittingResult(
            n_samples=n,
            rss_null_constant=rss_null,
            rss_candidate_exponential=rss_cand,
            bic_null_constant=float(bic_null),
            bic_candidate_exponential=float(bic_cand),
            delta_bic=delta_bic,
            bayes_factor=bayes_factor,
            decay_length_xi=decay_length,
            favors_exponential_decay=favors_exponential,
        )
