"""Pure physics service for vortex braiding, Berry phases, and non-Abelian anyonic holonomy."""

from dataclasses import dataclass

import numpy as np

from kitaev_mf.domain.lattice import make_honeycomb_2site
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.services.real_space import RealSpaceTorusCluster
from kitaev_mf.services.seeds import make_qsl_seed


@dataclass(frozen=True)
class BraidingHolonomyResult:
    """Berry phases and statistical holonomy accumulated along closed dual vortex loops."""

    L: int
    h: float
    gamma_enclosing: float
    gamma_non_enclosing: float
    delta_gamma: float
    delta_gamma_over_pi: float
    relative_error_from_pi: float


@dataclass(frozen=True)
class BraidingBayesianResult:
    """Bayesian model selection between trivial holonomy (0) and topological anyonic holonomy (pi)."""

    n_samples: int
    rss_null: float
    rss_candidate: float
    bic_null: float
    bic_candidate: float
    delta_bic: float
    bayes_factor: float
    mean_delta_gamma: float
    favors_pi_holonomy: bool


class VortexBraidingCalculator:
    """Evaluates the adiabatic Berry phase and statistical holonomy of vortex transport."""

    def __init__(self, L: int = 6, d_nn: float = 1.0) -> None:
        """Initialize honeycomb cluster and dual graph."""
        self.L = L
        self.d_nn = d_nn
        self.geom = RealSpaceTorusCluster(L=L, d_nn=d_nn)
        self.cell = make_honeycomb_2site(d_nn)
        self.qsl_seed = make_qsl_seed(self.cell)
        self.dual = self.geom.build_dual_graph()

    def find_crossing_bond(
        self,
        p_a: tuple[int, int],
        p_b: tuple[int, int],
    ) -> tuple[int, int, int]:
        """Find the physical bond crossing the dual edge between adjacent plaquettes p_a and p_b."""
        for nbr, b in self.dual[p_a]:
            if nbr == p_b:
                return b
        raise ValueError(f"No dual edge between {p_a} and {p_b}")

    def get_occupied_frame(self, A: np.ndarray) -> np.ndarray:
        """Extract occupied single-particle Bogoliubov frame U from Hermitian matrix H_c = i * A."""
        H_c = 1j * A
        vals, vecs = np.linalg.eigh(H_c)
        M = A.shape[0]
        # Choose positive-energy single-particle modes
        pos_indices = np.argsort(vals)[M // 2 :]
        return vecs[:, pos_indices]

    def compute_gauge_aligned_berry_phase(
        self,
        loop_plaquettes: list[tuple[int, int]],
        p_stat: tuple[int, int],
        enclosed_plaquettes: list[tuple[int, int]],
        params: ModelParameters,
        N_sub: int = 4,
    ) -> float:
        """Compute the gauge-aligned many-body Berry phase accumulated along a closed dual loop.

        Args:
            loop_plaquettes: Ordered list of plaquettes forming a closed path on the dual lattice.
            p_stat: Coordinate of the stationary reference vortex.
            enclosed_plaquettes: Plaquettes enclosed by the loop requiring gauge alignment.
            params: Model exchange parameters and magnetic field.
            N_sub: Number of adiabatic interpolation steps per bond flip.

        Returns:
            Geometric Berry phase in radians modulo 2*pi.
        """
        curr_cut = set(self.geom.generate_vortex_pair_cut(p_stat, loop_plaquettes[0]))
        frames: list[np.ndarray] = []

        for i in range(len(loop_plaquettes) - 1):
            p_curr = loop_plaquettes[i]
            p_next = loop_plaquettes[i + 1]
            b_flip = self.find_crossing_bond(p_curr, p_next)
            was_in_cut = b_flip in curr_cut

            for s in range(N_sub):
                lam = float(s) / float(N_sub)
                # Continuous interpolation between current sign and flipped sign
                disorder_val = np.cos(np.pi * lam) - 1.0

                A_step = self.geom.assemble_skew_hamiltonian(
                    params,
                    self.qsl_seed,
                    flipped_bonds=curr_cut,
                    bond_disorder={b_flip: disorder_val},
                )
                frames.append(self.get_occupied_frame(A_step))

            if was_in_cut:
                curr_cut.remove(b_flip)
            else:
                curr_cut.add(b_flip)

        # Final state closing the path
        A_final = self.geom.assemble_skew_hamiltonian(params, self.qsl_seed, flipped_bonds=curr_cut)
        u_final = self.get_occupied_frame(A_final)

        # Local gauge alignment operator G for enclosed plaquettes
        G_diag = np.ones(self.geom.num_majoranas, dtype=np.float64)
        for p_enc in enclosed_plaquettes:
            b_sites = self.geom.compute_plaquette_boundary_sites(p_enc)
            for s in b_sites:
                for mu in range(4):
                    G_diag[4 * s + mu] *= -1.0
        G = np.diag(G_diag)
        u_final_aligned = G @ u_final

        # Compute many-body Wilson loop: prod_k det(U_k^\dagger U_{k+1})
        total_log_det = 0.0 + 0.0j
        for k in range(len(frames) - 1):
            overlap = frames[k].conj().T @ frames[k + 1]
            s_val, log_det = np.linalg.slogdet(overlap)
            total_log_det += np.log(s_val) + log_det

        # Overlap closing step to final frame
        overlap_close = frames[-1].conj().T @ u_final
        s_c, log_c = np.linalg.slogdet(overlap_close)
        total_log_det += np.log(s_c) + log_c

        # Overlap closing to initial gauge-aligned frame
        overlap_gauge = u_final_aligned.conj().T @ frames[0]
        s_g, log_g = np.linalg.slogdet(overlap_gauge)
        total_log_det += np.log(s_g) + log_g

        gamma = -float(np.imag(total_log_det)) % (2.0 * np.pi)
        return gamma

    def evaluate_braiding_holonomy(
        self,
        params: ModelParameters,
        N_sub: int = 4,
    ) -> BraidingHolonomyResult:
        """Evaluate Berry phases for both enclosing and non-enclosing loops to compute holonomy."""
        L = self.L
        p_stat = (0, 0)

        # Loop enclosing (0, 0)
        loop_enc = [
            (1, 0),
            (0, 1),
            ((L - 1) % L, 1),
            ((L - 1) % L, 0),
            (0, (L - 1) % L),
            (1, (L - 1) % L),
            (1, 0),
        ]

        # Congruent 6-step loop enclosing empty plaquette (m, m)
        m = L // 2
        p_empty = (m, m)
        loop_non = [
            ((m + 1) % L, m),
            (m, (m + 1) % L),
            ((m - 1) % L, (m + 1) % L),
            ((m - 1) % L, m),
            (m, (m - 1) % L),
            ((m + 1) % L, (m - 1) % L),
            ((m + 1) % L, m),
        ]

        gamma_enc = self.compute_gauge_aligned_berry_phase(
            loop_enc, p_stat, [p_stat], params, N_sub=N_sub
        )
        gamma_non = self.compute_gauge_aligned_berry_phase(
            loop_non, p_stat, [p_empty], params, N_sub=N_sub
        )

        # Phase difference modulo 2*pi, centered on [0, 2*pi)
        diff = (gamma_enc - gamma_non) % (2.0 * np.pi)
        if diff > np.pi * 1.5:
            diff = abs(diff - 2.0 * np.pi)

        # Metric of interest is distance to pi
        rel_error = float(abs(diff - np.pi) / np.pi)
        h_mag = float(np.linalg.norm(params.h))

        return BraidingHolonomyResult(
            L=L,
            h=h_mag,
            gamma_enclosing=gamma_enc,
            gamma_non_enclosing=gamma_non,
            delta_gamma=diff,
            delta_gamma_over_pi=float(diff / np.pi),
            relative_error_from_pi=rel_error,
        )


class BraidingBayesianFitter:
    """Performs Bayesian model selection between trivial holonomy (Delta_gamma = 0) and anyonic holonomy (pi)."""

    def fit_and_compare(
        self,
        delta_gamma_vals: np.ndarray,
    ) -> BraidingBayesianResult:
        """Compare Null Model M0 (Delta_gamma = 0) vs Candidate Model M1 (Delta_gamma = pi).

        Args:
            delta_gamma_vals: Array of measured statistical holonomy differences Delta_gamma.

        Returns:
            BraidingBayesianResult containing RSS, BIC, and Bayes factor metrics.
        """
        n = len(delta_gamma_vals)
        if n < 1:
            raise ValueError("Need at least 1 measurement for model comparison")

        # Null Model M0: Delta_gamma = 0 (k = 0 free parameters)
        rss_null = float(np.sum((delta_gamma_vals - 0.0) ** 2))

        # Candidate Model M1: Delta_gamma = pi (k = 0 free parameters)
        rss_cand = float(np.sum((delta_gamma_vals - np.pi) ** 2))

        # BIC for zero-parameter models: BIC = n * ln(RSS / n)
        bic_null = float(n * np.log(max(rss_null / n, 1e-12)))
        bic_cand = float(n * np.log(max(rss_cand / n, 1e-12)))

        delta_bic = float(bic_null - bic_cand)
        bayes_factor = float(np.exp(np.clip(delta_bic / 2.0, -50.0, 50.0)))
        mean_diff = float(np.mean(delta_gamma_vals))

        favors_pi = bool(delta_bic > 10.0 and bayes_factor > 148.0)

        return BraidingBayesianResult(
            n_samples=n,
            rss_null=rss_null,
            rss_candidate=rss_cand,
            bic_null=bic_null,
            bic_candidate=bic_cand,
            delta_bic=delta_bic,
            bayes_factor=bayes_factor,
            mean_delta_gamma=mean_diff,
            favors_pi_holonomy=favors_pi,
        )
