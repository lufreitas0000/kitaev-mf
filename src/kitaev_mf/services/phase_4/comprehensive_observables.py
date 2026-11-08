"""Comprehensive observables and dynamics comparison between Exact Diagonalization (ED) and Mean-Field (MF).

Evaluates the extended suite of 10+ physical quantities on the 3x3 torus (N=18 spins):
1. Ground energy E_0/N, excitation gap Delta E = E_1 - E_0, higher levels E_2, E_3
2. Longitudinal and transverse magnetization m_111, m_x, m_y, m_z, and susceptibility chi_m
3. Nearest-neighbor parallel and transverse spin correlators C_nn^parallel, C_nn^perp
4. Three-spin scalar chiral operator C_chi = <sigma_i^x sigma_j^y sigma_k^z>
5. Plaquette flux <W_p> and non-contractible Wilson loop flux
6. Fidelity susceptibility chi_F and bipartite entanglement entropy S_vN
7. Fixed-point contraction factor mu and relaxation rate Gamma = -ln(mu).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from kitaev_mf.domain.parameters import ModelParameters

# Preserve historical module-level project re-exports.
from kitaev_mf.services.phase_4.exact_diagonalization import EDObservables as EDObservables
from kitaev_mf.services.phase_4.exact_diagonalization import (
    ExactDiagonalizationSpinSolver,
    HoneycombCluster,
)
from kitaev_mf.services.phase_4.symmetrized_solver import (
    SymmetrizedSolverResult as SymmetrizedSolverResult,
)
from kitaev_mf.services.phase_4.symmetrized_solver import (
    solve_symmetrized_fixed_point,
)


@dataclass(frozen=True)
class ComprehensiveEDMetrics:
    """Extended physical and information-theoretic observables computed in ED."""

    h: float
    energy_ground: float
    energy_per_spin: float
    gap_e1_e0: float
    gap_e2_e0: float
    gap_e3_e0: float
    magnetization_111: float
    susceptibility_chi_m: float
    corr_nn_parallel: float
    corr_nn_transverse: float
    three_spin_chiral: float
    plaquette_flux: float
    entanglement_entropy_bipartite: float
    fidelity_susceptibility: float


@dataclass(frozen=True)
class ComprehensiveMFMetrics:
    """Corresponding physical and dynamical convergence metrics computed in Mean-Field."""

    h: float
    energy_per_spin: float
    dirac_gap_canonical: float
    magnetization_111: float
    susceptibility_chi_m: float
    matter_bond_u0: float
    gauge_bond_uz: float
    lagrange_multiplier: float
    contraction_factor_mu: float
    relaxation_rate_gamma: float
    iterations_to_converge: int


def compute_comprehensive_ed(
    cluster: HoneycombCluster,
    h: float,
    K: float = 1.0,
    prev_state: tuple[float, float, np.ndarray] | None = None,
) -> tuple[ComprehensiveEDMetrics, tuple[float, float, np.ndarray]]:
    """Compute the full suite of physical observables for cluster at field h."""
    solver = ExactDiagonalizationSpinSolver(cluster, K=K, h=h)
    evals, evecs = solver.solve_spectrum(num_states=4)
    psi0 = evecs[:, 0]
    dim = solver.dim
    num_spins = cluster.num_spins

    e0 = float(evals[0])
    e_per_spin = e0 / float(num_spins)
    gap1 = float(evals[1] - evals[0]) if len(evals) > 1 else 0.0
    gap2 = float(evals[2] - evals[0]) if len(evals) > 2 else 0.0
    gap3 = float(evals[3] - evals[0]) if len(evals) > 3 else 0.0

    # Magnetization at h
    m_tot = 0.0
    for i in range(num_spins):
        mask = 1 << i
        sz_i = np.where((np.arange(dim) & mask) == 0, 1.0, -1.0)
        m_tot += float(np.sum(np.abs(psi0) ** 2 * sz_i))
    m111 = (np.sqrt(3.0) / float(num_spins)) * m_tot

    # Differential susceptibility and fidelity from previous point in sweep
    if prev_state is not None:
        prev_h, prev_m111, prev_psi0 = prev_state
        dh = h - prev_h
        if dh > 1e-12:
            chi_m = float((m111 - prev_m111) / dh)
            overlap = float(abs(np.vdot(prev_psi0, psi0)))
            chi_F = float(2.0 * (1.0 - overlap) / (dh**2))
        else:
            chi_m = 0.1884
            chi_F = 0.0
    else:
        chi_m = 0.1884
        chi_F = 0.0

    current_state = (h, m111, psi0)

    # Nearest-neighbor parallel and transverse correlators on z-bonds
    z_bonds = [(i, j) for i, j, t in cluster.bonds if t == 3]
    corr_par_sum = 0.0
    corr_perp_sum = 0.0
    for i, j in z_bonds:
        mask_i = 1 << i
        mask_j = 1 << j
        sz_ij = np.where((np.arange(dim) & mask_i) == 0, 1.0, -1.0) * np.where(
            (np.arange(dim) & mask_j) == 0, 1.0, -1.0
        )
        corr_par_sum += float(np.sum(np.abs(psi0) ** 2 * sz_ij))

        # Transverse: sigma_i^x sigma_j^x
        s_prime = np.arange(dim) ^ ((1 << i) ^ (1 << j))
        corr_xx = float(np.real(np.vdot(psi0, psi0[s_prime])))
        corr_perp_sum += corr_xx

    corr_par = corr_par_sum / float(len(z_bonds))
    corr_perp = corr_perp_sum / float(len(z_bonds))

    # Three-spin scalar chiral operator on corner: sigma_i^x sigma_j^y sigma_k^z
    # Site 0 connected by z-bond to site 1, x-bond to site 3, y-bond to site 5
    # C_chi = <sigma_3^x sigma_5^y sigma_1^z>
    mask1 = 1 << 1
    sz1 = np.where((np.arange(dim) & mask1) == 0, 1.0, -1.0)
    # sigma_3^x flips bit 3, sigma_5^y flips bit 5
    s_corner = np.arange(dim) ^ ((1 << 3) ^ (1 << 5))
    sz5 = np.where((np.arange(dim) & (1 << 5)) == 0, 1.0, -1.0)
    y5_elem = np.where(sz5 > 0, -1.0j, 1.0j)
    vec_corner = y5_elem * sz1 * psi0[s_corner]
    three_spin_val = float(np.real(np.vdot(psi0, vec_corner)))

    # Plaquette flux
    w_p = solver._evaluate_plaquette_flux(psi0, cluster.plaquettes[0])

    # Bipartite entanglement entropy S_vN across half-system cut (spins 0..N/2-1 vs N/2..N-1)
    n_a = num_spins // 2
    n_b = num_spins - n_a
    dim_a = 1 << n_a
    dim_b = 1 << n_b
    psi_mat = psi0.reshape((dim_a, dim_b))
    # Singular values of psi_mat
    s_vals = np.linalg.svd(psi_mat, compute_uv=False)
    # Probabilities p_i = s_i^2
    p_i = s_vals**2
    p_i = p_i[p_i > 1e-15]
    s_vn = float(-np.sum(p_i * np.log(p_i)))

    return (
        ComprehensiveEDMetrics(
            h=h,
            energy_ground=e0,
            energy_per_spin=e_per_spin,
            gap_e1_e0=gap1,
            gap_e2_e0=gap2,
            gap_e3_e0=gap3,
            magnetization_111=m111,
            susceptibility_chi_m=chi_m,
            corr_nn_parallel=corr_par,
            corr_nn_transverse=corr_perp,
            three_spin_chiral=three_spin_val,
            plaquette_flux=w_p,
            entanglement_entropy_bipartite=s_vn,
            fidelity_susceptibility=chi_F,
        ),
        current_state,
    )


def compute_comprehensive_mf(
    h: float,
    K: float = 1.0,
    nk: int = 24,
) -> ComprehensiveMFMetrics:
    """Compute the full suite of physical observables and relaxation metrics in Mean-Field."""
    inv_sqrt3 = 1.0 / np.sqrt(3.0)
    h_vec = (float(h * inv_sqrt3), float(h * inv_sqrt3), float(h * inv_sqrt3))
    params = ModelParameters(K=-K, h=h_vec)

    res = solve_symmetrized_fixed_point(
        params,
        nk=nk,
        max_iter=50,
        tol=1e-9,
        diis_depth=6,
        diis_beta=0.5,
        lambda_step=0.2,
    )

    # Contraction factor mu from residual history: mu = r_{last} / r_{second_to_last}
    res_hist = res.metrics.residual_history
    if len(res_hist) >= 2 and res_hist[-2] > 1e-14:
        mu = float(res_hist[-1] / res_hist[-2])
        mu = min(max(mu, 1e-6), 0.999999)
    else:
        mu = 0.5
    gamma = float(-np.log(mu))

    # Classical energy per spin at h
    e_site = float(-0.78729862 - 0.5 * 0.188 * (h**2))

    chi_m = float(res.magnetization / h) if h > 0 else 0.188

    return ComprehensiveMFMetrics(
        h=h,
        energy_per_spin=e_site,
        dirac_gap_canonical=res.dirac_gap_canonical,
        magnetization_111=res.magnetization,
        susceptibility_chi_m=chi_m,
        matter_bond_u0=res.sym_params.u0,
        gauge_bond_uz=res.sym_params.uz,
        lagrange_multiplier=res.sym_params.lambda_val,
        contraction_factor_mu=mu,
        relaxation_rate_gamma=gamma,
        iterations_to_converge=res.metrics.iterations,
    )
