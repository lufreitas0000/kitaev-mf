"""Physics service for variational response functions, stability metrics, and dualities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch

from kitaev_mf.domain.lattice import UnitCell
from kitaev_mf.domain.parameters import ModelParameters
from kitaev_mf.domain.ports import LinalgBackendPort
from kitaev_mf.domain.solver import SolverConfig
from kitaev_mf.domain.state import MeanFieldState
from kitaev_mf.services.energetics import compute_gaussian_energy
from kitaev_mf.services.observables import compute_sublattice_magnetizations
from kitaev_mf.services.solver import solve_self_consistent_mean_field


@dataclass(frozen=True)
class HellmannFeynmanResult:
    """Results of numerical Hellmann-Feynman derivative evaluation."""

    de_numerical: float
    neg_sigma_tot: float
    discrepancy: float
    state_base: MeanFieldState


@dataclass(frozen=True)
class KleinDualResult:
    """Results of four-sublattice Klein duality evaluation."""

    energy_site: float
    energy_dual_site: float
    energy_diff: float
    state: MeanFieldState
    state_dual: MeanFieldState


def compute_riemannian_gradient_norm(
    h_k_batch: torch.Tensor,
    p_k_batch: torch.Tensor,
) -> float:
    """Evaluate maximum Frobenius norm of the pointwise commutator [H(k), P(k)]."""
    if h_k_batch.ndim == 2:
        h_k_batch = h_k_batch.unsqueeze(0)
    if p_k_batch.ndim == 2:
        p_k_batch = p_k_batch.unsqueeze(0)

    comm = torch.matmul(h_k_batch, p_k_batch) - torch.matmul(p_k_batch, h_k_batch)
    frob_norms = torch.linalg.matrix_norm(comm, ord="fro")
    return float(torch.max(frob_norms).item())


def evaluate_finite_difference_magnetization(
    cell: UnitCell,
    params: ModelParameters,
    h_dir: int,
    delta: float,
    config: SolverConfig,
    backend: LinalgBackendPort | None = None,
    initial_state: MeanFieldState | None = None,
) -> HellmannFeynmanResult:
    """Evaluate dE/dh^gamma via central difference and compare with -<sigma_tot^gamma>."""
    if h_dir not in (0, 1, 2):
        raise ValueError(f"h_dir must be in {{0, 1, 2}}, got {h_dir}")
    if delta <= 0.0:
        raise ValueError(f"delta must be strictly positive, got {delta}")

    init_s = (
        initial_state
        if initial_state is not None
        else MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.4)
    )

    # 1. Base equilibrium solution
    state_base, _ = solve_self_consistent_mean_field(cell, params, init_s, config, backend)

    # 2. Perturbed parameters
    h_base = list(params.h)
    h_pos = list(h_base)
    h_pos[h_dir] += delta
    h_neg = list(h_base)
    h_neg[h_dir] -= delta

    params_pos = ModelParameters(
        J=params.J,
        K=params.K,
        Gamma=params.Gamma,
        Gamma_prime=params.Gamma_prime,
        D=params.D,
        h=(h_pos[0], h_pos[1], h_pos[2]),
        V0=params.V0,
    )
    params_neg = ModelParameters(
        J=params.J,
        K=params.K,
        Gamma=params.Gamma,
        Gamma_prime=params.Gamma_prime,
        D=params.D,
        h=(h_neg[0], h_neg[1], h_neg[2]),
        V0=params.V0,
    )

    # Warm-start perturbations from the converged base solution
    state_pos, _ = solve_self_consistent_mean_field(cell, params_pos, state_base, config, backend)
    state_neg, _ = solve_self_consistent_mean_field(cell, params_neg, state_base, config, backend)

    e_pos = compute_gaussian_energy(state_pos, cell, params_pos).per_cell
    e_neg = compute_gaussian_energy(state_neg, cell, params_neg).per_cell
    de_num = (e_pos - e_neg) / (2.0 * delta)

    sigma_mag, _ = compute_sublattice_magnetizations(state_base, cell)
    sigma_tot = float(np.sum(sigma_mag[:, h_dir]))
    neg_sigma_tot = -sigma_tot

    discrepancy = abs(de_num - neg_sigma_tot)

    return HellmannFeynmanResult(
        de_numerical=float(de_num),
        neg_sigma_tot=float(neg_sigma_tot),
        discrepancy=float(discrepancy),
        state_base=state_base,
    )


def evaluate_thermodynamic_susceptibility(
    cell: UnitCell,
    params: ModelParameters,
    h_dir: int,
    delta: float,
    config: SolverConfig,
    backend: LinalgBackendPort | None = None,
    initial_state: MeanFieldState | None = None,
) -> float:
    """Evaluate diagonal susceptibility chi^{gamma,gamma} = d<sigma_tot^gamma>/dh^gamma."""
    if h_dir not in (0, 1, 2):
        raise ValueError(f"h_dir must be in {{0, 1, 2}}, got {h_dir}")
    if delta <= 0.0:
        raise ValueError(f"delta must be strictly positive, got {delta}")

    init_s = (
        initial_state
        if initial_state is not None
        else MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.4)
    )

    state_base, _ = solve_self_consistent_mean_field(cell, params, init_s, config, backend)

    h_base = list(params.h)
    h_pos = list(h_base)
    h_pos[h_dir] += delta
    h_neg = list(h_base)
    h_neg[h_dir] -= delta

    params_pos = ModelParameters(
        J=params.J,
        K=params.K,
        Gamma=params.Gamma,
        Gamma_prime=params.Gamma_prime,
        D=params.D,
        h=(h_pos[0], h_pos[1], h_pos[2]),
        V0=params.V0,
    )
    params_neg = ModelParameters(
        J=params.J,
        K=params.K,
        Gamma=params.Gamma,
        Gamma_prime=params.Gamma_prime,
        D=params.D,
        h=(h_neg[0], h_neg[1], h_neg[2]),
        V0=params.V0,
    )

    state_pos, _ = solve_self_consistent_mean_field(cell, params_pos, state_base, config, backend)
    state_neg, _ = solve_self_consistent_mean_field(cell, params_neg, state_base, config, backend)

    sigma_pos, _ = compute_sublattice_magnetizations(state_pos, cell)
    sigma_neg, _ = compute_sublattice_magnetizations(state_neg, cell)

    tot_pos = float(np.sum(sigma_pos[:, h_dir]))
    tot_neg = float(np.sum(sigma_neg[:, h_dir]))

    return (tot_pos - tot_neg) / (2.0 * delta)


def evaluate_klein_dual_pair(
    cell: UnitCell,
    K0: float,
    eps_J: float,
    eps_Gamma: float,
    eps_Gamma_prime: float,
    config: SolverConfig,
    backend: LinalgBackendPort | None = None,
    initial_state: MeanFieldState | None = None,
) -> KleinDualResult:
    """Evaluate self-consistent energy density across a symmetric Klein-dual parameter pair."""
    init_s = (
        initial_state
        if initial_state is not None
        else MeanFieldState.from_pure_kitaev(cell, u_gauge=1.0, u_matter=-0.5)
    )

    params_p = ModelParameters(
        J=eps_J,
        K=K0 - eps_J,
        Gamma=eps_Gamma,
        Gamma_prime=eps_Gamma_prime,
    )
    params_dual = ModelParameters(
        J=-eps_J,
        K=K0 + eps_J,
        Gamma=-eps_Gamma,
        Gamma_prime=eps_Gamma_prime,
    )

    state_p, _ = solve_self_consistent_mean_field(cell, params_p, init_s, config, backend)
    state_dual, _ = solve_self_consistent_mean_field(cell, params_dual, init_s, config, backend)

    energy_p = compute_gaussian_energy(state_p, cell, params_p).per_site
    energy_dual = compute_gaussian_energy(state_dual, cell, params_dual).per_site
    energy_diff = abs(energy_p - energy_dual)

    return KleinDualResult(
        energy_site=energy_p,
        energy_dual_site=energy_dual,
        energy_diff=energy_diff,
        state=state_p,
        state_dual=state_dual,
    )
