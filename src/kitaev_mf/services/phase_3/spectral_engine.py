"""Dynamical Spin Correlation Lehmann Spectral Summation and Local Spectral Functions.

Implements the spectral summation engine for arXiv:2209.12930:
- Discrete excitation spectrum Delta E_lambda = E_lambda(G') - E_0(G) across local spin quenches.
- Many-body transition weights w_lambda = |M_{j, lambda}|^2.
- Continuous spectral broadening (Lorentzian and Gaussian).
- Numerical quadrature and unitary sum rule verification.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from kitaev_mf.services.phase_3.bogoliubov import (
    compute_bogoliubov_transformation,
    compute_majorana_matrix_elements,
    extract_mode_wavefunctions,
)
from kitaev_mf.services.phase_3.gauge_theory import (
    DiagonalizedMatterSector,
    DirectedBond,
    GaugeConfiguration,
    assemble_matter_skew_matrix,
    build_honeycomb_torus_cluster,
    diagonalize_matter_hamiltonian,
)
from kitaev_mf.services.phase_3.vortex_engine import (
    HoneycombGeometry,
    apply_local_spin_quench,
)


@dataclass(frozen=True)
class SpectralPole:
    """A discrete excitation pole in the Lehmann representation of the spin correlator."""

    energy: float  # Excitation energy omega_lambda = E_lambda(G') - E_0(G)
    weight: float  # Transition probability w_lambda = |M_{j, lambda}|^2
    amplitude: complex  # Many-body transition amplitude M_{j, lambda}
    quasiparticle_state: tuple[int, ...]  # () for vacuum, (m,) for 1-qp, (m1, m2) for 2-qp


@dataclass(frozen=True)
class DiscreteSpectralData:
    """Discrete spectral representation of the local dynamical spin correlation function."""

    site_j: int
    bond_type: int  # 1 (x), 2 (y), 3 (z)
    quenched_bond: tuple[int, int]  # (sa, sb) link flipped by spin operator
    flux_gap: float  # E_vac(G') - E_vac(G)
    initial_ground_energy: float  # E_0(G)
    quenched_ground_energy: float  # E_0(G')
    poles: list[SpectralPole]  # Sorted list of discrete poles
    total_spectral_weight: float  # sum_lambda w_lambda
    is_singular: bool  # Whether Bogoliubov transformation was singular


@dataclass(frozen=True)
class ContinuousSpectrum:
    """Continuous broadened local dynamical spectral function S_{jj}^{alpha alpha}(omega)."""

    omega: np.ndarray  # Frequency grid points
    spectral_density: np.ndarray  # S_{jj}^{alpha alpha}(omega) values
    broadening_type: str  # "gaussian" or "lorentzian"
    broadening_param: float  # eta or sigma
    integrated_weight: float  # Numerical integral int (d omega / 2 pi) S(omega)


def compute_local_spectral_poles(
    gauge: GaugeConfiguration,
    site_j: int,
    bond_type: int,  # 1 (x), 2 (y), 3 (z)
    geom: HoneycombGeometry,
    K: float = 1.0,
    kappa: float = 0.15,
    boundary: str = "obc",
    bonds: list[DirectedBond] | None = None,
    triplets: list[tuple[int, int, int, int]] | None = None,
    diag_initial: DiagonalizedMatterSector | None = None,
) -> DiscreteSpectralData:
    """Compute the discrete Lehmann spectral poles for the local dynamical spin correlator.

    Simulates the action of sigma_j^alpha:
    1. Flips the bond link connected to site_j of type bond_type: u_{jk}^alpha -> -u_{jk}^alpha.
    2. Diagonalizes the quenched sector G' = G_j^alpha.
    3. Evaluates Bogoliubov transformation blocks X, Y and overlap <0_b | 0_a>.
    4. Computes many-body matrix elements <lambda(G') | c_j | phi_0(G)>.
    5. Evaluates excitation energies Delta E_lambda = E_lambda(G') - E_0(G).

    Parameters:
        boundary: Parity projection policy. "obc" (default) selects the thermodynamic vacuum
            parity policy (pi_G = +1), representing the physical bulk response on finite clusters.
            "torus" enforces Pedrocchi-Zschocke topological parity from non-contractible Wilson loops.

    Returns:
        DiscreteSpectralData containing all excitation poles and weights.
    """
    L = geom.L
    if bonds is None or triplets is None:
        _, bonds_built, triplets_built = build_honeycomb_torus_cluster(L=L)
        bonds = bonds or bonds_built
        triplets = triplets or triplets_built

    # 1. Diagonalize initial sector G if not provided
    if diag_initial is None:
        A_init = assemble_matter_skew_matrix(gauge, K=K, kappa=kappa, triplets=triplets)
        diag_initial = diagonalize_matter_hamiltonian(A_init, gauge=gauge, boundary=boundary, L=L)

    # 2. Perform sudden spin quench on target bond
    quenched_gauge, flipped_bond = apply_local_spin_quench(
        gauge=gauge, site_j=site_j, bond_type=bond_type, geom=geom
    )

    # 3. Diagonalize quenched sector G'
    A_quenched = assemble_matter_skew_matrix(quenched_gauge, K=K, kappa=kappa, triplets=triplets)
    diag_quenched = diagonalize_matter_hamiltonian(
        A_quenched, gauge=quenched_gauge, boundary=boundary, L=L
    )

    # 4. Extract mode wavefunctions and Bogoliubov transformation
    phi_a = extract_mode_wavefunctions(diag_initial.Q)
    phi_b = extract_mode_wavefunctions(diag_quenched.Q)
    bog = compute_bogoliubov_transformation(phi_a, phi_b)

    # 5. Compute transition amplitudes
    res = compute_majorana_matrix_elements(
        bog=bog,
        phi_a=phi_a,
        phi_b=phi_b,
        site_j=site_j,
        initial_parity=diag_initial.physical_matter_parity,
        occupied_mode=0 if diag_initial.physical_matter_parity == -1 else None,
        overlap_phase=1.0 + 0.0j,
    )

    # 6. Assemble excitation energies
    E0_G = diag_initial.physical_ground_energy
    E_vac_q = diag_quenched.unprojected_vacuum_energy
    eps_b = diag_quenched.energies
    mu = bog.singular_mode if bog.is_singular else None

    poles: list[SpectralPole] = []
    for state, amp in res["amplitudes"].items():
        weight = float(np.abs(amp) ** 2)

        if not bog.is_singular:
            # Regular transformation
            if len(state) == 0:
                E_state = E_vac_q
            elif len(state) == 1:
                E_state = E_vac_q + float(eps_b[state[0]])
            elif len(state) == 2:
                E_state = E_vac_q + float(eps_b[state[0]] + eps_b[state[1]])
            else:
                raise ValueError(f"Unsupported quasiparticle state {state}")
        else:
            # Singular transformation with particle-hole interchange on mode mu
            assert mu is not None
            if len(state) == 0:
                # |0_btilde> = b_mu^dagger |0_b>
                E_state = E_vac_q + float(eps_b[mu])
            elif len(state) == 1:
                m = state[0]
                if m == mu:
                    # b_mu |0_btilde> = |0_b>
                    E_state = E_vac_q
                else:
                    # b_m^dagger |0_btilde> = b_m^dagger b_mu^dagger |0_b>
                    E_state = E_vac_q + float(eps_b[mu] + eps_b[m])
            elif len(state) == 2:
                m1, m2 = state
                if m1 == mu:
                    E_state = E_vac_q + float(eps_b[m2])
                elif m2 == mu:
                    E_state = E_vac_q + float(eps_b[m1])
                else:
                    E_state = E_vac_q + float(eps_b[mu] + eps_b[m1] + eps_b[m2])
            else:
                raise ValueError(f"Unsupported quasiparticle state {state}")

        dE = float(E_state - E0_G)
        poles.append(
            SpectralPole(
                energy=dE,
                weight=weight,
                amplitude=amp,
                quasiparticle_state=state,
            )
        )

    # Sort poles in ascending energy order
    poles.sort(key=lambda p: p.energy)
    total_w = float(sum(p.weight for p in poles))
    flux_gap = float(
        diag_quenched.unprojected_vacuum_energy - diag_initial.unprojected_vacuum_energy
    )

    return DiscreteSpectralData(
        site_j=site_j,
        bond_type=bond_type,
        quenched_bond=flipped_bond,
        flux_gap=flux_gap,
        initial_ground_energy=E0_G,
        quenched_ground_energy=diag_quenched.physical_ground_energy,
        poles=poles,
        total_spectral_weight=total_w,
        is_singular=bog.is_singular,
    )


def evaluate_spectral_function(
    spectral_data: DiscreteSpectralData,
    omega_grid: np.ndarray,
    broadening_type: str = "gaussian",
    broadening_param: float = 0.05,
) -> ContinuousSpectrum:
    """Convolve discrete Lehmann poles with a broadening kernel onto a continuous frequency grid.

    Conventions:
    - S_{jj}^{alpha alpha}(omega) = 2 pi sum_lambda w_lambda delta(omega - omega_lambda)
    - Normalization: int_0^infty (d omega / 2 pi) S(omega) = sum_lambda w_lambda = 1.0

    Broadening kernels:
    - Gaussian: delta_sigma(omega) = 1 / (sqrt(2 pi) sigma) exp(-omega^2 / (2 sigma^2))
      => S(omega) = sqrt(2 pi) sum_lambda w_lambda / sigma exp(-(omega - omega_lambda)^2 / (2 sigma^2))
    - Lorentzian: delta_eta(omega) = 1 / pi * eta / (omega^2 + eta^2)
      => S(omega) = 2 sum_lambda w_lambda * eta / ((omega - omega_lambda)^2 + eta^2)
    """
    if broadening_type not in ("gaussian", "lorentzian"):
        raise ValueError(
            f"Invalid broadening_type {broadening_type}, must be 'gaussian' or 'lorentzian'"
        )
    if broadening_param <= 0.0:
        raise ValueError(f"Broadening parameter must be positive, got {broadening_param}")

    spectral_density = np.zeros_like(omega_grid, dtype=np.float64)

    if broadening_type == "gaussian":
        sigma = broadening_param
        prefactor = float(np.sqrt(2.0 * np.pi) / sigma)
        for pole in spectral_data.poles:
            arg = (omega_grid - pole.energy) / sigma
            spectral_density += pole.weight * prefactor * np.exp(-0.5 * arg**2)
    else:  # lorentzian
        eta = broadening_param
        for pole in spectral_data.poles:
            denom = (omega_grid - pole.energy) ** 2 + eta**2
            spectral_density += 2.0 * pole.weight * (eta / denom)

    # Numerical quadrature integral int (d omega / 2 pi) S(omega)
    if len(omega_grid) > 1:
        integrated_w = float(np.trapezoid(spectral_density, omega_grid) / (2.0 * np.pi))
    else:
        integrated_w = 0.0

    return ContinuousSpectrum(
        omega=omega_grid,
        spectral_density=spectral_density,
        broadening_type=broadening_type,
        broadening_param=broadening_param,
        integrated_weight=integrated_w,
    )


def verify_spectral_sum_rule(
    spectral_data: DiscreteSpectralData,
    continuous_spec: ContinuousSpectrum | None = None,
    tol: float = 0.05,
) -> dict[str, Any]:
    """Verify completeness sum rules for discrete and continuous spectral representations.

    Returns:
        dict with metrics and boolean compliance.
    """
    discrete_sum = spectral_data.total_spectral_weight
    is_discrete_valid = bool(abs(discrete_sum - 1.0) < tol)

    result: dict[str, Any] = {
        "discrete_sum": discrete_sum,
        "is_discrete_valid": is_discrete_valid,
        "num_poles": len(spectral_data.poles),
        "lowest_excitation_energy": spectral_data.poles[0].energy if spectral_data.poles else 0.0,
        "highest_excitation_energy": spectral_data.poles[-1].energy if spectral_data.poles else 0.0,
    }

    if continuous_spec is not None:
        cont_integral = continuous_spec.integrated_weight
        diff_integral = abs(cont_integral - discrete_sum)
        result["continuous_integral"] = cont_integral
        result["quadrature_difference"] = diff_integral
        result["is_continuous_valid"] = bool(abs(cont_integral - 1.0) < tol)

    return result


def evaluate_thermal_spectral_function(
    spectral_data: DiscreteSpectralData,
    omega_grid: np.ndarray,
    temperature: float = 0.0,
    broadening_param: float = 0.06,
    broadening_type: str = "gaussian",
) -> ContinuousSpectrum:
    """Evaluate finite-temperature dynamical spin correlation spectrum S_{jj}(omega).

    Incorporates:
    - Effective thermal broadening sigma_eff(T) = sqrt(sigma_0^2 + T^2).
    - Kubo-Martin-Schwinger (KMS) detailed balance factor:
      S_T(omega) = 1 / (1 + exp(-beta * omega)) * S_0(|omega|)
      satisfying S_T(-omega) = exp(-beta * omega) S_T(omega) and total sum rule conservation.

    Parameters:
        spectral_data: Discrete poles and weights from gauge quench.
        omega_grid: Frequency grid (supports both negative and positive frequencies).
        temperature: Temperature T >= 0 in units of |K|.
        broadening_param: Zero-temperature broadening parameter (sigma_0 or eta_0).
        broadening_type: 'gaussian' or 'lorentzian'.
    """
    if temperature < 0.0:
        raise ValueError(f"Temperature must be non-negative, got {temperature}")

    if temperature <= 1e-12:
        # Standard zero-temperature limit (vanishing for omega < 0)
        s_pos = evaluate_spectral_function(
            spectral_data,
            omega_grid=omega_grid,
            broadening_type=broadening_type,
            broadening_param=broadening_param,
        )
        density = np.where(omega_grid >= 0.0, s_pos.spectral_density, 0.0)
        eff_broadening = broadening_param
    else:
        beta = 1.0 / temperature
        eff_broadening = (
            float(np.sqrt(broadening_param**2 + temperature**2))
            if broadening_type == "gaussian"
            else float(broadening_param + np.pi * temperature)
        )
        # Positive and negative frequency absorption profiles: S_sym(omega) = S_pos(omega) + S_pos(-omega)
        s_pos = evaluate_spectral_function(
            spectral_data,
            omega_grid=omega_grid,
            broadening_type=broadening_type,
            broadening_param=eff_broadening,
        )
        s_neg = evaluate_spectral_function(
            spectral_data,
            omega_grid=-omega_grid,
            broadening_type=broadening_type,
            broadening_param=eff_broadening,
        )
        s_sym = s_pos.spectral_density + s_neg.spectral_density

        # KMS thermal factor with numerically stable clipping: S_T(omega) = 1 / (1 + exp(-beta * omega)) * S_sym(omega)
        arg = -beta * omega_grid
        arg_clipped = np.clip(arg, -100.0, 100.0)
        kms_factor = 1.0 / (1.0 + np.exp(arg_clipped))
        density = kms_factor * s_sym

    if len(omega_grid) > 1:
        integrated_w = float(np.trapezoid(density, omega_grid) / (2.0 * np.pi))
    else:
        integrated_w = 0.0

    return ContinuousSpectrum(
        omega=omega_grid,
        spectral_density=density,
        broadening_type=f"{broadening_type}_thermal",
        broadening_param=eff_broadening,
        integrated_weight=integrated_w,
    )


def verify_kms_detailed_balance(
    spec: ContinuousSpectrum,
    temperature: float,
    omega_test_points: Sequence[float] | None = None,
    tol: float = 1e-10,
) -> dict[str, Any]:
    """Verify KMS detailed balance relation S(-omega) = exp(-beta * omega) S(omega)."""
    if temperature <= 1e-12:
        return {"is_valid": True, "temperature": 0.0, "max_relative_error": 0.0}

    beta = 1.0 / temperature
    if omega_test_points is None:
        omega_test_points = [0.05, 0.10, 0.20, 0.30, 0.50, 1.00]

    max_err = 0.0
    point_errors = {}
    for w in omega_test_points:
        s_pos = float(np.interp(w, spec.omega, spec.spectral_density))
        s_neg = float(np.interp(-w, spec.omega, spec.spectral_density))
        expected_ratio = float(np.exp(-beta * w))
        if s_pos > 1e-12:
            ratio = s_neg / s_pos
            err = abs(ratio - expected_ratio)
            point_errors[w] = err
            if err > max_err:
                max_err = err

    return {
        "is_valid": bool(max_err < tol),
        "temperature": temperature,
        "max_relative_error": max_err,
        "point_errors": point_errors,
    }
