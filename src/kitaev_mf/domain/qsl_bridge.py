"""Cross-QSL conventions-v2.0 SO(4) / W_p / (u_ij, v_ij) convention bridge.

Provides the opt-in algebraic and state-port bridge between `kitaev-mf`'s internal
4-Majorana convention ``c_j = (-c_j, b_j^x, b_j^y, b_j^z)^T`` (locked in
``src/kitaev_mf/domain/constants.py``) and the QSL-wide ``conventions-v2.0`` standard
(``research/qsl/contracts/qsl_conventions.json`` and ``ksl-sc-bilayer``) using
``gamma_j = (+c_j, b_j^x, b_j^y, b_j^z)^T`` and universal Cartesian-keyed generators
``SO4_L['x'|'y'|'z']``, ``SO4_R['x'|'y'|'z']``.

Key exact identities (verified on 16-state 2-site Fock and 4096-state 6-site plaquette spaces):
1. Basis parity involution: ``P_0 = diag(-1, 1, 1, 1)``, ``c_j = P_0 gamma_j``,
   preserving local Kitaev parity ``D_j = +c_j^0 c_j^1 c_j^2 c_j^3 = +b_j^x b_j^y b_j^z c_j = +1``.
2. Left <-> Right swap under ``P_0``:
   ``P_0 SO4_L[a] P_0 = SO4_R[a]``, ``P_0 SO4_R[a] P_0 = SO4_L[a]``,
   with ``+0.5 * M_MATRICES[idx] == SO4_R[ax]`` (physical spin in ``(-c, b)``),
   ``-0.5 * G_MATRICES[idx] == SO4_L[ax]`` (standard-sign parton gauge in ``(-c, b)``),
   and ``N_MATRICES[idx] == SO4_R[ax] - SO4_L[ax]``.
3. Lagrange multiplier bridge:
   ``+lambda_kmf^a G^a`` in ``H = (i/4) c^T A c`` equals ``+(-2 lambda_kmf^a) SO4_R[ax]``
   in ``H = (i/4) gamma^T A_v2 gamma``, i.e. ``lambda_v2 = -2 * lambda_kmf``.
4. Pure-Kitaev bond scalars ``(u_ij^alpha, v_ij^alpha)`` and plaquette Wilson loops ``W_p``:
   ``u_ij^alpha = <i b_i^alpha b_j^alpha> = U_b^{gamma, gamma}`` and
   ``v_ij^alpha = <i c_i c_j> = U_b^{0, 0}`` are invariant under ``P_0`` while mixed
   matter-gauge blocks ``U_b^{0, a}`` and ``U_b^{a, 0}`` negate.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import isfinite, sqrt
from types import MappingProxyType
from typing import Final, Literal

import numpy as np
from numpy.typing import NDArray

from kitaev_mf.domain.lattice import UnitCell
from kitaev_mf.domain.state import InhomogeneousMeanFieldState, MeanFieldState

MajoranaBasis = Literal["kitaev_mf", "conventions_v2"]

CONVENTIONS_VERSION: Final[str] = "conventions-v2.0"
CARTESIAN_AXES: Final[tuple[str, str, str]] = ("x", "y", "z")
AXIS_TO_BOND_TYPE: Final[Mapping[str, int]] = MappingProxyType({"x": 1, "y": 2, "z": 3})
BOND_TYPE_TO_AXIS: Final[Mapping[int, str]] = MappingProxyType({1: "x", 2: "y", 3: "z"})


def _make_readonly_array(values: NDArray[np.float64]) -> NDArray[np.float64]:
    arr = np.array(values, dtype=np.float64, copy=True)
    arr.setflags(write=False)
    return arr


def _build_p0_parity_bridge() -> NDArray[np.float64]:
    return _make_readonly_array(np.diag(np.array([-1.0, 1.0, 1.0, 1.0], dtype=np.float64)))


def _build_so4_lr_generators() -> tuple[
    Mapping[str, NDArray[np.float64]],
    Mapping[str, NDArray[np.float64]],
]:
    """Construct immutable universal 4x4 real skew-symmetric SO(4) generators (L^a, R^a).

    In the canonical quartet ``q = (q^0, q^x, q^y, q^z)^T``:
      (L^a)_{0,a} = -1/2, (L^a)_{a,0} = +1/2, (L^a)_{b,c} = -(1/2) eps_{abc},
      (R^a)_{0,a} = +1/2, (R^a)_{a,0} = -1/2, (R^a)_{b,c} = -(1/2) eps_{abc},
    satisfying ``[L^a, L^b] = eps_{abc} L^c``, ``[R^a, R^b] = eps_{abc} R^c``,
    and ``[L^a, R^b] = 0``.
    """
    left_mats: dict[str, NDArray[np.float64]] = {}
    right_mats: dict[str, NDArray[np.float64]] = {}

    for a_idx, ax in enumerate(CARTESIAN_AXES, start=1):
        mat_l = np.zeros((4, 4), dtype=np.float64)
        mat_r = np.zeros((4, 4), dtype=np.float64)

        mat_l[0, a_idx] = -0.5
        mat_l[a_idx, 0] = +0.5
        mat_r[0, a_idx] = +0.5
        mat_r[a_idx, 0] = -0.5

        for b_idx in (1, 2, 3):
            for c_idx in (1, 2, 3):
                eps = (b_idx - c_idx) * (c_idx - a_idx) * (a_idx - b_idx) // 2
                if eps != 0:
                    val = -0.5 * float(eps)
                    mat_l[b_idx, c_idx] = val
                    mat_r[b_idx, c_idx] = val

        mat_l.setflags(write=False)
        mat_r.setflags(write=False)
        left_mats[ax] = mat_l
        right_mats[ax] = mat_r

    return MappingProxyType(left_mats), MappingProxyType(right_mats)


P0_PARITY_BRIDGE: Final[NDArray[np.float64]] = _build_p0_parity_bridge()
SO4_L, SO4_R = _build_so4_lr_generators()
SO4_L_STACK: Final[NDArray[np.float64]] = _make_readonly_array(
    np.stack([SO4_L[ax] for ax in CARTESIAN_AXES], axis=0)
)
SO4_R_STACK: Final[NDArray[np.float64]] = _make_readonly_array(
    np.stack([SO4_R[ax] for ax in CARTESIAN_AXES], axis=0)
)


def _validate_finite_scalar(name: str, value: float) -> float:
    val = float(value)
    if not isfinite(val):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return val


def axis_to_bond_type(axis: str) -> int:
    """Map Cartesian axis label ``'x'|'y'|'z'`` to legacy integer ``bond_type`` in ``{1, 2, 3}``."""
    return AXIS_TO_BOND_TYPE[axis]


def bond_type_to_axis(bond_type: int) -> str:
    """Map legacy integer ``bond_type`` in ``{1, 2, 3}`` to Cartesian axis ``'x'|'y'|'z'``."""
    return BOND_TYPE_TO_AXIS[bond_type]


def kmf_spin_generator(axis: str) -> NDArray[np.float64]:
    """Return physical-spin generator matrix in ``kitaev-mf`` ``(-c, b)`` basis: ``+0.5 * M^a == R^alpha``."""
    return SO4_R[axis]


def kmf_gauge_generator(axis: str) -> NDArray[np.float64]:
    """Return standard-sign parton gauge generator matrix in ``(-c, b)`` basis: ``-0.5 * G^a == L^alpha``."""
    return SO4_L[axis]


def pauli_bilinear_matrix(
    axis: str,
    *,
    basis: MajoranaBasis = "kitaev_mf",
) -> NDArray[np.float64]:
    """Return 4x4 skew matrix ``N_basis^alpha`` such that ``sigma_j^alpha = (i/2) q_j^T N_basis^alpha q_j``.

    In ``kitaev_mf`` ``(-c, b)`` basis: ``SO4_R[axis] - SO4_L[axis] == +N_MATRICES[a]``.
    In ``conventions_v2`` ``(+c, b)`` basis: ``SO4_L[axis] - SO4_R[axis] == -N_MATRICES[a]``.
    """
    if basis == "kitaev_mf":
        return np.asarray(SO4_R[axis] - SO4_L[axis], dtype=np.float64)
    if basis == "conventions_v2":
        return np.asarray(SO4_L[axis] - SO4_R[axis], dtype=np.float64)
    raise ValueError(
        f"Unsupported Majorana basis {basis!r}; expected 'kitaev_mf' or 'conventions_v2'"
    )


def conjugate_site_block_p0(mat: NDArray[np.float64]) -> NDArray[np.float64]:
    """Apply the exact parity involution ``P_0 M P_0`` to any ``(..., 4, 4)`` real site/bond tensor.

    Negates the mixed matter-gauge row/column entries ``(0, 1..3)`` and ``(1..3, 0)``
    while preserving ``(0, 0)`` and the ``(1..3, 1..3)`` gauge subblock identically.
    """
    out: NDArray[np.float64] = np.array(mat, dtype=np.float64, copy=True)
    if out.ndim < 2 or out.shape[-2:] != (4, 4):
        raise ValueError(f"Expected array with trailing shape (4, 4), got {out.shape}")
    out[..., 0, 1:] *= -1.0
    out[..., 1:, 0] *= -1.0
    return out


def conjugate_multisite_p0(mat: NDArray[np.float64]) -> NDArray[np.float64]:
    """Apply ``P_0^{oplus N} M P_0^{oplus N}`` to any ``(4N, 4N)`` real Majorana matrix."""
    out: NDArray[np.float64] = np.array(mat, dtype=np.float64, copy=True)
    if out.ndim != 2 or out.shape[0] != out.shape[1] or out.shape[0] % 4 != 0:
        raise ValueError(f"Expected square (4N, 4N) matrix, got shape {out.shape}")
    out[0::4, :] *= -1.0
    out[:, 0::4] *= -1.0
    return out


def lambdas_kmf_to_v2(lambdas_kmf: NDArray[np.float64]) -> NDArray[np.float64]:
    """Convert ``kitaev-mf`` multipliers (``+lam_kmf^a G^a``) to ``conventions-v2.0`` (``+lam_v2^a R^a``)."""
    arr = np.asarray(lambdas_kmf, dtype=np.float64)
    if arr.ndim < 1 or arr.shape[-1] != 3:
        raise ValueError(f"Expected trailing axis length 3 for lambdas, got shape {arr.shape}")
    return np.asarray(-2.0 * arr, dtype=np.float64)


def lambdas_v2_to_kmf(lambdas_v2: NDArray[np.float64]) -> NDArray[np.float64]:
    """Convert ``conventions-v2.0`` multipliers (``+lam_v2^a R^a``) back to ``kitaev-mf`` (``+lam_kmf^a G^a``)."""
    arr = np.asarray(lambdas_v2, dtype=np.float64)
    if arr.ndim < 1 or arr.shape[-1] != 3:
        raise ValueError(f"Expected trailing axis length 3 for lambdas, got shape {arr.shape}")
    return np.asarray(-0.5 * arr, dtype=np.float64)


@dataclass(frozen=True)
class QSLMeanFieldPortState:
    """Translationally invariant mean-field state in canonical ``conventions-v2.0`` ``(+c, b)`` basis.

    Attributes:
        U_v2: Bond correlators ``<i gamma_i^mu gamma_j^nu>`` of shape ``(num_bonds, 4, 4)``.
        V_v2: On-site skew correlators ``<i gamma_s^mu gamma_s^nu>`` of shape ``(num_sub, 4, 4)``.
        lambdas_v2: Gauge Lagrange multipliers coupling via ``+sum_a lambda_v2^a G_s^{v2,a}``,
            shape ``(num_sub, 3)``.
        conventions_version: Locked port contract tag (default ``'conventions-v2.0'``).
    """

    U_v2: NDArray[np.float64]
    V_v2: NDArray[np.float64]
    lambdas_v2: NDArray[np.float64]
    conventions_version: str = CONVENTIONS_VERSION

    def __post_init__(self) -> None:
        if self.U_v2.ndim != 3 or self.U_v2.shape[1:] != (4, 4):
            raise ValueError(f"U_v2 must have shape (num_bonds, 4, 4), got {self.U_v2.shape}")
        if self.V_v2.ndim != 3 or self.V_v2.shape[1:] != (4, 4):
            raise ValueError(f"V_v2 must have shape (num_sublattices, 4, 4), got {self.V_v2.shape}")
        if self.lambdas_v2.ndim != 2 or self.lambdas_v2.shape[1] != 3:
            raise ValueError(
                f"lambdas_v2 must have shape (num_sublattices, 3), got {self.lambdas_v2.shape}"
            )
        if self.V_v2.shape[0] != self.lambdas_v2.shape[0]:
            raise ValueError(
                f"Sublattice mismatch: V_v2 ({self.V_v2.shape[0]}) != lambdas_v2 ({self.lambdas_v2.shape[0]})"
            )
        skew_err = float(np.max(np.abs(self.V_v2 + np.swapaxes(self.V_v2, -1, -2))))
        if skew_err > 1e-12:
            raise ValueError(f"V_v2 must be skew-symmetric, max error: {skew_err}")


@dataclass(frozen=True)
class QSLInhomogeneousPortState:
    """Real-space inhomogeneous mean-field state in canonical ``conventions-v2.0`` ``(+c, b)`` basis."""

    U_v2: NDArray[np.float64]
    V_v2: NDArray[np.float64]
    lambdas_v2: NDArray[np.float64]
    conventions_version: str = CONVENTIONS_VERSION

    def __post_init__(self) -> None:
        if self.U_v2.ndim != 4 or self.U_v2.shape[1:] != (3, 4, 4):
            raise ValueError(f"U_v2 must have shape (num_cells, 3, 4, 4), got {self.U_v2.shape}")
        if self.V_v2.ndim != 4 or self.V_v2.shape[1:] != (2, 4, 4):
            raise ValueError(f"V_v2 must have shape (num_cells, 2, 4, 4), got {self.V_v2.shape}")
        if self.lambdas_v2.ndim != 3 or self.lambdas_v2.shape[1:] != (2, 3):
            raise ValueError(
                f"lambdas_v2 must have shape (num_cells, 2, 3), got {self.lambdas_v2.shape}"
            )
        num_cells = self.U_v2.shape[0]
        if self.V_v2.shape[0] != num_cells or self.lambdas_v2.shape[0] != num_cells:
            raise ValueError(
                f"Cell count mismatch: U_v2 ({num_cells}), V_v2 ({self.V_v2.shape[0]}), "
                f"lambdas_v2 ({self.lambdas_v2.shape[0]})"
            )
        skew_err = float(np.max(np.abs(self.V_v2 + np.swapaxes(self.V_v2, -1, -2))))
        if skew_err > 1e-12:
            raise ValueError(f"V_v2 must be skew-symmetric, max error: {skew_err}")


@dataclass(frozen=True)
class SiteSU2Decomposition:
    """Single-site expectation values of the commuting ``SU(2)_L x SU(2)_R`` generators.

    Attributes:
        spin_s: Physical spin expectation ``<S_j> = <P_+ (sigma_j / 2) P_+>`` in ``R^3``.
        gauge_g: Standard-sign parton gauge residual ``<G_j^{v2}> = -<T_j^{kmf}>`` in ``R^3``
            (strictly zero on physical subspace ``D_j = +1``).
        pauli_sigma_unproj: Unprojected Kitaev bilinear ``<sigma_j> = 2(<S_j> - <G_j^{v2}>)`` in ``R^3``.
    """

    spin_s: NDArray[np.float64]
    gauge_g: NDArray[np.float64]
    pauli_sigma_unproj: NDArray[np.float64]


@dataclass(frozen=True)
class PureKitaevBondScalars:
    """Canonical ``conventions-v2.0`` pure-Kitaev directed bond scalars ``(u_ij^alpha, v_ij^alpha)``.

    Attributes:
        u_gauge: Directed ``A -> B`` gauge link expectations ``u_ij^alpha = <i b_i^alpha b_j^alpha>``.
        v_matter: Directed ``A -> B`` matter bond expectations ``v_ij^alpha = <i c_i c_j>``.
        bond_axes: Cartesian bond axis ``'x'|'y'|'z'`` per directed bond.
        off_diagonal_norm: Frobenius norm of all non-pure-Kitaev entries in ``U_b``
            (vanishes identically in pure ``h = 0`` Kitaev states).
    """

    u_gauge: NDArray[np.float64]
    v_matter: NDArray[np.float64]
    bond_axes: tuple[str, ...]
    off_diagonal_norm: float


def convert_state_kmf_to_v2(state: MeanFieldState) -> QSLMeanFieldPortState:
    """Losslessly map a ``kitaev-mf`` ``MeanFieldState`` into ``conventions-v2.0`` ``(+c, b)`` form."""
    return QSLMeanFieldPortState(
        U_v2=conjugate_site_block_p0(state.U),
        V_v2=conjugate_site_block_p0(state.V),
        lambdas_v2=lambdas_kmf_to_v2(state.lambdas),
    )


def convert_state_v2_to_kmf(port_state: QSLMeanFieldPortState) -> MeanFieldState:
    """Losslessly map a ``conventions-v2.0`` ``QSLMeanFieldPortState`` back to ``kitaev-mf``."""
    return MeanFieldState(
        U=conjugate_site_block_p0(port_state.U_v2),
        V=conjugate_site_block_p0(port_state.V_v2),
        lambdas=lambdas_v2_to_kmf(port_state.lambdas_v2),
    )


def convert_inhomogeneous_state_kmf_to_v2(
    state: InhomogeneousMeanFieldState,
) -> QSLInhomogeneousPortState:
    """Losslessly map an ``InhomogeneousMeanFieldState`` into ``conventions-v2.0`` ``(+c, b)`` form."""
    return QSLInhomogeneousPortState(
        U_v2=conjugate_site_block_p0(state.U),
        V_v2=conjugate_site_block_p0(state.V),
        lambdas_v2=lambdas_kmf_to_v2(state.lambdas),
    )


def convert_inhomogeneous_state_v2_to_kmf(
    port_state: QSLInhomogeneousPortState,
) -> InhomogeneousMeanFieldState:
    """Losslessly map a ``QSLInhomogeneousPortState`` back to ``InhomogeneousMeanFieldState``."""
    return InhomogeneousMeanFieldState(
        U=conjugate_site_block_p0(port_state.U_v2),
        V=conjugate_site_block_p0(port_state.V_v2),
        lambdas=lambdas_v2_to_kmf(port_state.lambdas_v2),
    )


def extract_site_su2_expectations(
    v_site: NDArray[np.float64],
    *,
    basis: MajoranaBasis = "kitaev_mf",
) -> SiteSU2Decomposition:
    """Extract basis-invariant ``<S_j>``, ``<G_j^{v2}>``, and ``<sigma_j>_unproj`` from a 4x4 skew covariance.

    Uses ``<O(M; q)> = -(1/4) Tr(M V)`` for ``V_{mu,nu} = <i q_mu q_nu>``.
    """
    v_arr = np.asarray(v_site, dtype=np.float64)
    if v_arr.shape != (4, 4):
        raise ValueError(f"Expected 4x4 on-site covariance matrix, got shape {v_arr.shape}")
    skew_err = float(np.max(np.abs(v_arr + v_arr.T)))
    if skew_err > 1e-12:
        raise ValueError(f"On-site covariance matrix must be skew-symmetric, max error: {skew_err}")

    if basis == "kitaev_mf":
        spin_mats, gauge_mats = SO4_R, SO4_L
    elif basis == "conventions_v2":
        spin_mats, gauge_mats = SO4_L, SO4_R
    else:
        raise ValueError(
            f"Unsupported Majorana basis {basis!r}; expected 'kitaev_mf' or 'conventions_v2'"
        )

    spin_s = np.array(
        [-0.25 * float(np.trace(spin_mats[ax] @ v_arr)) for ax in CARTESIAN_AXES],
        dtype=np.float64,
    )
    gauge_g = np.array(
        [-0.25 * float(np.trace(gauge_mats[ax] @ v_arr)) for ax in CARTESIAN_AXES],
        dtype=np.float64,
    )
    pauli_sigma = 2.0 * (spin_s - gauge_g)
    return SiteSU2Decomposition(
        spin_s=spin_s,
        gauge_g=gauge_g,
        pauli_sigma_unproj=pauli_sigma,
    )


def extract_pure_kitaev_bond_scalars(
    state: MeanFieldState | QSLMeanFieldPortState,
    cell: UnitCell,
) -> PureKitaevBondScalars:
    """Extract directed ``A -> B`` bond scalars ``u_ij^alpha`` and ``v_ij^alpha`` from a state.

    Because both ``(0, 0)`` and ``(gamma, gamma)`` entries are invariant under ``P_0``,
    the extracted scalars are identical across ``kitaev_mf`` and ``conventions_v2`` states.
    """
    u_tensor = state.U if isinstance(state, MeanFieldState) else state.U_v2
    num_bonds = len(cell.bonds)
    if u_tensor.shape[0] != num_bonds:
        raise ValueError(
            f"Bond count mismatch: state has {u_tensor.shape[0]} bonds, cell has {num_bonds}"
        )

    u_gauge = np.empty(num_bonds, dtype=np.float64)
    v_matter = np.empty(num_bonds, dtype=np.float64)
    bond_axes: list[str] = []
    residual_sq = 0.0

    for b_idx, bond in enumerate(cell.bonds):
        gamma = bond.bond_type
        u_b = u_tensor[b_idx]
        u_val = float(u_b[gamma, gamma])
        v_val = float(u_b[0, 0])
        u_gauge[b_idx] = u_val
        v_matter[b_idx] = v_val
        bond_axes.append(bond_type_to_axis(gamma))

        rem = u_b.copy()
        rem[0, 0] = 0.0
        rem[gamma, gamma] = 0.0
        residual_sq += float(np.sum(rem * rem))

    return PureKitaevBondScalars(
        u_gauge=u_gauge,
        v_matter=v_matter,
        bond_axes=tuple(bond_axes),
        off_diagonal_norm=float(sqrt(residual_sq)),
    )


def reconstruct_pure_kitaev_bond_tensors(
    u_gauge: NDArray[np.float64],
    v_matter: NDArray[np.float64],
    cell: UnitCell,
) -> NDArray[np.float64]:
    """Reconstruct pure-Kitaev diagonal ``(num_bonds, 4, 4)`` bond correlators from ``(u_ij, v_ij)``."""
    u_arr = np.asarray(u_gauge, dtype=np.float64)
    v_arr = np.asarray(v_matter, dtype=np.float64)
    num_bonds = len(cell.bonds)
    if u_arr.shape != (num_bonds,) or v_arr.shape != (num_bonds,):
        raise ValueError(
            f"Expected 1D arrays of length {num_bonds}, got u={u_arr.shape}, v={v_arr.shape}"
        )

    out = np.zeros((num_bonds, 4, 4), dtype=np.float64)
    for b_idx, bond in enumerate(cell.bonds):
        gamma = bond.bond_type
        out[b_idx, 0, 0] = v_arr[b_idx]
        out[b_idx, gamma, gamma] = u_arr[b_idx]
    return out


def compute_plaquette_wilson_loops_from_u(
    u_gauge: NDArray[np.float64],
    cell: UnitCell,
) -> NDArray[np.float64]:
    """Compute Kitaev plaquette Wilson loops ``W_p = prod_{b in partial p} u_b`` from directed ``A -> B`` links.

    Matches the exact 6-site projected spin plaquette identity
    ``P_+ W_p^{spin} P_+ = P_+ (prod_{b in partial p} i b_{A_b}^{alpha_b} b_{B_b}^{alpha_b}) P_+``
    with positive sign because the three ``A -> B`` and three ``B -> A`` hexagon steps yield
    ``(-i)^3 (+i)^3 = +1``.
    """
    u_arr = np.asarray(u_gauge, dtype=np.float64)
    num_bonds = len(cell.bonds)
    if u_arr.shape != (num_bonds,):
        raise ValueError(f"Expected u_gauge shape ({num_bonds},), got {u_arr.shape}")

    fluxes = np.ones(len(cell.plaquettes), dtype=np.float64)
    for p_idx, plaq in enumerate(cell.plaquettes):
        b_indices = (
            [b_idx for b_idx, _ in plaq.oriented_bonds]
            if plaq.oriented_bonds
            else list(plaq.bond_indices)
        )
        prod_val = 1.0
        for b_idx in b_indices:
            prod_val *= float(u_arr[b_idx])
        fluxes[p_idx] = prod_val
    return fluxes


def kitaev_coupling_pauli_to_bilayer(k_pauli: float) -> float:
    """Map external QSL ``+K_Pauli sum sigma.sigma`` (FM ``-1``) to ``ksl-sc-bilayer`` ``-K_bilayer`` (FM ``+1``)."""
    return -_validate_finite_scalar("k_pauli", k_pauli)


def kitaev_coupling_bilayer_to_pauli(k_bilayer: float) -> float:
    """Map ``ksl-sc-bilayer`` ``-K_bilayer sum sigma.sigma`` (FM ``+1``) to ``kitaev-mf`` ``+K_Pauli`` (FM ``-1``)."""
    return -_validate_finite_scalar("k_bilayer", k_bilayer)


def chiral_kappa_kwg_to_kmf_static(kappa_kwg: float) -> float:
    """Map ``kitaev-walk-geometry`` ``kappa_KWG`` to ``kitaev-mf`` static ``kappa_static = kappa_KWG / 2``."""
    return 0.5 * _validate_finite_scalar("kappa_kwg", kappa_kwg)


def chiral_kappa_kmf_static_to_kwg(kappa_static: float) -> float:
    """Map ``kitaev-mf`` static ``kappa_static`` to ``kitaev-walk-geometry`` ``kappa_KWG = 2 * kappa_static``."""
    return 2.0 * _validate_finite_scalar("kappa_static", kappa_static)


def chiral_kappa_kmf_static_to_kslsc(kappa_static: float) -> float:
    """Map ``kitaev-mf`` ``kappa_static`` to ``ksl-sc-bilayer`` ``kappa_KSLSC`` at matched full BdG gap.

    At ``kappa_KSLSC = kappa_static = kappa_KWG / 2``:
      ``Delta_can^MF = 12 sqrt(3) |kappa_static| = 2 Delta_KSL = Delta_K^KWG = 6 sqrt(3) |kappa_KWG|``.
    """
    return _validate_finite_scalar("kappa_static", kappa_static)


def chiral_full_gap_canonical(kappa_static: float) -> float:
    """Return canonical full particle-hole Dirac mass gap ``12 sqrt(3) |kappa_static|``."""
    return 12.0 * sqrt(3.0) * abs(_validate_finite_scalar("kappa_static", kappa_static))


class DefaultQSLConventionBridge:
    """Pure-domain implementation of ``QSLConventionBridgePort`` for ``conventions-v2.0`` interop."""

    @property
    def conventions_version(self) -> str:
        """Locked cross-QSL contract version identifier."""
        return CONVENTIONS_VERSION

    def to_v2_state(self, state: MeanFieldState) -> QSLMeanFieldPortState:
        """Convert ``MeanFieldState`` from ``(-c, b)`` to ``conventions-v2.0`` ``(+c, b)`` basis."""
        return convert_state_kmf_to_v2(state)

    def from_v2_state(self, port_state: QSLMeanFieldPortState) -> MeanFieldState:
        """Convert ``QSLMeanFieldPortState`` from ``(+c, b)`` back to ``kitaev-mf`` ``(-c, b)`` basis."""
        return convert_state_v2_to_kmf(port_state)

    def extract_bond_scalars(
        self,
        state: MeanFieldState | QSLMeanFieldPortState,
        cell: UnitCell,
    ) -> PureKitaevBondScalars:
        """Extract directed ``A -> B`` pure-Kitaev scalars ``(u_ij^alpha, v_ij^alpha)``."""
        return extract_pure_kitaev_bond_scalars(state, cell)

    def compute_plaquette_wilson_loops(
        self,
        u_gauge: NDArray[np.float64],
        cell: UnitCell,
    ) -> NDArray[np.float64]:
        """Compute plaquette Wilson loops ``W_p`` from directed gauge-link scalars ``u_ij^alpha``."""
        return compute_plaquette_wilson_loops_from_u(u_gauge, cell)

    def extract_site_su2(
        self,
        v_site: NDArray[np.float64],
        *,
        basis: MajoranaBasis = "kitaev_mf",
    ) -> SiteSU2Decomposition:
        """Extract basis-invariant ``(<S_j>, <G_j^{v2}>, <sigma_j>_unproj)`` from a 4x4 skew covariance."""
        return extract_site_su2_expectations(v_site, basis=basis)
