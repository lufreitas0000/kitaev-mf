---
id: RND-P1-M1.10-T1-K019A
title: "K-019-A — Second-order perturbation erratum and exact-track reuse"
phase: 1
milestone: M1.10
task: T1
status: draft
questions:
  - "Why does a vanishing first-order projected field term not imply a vanishing second-order effective Hamiltonian?"
  - "Which projector and resolvent domains are needed before the third-order chiral-mass calculation?"
depends_on:
  - M1.10-T1
contracts:
  - K-019
artifacts:
  - doc/rnd/phase_1/m1_10_t1_exact_perturbation_and_normalization_reconciliation.md
  - doc/rnd/exact-perturbation-theory/STATUS.md
  - doc/rnd/exact-perturbation-theory/notes/01_flux_sectors_and_spectrum.md
---

# K-019-A — Exact second-order erratum and prior-work inventory

**Disposition:** Mathematical correction proposed for independent review; not a new hypothesis, perturbative coefficient, or a change to the legacy K-019 historical status. **Owner assertion (2026-11-14):** the old implication is false, and the existing dedicated perturbation research must be reused before new work. Source snapshot: [kitaev-mf main e093e417](https://github.com/lufreitas0000/kitaev-mf/tree/e093e417494de7b7cb78ae415338189d628a8194).

## 1. Definite mathematical error

The [historical M1.10-T1 note, section 2](https://github.com/lufreitas0000/kitaev-mf/blob/e093e417494de7b7cb78ae415338189d628a8194/doc/rnd/phase_1/m1_10_t1_exact_perturbation_and_normalization_reconciliation.md#L43-L51) states that the flux selection rule $PVP=0$ implies $H_{\mathrm{eff}}^{(1)}=H_{\mathrm{eff}}^{(2)}=0$. This inference is **false**. It must not be used in subsequent derivations, regardless of the historical note's "Completed" metadata.

For a finite-dimensional self-adjoint $H_0$, let $P$ be its entire $E_0$ ground eigenspace, $Q=I-P$, and assume the reduced gap $QH_0Q-E_0 Q$ is strictly positive on $\operatorname{Ran}Q$. At an energy-independent Rayleigh–Schrödinger expansion around $E_0$,

$$
H_{\mathrm{eff}}^{(1)}=PVP,\qquad
H_{\mathrm{eff}}^{(2)}
=PVQ(E_0-QH_0Q)^{-1}QVP
=-PVQ(QH_0Q-E_0)^{-1}QVP.
$$

Hence for any $|\psi\rangle\in\operatorname{Ran}P$,

$$
\langle\psi|H_{\mathrm{eff}}^{(2)}|\psi\rangle
=-\left\|(QH_0Q-E_0)^{-1/2}QV|\psi\rangle\right\|^2\le0.
$$

The expectation is strictly negative whenever $QV|\psi\rangle\ne0$. Thus **$PVP=0$ does not imply $H_{\mathrm{eff}}^{(2)}=0$**. This is a mathematical statement about the stipulated finite-dimensional gapped block, not a claim that all second-order operators in the Kitaev model have been computed.

An explicit two-level counterexample suffices: take

$$
H_0=\begin{pmatrix}0&0\\0&\Delta\end{pmatrix},\quad
V=v\begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
P=\begin{pmatrix}1&0\\0&0\end{pmatrix},\qquad
\Delta>0,\ v\ne0.
$$

Then $PVP=0$ while $H_{\mathrm{eff}}^{(2)}=-(v^2/\Delta)P\ne0$.

## 2. What still requires a Kitaev-specific derivation

- **Projector choice:** The local physical gauge projector $P_G$, the full physical flux-free sector $P_{\mathrm{ff}}$, and the ground-energy eigenspace projector $P_{E_0}$ are different objects. The preceding simple reduced-gap formula is not automatically a valid energy-independent operator across the entire flux-free *matter continuum*. For an arbitrary retained block $P$, use an energy-dependent Feshbach map only for $E\in\rho(QH_0Q)$.
- **Order and symmetry:** A first-order flux-selection rule may eliminate $P_{\mathrm{ff}}VP_{\mathrm{ff}}$. It does not eliminate even-order virtual returns. A time-reversal-even second-order energy shift or other term need not generate a time-reversal-odd chiral mass. The onset of the targeted mass at third order requires operator-specific selection rules; it cannot follow from $PVP=0$ alone.
- **Third order:** The three-spin/next-nearest-neighbor matter hopping requires a consistent local physical constraint, ordered gauge links, site orientation and projected Hilbert space. A fixed-gauge algebraic identity does not compute its resolvent-weighted coefficient.
- **Normalizations:** The spectrum of the production $iA/4$ matrix is scaled relative to the canonical $iA$ generator, but the two-band separation, positive canonical quasiparticle energy and projected physical-spin gap must remain distinct. No numeric coefficient is certified here.

The [merged Wave 2 pilot, sections 4–5](https://github.com/lufreitas0000/kitaev-mf/blob/e093e417494de7b7cb78ae415338189d628a8194/doc/rnd/phase_1/m1_10_t1_wave2_cubic_gap_evidence_pilot.md#L88-L170) already records these qualifications. They are not new discoveries claimed by this erratum.

## 3. Reuse the existing perturbation branch instead of restarting

The dedicated branch [research/pure-kitaev-exact-perturbation-20261102](https://github.com/lufreitas0000/kitaev-mf/tree/research/pure-kitaev-exact-perturbation-20261102) is, at this snapshot, an **ancestor of main**: behind by seven commits, ahead by zero, with no branch-exclusive file changes. Its research is retained on main:

- [STATUS.md](https://github.com/lufreitas0000/kitaev-mf/blob/e093e417494de7b7cb78ae415338189d628a8194/doc/rnd/exact-perturbation-theory/STATUS.md): M0 spectrum/energy benchmarks and pending projector, parity, flux and resolvent tasks; explicitly states that no perturbative coefficient has been implemented.
- [01_flux_sectors_and_spectrum.md](https://github.com/lufreitas0000/kitaev-mf/blob/e093e417494de7b7cb78ae415338189d628a8194/doc/rnd/exact-perturbation-theory/notes/01_flux_sectors_and_spectrum.md): fixed-gauge singular-value normalization, flux-free spectrum and matter parity caveats.
- [02_normalization_and_bloch_conventions.md](https://github.com/lufreitas0000/kitaev-mf/blob/e093e417494de7b7cb78ae415338189d628a8194/doc/rnd/exact-perturbation-theory/notes/02_normalization_and_bloch_conventions.md): coupling and Bloch generator normalization, including the factor-of-four matrix-level distinction.
- [M1.1 physical projector](https://github.com/lufreitas0000/kitaev-mf/blob/e093e417494de7b7cb78ae415338189d628a8194/doc/rnd/phase_1/m1_1_t3_physical_projection_pauli_algebra.md): explicit physical-sector convention and Pauli algebra.

This inventory prevents duplicate tasks but does **not** promote the incomplete exact spin perturbation track to completed analytical status.

## 4. Handoff and scope boundary

**Proved here:** the finite-dimensional implication is false; the counterexample and conditional quadratic-form identity follow directly from the stated algebra.

**Not proved here:** zero/nonzero character of each actual Kitaev second-order operator; coefficient of the third-order chiral mass; full physical gauge projection on a torus; finite-field KKT dressing; identity of a mean-field fitted gap with a physical-spin spectral threshold.

**Proposed next independent audit, only if needed:** reconcile the existing local parity conventions, oriented two-link three-spin identity, and the choice of $P/Q$ and target spectral gap, using the files above as the starting corpus. No repetitive implementation or new symbolic/numerical tests are required merely to refute the false implication.

**STOP:** Mathematical erratum submitted for review; do not rewrite the legacy K-019 source, change statuses, implement tests or infer a completed continuum coefficient without a separate reviewed task.
