# Note 00f: M1.3 Project-State Synchronization and K-006 Provenance

## 1. Purpose

This note records the M1.3 synchronization decision after completion of the M1.1 foundations batch. It is complementary to Notes 00a--00e; those existing R&D notes are not modified by this batch.

M1.1 has been locally validated by the user on WSL after the latest development-branch updates. The validation applies to the four M1.1 test modules and the associated Ruff and Mypy checks that were requested. This is a user-reported result; no claim is made that the assistant executed those checks.

## 2. K-006 statement and assumptions

K-006 concerns the zero-temperature, zero-flux, isotropic ferromagnetic pure-Kitaev reference state on the two-sublattice primitive honeycomb cell. The assumptions are

- $K<0$;
- $J=\Gamma=\Gamma'=D=\mathbf h=\boldsymbol\lambda=0$;
- directed gauge-link correlators $U_b^{\gamma\gamma}=+1$;
- translationally invariant matter bond correlator $U_b^{00}=u_{\mathrm{matter}}$;
- the negative-energy matter band is occupied, with particle-hole-symmetric half occupation assigned to exact zero modes.

The analytical result in Note 03, Section 3 is

$$
U_b^{00}=-\bar U,
\qquad
\bar U=\int_{[0,1)^2}\!dk_1\,dk_2\,
\frac{1+\cos(2\pi k_1)+\cos(2\pi k_2)}
{\sqrt{3+2\cos(2\pi k_1)+2\cos(2\pi k_2)+2\cos(2\pi(k_1-k_2))}},
$$

with the reported thermodynamic-limit quadrature value $\bar U\approx0.5248645$. The same derivation gives $V_s=0$, vanishing local gauge residuals, and $U_b^{\gamma\gamma}=+1$ for all three bond directions in the stated convention.

The integral expression and its sign follow from the negative-energy eigenprojector of the two-band matter block and the Brillouin-zone correlator formula. Equality of the three bond values follows from the $C_3$-related relabeling of the Brillouin-zone torus. The decimal value is a numerical quadrature result, not an exact closed-form constant.

## 3. Implementation and evidence map

| Layer | Repository evidence | What it checks |
|---|---|---|
| Analytical dispersion | `doc/rnd/note_02_momentum_bloch_formalism.md`, Section 4 | Pure-Kitaev matter and flat-gauge spectra |
| Spectral projector | `src/kitaev_mf/services/correlators.py::compute_spectral_projector` | Negative-energy projector with half-weight zero modes |
| BZ correlators | `src/kitaev_mf/services/correlators.py::compute_correlators_from_bloch` | Fourier-phase integration into real-space $U,V$ |
| K-006 derivation | `doc/rnd/note_03_correlators_saddle_point.md`, Section 3 | Analytical integral, sign, and gauge-sector statements |
| Regression test | `test/unit/test_correlators.py::test_pure_kitaev_correlator_invariant_k006` | $36\times36$ mesh check of $U_b^{00}$, $U_b^{\gamma\gamma}$, $V$, and gauge residuals |

The regression test uses absolute tolerance $3\times10^{-3}$ for $U_b^{00}$ relative to the quoted quadrature value, $10^{-12}$ for the gauge-bond correlator, and $10^{-14}$ for $V$ and the gauge residuals. The matter-correlator tolerance is deliberately a finite-mesh regression tolerance; it does not numerically certify the seven displayed decimal places.

The quoted $0.5248645$ value is present in Note 03, but a standalone high-precision quadrature script and its convergence study are not currently identified among the committed artifacts. This is a provenance limitation, not a failure of the analytical integral formula. If the decimal benchmark is to become a strict numerical reference, a later dedicated quadrature artifact should document mesh refinement or an adaptive integration error estimate.

## 4. Status classification

K-006 is recorded as:

- **Analytical formula:** derived in Note 03 under the assumptions above.
- **Numerical decimal:** reported high-precision quadrature value in Note 03; independent reproducibility artifact not yet committed.
- **Regression test:** implemented in `test_correlators.py`; its finite-mesh tolerance checks the expected value at coarse-grid accuracy.
- **Local validation:** pre-existing test status is retained; this M1.3 synchronization does not claim a new WSL run of the full repository suite.

This classification avoids conflating an analytical integral, a quoted quadrature value, and a finite-mesh regression test.

## 5. M1.3 completion and next milestone

M1.3 synchronizes the machine-readable project state with the completed M1.1 and M1.2 work, records K-006 provenance and evidence quality, and identifies M1.4.1 as the next batch.

M1.4.1 is a derivation/specification batch for observables and energetics. It should define the physical observables, normalization conventions, energy functional, and analytic limiting cases before implementation changes are made. It must not silently promote a finite-mesh numerical benchmark to an exact identity.
