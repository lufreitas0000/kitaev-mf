# Foundational gauge audit and common QSL notation — research decision record

**AI contribution:** ChatGPT (OpenAI GPT-6.2), 2026-11-25, analytical derivations and static GitHub source inspection, following discussion with the project owner. Author attribution identifies the drafting model, not independent scientific approval. **Status:** proposed, independently unverified. No numerical tests executed and no claim of literature novelty. The owner selected the bilayer notation as the common target; implementation and scientific acceptance remain gated.

**Sources:** [kitaev-mf baseline 2971ad21ce72](https://github.com/lufreitas0000/kitaev-mf/tree/2971ad21ce7207bcd1be5a4e44456be9a9d28886), [bilayer conventions-v2.0 at c0554a9](https://github.com/lufreitas0000/ksl-sc-bilayer/blob/c0554a9a18359a44a7bcea37a35c9078f7408917/CONVENTIONS.md), [bilayer foundational derivation at ae5e68e4](https://github.com/lufreitas0000/ksl-sc-bilayer/blob/ae5e68e4bbb7386315d8e7d353e6bb19d7c9f9df/doc/rnd/mf_analytic/gauge_constraint_physical_projection_gaussian_obstruction.md). Authoritative work packets: [notation #25](https://github.com/lufreitas0000/kitaev-mf/issues/25), [methodology #27](https://github.com/lufreitas0000/kitaev-mf/issues/27), [bilayer #17](https://github.com/lufreitas0000/ksl-sc-bilayer/issues/17).

## Exact dictionary (proposed regression contract)

Write $\gamma=(c,b^x,b^y,b^z)^T$, $c^{\rm old}=P_0\gamma$, $P_0=\operatorname{diag}(-1,1,1,1)$. Comparing explicit matrices in the pinned [constants.py L19–110](https://github.com/lufreitas0000/kitaev-mf/blob/2971ad21ce7207bcd1be5a4e44456be9a9d28886/src/kitaev_mf/domain/constants.py#L19-L110) with the bilayer convention gives
$$
P_0M_{\rm old}^\alpha P_0/2=L^\alpha,\quad
P_0G_{\rm old}^\alpha P_0/2=-R^\alpha,\quad
P_0N_{\rm old}^\alpha P_0=L^\alpha-R^\alpha.
$$
In particular, $P_0 L P_0=R$, $P_0 R P_0=L$. The old spin $S=i(c^{\rm old})^TM_{\rm old}c^{\rm old}/8$ becomes $i\gamma^TL\gamma/4$. The old gauge $T=i(c^{\rm old})^TG_{\rm old}c^{\rm old}/8$ becomes $-G$ with $G=i\gamma^TR\gamma/4$. With the old multiplier $2\lambda_{\rm old}\cdot T$, the new multiplier is $\lambda=-2\lambda_{\rm old}$. This is operator normalization, not a claim about raw legacy Bloch matrix units.

$$
D_{\rm old}=c^0c^1c^2c^3=-cb^xb^yb^z=b^xb^yb^zc=D.
$$

Both physical sectors are $D=+1$. The old Pfaffian expression is $\langle D\rangle=-\operatorname{Pf}(B_{\rm old})$, the new is $+\operatorname{Pf}(B)$, since $B=P_0 B_{\rm old}P_0$ and $\det P_0=-1$. The signs cancel exactly. For many sites transform all covariance and Hamiltonian blocks, and track the orientation of the full basis and the parity of the Thouless reference. In particular, a sign-reversing coordinate transform need not be an even physical gauge rotation. Retain all sector and Fock-basis conventions explicitly.

The exact physical Pauli extension used by old $N$ maps to $2(S-G)$, not $2S$ as an enlarged-space identity. For pure notation migration preserve $2(S-G)$ and all model units; the choice of extension is a separate scientific task.

## Methodological results and required tests

For Hermitian constraints, $\langle\sum_\alpha(G^\alpha)^2\rangle=\sum_\alpha\|G^\alpha\psi\|^2$ for pure states and $\sum_\alpha\|G^\alpha\sqrt\rho\|_{\rm HS}^2$ for mixed states. Zero is equivalent to physical support. The old local Casimir is already invariant K-002. With $\mathbf G^2=3(1-D)/8$, Gaussian Wick factorization gives $\langle D\rangle=\operatorname{Pf}(B)=4(|\langle\mathbf S\rangle|^2-|\langle\mathbf G\rangle|^2)$ in the new convention. Thus exact local physicality of a Gaussian spin requires maximal polarization and a pure local spin marginal. A locally unpolarized gauge-neutral Gaussian has leakage $1/2$. This is a structural property of Gaussian parents, not a universal relative error of observables.

A proposed explicit counterexample to a physical upper bound is a ferromagnetic Heisenberg dimer with $H=-J\sum_\alpha\sigma_1^\alpha\sigma_2^\alpha$, $J>0$. Its physical minimum is $-J$. A Gaussian covariance with zero onsite blocks and cross block $Q=\operatorname{diag}(-1,1,1,1)$ is pure, gauge-neutral and gives each unprojected legacy Pauli bond correlation $-Q_{00}Q_{\alpha\alpha}=1$. Its extension energy is $-3J$ and joint physical weight is zero (opposite local parities). Two copies permit even global Dirac parity. Verify by independent Fock matrices. Physical Rayleigh–Ritz applies to nonzero normalized projected states, not generally unprojected candidates.

The earlier [quasifree note L118–173](https://github.com/lufreitas0000/kitaev-mf/blob/2971ad21ce7207bcd1be5a4e44456be9a9d28886/doc/rnd/phase_1/m1_5_t3_quasifree_states_wick_variational.md#L118-L173) requires scope corrections: the sign of parity must agree with the actual basis; commutator stationarity is weaker than ground-state spectral selection; an affine tangent is not automatically a lower bound for a nonconvex Wick functional; a physical upper bound is not guaranteed. Preserve exact local operator algebra and numerical solver benchmarks with their stated auxiliary scope.

Proposed experiment: compare $X_i^\alpha(\eta)=2(S_i^\alpha-\eta G_i^\alpha)$ for $\eta=-1,0,1$. Their projected spin Hamiltonians agree; Gaussian parent states may differ. Use existing [projection engine](https://github.com/lufreitas0000/kitaev-mf/blob/2971ad21ce7207bcd1be5a4e44456be9a9d28886/src/kitaev_mf/services/gutzwiller_projection.py) and [HC-02 runner](https://github.com/lufreitas0000/kitaev-mf/blob/2971ad21ce7207bcd1be5a4e44456be9a9d28886/scripts/wave1/run_hc02_variational_sandwich.py#L55-L173) as baseline. Retain all candidates before projection; match finite geometry, parent optimization, parity, twists and physical Hamiltonian. Compare KKT on/off, PAV, VAP and independently sourced ED/SDP. Separate absolute physical projection weight from relative Thouless amplitude norm. Distinguish exact energy from VMC estimates and chart singularity from genuinely vanishing projection. Reassess HC-01 auxiliary gap, HC-02 projected phase selection, HC-03 visons, and thermal entropy before physical claims. No exact microscopic Hamiltonian change inferred.

## Proposed execution dependency

1. N0/F0: independent algebra, signs, exact Fock tests and Gaussian obstruction.
2. N1–N4: notation migration as a behavior-preserving coordinate change, regression tests and independent review; this may run alongside independent oracle development.
3. F1/F2: extensions and projection engine verification.
4. F3/F4: matched numerical comparisons and projected variation.
5. F5–F7: revise scientific interpretation, response, defects and finite-temperature methodology only after reviewed evidence.

Issue #25 owns the notation checklist; #27 owns the methodology checklist. The bilayer #17 owns its spin–electron-specific tests. Do not copy ED/SDP production solvers into this repository (ADR-003); import verified data. No threshold changes or automatic promotion of hypotheses. No new simulation evidence in this record.
