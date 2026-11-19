# Research idea backlog — kitaev-mf

**Status:** non-committal idea register, 2026-11-14 source reconciliation. This is **not** the ordered milestone/product backlog: that authority remains in [doc/MILESTONE_ROADMAP.md](doc/MILESTONE_ROADMAP.md). Active tasks belong to [PROJECT_STATE.yaml](PROJECT_STATE.yaml) and the [sprint brief](CURRENT_SPRINT.md); discoveries and negative results belong to [LOG_DISCOVERY_AND_PIVOTS.md](LOG_DISCOVERY_AND_PIVOTS.md). A GitHub Issue is needed only for an independently reviewable deliverable, discrepancy or cross-project handoff.

The items below capture **uncommitted scientific methods and alternatives** from the 2026-11-14 user-supplied `maximalist_future_tasks_backlog.md` (§2.2–2.5) and related open research. The original document is an historical source, not a verified statement that each suggestion is feasible, correct or approved. Prior milestone or registered-HC status is never superseded by an IDEA label. IDs remain stable when promoted, rejected or superseded.

## IDEA-001 — Analytic weak-field scMFT mass coefficient

Status: investigating  
Domain: mathematical physics / Gaussian mean field  
Source: user backlog §2.2 (K-010/K-019 P1–P6); [Issue #6](https://github.com/lufreitas0000/kitaev-mf/issues/6); [noncanonical Wave 2 pilot](doc/rnd/phase_1/m1_10_t1_wave2_cubic_gap_evidence_pilot.md).  
Question: derive the order-h³ coefficient in the fully self-consistent 8×8 matter/gauge and KKT equations, including field-dependent Fock vertices, rather than equating the fit `c_MF=2.5921±0.0034` with a theorem.  
Evidence now: numerical finite-grid power-law fit; a proposed analytic derivation is *not established* by that fit.  
First feasibility test: derive a transparent weak-field Schur/resolvent expansion in one frozen convention, explicitly bounding infrared and finite-mesh terms.  
Prerequisites: K-019-A projector/operator consistency; review the order of limits and source normalization.

## IDEA-002 — Exact spin perturbation coefficient with continuum resolvent

Status: investigating  
Domain: exact Kitaev perturbation  
Source: user backlog §2.2 and §2.3; [Issue #12](https://github.com/lufreitas0000/kitaev-mf/issues/12); [K-019-A #10](https://github.com/lufreitas0000/kitaev-mf/issues/10).  
Question: compute the *physical* third-order spin effective operator using its two-vison and itinerant-Majorana resolvent, not an uncontrolled scalar `1/Delta_flux` substitution.  
Evidence now: formal finite-volume Feshbach identities under their resolvent-domain assumptions; thermodynamic uniformity, coefficient and projected physical spectrum remain open.  
First feasibility test: independently verify finite-volume parity and an explicit matrix-element oracle with the exact `P,Q` spaces.  
Owner choices: physical sector, boundary/twist and intended energy window.

## IDEA-003 — Symmetry-fixed covariance geometry and bounded multipliers

Status: investigating  
Domain: representation theory / nonlinear KKT  
Source: user backlog §2.2 L2; [Issue #8](https://github.com/lufreitas0000/kitaev-mf/issues/8).  
Question: distinguish 86 stored array slots from the admissible covariance manifold, true C3/sublattice fixed subspace and multipliers. Is any proposed 5-parameter slice invariant under the complete fixed-point map?  
Evidence now: 86 storage slots and 66 independent real slots after imposing only onsite skew-symmetry; neither is a physical dimension or proof of five-dimensional closure.  
First feasibility test: derive the action on directed U, onsite V and KKT variables symbolically before exploratory T1–T7 numerics.  
Owner choices: meaning of restricted solve (symmetry-closed, constrained, or truncated), parameter norm/metric and energy comparison.

## IDEA-004 — Certified stationary-state enumeration

Status: deferred  
Domain: validated numerics / nonlinear geometry  
Source: user backlog §2.2 L4.  
Question: when does interval Newton or Krawczyk with **validated Brillouin quadrature** certify all roots of a compact admissible reduced family?  
Evidence now: exploratory proposal only. A degree/index sum alone cannot certify completeness; singular KKT charts and unbounded multipliers require treatment.  
First feasibility test: independently prove compactness and interval enclosure on a genuinely symmetry-closed finite-dimensional model.  
Trigger: successful IDEA-003 and approved domain/precision contract.

## IDEA-005 — Projected non-Gaussian diagnostics and relative entropy

Status: deferred  
Domain: quantum information / predictive error  
Source: user backlog §2.4, BL-07 and H1/H4; [M2.1 roadmap](doc/MILESTONE_ROADMAP.md).  
Question: can a *defined* relative entropy or physical projection distance be related to controlled observable errors against finite ED or another independent oracle?  
Evidence now: proposed diagnostic; relative entropy alone is not automatically an a posteriori bound on arbitrary observables or on thermodynamic predictions.  
First feasibility test: specify density matrices on the *same physical spin Hilbert space*, their supports and a correct inequality applicable to bounded observables.  
Trigger: independently accepted physical projector and matched finite-cluster comparison.

## IDEA-006 — Infrared-aware Jacobian preconditioning at Dirac nodes

Status: deferred  
Domain: numerical analysis  
Source: user backlog §2.4 BL-05.  
Question: can a regularized spectral preconditioner improve JFNK stability near gapless points without changing the actual fixed-point solution?  
Evidence now: proposed algorithm, not a proven speedup or convergence guarantee.  
First feasibility test: define regularized operator, invariants and cap-independent work accounting; compare to an unpreconditioned exact same starting state.  
Trigger: solver evaluation-cost protocol approved under [Issue #7](https://github.com/lufreitas0000/kitaev-mf/issues/7).

## IDEA-007 — Explicit SO(2) / bond-permutation coordinate intertwiners

Status: deferred  
Domain: symmetry / model implementation  
Source: user backlog §2.4 BL-04.  
Question: reconcile `multisublattice_saddle_survey.md` and `domain/lattice.py` with explicit coordinate shears, the S3 bond action, and transformed physical operators.  
Evidence now: possible representation/coordinate mismatch, not a verified symmetry defect.  
First feasibility test: derive an operator-level intertwiner on a finite cell and independently test matrix equality.

## IDEA-008 — Incommensurate spiral pitch optimization

Status: deferred  
Domain: extended J–K–Gamma model phases  
Source: user backlog §2.4 BL-03; roadmap Phase 5.  
Question: how should a continuous pitch vector `Q` in the BZ enter an admissible and properly normalized multicomponent ansatz without presupposing a periodic supercell?  
Evidence now: exploratory method, outside the active pure ferromagnetic Kitaev [111] model.  
First feasibility test: define the reference Hamiltonian and an exactly tractable limit before any numerical optimization.  
Trigger: explicit model-family priority selection.

## IDEA-009 — Per-flavor sign symmetry versus KKT constraints

Status: investigating  
Domain: gauge and discrete symmetry  
Source: user backlog §2.2 Appendix A A.5; [Issue #8](https://github.com/lufreitas0000/kitaev-mf/issues/8).  
Question: do proposed independent Majorana-flavor sign operations preserve the full mean-field equations, field coupling and gauge constraints at h=0, and what breaks at h≠0?  
Evidence now: an open test of a proposed symmetry, not a proven quotient.  
First feasibility test: compute transformations of every exchange, onsite, field and multiplier term before reducing coordinates.

## IDEA-010 — Competing finite-field vortex constraints

Status: investigating  
Domain: vortex / approximate flux formulations  
Source: user backlog §2.2 L3; [Issue #15](https://github.com/lufreitas0000/kitaev-mf/issues/15).  
Question: compare fixed-background link sectors, constrained Gaussian flux expectations and initially seeded, then relaxed, defects **without** assuming conserved physical plaquette flux at nonzero field.  
Evidence now: unapproved alternative definitions; the old fixed-background number `0.167911|K|` is not a parity-projected relaxed physical vortex ground-state energy.  
First feasibility test: prove feasible constraint families and derive a consistent energy functional, stationarity and second variation.  
Owner gate: R3 physics choices, thresholds and comparison sector.

## IDEA-011 — Spin-channel choice under physical projection

Status: deferred  
Domain: gauge projection / Gaussian spin observables  
Source: user backlog §2.4 H4; roadmap M2.1.  
Question: how should normal `N^a` versus chiral `M^a` spin channels be compared after physical projection and under the exact same physical variational functional?  
Evidence now: interpretation ambiguity; no rigorous physical-spin energy upper bound from unprojected scMFT follows automatically.  
First feasibility test: explicit two-site projected matrix elements and independent small spin-ED oracle.  
Trigger: [Issues #10](https://github.com/lufreitas0000/kitaev-mf/issues/10), [#12](https://github.com/lufreitas0000/kitaev-mf/issues/12), [#14](https://github.com/lufreitas0000/kitaev-mf/issues/14).

## IDEA-012 — Alternative exact–MF dynamical correlator approximations

Status: deferred  
Domain: dynamical spin correlations  
Source: user backlog §2.5 Phase 3; roadmap M3.0–M3.5.  
Question: which approximation to the gauge-quench response (Gaussian bubble, vertex-corrected TD-HFB, RPA, or an exact finite-volume Lehmann oracle) preserves each intended physical sum rule?  
Evidence now: method choices, **not** completed M3 physics.  
First feasibility test: precise finite-volume operator and frequency-Fourier normalization plus a spin-ED sum-rule check.  
Trigger: the owner selects Phase 3 and its benchmark scope.

## IDEA-013 — Multi-parameter continuation holonomy and real-space fixed-point relaxation visualization

Status: deferred  
Domain: continuation geometry / real-space fixed-point diagnostics  
Source: legacy harvest chats [`47bd81c0-b1ab-4562-be97-4944fb2ab1a4`](conversation://47bd81c0-b1ab-4562-be97-4944fb2ab1a4) (Prompt 6) and [`13459723-61fc-426c-b7ea-41c031eccff3`](conversation://13459723-61fc-426c-b7ea-41c031eccff3) (Prompt 3).  
Question: (1) Does closed-loop multi-parameter continuation around spinodal cusps or a finite-$T$ critical end point $(h_{\mathrm{CEP}}, T_{\mathrm{CEP}})$ (e.g., tilted field angles $(\theta, \phi)$ or $(J, K, \Gamma, h)$) exhibit branch-switching holonomy beyond 1D $\mathbf{h} \parallel [111]$ forward/backward hysteresis? (2) Can real-space fixed-point map iteration visualizer movies on $L=32$ clusters (encoding diagonal bond amplitudes $U_b$ in edge width/color, plaquette fluxes $W_p \in [-1,+1]$ in hexagon shading, and comparing $L \bmod 2$ vs $L \bmod 3$) diagnose domain-wall vs metastable convergence basins?  
Evidence now: exploratory proposal only; no holomorphic/cusp continuation theorem or certified M2.3 dataset exists yet.  
First feasibility test: trace one small closed loop in $(h, \theta)$ on $L=6$ / $L=12$ using existing fixed-point/Newton solvers and prototype a static $L=32$ relaxation frame series.  
Trigger: M2.3 spinodal/CEP extension or Phase 5 / `NP-007` visualizer milestone.

## Non-IDEA items retained in their authoritative records

- M2.1–M2.4, M3.0–M3.5 and future model phases are already **roadmap milestones**; do not duplicate them as new Issues or active tasks here.
- K-010/011/019/020 are **historical recorded claims/evidence** requiring source-aware correction if disputed, not permission to rewrite their original status.
- L1 paired benchmark repair is tracked by existing Issue #7 and the [merged, unregistered W1 protocol PR #16](https://github.com/lufreitas0000/kitaev-mf/pull/16) (reviewed for documentary scope only); L3 fixed vortex constraint by Issue #15. Retain stable source links rather than duplicating their executable task list.
- Multi-method ED, iPEPS, SDP, HTE and tensor-network implementation belongs primarily in the independently owned `kitaev-numerics` roadmap/backlog unless a concrete cross-project interface is requested.
- Governance implementation waves, website and Lean/monograph work are owned by their respective repositories. The source-to-destination crosswalk for this migration is [doc/PLANNING_RECORDS_MIGRATION.md](doc/PLANNING_RECORDS_MIGRATION.md).

## Triage policy

Ideas are reviewed when scientific prerequisites or user priorities change, not automatically every sprint. `promoted` requires a link to an approved owning milestone or bounded Issue; `withdrawn` retains the reason and old evidence. No IDEA may be silently represented as a confirmed analytical theorem, numeric convergence guarantee or user-approved acceptance threshold.
