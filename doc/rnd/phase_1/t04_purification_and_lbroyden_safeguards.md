---
id: T04
title: "Purification Warm-Start and Limited-Memory Broyden Safeguards Verification"
status: Completed
date: "2026-11-07"
questions:
  - "How can purify_newton_schulz warm-start be structured so it does not freeze at Y_0 = Omega_{k-1} or discard Hamiltonian updates?"
  - "Does setting default alpha=0.5 and enforcing tol_singular scaling in LimitedMemoryBroydenSolver stabilize the orthogonal complement across a 30-trial perturbed benchmark (Contract H9)?"
depends_on:
  - T03
  - M1.9-T1
  - M1.9-T3
contracts:
  - K-013
  - K-018
  - K-021
  - Contract H9
artifacts:
  - src/kitaev_mf/services/purification.py
  - src/kitaev_mf/services/fixed_point.py
  - test/unit/test_profiler_and_lbroyden.py
  - test/unit/test_real_space_and_purification.py
---

# T04: Purification Warm-Start and Limited-Memory Broyden Safeguards Verification

## 1. Scope & Adversarial Audit Findings (`[B4]`, `[B6]`)

During the Wave 3/4 adversarial audit, two core numerical defects were identified in [`src/kitaev_mf/services/purification.py`](file:///home/lucas/Projects/kitaev-mf/src/kitaev_mf/services/purification.py) and [`src/kitaev_mf/services/fixed_point.py`](file:///home/lucas/Projects/kitaev-mf/src/kitaev_mf/services/fixed_point.py):

1. **Defect `[B4]` (Broken `warm_start` in `purify_newton_schulz`):**
   - In `purify_newton_schulz`, passing `warm_start` directly set `Y = 0.5 * (warm_start - warm_start.T)` and immediately executed the Newton-Schulz cubic map $Y_{k+1} = \frac{1}{2}Y_k(3I + Y_k^2)$.
   - Because any converged covariance matrix satisfies $Y^2 = -I$, this yielded $Y_1 = \frac{1}{2}Y(3I - I) = Y$, instantly terminating the iteration without ever touching the new Hamiltonian matrix $A$!
   - To make `warm_start` mathematically valid, the initial iterate must incorporate the updated Hamiltonian $A$. As derived in `doc/rnd/phase_1/m1_9_t1_self_consistent_purification_warm_starts.md`, the first-order commutator gradient correction is:
     $$Y_0 = \Omega_0 - \tau \left( A + \Omega_0 A \Omega_0 \right),$$
     where $A + \Omega_0 A \Omega_0 = [ [A, \Omega_0], \Omega_0 ]$ is the off-diagonal particle-hole generator that breaks idempotency in the direction of the updated Hamiltonian gradient. With $\tau \le \frac{1}{2 \|A\|_2}$, $Y_0$ stays within the cubic convergence basin while reducing Newton-Schulz iterations from 11--12 to 4--6 steps. If $A$ is unchanged, it preserves $\Omega_0$; if $A$ has changed, it drives convergence toward the new ground-state projector.
2. **Defect `[B6]` (`LimitedMemoryBroydenSolver` Safeguards & `Contract H9`):**
   - In `LimitedMemoryBroydenSolver.__init__`, `alpha` defaulted to `1.0` instead of `0.5`, contradicting the damping protocol established in `K-021` and `doc/rnd/phase_1/m1_9_t3_lbroyden_damping_and_multisecant_reconciliation.md`.
   - In `step()`, the denominator check `abs(denom) > 1e-14` was hardcoded instead of checking `abs(denom) > self.tol_singular * max(norm(s_j) * norm(hy_j), 1e-14)`, allowing near-singular secant updates to corrupt the inverse Jacobian recursion.
   - We set default `alpha = 0.5`, enforce relative scale-invariant denominator gating, and verify stability across a 30-trial perturbed linear and mean-field benchmark (`Contract H9`).

---

## 2. Verification Plan
- Update `purify_newton_schulz` in `src/kitaev_mf/services/purification.py`.
- Update `LimitedMemoryBroydenSolver` in `src/kitaev_mf/services/fixed_point.py`.
- Add comprehensive unit tests in `test/unit/test_profiler_and_lbroyden.py` and `test/unit/test_real_space_and_purification.py`.
- Run `make format && make verify`.
- Deploy adversarial referee subagent `c5017fbe` for gate verification.
