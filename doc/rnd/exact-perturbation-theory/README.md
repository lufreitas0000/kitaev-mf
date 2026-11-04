# Exact perturbation theory of the pure Kitaev model

This is an exploratory research track for deriving and validating perturbation
theory about the isotropic, vortex-free Kitaev spin liquid. It is isolated from
the existing Majorana mean-field implementation.

## Scope

- Establish conventions for the pure honeycomb Kitaev Hamiltonian.
- Compute matter spectra and flux-sector energies from quadratic Majorana problems.
- Enforce the physical gauge constraint and matter-parity selection rule.
- Define the low-energy external-state window and flux-sector resolvent domain.
- Use small-cluster spin-space exact diagonalization as an independent check.

Analytical derivations are primary. Numerical calculations provide validation
and finite-size diagnostics.

## Planning

- `ROADMAP.md`: phased research plan and completion criteria.
- `BACKLOG.md`: prioritized implementation and derivation tasks.

## Structure

- `STATUS.md`: scope, milestones, and unresolved points.
- `notes/`: mathematical and physical research notes.
- `code/README.md`: implementation boundary and numerical contracts.
- `results/README.md`: generated-data policy.

No production Python modules are introduced at this stage. Their design follows
after the mathematical contracts and tests are agreed.

## Reference

A. Kitaev, *Anyons in an exactly solved model and beyond*,
Annals of Physics **321**, 2–111 (2006), arXiv:cond-mat/0506438.

Numerical vortex-energy values quoted in these notes are finite-size estimates
reported in Appendix A, not exact thermodynamic constants.
