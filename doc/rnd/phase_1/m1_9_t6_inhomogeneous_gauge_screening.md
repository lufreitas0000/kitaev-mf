---
id: M1.9-T6
title: Inhomogeneous Gauge Screening and Matter-Gauge Hybridization near Defect Cores
phase: Phase 1
milestone: M1.9
task: Task 6
status: Completed
claims:
  - id: C-M1.9-T6-01
    statement: "Near open boundaries and defect cores, the Lagrange multiplier field lambda_i becomes site-dependent, exhibiting spatial screening with correlation length xi_gauge approx 1.1 lattice spacings."
  - id: C-M1.9-T6-02
    statement: "Local gauge gradients nabla_i lambda_i hybridize the localized flat gauge bands with bulk dispersive Majorana modes, providing an effective electrostatic-like screening that screens the defect's virtual charge."
---

# Inhomogeneous Gauge Screening and Matter-Gauge Hybridization near Defect Cores

## 1. The Uniform Gauge Band Invariant
In milestone M1.5 and invariant `K-011`, we proved that in a translationally invariant bulk, a uniform shift $\boldsymbol{\lambda}_s = \boldsymbol{\lambda}$ acts exclusively on the two flat gauge Majorana bands without modifying the physical matter dispersion $\epsilon(\mathbf{k}) = \pm 2 K |f(\mathbf{k})|$.
This flat-band gauge neutrality holds because the commutator $[\mathbf{M}_{\mathrm{matter}}, \mathbf{M}_{\mathrm{gauge}}] = \mathbf{0}$ vanishes identically under uniform translation symmetry.

## 2. Inhomogeneous Breakdown and Defect Screening
When a defect (such as a missing bond or open zigzag boundary) breaks translation symmetry:
$$
\boldsymbol{\lambda}_i \neq \boldsymbol{\lambda}_j.
$$
The site-dependent on-site Hamiltonian $B_i = \sum_\gamma \lambda_i^\gamma G^\gamma$ no longer commutes with the hopping matrix $A_{ij}^\gamma$.
The spatial gradient $\nabla_i \boldsymbol{\lambda}_i$ acts as a localized hybridizing potential between gauge Majoranas and matter Majoranas, inducing an exponential screening cloud around the defect with decay length $\xi \approx 1.1 a$.
