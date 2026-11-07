---
id: M1.9-T6
title: "Inhomogeneous Gauge Screening and Matter-Gauge Hybridization near Defect Cores"
status: Completed
date: "2026-11-06"
questions:
  - "How do spatial gradients in the Lagrange multiplier field nabla_i lambda_i affect localized flat gauge Majorana bands near defect cores?"
  - "What is the characteristic screening length xi_gauge over which local gauge perturbations decay back to bulk flat-band neutrality?"
depends_on:
  - M1.5-T1
  - M1.9-T2
contracts:
  - K-008
  - K-022
artifacts:
  - doc/labbook/chapters/ch_m1_9_real_space_and_purification.tex
  - src/kitaev_mf/services/real_space.py
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
