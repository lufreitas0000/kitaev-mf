# W001-T2: Earlier perturbation work and preservation chronology

Date: 2026-11-07. Source commits and full Git blobs are in
`w001_source_inventory.json`; this report does not certify scientific correctness.

## Corrected question

Lu Freitas clarified that the earlier perturbation work was **ignored**, not
believed altered. The earlier T1 wording about possible damaged findings is
superseded. The task is to make the original contracts operative in subsequent
review, not restore allegedly lost mathematics.

## Content comparison

Baseline: earlier research ref at `13db8f2cc50c43cbe893e251bd6e09ea29456820`.
Compared against Phase 2 `8b242a5`, Phase 3 `21a27d8`, and Phase 4 `1d82e6a`.

| Original note | Baseline blob | Phase 2 | Phase 3 | Phase 4 |
|---|---|---|---|---|
| `00_model_and_conventions.md` | `2255bd043d52` | identical | identical | identical |
| `01_flux_sectors_and_spectrum.md` | `d8363aa68945` | identical | identical | identical |
| `02_normalization_and_bloch_conventions.md` | `c9338ad54904` | identical | identical | identical |
| `02_vortex_energies.md` | `bee469983f88` | identical | identical | identical |
| `03_low_energy_projector_and_feshbach.md` | `f97cf42527ec` | identical | identical | identical |
| `04_numerical_validation.md` | `5bde0af6152e` | identical | identical | identical |
| `05_finite_torus_conventions.md` | `f376b18f1e46` | identical | identical | identical |

All original notes listed above remain byte-identical at those compared refs.
This bounds the observation to these paths and commits; it does not assess every
historical artifact or all mathematical claims in the repository.

## Chronology and neglected constraints

1. `a8847ea`: initial exact-Kitaev framework. Subsequent original derivations
   include `1c4e2a8` (dispersion), `deb52c8` (vacuum energy), and `22b4b24`
   (one-particle versus many-body widths).
2. `f79ae1a` and later November 4 notes/backlog changes specify projection,
   physical thresholds, ordering, and resolvent-domain prerequisites. The earlier
   STATUS.md records open questions rather than completed perturbative coefficients.
3. `4127d0e` adds T12-T15 and exact_matter/feshbach/exact_perturbation services,
   and changes many BACKLOG checkboxes to completed. These additions coexist with
   the earlier notes; their presence does not mean the earlier caveats were met.
4. `7c75e70` integrates that research history into the later source lineage.
5. The historical Phase 4 study then presents an MF/ED prefactor explanation
   disapproved by the PI. Its preservation does not reinstate acceptance.

The earlier `notes/03_low_energy_projector_and_feshbach.md` explicitly distinguishes
all physical flux-free states from an external matter-energy window, requires
matched physical Hilbert spaces and Wilson-sector conventions, and says that an
unperturbed flux gap alone does not control the perturbed eliminated block.
It asks for an actual sector/local/linked-cluster argument when a global norm
estimate is insufficient. These are prerequisites, not an explanation already
provided by invoking locality in prose.

Later T14 asserts local-ratio control without establishing the required uniform
linked-cluster bound. T15 explicitly uses a static two-vison resolvent
approximation while calling its coefficient exact, and claims an error estimate
with an unspecified constant. The scalar-gap routine implements that assumed
formula, not a general operator-valued, all-orders expansion. These differences
require W001-T4 review; no substitute coefficient is accepted in this task.
The original notes remain the source of the stricter questions and contracts.

## Preservation

Five original untracked Phase 3 files were captured in a separate snapshot commit
whose parent is the recorded Phase 3 source tip. It changes only those five paths.
Their source worktree, original index, and branch HEAD were left untouched.
`archive/w001/phase3-untracked-inputs` points to that commit. Source-tip archive
tags preserve policy baseline, earlier perturbation ref, and the three research
inputs without opening additional active research branches. Full refs and hashes
are recorded in the inventory. No archive tag is a release or validity certificate.

## T3 handoff

- Integrate the earlier contracts alongside later results and use the acceptance
  register to qualify unsupported completion labels. Do not overwrite original notes.
- Distinguish retained reference identity from scientific acceptance.
- Reconcile local-main/remote-main drift before main-targeted PRs.
- Create exploration and update CI branch/PR filters; current filters do not yet
  cover research/, promote/, or exploration consistently.
- Keep disapproved Phase 4 sources discoverable as archival/exploratory material,
  never as authoritative input to an accepted prefactor explanation.

## Verification scope

Git blob identity checks establish preservation, not physical validity. Snapshot
file hashes and changed-path checks are required before push. No new physics
calculations, numerical tests, or document publication were run for this audit.
