# Rule 09: Research Intake, Classification, Acceptance, and Publication

## Before implementation

Every new task must record, before adding a module/chapter or running its research:

1. Guiding question and why it matters.
2. Prerequisites and exact source/claim dependencies, including unresolved ones.
3. Mathematical/engineering contract: objects, parameter domain, conventions,
   assumptions, boundaries, approximation, and failure conditions.
4. Intended evidence, independent checks, verification commands, and acceptance
   criteria. Mark statistical criteria not applicable when inappropriate.
5. Intended document destination: main narrative topic, wave appendix, manuscript
   candidate, or exploration only. A destination is planned, not an acceptance claim.
6. Target files, integration target, wave/task ID, tags, and stopping condition.

Record research intake in YAML-frontmatter notes; routine policy work can use the
wave plan. Missing prerequisites block dependent acceptance, not unrelated tasks.

## Classification and growth

Use multiple tags; the review order is not a mutually exclusive ontology. Initial
subject tags: `conventions`, `physical-state-space`, `mean-field`, `numerics`,
`bulk`, `defects`, `exact-benchmarks`, `perturbation`, `thermodynamics`,
`static-correlations`, `dynamics`. Record model scope separately (couplings, field,
lattice), method, epistemic kind, wave, claim IDs, and document destinations.

Example: a thermal vortex calculation can carry `defects`, `thermodynamics`,
`mean-field`, and `numerics`; its model parameters state its actual domain.

Use one canonical file with links from multiple indexes. Add tags when needed,
with a definition and aliases for replaced terms. Do not move files whenever a
tag changes. Stable claim IDs and source locators survive chapter reordering.
Do not impose every field on every old file; backfill metadata during review.

## Scientific selection

Review one claim with its evidence and dependencies. Track scientific status
separately from implementation completeness and project priority. Decisions are
`unreviewed`, `supported-within-scope`, `needs-revision`, `disapproved-by-PI`,
`refuted`, or `deferred`; disapproval is not an analytical disproof.

Source/implementation existence, passing tests, BIC, and referee PASS are not
sufficient for a general physical claim. Require matched observables/conventions
and independent checks. Do not infer full-Brillouin-zone self-consistency from
Dirac-point evaluation. Do not replace operator resolvents by a scalar gap and
call the resulting coefficient exact without a justified reduction. Never infer
all-orders control from a third-order model or a single finite-cluster splitting.

The recent historical Phase 4 calculations and interpretations are disapproved
by the PI. Preserve, inventory, and audit them separately from older perturbation
work. No automatic restoration of acceptance, deletion, or overwriting older work.

## Labbook and publication ownership

Maintain one canonical labbook source root. Main chapters follow logical
questions; technical appendices follow historical waves/milestones and carry
current evidence status. Proofs and derivations live in appendices with precise
cross-references; main chapters retain enough explanation to be intelligible.
Chapters evolve with accepted knowledge and new model families; they are not a
fixed nine-chapter template or a copy of directory names. State the scope of each
result explicitly; new couplings/lattices require new contracts and comparisons.

A task/wave may produce standalone incremental LaTeX notes/PDFs including only
its contribution and required context. They are review artifacts, not competing
labbooks, and must identify their base, status, and eventual canonical destination.
Previewing the canonical book on a branch does not create a new official edition.

Use the neutral manuscript root `doc/manuscript/manuscript_kitaev_mf.tex` when the
publication migration task runs; remove SciPost branding from live entrypoints
and adopt PRL/REVTeX formatting then. Do not relabel an old PDF as a new PRL build.
Keep historical records intact; publication claim selection is a separate gate.
