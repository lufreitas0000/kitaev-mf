# Provisional research integration — W001

This provisional integration candidate collects source history; collection is not scientific
acceptance. `doc/audits/w001_acceptance.yaml` owns current review decisions.
Historical "complete", "certified", "exact", or "PASS" labels in imported notes,
tests, datasets, and changelogs are original assertions pending scoped review.

| Area | Current disposition |
|---|---|
| Original exact-perturbation notes | Preserved unchanged; their projection/resolvent contracts guide review |
| Later exact-perturbation additions and K-019 reconciliation | Needs review; scalar-gap formula is not an accepted exact all-orders solution |
| Historical Phase 2 sources | Unreviewed physical claims; software tests have bounded meaning |
| Historical Phase 3 sources and five preserved inputs | Unreviewed reproduction/physical claims; imported for review |
| Historical Phase 4 calculations, explanation, and standalone notes | Disapproved by Lu Freitas; preserved, not an accepted result |

The earlier perturbation work was ignored, not presumed altered. See
`doc/audits/w001_perturbation_recovery.md` for byte-identity evidence and caveats.

No disapproved Phase 4 chapter is included in the canonical labbook root. Its
standalone TeX source remains a historical review artifact, explicitly marked
with disapproval. The old standalone PDF remains in its archive tag (and a local
ignored preservation copy), not among active deliverables. Existing canonical
PDFs are inherited snapshots until rebuilt with a source receipt; do not treat
their presence as a fresh publication or claim review.

Ordinary pushes/PRs build previews and run checks. Pages/release publication is
separate: a manual main workflow dispatch with publish=true is required by the
candidate workflows. Remote main has not yet received these workflow changes.
