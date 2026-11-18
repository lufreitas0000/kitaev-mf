# External Reference Benchmark Data (`data/external/`)

Producer–Consumer boundary contract under `ADR-003` separating independent physical/mathematical oracle software from `kitaev-mf` 4-Majorana Gaussian `scMFT`:

| Subdirectory | Upstream Producer Repository | Exported Reference Data |
| :--- | :--- | :--- |
| `kitaev_numerics/` | [`kitaev-numerics`](https://github.com/lufreitas0000/kitaev-numerics) | Spin-space ED ($N \le 32$), KPM/OFTLM thermodynamics, Clarabel 2-RDM SDP lower bounds, HTE, and iPEPS/VUMPS curves |
| `kitaev_walk_geometry/` | [`kitaev-walk-geometry`](https://github.com/lufreitas0000/kitaev-walk-geometry) | Exact Majorana random-walk moments, Cauchy PV resolvents, and single-defect LDOS tables |
| `historical_unverified/` | Legacy pre-split pilot runs (2026-11-07) | Unverified Spin-ED and 2-RDM SDP pilot `.dat` files retained until replaced by `kitaev-numerics` Sprint 1 exports |

## Mandatory Companion Metadata (`<stem>.meta.json`)

Every `.dat` file imported into `kitaev_numerics/` or `kitaev_walk_geometry/` must include a companion `<stem>.meta.json` checked by `scripts/sync_external_data.py --verify`:

```json
{
  "source_repo": "kitaev-numerics",
  "git_sha": "<40-char commit SHA>",
  "sha256": "<sha256 of sibling .dat file>",
  "unit_convention": "pauli_K_minus_1_h111",
  "geometry": "honeycomb_torus_3x3_N18",
  "parameters": {"K": -1.0, "field_axis": [1, 1, 1]}
}
```
