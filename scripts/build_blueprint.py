#!/usr/bin/env python3
"""Universal R&D Blueprint & Epistemic Graph Generator for GitHub Pages.

Compiles repository metadata into a self-contained static HTML & Markdown portal:
- Reads PROJECT_STATE.yaml (milestones, active tasks, verified invariants)
- Reads YAML frontmatter from doc/rnd/**/*.md (questions, contracts, claims, dependencies)
- Reads LOG_DISCOVERY_AND_PIVOTS.md (non-linear anomalies, dead ends, and pivots)
- Reads provenance headers from data/processed/benchmarks/*.dat
- Renders an interactive Mermaid Epistemic DAG and warm-parchment (#FAF6EE) HTML portal.
"""

from __future__ import annotations

import argparse
import datetime
import html
import re
import shutil
from pathlib import Path
from typing import Any

import yaml


def parse_frontmatter(md_path: Path) -> dict[str, Any]:
    try:
        text = md_path.read_text(encoding="utf-8")
    except Exception:
        return {}
    match = re.match(r"^---\s*\n(.*?)\n---\s*\n", text, flags=re.DOTALL)
    if not match:
        return {}
    try:
        return yaml.safe_load(match.group(1)) or {}
    except Exception:
        return {}


def parse_pivot_log(log_path: Path) -> list[dict[str, str]]:
    if not log_path.is_file():
        return []
    rows: list[dict[str, str]] = []
    for line in log_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line.startswith("|") or "---" in line:
            continue
        cols = [c.strip() for c in line.strip("|").split("|")]
        if len(cols) >= 5 and cols[0].lower() != "date":
            rows.append(
                {
                    "date": cols[0],
                    "context": cols[1],
                    "obstruction": cols[2],
                    "resolution": cols[3],
                    "status": cols[4],
                }
            )
    return rows


def collect_dat_benchmarks(root: Path) -> list[dict[str, str]]:
    benchmarks: list[dict[str, str]] = []
    for dat_path in sorted(root.glob("data/processed/benchmarks/**/*.dat")):
        headers: list[str] = []
        try:
            for line in dat_path.read_text(encoding="utf-8").splitlines()[:10]:
                if line.startswith("#"):
                    headers.append(line.lstrip("# ").strip())
        except Exception:
            pass
        benchmarks.append(
            {
                "path": dat_path.relative_to(root).as_posix(),
                "summary": " | ".join(headers[:3]) if headers else "ASCII benchmark table",
            }
        )
    return benchmarks


def sanitize_mermaid_id(raw: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_]", "_", raw)


def build_mermaid_dag(
    state: dict[str, Any],
    rnd_notes: list[dict[str, Any]],
    pivots: list[dict[str, str]],
) -> str:
    lines: list[str] = [
        "flowchart LR",
        "  classDef milestone fill:#1A365D,color:#FAF6EE,stroke:#1A365D,stroke-width:2px;",
        "  classDef note fill:#FAF6EE,color:#2B2B2B,stroke:#1A365D,stroke-width:1.5px;",
        "  classDef invariant fill:#E6F4EA,color:#137333,stroke:#137333,stroke-width:1.5px;",
        "  classDef pivot fill:#FCE8E6,color:#C5221F,stroke:#C5221F,stroke-width:1.5px,stroke-dasharray: 4 4;",
    ]

    milestones = list(state.get("milestone_history") or [])
    active_ms = state.get("active_milestone")
    if active_ms and active_ms.get("id"):
        milestones.append(
            {
                "id": active_ms.get("id"),
                "name": active_ms.get("description", "Active Milestone"),
                "status": active_ms.get("status", "in_progress"),
            }
        )

    for ms in milestones:
        ms_id = str(ms.get("id", "M"))
        node_id = sanitize_mermaid_id(f"MS_{ms_id}")
        status = ms.get("status", "")
        label = f"{ms_id} ({status})".replace('"', "'")
        lines.append(f'  {node_id}["{label}"]:::milestone')

    for note in rnd_notes:
        nid = str(note.get("id") or note.get("file", "note"))
        m_nid = sanitize_mermaid_id(f"N_{nid}")
        title = str(note.get("title") or note.get("file", "")).replace('"', "'")
        if len(title) > 42:
            title = title[:39] + "..."
        lines.append(f'  {m_nid}["{nid}: {title}"]:::note')
        ms_ref = str(note.get("milestone", ""))
        if ms_ref:
            ms_node = sanitize_mermaid_id(f"MS_{ms_ref}")
            lines.append(f"  {ms_node} --> {m_nid}")

        for claim in note.get("claims") or []:
            c_node = sanitize_mermaid_id(f"CLM_{claim}")
            lines.append(f'  {c_node}["Claim {claim}"]:::invariant')
            lines.append(f"  {m_nid} --> {c_node}")

    for inv in state.get("verified_physics_invariants") or []:
        inv_id = str(inv.get("id", "INV"))
        inv_node = sanitize_mermaid_id(f"INV_{inv_id}")
        inv_name = str(inv.get("name", "")).replace('"', "'")
        if len(inv_name) > 38:
            inv_name = inv_name[:35] + "..."
        lines.append(f'  {inv_node}["{inv_id}: {inv_name}"]:::invariant')

    for idx, piv in enumerate(pivots[:8], start=1):
        p_node = f"PIV_{idx}"
        ctx = piv["context"].replace('"', "'")[:30]
        lines.append(f'  {p_node}["Pivot: {ctx}"]:::pivot')

    return "\n".join(lines)


def generate_blueprint(root: Path, output_dir: Path) -> None:
    state_path = root / "PROJECT_STATE.yaml"
    state: dict[str, Any] = {}
    if state_path.is_file():
        state = yaml.safe_load(state_path.read_text(encoding="utf-8")) or {}

    rnd_notes: list[dict[str, Any]] = []
    for p in sorted(root.glob("doc/rnd/**/*.md")):
        meta = parse_frontmatter(p)
        meta["file"] = p.relative_to(root).as_posix()
        if "id" not in meta:
            meta["id"] = p.stem
        rnd_notes.append(meta)

    pivots = parse_pivot_log(root / "LOG_DISCOVERY_AND_PIVOTS.md")
    benchmarks = collect_dat_benchmarks(root)
    mermaid_code = build_mermaid_dag(state, rnd_notes, pivots)

    project_name = str(state.get("project_name") or root.name)
    active_phase = str(state.get("active_phase") or "R&D Active")
    current_ms = str(state.get("current_milestone") or "N/A")
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Copy compiled PDFs and Hypothesis Contracts to output_dir so they are hosted directly on GitHub Pages
    labbook_src = root / "doc/labbook/labbook_kitaev_mf.pdf"
    manuscript_src = root / "doc/manuscript/scipost_kitaev_mf.pdf"
    contracts_src = root / "doc/HYPOTHESIS_CONTRACTS.md"
    if labbook_src.is_file():
        shutil.copy2(labbook_src, output_dir / "kitaev-mf-labbook.pdf")
        shutil.copy2(labbook_src, output_dir / "labbook_kitaev_mf.pdf")
    if manuscript_src.is_file():
        shutil.copy2(manuscript_src, output_dir / "kitaev-mf-manuscript.pdf")
        shutil.copy2(manuscript_src, output_dir / "scipost_kitaev_mf.pdf")
    if contracts_src.is_file():
        shutil.copy2(contracts_src, output_dir / "HYPOTHESIS_CONTRACTS.md")

    # 1. Write Markdown summary
    md_out = output_dir / "BLUEPRINT.md"
    md_lines = [
        f"# R&D Epistemic Blueprint: `{project_name}`",
        f"- **Active Phase:** {active_phase}",
        f"- **Current Milestone:** `{current_ms}`",
        f"- **Generated:** `{now}`\n",
        "## 1. Non-Linear Epistemic Graph (Milestones, Notes, Claims & Pivots)",
        f"```mermaid\n{mermaid_code}\n```\n",
    ]
    md_out.write_text("\n".join(md_lines), encoding="utf-8")

    # 2. Write Self-Contained Warm-Parchment HTML Portal
    invariants_rows = ""
    for inv in state.get("verified_physics_invariants") or []:
        invariants_rows += (
            f"<tr><td><code>{html.escape(str(inv.get('id', '')))}</code></td>"
            f"<td>{html.escape(str(inv.get('name', '')))}</td>"
            f"<td><code>{html.escape(str(inv.get('equation', '')))}</code></td>"
            f"<td><span class='badge'>{html.escape(str(inv.get('status', '')))}</span></td></tr>"
        )

    notes_rows = ""
    for note in rnd_notes:
        contracts = ", ".join(str(c) for c in (note.get("contracts") or []))
        claims = ", ".join(str(c) for c in (note.get("claims") or []))
        questions = "; ".join(str(q) for q in (note.get("questions") or []))
        notes_rows += (
            f"<tr><td><code>{html.escape(str(note.get('id', '')))}</code></td>"
            f"<td><code>{html.escape(str(note.get('milestone', '-')))}</code></td>"
            f"<td>{html.escape(str(note.get('title', note.get('file', ''))))}"
            + (f"<br/><small><em>Q: {html.escape(questions)}</em></small>" if questions else "")
            + "</td>"
            f"<td><code>{html.escape(contracts or '-')}</code></td>"
            f"<td><code>{html.escape(claims or '-')}</code></td>"
            f"<td><code>{html.escape(str(note.get('file', '')))}</code></td></tr>"
        )

    pivot_rows = ""
    for piv in pivots:
        pivot_rows += (
            f"<tr><td>{html.escape(piv['date'])}</td>"
            f"<td><code>{html.escape(piv['context'])}</code></td>"
            f"<td>{html.escape(piv['obstruction'])}</td>"
            f"<td>{html.escape(piv['resolution'])}</td>"
            f"<td><span class='badge'>{html.escape(piv['status'])}</span></td></tr>"
        )

    dat_rows = ""
    for b in benchmarks:
        dat_rows += (
            f"<tr><td><code>{html.escape(b['path'])}</code></td>"
            f"<td>{html.escape(b['summary'])}</td></tr>"
        )

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8"/>
  <title>{html.escape(project_name)} — R&D Epistemic Blueprint</title>
  <script type="module">
    import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';
    mermaid.initialize({{ startOnLoad: true, theme: 'neutral' }});
  </script>
  <style>
    :root {{
      --bg: #FAF6EE;
      --text: #2B2B2B;
      --accent: #1A365D;
      --rule: #D1C7B7;
      --card: #F4EFE4;
    }}
    body {{
      background-color: var(--bg);
      color: var(--text);
      font-family: "Georgia", "Times New Roman", serif;
      max-width: 1180px;
      margin: 0 auto;
      padding: 2rem;
      line-height: 1.6;
    }}
    h1, h2, h3 {{ color: var(--accent); border-bottom: 1px solid var(--rule); padding-bottom: 0.35rem; }}
    .meta {{ background: var(--card); padding: 1rem 1.25rem; border-left: 4px solid var(--accent); margin-bottom: 1.5rem; }}
    table {{ width: 100%; border-collapse: collapse; margin: 1rem 0 2rem 0; font-size: 0.95rem; }}
    th, td {{ border-bottom: 1px solid var(--rule); padding: 0.6rem 0.75rem; text-align: left; vertical-align: top; }}
    th {{ color: var(--accent); font-weight: bold; border-top: 2px solid var(--accent); border-bottom: 2px solid var(--accent); }}
    code {{ font-family: "JetBrains Mono", "Courier New", monospace; font-size: 0.88em; background: #EAE3D2; padding: 0.1rem 0.3rem; border-radius: 3px; }}
    .badge {{ background: var(--accent); color: var(--bg); padding: 0.15rem 0.45rem; border-radius: 3px; font-size: 0.8rem; font-family: sans-serif; }}
    .graph-box {{ background: #FFFFFF; border: 1px solid var(--rule); padding: 1rem; border-radius: 6px; overflow-x: auto; }}
    .pdf-cards {{ display: flex; gap: 1.25rem; flex-wrap: wrap; margin: 1.25rem 0 2rem 0; }}
    .pdf-card {{ flex: 1 1 320px; background: var(--card); border: 2px solid var(--accent); border-radius: 6px; padding: 1rem 1.25rem; }}
    .pdf-card h3 {{ margin-top: 0; border-bottom: none; padding-bottom: 0; }}
    .btn {{ display: inline-block; background: var(--accent); color: var(--bg); text-decoration: none; padding: 0.45rem 0.9rem; border-radius: 4px; font-family: sans-serif; font-weight: bold; font-size: 0.9rem; margin-right: 0.5rem; margin-top: 0.5rem; }}
    .btn-secondary {{ background: transparent; color: var(--accent); border: 1px solid var(--accent); }}
  </style>
</head>
<body>
  <h1>{html.escape(project_name)}: R&amp;D Epistemic Blueprint &amp; Compiled Monographs</h1>
  <div class="meta">
    <strong>Active Phase:</strong> {html.escape(active_phase)} &nbsp;|&nbsp;
    <strong>Current Milestone:</strong> <code>{html.escape(current_ms)}</code> &nbsp;|&nbsp;
    <strong>Updated:</strong> <code>{html.escape(now)}</code>
  </div>

  <h2>0. Compiled PDF Monographs &amp; Key Deliverables</h2>
  <div class="pdf-cards">
    <div class="pdf-card">
      <h3>Living Lab Book Monograph (PDF)</h3>
      <p>Comprehensive <strong>74-page</strong> R&amp;D monograph (11 chapters, full proofs, Grassmannian schematics, and PGFPlots benchmarks).</p>
      <a class="btn" href="kitaev-mf-labbook.pdf" target="_blank">View / Download PDF (GitHub Pages)</a>
      <a class="btn btn-secondary" href="https://github.com/lufreitas0000/kitaev-mf/releases/download/latest/kitaev-mf-labbook.pdf" target="_blank">Release PDF</a>
    </div>
    <div class="pdf-card">
      <h3>SciPost Executive Monograph (PDF)</h3>
      <p>Publication-grade <strong>9-page</strong> executive paper synthesizing the core analytical and numerical breakthroughs of Phase 1.</p>
      <a class="btn" href="kitaev-mf-manuscript.pdf" target="_blank">View / Download PDF (GitHub Pages)</a>
      <a class="btn btn-secondary" href="https://github.com/lufreitas0000/kitaev-mf/releases/download/latest/kitaev-mf-manuscript.pdf" target="_blank">Release PDF</a>
    </div>
    <div class="pdf-card">
      <h3>Hypothesis Contracts &amp; Repository</h3>
      <p>Centralized registry of all <strong>22 Hypothesis Contracts</strong> (passed, refuted, and Phase 2 gates) and source code.</p>
      <a class="btn" href="https://github.com/lufreitas0000/kitaev-mf/blob/main/doc/HYPOTHESIS_CONTRACTS.md" target="_blank">Hypothesis Contracts</a>
      <a class="btn btn-secondary" href="https://github.com/lufreitas0000/kitaev-mf" target="_blank">GitHub Repository</a>
    </div>
  </div>

  <h2>1. Epistemic Dependency Graph (Questions &rarr; R&amp;D Notes &rarr; Claims &rarr; Pivots)</h2>
  <div class="graph-box">
    <pre class="mermaid">
{html.escape(mermaid_code)}
    </pre>
  </div>

  <h2>2. Verified Physics &amp; Mathematical Invariants (<code>PROJECT_STATE.yaml</code>)</h2>
  <table>
    <thead><tr><th>ID</th><th>Invariant Name</th><th>Defining Identity / Equation</th><th>Status</th></tr></thead>
    <tbody>{invariants_rows or "<tr><td colspan='4'>No invariants recorded yet.</td></tr>"}</tbody>
  </table>

  <h2>3. Front 1 — Agentic R&amp;D Notes (Context of Discovery: <code>doc/rnd/</code>)</h2>
  <table>
    <thead><tr><th>Note ID</th><th>Milestone</th><th>Title &amp; Guiding Question</th><th>Contracts</th><th>Claims</th><th>Path</th></tr></thead>
    <tbody>{notes_rows or "<tr><td colspan='6'>No R&D notes indexed yet.</td></tr>"}</tbody>
  </table>

  <h2>4. Non-Linear Discovery &amp; Pivot Chronicle (<code>LOG_DISCOVERY_AND_PIVOTS.md</code>)</h2>
  <table>
    <thead><tr><th>Date</th><th>Module / Context</th><th>Obstruction / Anomaly</th><th>Pivot / Resolution</th><th>Status</th></tr></thead>
    <tbody>{pivot_rows or "<tr><td colspan='5'>No obstructions or pivots logged yet.</td></tr>"}</tbody>
  </table>

  <h2>5. Reproducible ASCII Benchmark Data (<code>data/processed/benchmarks/*.dat</code>)</h2>
  <table>
    <thead><tr><th>Benchmark File</th><th>Provenance Header Summary</th></tr></thead>
    <tbody>{dat_rows or "<tr><td colspan='2'>No .dat benchmarks exported yet.</td></tr>"}</tbody>
  </table>
</body>
</html>
"""
    html_out = output_dir / "index.html"
    html_out.write_text(html_content, encoding="utf-8")
    print(
        f"✓ Built R&D Blueprint Portal at {html_out} "
        f"({len(rnd_notes)} R&D notes, {len(state.get('verified_physics_invariants') or [])} invariants, "
        f"{len(pivots)} pivots, {len(benchmarks)} .dat tables)"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Build R&D Blueprint static portal.")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Repository root directory",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for HTML/MD blueprint (default: <root>/build/blueprint)",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    out_dir = args.output_dir.resolve() if args.output_dir else (root / "build" / "blueprint")
    generate_blueprint(root, out_dir)


if __name__ == "__main__":
    main()
