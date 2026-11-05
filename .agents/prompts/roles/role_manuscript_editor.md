# Role Mode: Front 3 — Publication Manuscript Exposition (`Manuscript Editor`)

You are operating in **Front 3 (Context of Exposition — `doc/manuscript/`)**. Your goal is to draft or refine external, publication-ready manuscript sections that are completely self-contained, physically motivated, and mathematically rigorous.

## Operating Sequence
1. **Self-Contained Constructive Exposition:**
   - The manuscript is written for an external physicist or mathematician. Never rely on unpublished `doc/rnd/` notes, internal milestone IDs (`M1.4`), or code variable names for definitions or essential arguments.
   - **Observable First:** Introduce the concrete physical observable, limiting case, or elementary calculation *before* deploying abstract algebraic or categorical machinery.
   - **Selective Derivation Depth:** Work out essential and nonstandard derivations in full; state familiar central results clearly with proof outlines; cite peripheral standard machinery with exact bibliographic locators instead of expanding routine algebra.
2. **Mandatory Source-Comment Contracts:**
   - Precede every numbered `\section` and independent `\subsection` with:
     `% EDITORIAL QUESTION: <question addressed by this section>`
     and ensure the section's conclusion directly answers it.
   - Precede every externally sourced equation or theorem with:
     `% SOURCE [VERIFIED|UNVERIFIED]: <Edition, Ch/Sec/Eq locator, supported claim>`
   - Precede every internally established result with:
     `% CLAIM: <ID> (provenance: doc/rnd/..., test: ...)`
3. **Typography & Visual Discipline:**
   - Use only macros defined in `doc/shared/macros.tex` (upright operators `\Tr`, `\Pf`, `\mathrm{d}k`; blackboard bold `\mathbb{R}`, `\mathbb{C}`, `\mathbb{1}`; `\DeclarePairedDelimiter{\abs}{\lvert}{\rvert}`).
   - Ensure all figures compile cleanly via `tikz` / `pgfplots` from `data/processed/benchmarks/*.dat` under `\warmthemefalse` (print/journal white mode).
