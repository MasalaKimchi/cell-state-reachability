# CONSOLIDATION.md — Final Submission Cleanup

*A complete record of what was merged, moved, or removed to bring
`cell-state-reachability` to its final, navigable submission state, and why.
Every removal is recoverable from the system Trash and from git history.*

---

## Why

The repository had accumulated the by-products of a week-long research push:
eleven loose narrative markdowns in `docs/`, duplicate figure copies across
three directories, throwaway helper scripts alongside the real analysis
drivers, and interpretation-only result sidecars. The science was sound but the
structure was hard to follow. This pass keeps **only final files, code, and
modifiable PDFs** — research-direction progressions are intentionally out of
scope — and adds a clear navigation front door.

The work was executed by four parallel sub-agents (Dossier, Results audit,
Code & manuscript hygiene, Front door), each operating on an isolated copy of
its inputs; **no sub-agent deleted or moved anything on disk**. All removals
were merged into a single, user-approved Trash operation applied by the
coordinator, after a content-preservation audit and a repo-wide
reference-repointing pass.

---

## 1 · New deliverables (created)

| File | What it is |
|---|---|
| `SUBMISSION.md` | Top-level landing page — 2-minute elevator pitch (Problem → Method → Key Result → Impact) with a "Start here" nav table. Every headline number verified verbatim against `manuscript/manuscript_facts.json`. |
| `docs/Technical_Dossier.pdf` | **162-page** merged dossier — all technical write-ups in one modifiable PDF (pandoc → weasyprint). Clickable TOC, 15 figures, 40 tables. |
| `docs/Technical_Dossier.md` | Editable markdown source of the dossier (self-contained: all figures reference `figures/*.png` relatively). |
| `results/README.md` | Catalog of all 61 retained result files (column-level descriptions + a "Removed" table documenting the 4 deletions). |
| `README.md` | **Rewritten** as the repository map / navigation front door (was a 19.8 KB mixed framing+rubric doc; now a 10.6 KB clean map). |

`app/index.html` was polished in place (stale footer that named four
already-superseded docs repointed to the dossier); its 7 explorer links and all
tile numbers were verified. `DEMO_VIDEO_SCRIPT.md` was verified against the
facts sheet with **zero hard-number drift** and left unchanged.

---

## 2 · The Technical Dossier merge

Nine narrative documents were merged **verbatim** (every number, gene name,
table row, and citation preserved) into one dossier:

| Source doc | Became |
|---|---|
| `RESULTS.md` | Part 1 · Results |
| `NOVELTY.md` | Part 2 · Novelty, Impact & Field Positioning |
| `RELATED_WORK.md` | Part 3 · Related Work (91-method survey) |
| `CAUSAL.md` | Part 4 · Trust & Causal Inference |
| `GENERALIZABILITY_SURVEY.md` | Appendix A · Generalizability |
| `CROSS_CELLTYPE_TRANSFER.md` | Appendix B · Cross-Cell-Type Transfer |
| `REINFORCEMENT_RESULTS.md` | Appendix C · Reinforcement Analyses |
| `FORMULATION_ASSESSMENT.md` | Appendix D · Formulation Assessment |
| `REVIEWER_RESPONSE.md` | Appendix E · Response to Reviewer 2 |

`REVIEWER_RESPONSE.md` was added as Appendix E after a content-preservation
audit found its four resolved reviewer verdicts (G4/D4/T2/T3), reframes, and
ready-to-paste Limitations paragraph were **unique final content** not present
elsewhere — it is the write-up for the retained `results/reviewer2_*.csv`
files. During the merge, 44 intra-document cross-references were repointed to
Part/Appendix references and 3 prior-merge bookkeeping blurbs were de-duplicated
(all provenance numbers kept).

---

## 3 · Moved to Trash (34 items)

All items below were moved to the system Trash under a single user approval.
**Recoverable from Trash; also in git history.**

### Documents merged into the dossier (9) + process docs (2)
The 9 merged docs listed in §2 — the 8 from the merge bundle plus
`docs/REVIEWER_RESPONSE.md`, which was added as **Appendix E** after the
content-preservation audit. All 9 were Trashed *because* their content is
preserved verbatim in `docs/Technical_Dossier.md/.pdf`. Plus two
process/progression documents with no final content:
- `docs/ROADMAP.md` — 3-day hackathon build plan (progression; out of scope). Its plain-language framing already lives in the new README.
- `docs/FIGURE_AUDIT.md` — figure-standardization process log; zero references, no final content.

### Results files (4)
| File | Category | Reason |
|---|---|---|
| `results/a6_regulator_index.json` | unreferenced stray | Regenerable matrix-row lookup; no scientific result; zero-ref. |
| `results/reviewer2_deg_survival_meta.json` | unreferenced stray | Interpretation-only sidecar; all numbers in its referenced CSV twin. |
| `results/reviewer2_stim48_confound_meta.json` | unreferenced stray | Interpretation-only sidecar; all numbers in its referenced CSV twin. |
| `results/verification_dashboard.png` | superseded | Oldest file (Jul 7), zero-ref; panels now in `docs/figures/` + `reviewer2_trust_checks.png`. |

*(Of 15 zero-reference results files flagged in the blueprint, only these 4 were
removed — 11 were kept as unique load-bearing content: `benchmark_comparison.csv`,
`design_cards.pdf`, `norman_table2/3/4`, `heldout_modality_per_seed.csv`,
`figure_facts.json`, `a5_control_panel.json`, and three methods sidecars carrying
unique constants. A zero ref-count was treated as a signal, not an auto-removal.
`a2_expected_output_schema.json` was corrected from remove→keep after it was found
referenced 4×.)*

### Throwaway code (6) — all verified zero live references
| File | Reason |
|---|---|
| `notebooks/_nb03_script.py` | Flattened `.py` mirror of notebook 03; figures saved by the notebook itself. |
| `notebooks/_nb04_script.py` | Flattened `.py` mirror of notebook 04; figures saved by the notebook itself. |
| `notebooks/_nb07_verify.py` | Print-only cross-cell-type verification scratch; subsumed by notebook 07. |
| `notebooks/_run_nb_inprocess.py` | One-off in-process notebook runner; not in the reproduce path. |
| `scripts/_split_stability_worker.py` | One-off ProcessPool worker; output consumed by nothing. |
| `scripts/build_nbB.py` | One-off builder that *generated* notebook 09 (now a retained deliverable). |

### Manuscript stray (1)
- `manuscript/references.bib.bak` — editor backup of `references.bib`.

### Figure de-duplication (7)
- `notebooks/figures/` — whole directory (46 files, 100% gitignored, regenerable notebook outputs; 18 verified byte-identical to `docs/figures/`, the rest unique regenerable figures).
- 6 untracked `docs/figures/*.pdf` strays (`fig1_reachability_spectrum`, `fig2_decomposition_certificate`, `fig_certificate_split_stability`, `nb04_fig1_design_summary`, `nb04_fig2_modality_triage`, `nb04_fig3_reliability`) — each a PDF twin of a committed `.png` the docs actually embed.

### Build cruft (5)
`__pycache__/`, `notebooks/__pycache__/`, `handoff/` (untracked scratch),
`analysis_cache/formulation_probe.npz` (untracked stray), `.DS_Store`.

---

## 4 · Reference repointing (no dangling links)

Because 14 retained files referenced the removed docs, every reference was
repointed to the corresponding dossier Part/Appendix **before** removal, so
nothing dangles:

- **Text/docs:** `CLAUDE.md`, `data/README.md`, `notebooks/README.md`, `results/README.md`, `scripts/a2_conditional_reachability_scaffold.py`, `scripts/run_iv_compliance.py`
- **App HTML:** `app/pharma_capability.html`, `app/reachability_explorer.html`
- **Notebooks (raw-text edits, JSON validity + minimal diffs preserved):** `02`, `03`, `04`, `05`, `07`, `09`
- **Provenance strings:** `manuscript/manuscript_facts.json`

A repo-wide sweep confirms **zero remaining references** to any removed
document. `README.md` was rewritten wholesale, so its 43 old references were
replaced by the new repo map.

---

## 5 · Verified and kept (integrity checks)

- **Manuscript is complete and self-contained** for a `pdflatex + bibtex` build: `main.tex` + all 7 `sections/` inputs, `references.bib`, `main.bbl` (66 `\bibitem`s), and all 7 style/bst files (`icml2026.sty/.bst`, `algorithm`, `algorithmic`, `fancyhdr`, `forloop`, `multido`). All 14 figures (`fig1–5`, `figS1–9`; png+pdf) present. `limitations_and_reinforcement_plan.tex` is a standalone companion. (Full report: sub-agent artifact `manuscript_integrity.json`.)
- **`table*` vs `norman_table*` — BOTH kept as canonical.** They are *not* duplicates: `table1–5` are the headline result on primary human CD4⁺ T cells (Zhu 2025 CRISPRi, target toward_Th1); `norman_table1–5` are the second-dataset generalization on Norman 2019 K562 CRISPRa (target held-out CEBPA). Different datasets, assays, cell types, and schemas.
- **Figure source of truth:** `docs/figures/*.png` (13 committed PNGs) for the repo-facing dossier; `manuscript/figures/*` (fig1–5, figS1–9, png+pdf) for the paper.
- **Notebooks** all re-validated as valid `nbformat 4` JSON after repointing, with minimal (reference-line-only) diffs.

---

## 6 · Recovery

Nothing is permanently lost. Every Trashed item is in the system Trash and in
git history (`git show HEAD:<path>`, or `git checkout HEAD -- <path>`). The 9
merged documents survive verbatim inside `docs/Technical_Dossier.md`.

---

*Generated during final submission consolidation, 2026-07-11.*
