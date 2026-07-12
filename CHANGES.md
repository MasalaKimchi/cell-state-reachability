# Changes — final-submission polish

A pass over the repo for the final submission: five technical corrections (A–E), removal of
process-documentation markdown, and consolidation of the seven interactive explorers into one guided
walkthrough. Both PDFs were rebuilt from corrected source. The working tree is left uncommitted for
review; nothing is committed and no history is rewritten.

Baseline: `6b060b2` ("Clean-up"). Net: 21 text files changed (+2143 / −411), 2 files deleted, 7 files
moved, 1 file added, 2 PDFs rebuilt.

---

## A — Dead reproduction command (Technical Dossier)

`docs/Technical_Dossier.md` §7.5 gave a reproduction command that pointed at a script
(`_nb04_script.py`) that does not exist in the repo.

- Replaced `cd notebooks && python _nb04_script.py` with `bash reproduce.sh` (pytest + self-test) plus
  `jupyter nbconvert --execute notebooks/04_experimental_design_toolkit.ipynb` to regenerate the design
  cards, and reworded the provenance sentence to point at the notebook rather than the missing script.
- Added a self-consistency clause: the dossier separately notes that `nbconvert --execute` cannot run
  inside the build sandbox (the kernel's TCP-socket bind is blocked), so the §7.5 command is annotated
  "needs a local Jupyter kernel" and cross-references that note. The command works on a normal machine;
  the annotation removes the apparent contradiction inside the document.

## B — "71% of the ceiling" identity disclosed (manuscript)

The √f_LOF "ceiling" and the reported 71% figure are not an independent validation — they are the same
quantity. `sections/30_results.tex` now states this outright: *"This bound is not an independent check:
√f_LOF coincides with the in-sample reachable cosine (0.627 here, to 10⁻⁴), so '71% of the ceiling' is
the same statement as the held-out cosine recovering 71% of the in-sample fit."* A matching parenthetical
was added to `sections/50_discussion.tex`.

## C — Missing foundational citations (manuscript)

Added two references the method's mathematics rests on:

- **Kuhn & Tucker (1951)** — cited alongside Farkas at the certificate's introduction
  (`sections/10_introduction.tex`).
- **Moreau (1962)** — cited at the signed-decomposition definition (`sections/20_methods.tex`), which is
  an orthogonal split along mutually polar cones in the sense of Moreau.

Both added to `references.bib`; `main.bbl` regenerated.

## D — `Mejia` reference reconciled to the verified preprint (manuscript + dossier)

The shipped entry cited "Needles in the Haystack" as an ICML 2026 paper, which resolves to no record.
Reconciled to the verified arXiv preprint the authors actually posted:

- **Mejia et al., *Diversity by Design: Addressing Mode Collapse Improves scRNA-seq Perturbation
  Modeling on Well-Calibrated Metrics*, arXiv:2506.22641 (2025).**
- Manuscript: rekeyed `Mejia2026` → `Mejia2025`, rewrote the bib entry, fixed the three in-text `\citet`
  sites and two stale mentions in `manuscript/README.md` and `manuscript_facts.json`.
- Dossier: retitled all seven mentions (survey entry, trust table, prose), dropped an unverifiable second
  bioRxiv DOI in favor of the single verified arXiv record.

## E — Polish (manuscript + SUBMISSION.md)

- **Attrition statistics sourced** — `SUBMISSION.md` now carries a footnote citing the sources behind
  the "9 in 10 fail" / "Phase II is the lowest transition" figures (Thomas et al. 2016; Wong, Siah & Lo
  2019; Harrison 2016).
- **Proxy caveat** — the Limitations paragraph in `sections/50_discussion.tex` now notes that matching a
  target in expression space is not the same as achieving the corresponding cellular phenotype.
- **Fig 1C caption softened** — `sections/30_results.tex` now frames the realizing intervention by gene
  role (activation for effectors, de-repression for negative regulators such as IKZF3/CBLB) rather than
  overclaiming a single mechanism.
- **Base-rate citation** — added Wong, Siah & Lo (2019) to the base-rate sentence in
  `sections/40_related_work.tex`.

---

## Documentation cleanup

Removed two files that documented the *process* of building the repo rather than reporting its results
— of interest to the authors, not to readers:

- `CONSOLIDATION.md` (deleted) — a log of a file-consolidation pass.
- `REPRODUCIBILITY_CHECK.md` (deleted) — a self-audit worksheet.

Both are git-tracked, so they remain fully recoverable from history; they were moved to the system Trash
rather than force-removed. No file references either one.

Tightened three reader-facing markdown files that had drifted toward narrating their own edit history:

- `results/README.md` — removed a "consolidation pass" tagline and a "removed files" table; kept the
  exact canonical-file counts.
- `CLAUDE.md` — rewrote a backward-looking figure-audit log into a forward guardrail (the two
  decomposition triples — signed 39/25/35 and one-sided 39/31/30 — are both canonical and must not be
  conflated).
- `SUBMISSION.md` — see E above.

---

## Interactive explorers → one guided walkthrough

The seven explorers were strong individually but arrived as a flat list with no connecting story.
Consolidated into a single scrollytelling page while keeping each explorer fully intact.

- **`app/index.html`** is now a six-chapter narrative — **Problem → Reframe → Verdict → Trust → Impact
  → Novelty** — with interpretive prose ("what you're looking at" / "what it means") framing each step.
  Each explorer is embedded **live** via an inline `<iframe srcdoc>`, so the page is one self-contained
  file that works on double-click (no server, no external dependencies) and preserves every explorer
  bit-for-bit. Sticky chapter nav, a hero with the four headline numbers, and a footer linking the
  paper, dossier, and reproduce script.
- **`app/explorers/`** — the seven explorers moved here (via `git mv`, history preserved) and remain
  independently openable; the narrative links each as "open standalone ↗". These are what you
  screen-record for the demo.
- **`app/_build_index.py`** (new) regenerates `index.html` from the explorers plus the prose — re-run it
  after editing any explorer.
- Path references updated in `app/DEPLOY.md`, `DEMO_VIDEO_SCRIPT.md`, and `README.md`.

### Null-model reconciliation

The demo script had flagged that the flagship explorer and Fig 1C report the null at different z-values
(≈45 vs ≈24). The cosine estimates describe the same result, but the original explanation of the z gap
was incorrect. The explorer's `z = 45` (cosine 0.446) comes from only eight shuffles; the manuscript's
`z ≈ 24` (cosine 0.448) uses 60. A z-score does not lawfully scale with shuffle count: eight draws give
an unstable estimate of the null standard deviation. The current explorer labels the eight-shuffle z as
descriptive and uses the 60-shuffle result (above all 60 controls; empirical p=1/61) for inference.

---

## Rebuilt artifacts

- **`manuscript/main.pdf`** — rebuilt from corrected source (bundled TinyTeX; pdflatex → bibtex →
  Tectonic/BibTeX, zero undefined references). 25 pages. Citation graph is a clean 1:1 (69 cite keys, 69
  bibitems, none undefined, none unused).
- **`docs/Technical_Dossier.pdf`** — rebuilt from corrected `.md` (pandoc → print-quality browser PDF). 101 pages after tighter, more readable typesetting.
  Protected one-sided decomposition numbers (0.31 / 0.30, "heuristic upper bound", 39/25/35) verified
  intact; digits render correctly.

## Known limitation

`app/previews/index_preview.png` is a static thumbnail of the *previous* hub design and is now stale.
It could not be regenerated here — headless browsers cannot launch in this environment (an OS-level
sandbox restriction). It is not shown as an image anywhere (only listed in a directory tree), so it is
not reader-visible; regenerate it from `app/index.html` in any browser before publishing if you want the
thumbnail current.
