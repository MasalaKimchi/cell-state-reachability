# Cell-State Reachability Oracle — Manuscript Source

Preprint manuscript + standalone Limitations & Reinforcement Plan for the
convex-cone cell-state reachability oracle (Built with Claude, Life Sciences).

## Build
Requires a TeX distribution with pdflatex + bibtex (natbib, authblk, subcaption,
booktabs, tabularx, placeins). On this project a TinyTeX bundle was used.

```
pdflatex main && bibtex main && pdflatex main && pdflatex main
pdflatex limitations_and_reinforcement_plan   # (x2)
```

## Layout
- `main.tex`                — preprint driver (preamble + \input of sections/)
- `sections/`               — 00_abstract, 10_introduction, 20_methods,
                              30_results, 40_related_work, 50_discussion, 90_supplement
- `references.bib`          — 103 verified BibTeX entries (101 prior + Replogle2022 and
                              Mejia2026 added for the reinforcement threads. Replogle2022
                              verified by DOI (OpenAlex) and PMID (PubMed E-utilities);
                              Mejia2026 is the ICML 2026 poster "Needles in the Haystack"
                              (OpenReview XsrXLPxBJw), whose expanded preprint DOI
                              10.1101/2025.10.20.683304 is cited in the entry note; 0 fabricated)
- `figures/`                — fig1-5 (main) + figS1-9 (supp), each as 300 dpi PNG + vector PDF
                              (figS5-9 added: DEG-weighting, calibration, cross-cell-type
                              transfer, certificate split-stability, certificate cross-reference)
- `manuscript_facts.json`   — locked single-source-of-truth facts sheet (every number in the
                              text is verbatim from here)
- `limitations_and_reinforcement_plan.tex` — standalone self-critique (also embedded as the
                              paper's Discussion)
- `editorial_verdict.json`  — paper-narrative handling-editor arc assessment
- `verification_log.csv`    — per-citation verification record (in the citations artifacts)

## Provenance
Data: Zhu et al. 2025 (CD4+ T-cell CRISPRi Perturb-seq, CZI VCP) and Norman et al. 2019
(K562 CRISPRa, GSE133344). Method: reachability.py (convex-cone NNLS + Farkas/KKT certificate).
