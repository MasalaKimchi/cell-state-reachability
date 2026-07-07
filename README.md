# Counterfactual Biology Explorer

**Built with Claude — Life Sciences Hackathon (Research Track)**

*What is the minimal set of gene perturbations that shifts a primary human CD4+ T cell from one transcriptional state toward a target "healthier" state — and how confident should we be?*

---

## TL;DR

Rather than reporting *what is different* between two cell states (differential
expression), this project generates *interpretable, falsifiable counterfactual
hypotheses* about **what changes might move a cell toward a target state**. It does
this by treating each genome-scale CRISPRi perturbation as a measured causal
"effect vector" and finding the **smallest set of perturbations** whose combined
effect best matches a desired transcriptomic shift. Every hypothesis ships with an
explicit **confidence score** built from the dataset's own reproducibility metrics,
held-out validation, and orthogonal literature/database evidence.

## Why this framing (and an honest scope note)

The hackathon abstract originally targeted Acute Myeloid Leukemia (AML) using the
Marson/Pritchard Perturb-seq dataset. **That dataset is a genome-scale CRISPRi
Perturb-seq screen in primary human CD4+ T cells — not AML/HSPC data.** We
therefore reframe the identical *method* onto a state transition that this dataset
can actually support:

- **Primary axis:** Th1 ↔ Th2 helper-T polarization balance
- **Secondary axis:** CD4+ T-cell aging signature (aged → young-like)

Both target-state signatures are provided directly by the dataset authors. The
polarization signature ships from **two independent source contrasts** (Ota 2021,
Höllbacher 2021), which we use for a cross-source robustness check in lieu of an
arrayed wet-lab validation table (not present in our local data). See
[`ROADMAP.md`](./ROADMAP.md) for the full rationale and the analysis catalog.

## Data

Marson & Pritchard labs, *Genome-scale perturb-seq in primary human CD4+ T cells*
(Zhu et al., bioRxiv 2025). Hosted on the CZI Virtual Cells Platform under an MIT
license.

- Dataset card: https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- Preprint: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1
- Analysis repo & supplementary tables: https://github.com/emdann/GWT_perturbseq_analysis_2025

**Feasibility on a laptop (CPU-only):** the raw dataset is ~22M cells and is *not*
tractable on a laptop. We build on the authors' **precomputed derived artifacts**,
and the work is split into two tiers by what's needed.

**Tier 1 — already local (checked into `data/`, ~36 MB of CSVs):** enough for the
graded core — target signatures, per-perturbation QC/effect summaries, guide QC,
off-target design, autoimmune-disease enrichment, and donor metadata.

| Local file | Rows | What it is |
|---|---|---|
| `DE_stats.suppl_table.csv` | 33,983 | Per **perturbation × condition** summary (DE-gene counts, on/off-target flags, cross-donor & cross-guide reproducibility). *Summary, not the gene-level matrix.* |
| `Th2_Th1_polarization_signature…csv` | 37,288 | Th2→Th1 target signature, two source contrasts |
| `CD4T_aging_signature…csv` | 10,000 | Aged→young target signature |
| `guide_kd_efficiency.suppl_table.csv` | 73,765 | Per-guide knockdown QC |
| `sgrna_library_metadata.suppl_table.csv` | 31,110 | Guide design / off-target annotation |
| `cluster_autoimmune_enrichment…csv` | 5,236 | Perturbation-cluster × 17 autoimmune diseases |
| `sample_metadata.suppl_table.csv` | 12 | 4 donors × 3 conditions, demographics |

**Tier 2 — one download away (not local yet):** the full gene-level effect matrix,
needed only for the sparse reconstruction solver.

| Artifact | Shape | What it is |
|---|---|---|
| `GWCD4i.DE_stats.h5ad` | 33,983 pert×cond × 10,282 genes | Per-perturbation effect matrix (logFC, z-score, p) — **input to the reconstruction solver** |
| `GWCD4i.pseudobulk_merged.h5ad` | guide×donor×cond × 18,129 genes | Pseudobulk expression profiles |

Working from the ~34k × 10k matrix (subset to significant on-target perturbations +
HVGs) is comfortably CPU-tractable once fetched. See [`data/README.md`](./data/README.md).

## Method (baseline-first, honestly benchmarked)

1. **Perturbation dictionary** `E ∈ R^{P×G}`: each row is one perturbation's measured
   causal effect on the transcriptome (z-scored logFC), per stimulation condition.
2. **Target direction** `d ∈ R^G`: the desired transcriptomic shift (e.g. the Th2→Th1
   signature, or the reverse-aging vector).
3. **Minimal set** = sparse selection (LASSO / Orthogonal Matching Pursuit / greedy
   forward selection) of perturbations whose weighted sum best reconstructs `d`.
   Sparsity is the operationalization of "minimal."
4. **Confidence** = combination of (a) dataset-native reproducibility (cross-guide and
   cross-donor correlation, on-target knockdown significance, off-target flags),
   (b) bootstrap **stability selection** frequency, (c) held-out-donor generalization,
   and (d) orthogonal **literature/Open Targets** evidence for the nominated gene.

A deliberately simple **linear/additive baseline is the primary model, not a
fallback.** A 2025 *Nature Methods* benchmark (Ahlmann-Eltze, Huber & Anders) shows
current deep-learning perturbation predictors do **not** yet beat simple linear
baselines — so any DL component here is an explicitly-optional comparison, never an
unbenchmarked claim.

## Honest limitations (read before trusting any output)

- **CRISPRi is loss-of-function only.** We can directly nominate *knockdowns*.
  Gain-of-function hypotheses are extrapolations that this assay cannot test.
- **Additivity is an assumption.** The screen perturbs single genes; any multi-gene
  "minimal set" is an *untested combinatorial extrapolation* that ignores epistasis.
- **Transcriptome ≠ phenotype.** Matching a transcriptional signature does not prove
  functional rescue.
- **The target state is a proxy.** "Healthier" is operationalized as a transcriptomic
  signature, not a clinical outcome.
- Outputs are **ranked, falsifiable hypotheses for future experimental validation** —
  not conclusions.

## Repository layout

```
counterfactual-biology-explorer/
├── README.md
├── ROADMAP.md            # 7-day plan, rationale, risks
├── LICENSE               # MIT
├── requirements.txt
├── environment.yml
├── data/                 # Tier-1 supplementary CSVs are tracked; large h5ad stays gitignored (data/README.md)
├── src/
│   ├── data_loader.py    # load DE_stats / pseudobulk, subset, QC
│   ├── target_states.py  # build Th1/Th2 & aging target vectors
│   ├── counterfactual.py # sparse minimal-set solver
│   ├── confidence.py     # reproducibility + stability + held-out scoring
│   ├── evidence.py       # PubMed / Open Targets literature support
│   └── pathways.py       # gene-set / pathway interpretation
├── app/
│   └── explorer.py       # lightweight Streamlit hypothesis explorer
├── notebooks/            # exploratory analysis (see notebooks/README.md)
└── tests/
    └── test_counterfactual.py
```

## Quickstart

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# See data/README.md to fetch the DE_stats artifact, then:
python -m src.data_loader --check
```

## License

Code: MIT. Data: MIT (CZI Virtual Cells Platform). Please cite Zhu et al. 2025.
