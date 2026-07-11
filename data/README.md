# The data

**Two tiers, both now local.** The no-login supplementary CSVs (this folder) are the
**Tier-1 warm-up** — they need no CZI / Synapse / Wiley account. The **headline method
runs on Tier-2**: the gene-level `GWCD4i.DE_stats.h5ad` effect matrix, which is **present
locally at `data/GWCD4i.DE_stats.h5ad` (16.8 GB)**. Every result in the Technical Dossier (Part 1 — Results)
(the headline verdict, the modality triage, and the generalizability transfers) is computed on it. The raw ~22M-cell dataset is
**not** needed (hundreds of GB — only this derived matrix). The fetch instructions below
(Options A/B) are how to re-obtain both tiers on a fresh clone; you do not need to run
them if the files are already in `data/`.

> Note (`.gitignore`): `data/*.csv` is currently ignored, so the 7 CSVs are present
> locally but **not tracked** by git. To ship a self-contained repo, force-add them:
> `git add -f data/*.suppl_table.csv`.

## What is local

Tier-1 **supplementary CSV tables (~36 MB total)**, present in this folder (gitignored;
`git add -f` to ship). These drive the Tier-1 warm-up and the disease-genetics layer:

| File | Rows | What it is |
|---|---|---|
| `DE_stats.suppl_table.csv` | 33,983 | Per **perturbation × culture-condition** summary: DE-gene counts, on-target knockdown effect/significance, off-target flag, cross-donor & cross-guide reproducibility. **Not** the full gene-level effect matrix. |
| `Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv` | 37,288 | Th2-vs-Th1 target signature (per-gene logFC/zscore) from **two** source contrasts (Ota 2021, Höllbacher 2021). |
| `CD4T_aging_signature_DE_results_full.suppl_table.csv` | 10,000 | CD4+ T-cell aging signature (per-gene logFC/zscore), Yaza 2022 discovery contrast. |
| `guide_kd_efficiency.suppl_table.csv` | 73,765 | Per-guide × condition knockdown QC (guide vs NTC expression, t-stat, `signif_knockdown`). |
| `sgrna_library_metadata.suppl_table.csv` | 26,504 | Per-guide design + off-target annotation (TSS distance, nearby/non-target genes, alternate alignments). *(26,504 logical rows; `wc -l` reports ~31,110 because some annotation fields contain embedded newlines.)* |
| `cluster_autoimmune_enrichment_results.suppl_table.csv` | 5,236 | Perturbation-cluster × **autoimmune-disease** GWAS-gene enrichment (odds ratio, FDR, intersecting genes) across 17 diseases and 4 gene sets. |
| `sample_metadata.suppl_table.csv` | 12 | Sample sheet: 4 donors × 3 conditions, donor demographics (age, sex, ethnicity). |

**Tier-2 (local, the headline core):** `GWCD4i.DE_stats.h5ad` — the
gene×perturbation *effect matrix*, present at `data/GWCD4i.DE_stats.h5ad` (16.8 GB,
gitignored). This is the reachability dictionary `E`: the headline Th2→Th1 verdict, the
12-cell atlas, and the modality triage all run on it. It is read selectively (one layer,
`float32`, subset to significant on-target perturbations + HVGs) and cached to
`atlas_work/inputs.npz` — never loaded whole on an 18 GB-RAM laptop. To re-obtain it on a
fresh clone, use Option A below.

## Option A — CZI Virtual Cells Platform CLI (recommended)

1. Register (free) at https://virtualcellmodels.cziscience.com/
2. Install the CLI: see https://chanzuckerberg.github.io/vcp-cli/usage/data.html
3. Search & download:

```bash
vcp data search "Primary Human CD4+ T Cell Perturb-seq" --exact
# then download the artifact(s) you need, e.g. GWCD4i.DE_stats.h5ad
```

Core artifact for this project: **`GWCD4i.DE_stats.h5ad`**
(33,983 perturbation×condition rows × 10,282 genes; layers: `log_fc`, `zscore`,
`p_value`, `adj_p_value`, `baseMean`, `lfcSE`).

## Option B — No-auth fallback (supplementary tables on GitHub)

Open, MIT-licensed supplementary tables (enough to prototype the pipeline):
https://github.com/emdann/GWT_perturbseq_analysis_2025/tree/master/metadata/suppl_tables

The files listed in the table above are the ones we already pulled from here. If
you re-clone the analysis repo, also look for an **arrayed CRISPRi validation
table** (e.g. a `Th1Th2_validation*` file) if the authors publish one — we do
**not** currently have it locally, so our polarization ground-truth check instead
relies on internal reproducibility (`guide_kd_efficiency` + `DE_stats`
cross-donor/cross-guide columns) and orthogonal literature/Open Targets evidence.

## Memory tips (CPU-only laptop)

- Load a single layer (e.g. `zscore`) as `float32`, not the full AnnData with all layers.
- Subset to perturbations with `keep_test_genes == True` and the top ~2,000 HVGs.
- Cache the reduced matrix to `.npz` so you never reload the h5ad during iteration.
