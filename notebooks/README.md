# Notebooks

Exploratory, day-by-day. Keep heavy exploration here; promote stable logic into `src/`.

Suggested sequence (mirrors `ROADMAP.md`):

1. `01_data_load_qc.ipynb` — load the local Tier-1 CSVs, assert shapes, index guide
   QC, sanity-check a known regulator (Day 1). *(Tier-2 h5ad optional, if fetched.)*
2. `02_target_states.ipynb` — build & visualize Th1/Th2 and aging target vectors,
   plus the Ota-vs-Höllbacher concordance core (Day 2).
3. `03_counterfactual.ipynb` — Tier-1 directional ranking; Tier-2 greedy/OMP/lasso
   with k-vs-alignment curves (Day 3).
4. `04_confidence_benchmark.ipynb` — reproducibility, off-target audit, stability
   selection, leave-one-donor-out (n=4), random null, linear baseline (Day 4).
5. `05_disease_evidence_pathways.ipynb` — autoimmune-enrichment linkage +
   PubMed/Open Targets/Consensus evidence + enrichment (Day 5).

Kept out of version control by default (see `.gitignore` if you add outputs).
