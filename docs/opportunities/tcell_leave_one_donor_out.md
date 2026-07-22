# A2 + E3 — CD4⁺ T-cell genome-scale Perturb-seq: acquisition + leakage-free leave-one-donor-out

**Track status: DONE (real data acquired; leakage-free LODO executed).**
This was the designated *honest-risk* track — the target dataset is a ~22-million-cell
screen released in Dec 2025 and could have been gated, unreleased, or too large to
reach. It turned out to be publicly streamable, so the acquisition **succeeded** and a
genuine leakage-free donor-holdout was run on real per-donor data.

Everything below is labelled **COMPUTED** (from real repo/remote data), **DESIGN-ONLY**
(harness/protocol), or **CONTEXT**.

---

## 1. Acquisition attempt — result: SUCCEEDED

**Target** — *Genome-scale perturb-seq in primary human CD4⁺ T cells maps
context-specific regulators of T cell programs and human immune traits*
(Zhu R., Dann E., et al., 2025; **doi 10.64898/2025.12.23.696273**). ~22 M cells,
4 donors, rest + stimulation, probe-based Perturb-seq.

**What was located (COMPUTED / verified this session):**

| Item | Finding |
|---|---|
| Article full text | Fetched via Unpaywall (green OA). Confirms 22 M cells, 4 donors, rest/stim. |
| Analysis code | GitHub `emdann/GWT_perturbseq_analysis_2025` (public). |
| Processed data | **Public + anonymous** on the CZI Virtual Cells Platform S3 bucket `s3://genome-scale-tcell-perturb-seq/marson2025_data/`. |
| Raw sequencing | "will be made available" via **SRA SRP643211 / GEO GSE314342** (not needed here). |

**File inventory (COMPUTED from anonymous S3 REST listing, KeyCount=32):**

| File | Size | Usable for leakage-free LODO? |
|---|---|---|
| `D*_*.assigned_guide.h5ad` (12 cell-level files) | **119–173 GB each** | No — individually prohibitive. |
| `GWCD4i.pseudobulk_merged.h5ad` | **44.57 GB** | **Yes** — aggregated by guide × **donor** × condition → retains per-donor resolution. |
| `GWCD4i.DE_stats.h5ad` | 16.79 GB | Donor-pooled DE (`~ log10_n_cells + donor_id + target`) → leaks; not used to build atoms. |
| `GWCD4i.DE_stats.by_donors.h5mu` | 16.87 GB | Two-donor **group** summaries → cannot isolate one donor. |
| `suppl_tables/DE_stats.suppl_table.csv` | 4.8 MB | Used only to **select** the panel (not to build atoms). |

**Constraints handled honestly:** local machine has ~1–2 GiB free RAM and ~34 GiB disk;
the pseudobulk file (44.57 GB) exceeds disk and RAM, and the cell-level files
(119–173 GB) exceed both by far. The full matrix was therefore **never downloaded or
loaded whole**. Instead the file was opened lazily over **HTTP range requests**
(`HTTPRangeFile` → `h5py` `driver='fileobj'`); only the planned pseudobulk rows were
streamed. `aws`/`boto3`/`s3fs` were unavailable and the `ros3` HDF5 driver failed TLS
through the sandbox proxy — the range-request path is the workaround, and it is
reproducible (`extract_tcell_lodo_substrate.py`).

**Bytes actually transferred: ~882 MB** to extract the full working substrate — 0.02×
the size of the pseudobulk file (44.57 GB) and under 0.01× any single cell-level file.

---

## 2. What kind of screen this is (decisive for the analysis)

**COMPUTED.** This is a **genome-scale *singles* screen**: each cell carries **one**
guide (`guide_type ∈ {targeting, non-targeting}`, one `perturbed_gene_name` per
pseudobulk). There are **no measured multi-gene combinations** in the data.

CombiCone's headline product — *certified combinatorial emergence* — requires a
**measured combination** to test against the single-gene cone. This screen does not
contain any, so **combination-emergence certification is not applicable here** and is
**not** claimed. The well-posed leakage-free question that this data *can* answer is the
one the shipped `results/donor_pair_transfer.json` explicitly lists as **out of scope**:
**state reachability under a leakage-free donor holdout** — i.e. does a held-out donor's
single-gene effect lie in the non-negative cone of the *other* donors' effects?
That is exactly what was run.

---

## 3. Leakage-free LODO — method (COMPUTED)

**Substrate build (`extract_tcell_lodo_substrate.py`).** Rest condition. Panel = **40**
perturbations, chosen as the top on-target-significant genes by number of downstream
trans-effects (from the pooled DE summary), **each with effective guides in all 4
donors** (9,849 genes qualify on the all-4-donor criterion; 40 taken as a
RAM/network-bounded panel). For **each donor independently**: every pseudobulk
normalised to log1p-CP10K, guide-averaged within (donor, gene), then the donor's **own**
non-targeting (NTC) mean subtracted → that donor's single-gene effect atoms. **No
cross-donor pooling at any step** — this disjointness is what makes the holdout
leakage-free. Result: 40 atoms × 18,129 genes × 4 donors, **full 40/40 coverage in
every donor**; NTC baseline norms ≈ 69.8 across donors (consistent).

**Folds (`leave_one_donor_out.py`).** For each donor *d*: build the non-negative conic
hull from the **other three** donors' atoms (`reachability.project_cone`, NNLS +
Farkas/KKT separator); project donor *d*'s atoms onto it; record unreachable-fraction,
cosine, KKT violation. Two settings: **self-gene excluded** (drop the same-gene atom
from the cone → pure trans-state geometry) and **self-gene included** (leakage-free
cross-donor reproducibility of the same perturbation). Trans-gene axis = top-2000 genes
by mean |effect|, on-target panel genes removed. KKT violation ≤ 2.2 × 10⁻¹⁵ across all
8 folds (both settings; max 2.11 × 10⁻¹⁵ at self-gene-included D3) — far below the
1 × 10⁻⁸ certification tolerance (numerically certified).

---

## 4. Leakage-free LODO — results (COMPUTED, real CD4⁺ T-cell data)

**Per-fold (self-gene included / excluded), median over 40 atoms:**

| Hold-out | Reference | Median unreachable frac (incl / excl) | Median cosine (incl) |
|---|---|---|---|
| D1 | D2,D3,D4 | 0.628 / 0.738 | 0.778 |
| D2 | D1,D3,D4 | 0.686 / 0.744 | 0.728 |
| D3 | D1,D2,D4 | 0.586 / 0.671 | 0.810 |
| D4 | D1,D2,D3 | 0.584 / 0.727 | 0.812 |

**Cross-donor transfer is real.** Held-out-donor cosine to the other-donor cone
(self-included) has median **0.78**, versus a permutation null (shuffled cone) median of
**0.04**. A held-out donor's single-gene effects are largely — but not entirely —
reachable from the other donors' cone; the residual (~0.6 unreachable fraction) is the
donor-specific component that a leakage-free holdout is designed to expose.

**The leakage contrast (headline).** When the held-out donor's **own** same-gene atom is
(improperly) added to the reference cone, the median unreachable fraction collapses from
**0.62 → 0.00** across all 160 atom-folds — every target becomes trivially "reachable".
This is the concrete inflation that a non-leakage-free donor analysis would produce, and
the direct, quantitative reason leakage-free holdout matters.

Full numbers: `results/tcell_lodo_result.json`. Figure: `results/tcell_lodo.png`.

---

## 5. How this differs from `results/donor_pair_transfer.json` (CONTEXT)

`donor_pair_transfer.json` — its own `claim_ceiling` (verbatim): *"fixed-four-donor,
published-eligibility robustness sensitivity; **not** leakage-free donor holdout, donor-
population inference, predictive utility, or **state reachability**"*. It operates on the
released **two-donor-group** summaries (6 modalities = donor pairs), so a single donor
**cannot** be isolated from the reference — hence it is explicitly not a leakage-free
holdout.

**This track closes exactly that gap.** It uses **individual per-donor** pseudobulk
effect vectors and holds one donor **fully out** of the reference cone, and it scores
**state reachability** — both items on that claim_ceiling's exclusion list. The two are
complementary, not redundant: the shipped result is a robustness sensitivity on
published summaries; this is a leakage-free donor-holdout reachability test on
individual-donor data.

---

## 6. Deliverables

| File | Type | Note |
|---|---|---|
| `leave_one_donor_out.py` | **DESIGN + COMPUTED** | Parametrised LODO harness (donor key). Two paths: singles cone-transfer (used here) and combination-certification (for any *combinatorial* multi-donor screen — untested here as this screen has no combos). Includes `build_per_donor_atoms_from_pseudobulk()` and a CLI (`--substrate … --held-out-donor …`). |
| `extract_tcell_lodo_substrate.py` | **COMPUTED** | Streams only the planned pseudobulk rows over HTTP range; builds per-donor atoms; reproducible. |
| `results/tcell_lodo_substrate.npz` | **COMPUTED** | 40 × 18,129 × 4-donor atom matrices + coverage + NTC baselines. |
| `results/tcell_lodo_result.json` | **COMPUTED** | All folds, both settings, permutation null, leakage contrast, full provenance + claim_boundary. |
| `results/tcell_lodo.png` | **COMPUTED** | 3-panel figure (leakage contrast · transfer-vs-null · per-fold stability). |
| `docs/tcell_leave_one_donor_out.md` | this note | |

---

## 7. Claim boundary (binding)

- **COMPUTED from real per-donor CD4⁺ T-cell pseudobulk** (Rest, 40-gene panel,
  4 donors). Leakage-free donor holdout via single-gene cone transfer, demonstrated end
  to end on real data.
- **NOT combinatorial emergence certification** — this screen has no measured
  combinations; the harness's combo path is provided but untested on this dataset.
- **NOT genome-wide** — a 40-gene panel bounded by the ~1–2 GiB RAM / streaming budget,
  not a screen-wide result. 9,849 genes satisfy the all-4-donor coverage criterion and
  are reachable by the same harness with more compute (or a machine co-located with the
  S3 bucket).
- **Rest condition only** — Stim8hr / Stim48hr are equally available via the same script
  (change `condition`), not run here.
- Reachability / unreachable-fraction are **model-relative** under the chosen metric,
  **not** claims of biological impossibility (the repo's own `_SCOPE`).
