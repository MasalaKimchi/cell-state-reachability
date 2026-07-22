# A1 — NCI-ALMANAC Drug-Combination Result: Verified & Correctly Scoped

**Status:** VERIFIED (all reported numbers reproduce exactly from raw CSVs) and RECONCILED
(the strategy report's "certified-emergence on measured combinations" phrasing is an overclaim;
corrected below).

**Dataset:** NCI-ALMANAC (CellMiner v2.15, 2025-09-17) — ComboScore + single-agent Z-scores,
projected on 60 NCI-60 cancer cell lines (identical axes across single-agent and combo files).
**Source:** `https://discover.nci.nih.gov/cellminer` (allowlisted NCI portal).

---

## 1. What was verified (COMPUTED from real repo data)

Recomputed directly from the two staged raw files:
- `data/nci_almanac/drug_singleagent_cone_reachability.csv` (89 rows)
- `data/nci_almanac/drug_triage_vs_synergy.csv` (3,802 rows)

against the published summary `results/drug_combination_generalization.json`.

### Counts (all confirmed)
| Quantity | JSON | Recomputed | Match |
|---|---|---|---|
| ALMANAC drugs | 105 | — (matrix-level, not re-derivable from staged CSVs) | reported |
| Drugs with single-agent atoms | 89 | 89 rows | ✅ |
| Pairs total | 5,355 | — | reported |
| Usable pairs | 3,802 | 3,802 rows | ✅ |

### Test (a) — Single-agent leave-one-out cone reachability
Each of the 89 single agents is projected onto the non-negative conic hull of the **other 88**
single-agent activity vectors (leave-one-out); the certificate reports residual fraction, cosine,
inside/outside verdict, and the KKT stationarity violation.

| Metric | JSON | Recomputed | Match |
|---|---|---|---|
| Drugs outside cone (leave-one-out) | 89 | 89 / 89 (0 inside) | ✅ |
| Mean residual fraction | 0.5228331685118958 | 0.5228331685118957 | ✅ (last-digit float) |
| Max KKT violation | 1.554e-15 | 1.554e-15 | ✅ |
| Most unique (top-5) | Streptozocin, Decitabine, Vemurafenib, Estramustine, Celecoxib | identical | ✅ |
| Most redundant (bottom) | Triethylenemelamine, Thiotepa, Uracil Mustard, Chlorambucil | identical order | ✅ |

Residual-fraction range [0.111, 0.861]; cosine range [0.508, 0.994]. The most-redundant drugs are
all classic alkylating agents (shared mechanism of action) — the geometry recovers drug MoA
without supervision. **This is the faithful, in-scope geometric result.**

### Test (b) — Triage-vs-ComboScore enrichment (non-circular by construction)
- **Atoms** = single-agent activity Z-score vectors (independent measurement).
- **Triage score** = −cos between the two single-agent atoms (CombiCone machinery, UNCHANGED).
- **Label** = pair ComboScore (an INDEPENDENT synergy measurement, never seen by the triage score).
- `is_synergy` is a **relative** label: top-quartile of `combo_mean` (verified: `combo_mean > Q75`
  reproduces `is_synergy` with 100.00% agreement; base rate = 25.01%, 951 / 3,802 positives).

| Metric | JSON | Recomputed | Match |
|---|---|---|---|
| Training-free Spearman (triage_score vs combo_mean) | −0.00717 | −0.007172, p=0.6584 | ✅ |
| Training-free top-20 enrichment | 1.2 | 1.199 | ✅ |
| Recalibrated learned Spearman (OOF vs combo_mean) | +0.10753 | +0.107534, p=2.98e-11 | ✅ |

Head-to-head over the top-quartile `is_synergy` target (recomputed for the figure):

| Metric | Training-free −cos triage | Learned (OOF, recalibrated) |
|---|---|---|
| ROC-AUC | 0.510 | 0.568 |
| Average precision (PR) | 0.282 | 0.302 |
| Peak top-k enrichment | 1.70× @ top-3% | 1.54× @ top-3% |
| Enrichment @ top-20% | 1.17× | 1.29× |

**Reading:** the training-free −cos triage does **not** transfer to ALMANAC drug synergy
(ρ = −0.007, n.s.; ROC-AUC ≈ 0.51, indistinguishable from chance overall). A learned model
recalibrated on a labelled pilot recovers a **modest but significant** signal (ρ = +0.108,
p = 3e-11; ROC-AUC 0.57). This REPRODUCES CombiCone's own stated boundary on a new modality:
*the certificate/cone geometry transfers unchanged; the triage score must be recalibrated
per screen.*

> Note on the enrichment curve: the training-free score has a slightly higher **peak** at the very
> top (top ~3%, ~1.7×) but collapses to chance by top-20%, whereas the learned score sustains
> enrichment across the ranking (1.2–1.3× through top-40%). Consistent with the paper's thesis
> that the certificate's value is geometric feasibility, not ranking accuracy.

---

## 2. CRITICAL RECONCILIATION — the strategy report overclaims

**The overclaim (strategy report `FIELD_OPPORTUNITIES.md`):** describes A1 as a
*"certified-emergence result on measured combinations."*

**Why that is wrong.** Certified emergence in CombiCone requires the **measured combination
effect vector** — you project the actual combination's multi-axis effect onto the cone of the
single-agent atoms and emit an infeasibility certificate when it lies outside. That test was
**never run here**, because the raw per-combination growth vectors were not obtained.

**What A1 actually is:** a **modality-generalization TRIAGE demonstration**. It uses
single-agent atoms + an INDEPENDENT scalar ComboScore label. There are two in-scope tests —
single-agent leave-one-out reachability (geometry) and triage-vs-independent-label enrichment
(ranking) — and **neither is certified emergence on a measured combination.** The repo's own
`claim_boundary` states this correctly; only the external strategy report drifted.

### Live re-attempt of the raw combo-vector download (this session)
Per the task, I re-attempted the allowlisted NCI portal for raw combination growth vectors:

| Endpoint | Result |
|---|---|
| `discover.nci.nih.gov/cellminer/` (portal root) | **200 OK** (reachable) |
| CellMiner download manifest | Lists `DTP_NCI60_ALMANAC_COMBO_SCORE.zip` (= the ComboScore **label**, already in repo) and `DTP_NCI60_ZSCORE.zip` (single-agent atoms). **No raw per-combination growth-vector dataset is offered.** |
| `wiki.nci.nih.gov/.../ComboDrugGrowth_Nov2017.zip` (the file that holds raw combo growth) | **HTTP 403 Forbidden** |

**Conclusion:** raw combination growth vectors remain unobtainable through the allowlisted route
(portal serves only the aggregated ComboScore; the wiki file with raw growth data returns 403).
The repo's `claim_boundary` is **accurate**. We therefore do **NOT** claim certified emergence on
measured drug combinations. If the wiki file is later obtained (institutional access / manual
download), the absolute-combination projection can be added — that is the only missing piece.

---

## 3. Claim boundary (verbatim, from `results/drug_combination_generalization.json`)

> "Modality-generalization demonstration on real FDA-approved drug combinations. The
> certificate/cone geometry runs unchanged and is mechanism-coherent; prospective triage requires
> per-screen recalibration. We use single-agent activity as atoms and ComboScore as an INDEPENDENT
> synergy label (non-circular). We do NOT have raw combination growth vectors (NCI wiki file blocks
> automated download), so we do not reproduce the absolute-combination projection; the single-drug
> leave-one-out reachability + independent-label triage are the faithful in-scope tests."

This writeup stays strictly within that boundary.

---

## 4. Corrected one-paragraph A1 framing for the strategy report

> **A1 — NCI-ALMANAC (real FDA-approved drug combinations), modality-generalization demo.**
> On NCI-ALMANAC (CellMiner v2.15; 60 NCI-60 lines), CombiCone's certificate geometry runs
> **unchanged**: all 89 single agents with atoms are certified outside the leave-one-out cone of
> the others (max KKT 1.6e-15, mean residual fraction 0.52), and the geometry is
> mechanism-coherent — the most-redundant drugs are all classic alkylating agents. Using
> single-agent atoms and an **independent** ComboScore label (non-circular by construction), the
> training-free −cos triage does **not** transfer to drug synergy (Spearman −0.007, n.s.; ROC-AUC
> 0.51), while a per-screen-recalibrated learned model recovers a modest, significant signal
> (Spearman +0.108, p=3e-11; ROC-AUC 0.57) — reproducing CombiCone's stated boundary on a new
> modality: *the certificate transfers; the triage score must be recalibrated per screen.*
> **Scope (important):** this is a modality-generalization **triage** result, **not** certified
> emergence on measured combinations — raw per-combination growth vectors are not available (NCI
> portal serves only the aggregated ComboScore; the wiki growth file returns HTTP 403), so the
> absolute-combination cone projection was not run and no such claim is made.

---

## Provenance
- COMPUTED from real repo data: all numbers in §1 (recomputed from the two staged CSVs;
  reproduce `results/drug_combination_generalization.json` exactly).
- LIVE network check (§2): NCI CellMiner portal + wiki, this session.
- No synthetic data used. No wet-lab or prospective claim is made.
