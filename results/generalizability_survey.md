# Generalizability and application map — convex-cone reachability for cell-state engineering

*Companion to `dataset_catalog.csv` (13 candidate datasets) and `tractability_grounding.csv` (6 targets). Every accession below was retrieved live from NCBI GEO this session; no accession is asserted from memory.*

## (a) Is the method confined to the CD4+ T-cell CRISPRi dataset?

No. The method's confinement is decided by its **input contract**, not by the dataset it was demonstrated on. The algorithm consumes exactly two objects:

1. a matrix **E** (P perturbations x G genes) of *measured* perturbation effect vectors, and
2. a target direction **d** in the same G-dimensional gene space.

It then asks a feasibility question — is `d` inside the convex cone `{ E^T w : w >= 0 }` — and returns a null-calibrated verdict, a minimal recipe if reachable, or a Farkas infeasibility certificate (the specific genes `d` wants UP that no non-negative combination can deliver) if outside. This is plain NNLS plus linear algebra: CPU-only, no training, no GPU, no dependence on assay chemistry. Any experiment that yields (effect matrix, target direction) in a shared gene space is admissible. The repo's Th2->Th1 demo is one instance of that contract; the parallel Norman/Weissman 2019 K562 CRISPRa demo (GSE133344) is a second instance of the *same code, unchanged*.

What *is* dataset-dependent is the validity of three modelling assumptions. The method is dataset-agnostic within the contract; its **interpretation** is only as good as where these hold:

- **Non-negativity (`w >= 0`).** Encodes "a perturbation can be applied or not, but not negatively applied." This matches loss-of-function screens (CRISPRi/KO: you can knock a gene down, not "anti-knock-down" it) and one-directional dictionaries generally (a drug at a dose, a cytokine added). It is *relaxed* by the signed extension, which splits the target norm into a loss-of-function-reachable part, an activation-only part, and neither — so the same machinery covers gain-of-function bases (ORF overexpression, CRISPRa) by treating activation as a distinct sign.
- **Additivity (co-perturbation = sum of single effects).** A modelling prior, not a law. It is *directly testable* on any combinatorial dataset: fit the observed double-perturbation vector against the cone of its singletons and read the residual. Datasets with paired perturbations (GSE146194, GSE90546) are exactly the substrate for validating or bounding this assumption.
- **Linear / shared gene space.** Effect vectors and target must live in one comparably-normalized coordinate system. Within an assay this holds by construction; *across* assays (e.g. CRISPRi effect vectors vs L1000 drug signatures) it requires mapping to a common gene set and comparable units (log-fold-change space). This is a preprocessing obligation, not a property of the method.

## (b) Candidate-dataset catalog

Thirteen public datasets spanning five perturbation modalities. All 13 are GEO series hosted on `ftp.ncbi.nlm.nih.gov`, which is on this sandbox's network allowlist; three representative accessions were HTTP-HEAD probed this session (GSE314342, GSE92872, and the L1000 series GSE92742) and each returned HTTP 200, so every catalog series is reachable over GEO FTP-over-HTTPS via that same server. By contrast, portals that host processed matrices off-GEO — CELLxGENE downloads, figshare, Single Cell Portal, clue.io, Zenodo — all failed at the socket (HTTP 000) when probed and are marked accordingly in the CSV. Full per-row reachability questions are in `dataset_catalog.csv`.

| Accession | Modality | Cell system | Organism | ~Perturbations | Combo |
|---|---|---|---|---|---|
| GSE314342 | Perturb-seq (CRISPRi knockdown) | Primary human CD4+ T cells | H. sapiens | ~genome-scale (thousands of genes) | no |
| GSE146194 | Perturb-seq (CRISPRi + CRISPRa, single & paired guides) | K562 (CML) | H. sapiens | ~100s (singles + pairs) | yes |
| GSE90546 | Perturb-seq (CRISPRi knockdown) | K562 (CML) | H. sapiens | ~90 (UPR branch + epistasis pairs) | yes |
| GSE90063 | Perturb-seq (Cas9 knockout) | BMDCs + K562 | M. musculus; H. sapiens | ~24 TFs | no |
| GSE153056 | CRISPR knockout + surface-protein readout (ECCITE-seq) | THP-1 (monocytic) | H. sapiens | ~26 immune-checkpoint regulators | no |
| GSE139944 | Small-molecule / chemical perturbation | A549, K562, MCF7 | H. sapiens; M. musculus | 188 compounds (multi-dose) | no |
| GSE92742 | Small-molecule + shRNA/ORF L1000 signatures | ~70+ cell lines | H. sapiens | ~20,000+ perturbagens (compound + genetic) | no |
| GSE70138 | Small-molecule L1000 signatures | ~cell-line panel | H. sapiens | ~5,000+ compounds | no |
| GSE216463 | ORF overexpression (gain-of-function, non-CRISPR) | Human embryonic stem cells (hESC) | H. sapiens | ~1,800 TFs (3,548 ORF isoforms) | no |
| GSE202186 | Cytokine / ligand perturbation (non-CRISPR) | In vivo lymph-node immune cells | M. musculus | 86 cytokines | no |
| GSE92872 | Perturb-seq (Cas9 knockout, CROP-seq) | Jurkat (CROP-seq) + HEK293T/3T3 (species-mixing control) | H. sapiens; M. musculus | ~29-89 TCR-pathway guides | no |
| GSE124703 | Perturb-seq (CRISPRi knockdown) | Human iPSC-derived neurons | H. sapiens | ~hundreds (survival + neuronal genes) | no |
| GSE278572 | CRISPR + surface-protein readout (Perturb-CITE-seq) | Primary human CD4+ Tregs / Teffs | H. sapiens | ~tens (regulatory-circuit genes) | no |

**What each modality family unlocks:**

- **Genome-scale / pooled CRISPR Perturb-seq (GSE314342, GSE90063, GSE92872, GSE124703, GSE153056, GSE278572).** The archetypal loss-of-function basis. Establishes cross-cell-type generality: the same non-negative-knockdown cone is built in CD4+ T cells, dendritic cells, monocytes, iPSC-derived neurons, Tregs/Teffs, and a Jurkat T-cell line. Each asks whether a stimulus- or disease-defined state is reachable by knockdown alone, and — when outside — hands back an activation certificate.
- **Combinatorial CRISPR (GSE146194, GSE90546).** Carry paired perturbations, so they are the substrate for *testing the additivity assumption* rather than assuming it: does the measured double lie in the cone of its singletons?
- **Chemical / small-molecule (GSE139944 sci-Plex; GSE92742 / GSE70138 LINCS L1000).** Replace the genetic basis with a drug basis. This is the cross-modality bridge to therapeutics (section c).
- **Gain-of-function overexpression (GSE216463, Joung TF atlas).** A ~1,800-TF *activation* dictionary. Directly instantiates the sign that CRISPRi cannot reach, so an "activation-only" residual from a knockdown certificate becomes a concrete reachability query on this basis.
- **Cytokine / ligand (GSE202186, Immune Dictionary).** A non-genetic, physiological perturbation basis (86 cytokines in vivo) — reachability of an immune-cell state by a ligand cocktail rather than an edit.

## (c) Cross-modality mapping to drug / small-molecule signatures (LINCS L1000 / CMap)

The input contract does not care whether a row of E came from a genetic edit or a compound. Substituting the perturbation basis with **drug effect vectors** (L1000 landmark-gene signatures) turns the same feasibility question into two therapeutic questions:

1. **Drug-combination design.** With E = drug signatures and d = a desired transcriptional state, "is `d` in the cone?" returns the *minimal non-negative mixture of compounds* whose additive signature reaches the target — a combination-design recipe, with the null-calibrated verdict guarding against over-fitting a target into a high-dimensional drug space.
2. **Disease-signature reversal.** Set `d = -(disease signature)` (the reverse of a case-vs-control shift). Reachability then asks whether some non-negative drug mixture *reverses* the disease signature (the Connectivity Map premise, generalized from single-drug to combinations). When it is **outside**, the Farkas certificate names the genes the reversal requires UP that no available drug delivers — i.e. the residual that needs a genetic or gain-of-function modality.

**Data-access caveat.** L1000 processed signatures (Level 5 GCTx) are large and primarily distributed via clue.io, which is not reachable here; the raw GEO series (GSE92742, GSE70138) *are* reachable over GEO FTP but require the L1000 processing stack to assemble signatures. So the L1000 mapping is architecturally in-contract and the accessions are live, but building the signature matrix is a data-engineering step outside this sandbox's network reach.

## (d) Four-direction application map

Each direction is an (input basis -> target direction -> decision the verdict supports) triple.

**1. Other Perturb-seq / CRISPR screens (cross-cell-type generality).**
- *Input:* measured single-gene knockdown/knockout effect matrix from any pooled screen (GSE90063, GSE92872, GSE124703, GSE153056, GSE278572).
- *Target:* a stimulus- or polarization-defined transcriptional shift in that cell system (e.g. LPS-stimulated DC state; TCR-activated T-cell state).
- *Decision supported:* whether the state is achievable by knockdown alone in that cell type, the minimal guide set if so, and — if outside — which genes must instead be activated. Demonstrates the method is not tied to one cell type or one screen chemistry.

**2. Drug / small-molecule signatures (combination design + signature reversal).**
- *Input:* drug effect vectors (sci-Plex GSE139944; L1000 GSE92742/GSE70138).
- *Target:* a desired cell state, or the reverse of a disease signature.
- *Decision supported:* the minimal drug combination that reaches the state, or the verdict that reversal is outside the drug cone plus the certificate genes that need a non-drug modality. Turns a ranking problem into a feasibility-plus-recipe problem.

**3. Cell reprogramming / differentiation / fate control.**
- *Input:* a gain-of-function basis — TF overexpression (Joung TF atlas GSE216463) and/or a cytokine basis (GSE202186).
- *Target:* a developmental or therapeutic cell state defined from a reference atlas (e.g. the CELLxGENE-hosted dopaminergic-neuron signature, Kamath et al. 2022, collection b0f0b447 on CELLxGENE — metadata reachable via connector; matrix download blocked).
- *Decision supported:* which non-negative combination of TF (or cytokine) perturbations reaches the target fate — a directed-differentiation recipe on an activation basis, and the certificate for fates that no measured activator spans.

**4. Disease reversal / target discovery.**
- *Input:* any measured perturbation basis in the disease-relevant cell type (genetic or drug).
- *Target:* `-(disease signature)`.
- *Decision supported:* reachable -> a combinable, prioritizable intervention set; outside -> a Farkas certificate of genes that must be moved UP, which is then filtered by therapeutic tractability (section e) to decide the intervention *modality* — repurposing vs cell-engineering.

## (e) Therapeutic tractability grounding — from verdict to modality

A certificate says *which genes* a target demands be activated (or which the recipe knocks down). It does not say whether those genes are *druggable*. Grounding the certificate against Open Targets tractability (live GraphQL, this session) converts the reachability verdict into a modality decision. Full table in `tractability_grounding.csv`.

| Gene | Role in a certificate | SM/AB tractable? | # drugs/candidates | Modality implied |
|---|---|---|---|---|
| TBX21 | Th1 master TF — ACTIVATE | no (only PROTAC/degrader annotations) | 0 | cell-engineering / CRISPRa |
| STAT4 | Th1 signalling TF — ACTIVATE | no | 0 | cell-engineering / CRISPRa |
| FOXP3 | Treg master TF — ACTIVATE | no | 0 | cell-engineering / CRISPRa |
| GATA3 | Th2 master TF — (knockdown for Th1) | no | 0 | CRISPRi achievable; activation not druggable |
| JAK1 | signalling kinase — INHIBIT | yes (Approved Drug, High-Quality Ligand, Druggable Family) | 25 | drug repurposing (ChEMBL CHEMBL2835: 17 curated inhibitor mechanisms) |
| DGKA | immuno-metabolic enzyme — INHIBIT | high-quality ligand (SM), 0 approved | 0 | medicinal chemistry (chemically tractable enzyme) |

The split is the bridge from verdict to intervention:

- **Certificates that demand transcription-factor ACTIVATION are largely outside small-molecule space.** All four master TFs (TBX21, STAT4, FOXP3, GATA3) return **zero small-molecule and zero antibody tractability and zero drugs/clinical candidates**. The only positive tractability buckets for them are PROTAC/degrader annotations (`PR:`) — a *loss-of-function* modality, i.e. the wrong direction for a certificate that requires the gene UP. So a "must-activate-TBX21" certificate is a **cell-engineering / CRISPRa problem**, not a repurposing one. This is exactly why the Th2->Th1 axis is a cell-engineering target: the genes the certificate demands activated cannot be turned on with a pill.
- **Certificates that demand enzyme INHIBITION map onto existing drug modalities.** JAK1 is fully small-molecule tractable (approved drugs, high-quality ligands, druggable family), with 25 drugs/candidates in Open Targets and 17 curated inhibitor mechanisms in ChEMBL (target CHEMBL2835). DGKA is chemically tractable (high-quality ligand) though without an approved drug. An inhibition certificate over such genes is a **repurposing / medicinal-chemistry problem**.

**Consequence.** The reachability verdict plus the certificate plus tractability annotation together route each engineering goal to the correct intervention class: a knockdown recipe over CRISPRi-achievable genes, a drug combination over SM-tractable genes, or a CRISPRa / synthetic-biology programme for the TF-activation residual that no measured non-negative combination and no existing drug can deliver.

---
*Reachability probes (this session): `ftp.ncbi.nlm.nih.gov` returned HTTP 200 for the three representative accessions HTTP-HEAD probed (GSE314342, GSE92872, GSE92742); all 13 catalog series are GEO deposits on that same allowlisted host and are marked reachable on that basis. Off-GEO hosts probed — CELLxGENE download host, figshare, Single Cell Portal, clue.io, Zenodo — returned HTTP 000 (blocked). Replogle et al. 2022 genome-scale Perturb-seq (RPE1/K562) is **not deposited in GEO** — GEO returned zero series for its title and essential-gene queries — and its processed data is distributed via figshare (blocked here); genome-scale/combinatorial CRISPR is therefore represented in the catalog by Replogle et al. 2020 (GSE146194), which is GEO-hosted and reachable.*
