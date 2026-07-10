# Design card - toward_older (Rest)

**Verdict:** weakly reachable  |  **Confidence:** medium

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample, full panel) | 0.596 |
| gene-panel bootstrap reach (85%, N=12) | mean 0.608, 95% CI 0.603–0.612 |
| reach cosine (held-out, honest) | 0.390 |
| in-sample overstatement | 0.206 |
| held-out significance (z) | 14.3 |
| knockdown-reachable (LOF) | 35.6% |
| activation-only (GOF, CRISPRa hypothesis) | 24.9% |
| neither direction | 39.6% |
| optimal library size (knee) | k = 4 |

> 41% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 2/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress - CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | KDM6A | 0.099 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 2 | AGAP3 | 0.086 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 3 | MAP3K4 | 0.086 | CRISPRi or SM inhibitor (ligandable: High-Quality Ligand) |
| 4 | PRPF38B | 0.082 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 5 | TYW5 | 0.082 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 6 | MOSPD2 | 0.082 | CRISPRi or antibody blockade (surface protein) |
| 7 | SAC3D1 | 0.081 | CRISPRi/shRNA (no SM tractability) |
| 8 | UBA2 | 0.080 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |

## Activation recipe (induce - CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | TMEM265 | 0.075 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 2 | H2AC11 | 0.072 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 3 | HSPBP1 | 0.071 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 4 | SEPTIN7 | 0.070 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 5 | COPS7B | 0.066 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 6 | CCND3 | 0.065 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 7 | STEEP1 | 0.064 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 8 | AGL | 0.063 | CRISPRa or ORF overexpression (activation, not SM-druggable) |

## Optimal next-screen library (knee at k=4)

`C1D > SMG1 > SUPT7L > TARS2`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*