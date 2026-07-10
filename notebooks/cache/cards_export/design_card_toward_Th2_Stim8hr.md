# Design card - toward_Th2 (Stim8hr)

**Verdict:** weakly reachable  |  **Confidence:** medium

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample, full panel) | 0.537 |
| reach cosine (held-out, honest) | 0.304 |
| in-sample overstatement | 0.233 |
| held-out significance (z) | 20.4 |
| knockdown-reachable (LOF) | 28.8% |
| activation-only (GOF, CRISPRa hypothesis) | 24.9% |
| neither direction | 46.2% |
| optimal library size (knee) | k = 5 |

> 46% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 2/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress - CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | IFNGR1 | 0.218 | CRISPRi (SM inhibitor available: 1 drugs) |
| 2 | ZBTB49 | 0.171 | CRISPRi/shRNA (no SM tractability) |
| 3 | DZIP1 | 0.139 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 4 | GGA1 | 0.134 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 5 | PACS1 | 0.130 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 6 | DENND4B | 0.129 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 7 | BMAL1 | 0.129 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 8 | JAK2 | 0.129 | CRISPRi (SM inhibitor available: 31 drugs) |

## Activation recipe (induce - CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | AMPD3 | 0.127 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 2 | SEPTIN2 | 0.126 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 3 | BST2 | 0.120 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 4 | MAGEH1 | 0.118 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 5 | MAL | 0.115 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 6 | INF2 | 0.115 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 7 | TMEM184C | 0.113 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 8 | CHRAC1 | 0.108 | CRISPRa or ORF overexpression (activation, not SM-druggable) |

## Optimal next-screen library (knee at k=5)

`DENND4B > JAK2 > IL7R > IRF9 > MAPK14`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*