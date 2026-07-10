# Design card - toward_younger (Stim48hr)

**Verdict:** weakly reachable  |  **Confidence:** low

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample, full panel) | 0.597 |
| reach cosine (held-out, honest) | 0.363 |
| in-sample overstatement | 0.234 |
| held-out significance (z) | 6.4 |
| knockdown-reachable (LOF) | 35.6% |
| activation-only (GOF, CRISPRa hypothesis) | 23.2% |
| neither direction | 41.2% |
| optimal library size (knee) | k = 4 |

> 39% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 0/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress - CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | PLIN3 | 0.094 | CRISPRi or antibody blockade (surface protein) |
| 2 | WDR62 | 0.093 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 3 | CTNNA1 | 0.087 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 4 | NXPE3 | 0.087 | CRISPRi or antibody blockade (surface protein) |
| 5 | ZC3H12D | 0.079 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 6 | COX7A2L | 0.079 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 7 | NFATC2IP | 0.075 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 8 | UBASH3B | 0.067 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |

## Activation recipe (induce - CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | ALDH3A2 | 0.083 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 2 | CD2AP | 0.078 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 3 | TET2 | 0.076 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 4 | CENPS | 0.075 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 5 | CNKSR2 | 0.070 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 6 | KPNA3 | 0.066 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 7 | ACADVL | 0.066 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 8 | RASGRP1 | 0.062 | CRISPRa or SM agonist (ligandable: High-Quality Ligand) |

## Optimal next-screen library (knee at k=4)

`SLC9A1 > ATXN2 > MFSD10 > PLIN3`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*