# Design card - toward_older (Stim48hr)

**Verdict:** weakly reachable  |  **Confidence:** low

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample, full panel) | 0.566 |
| reach cosine (held-out, honest) | 0.331 |
| in-sample overstatement | 0.235 |
| held-out significance (z) | 6.9 |
| knockdown-reachable (LOF) | 32.1% |
| activation-only (GOF, CRISPRa hypothesis) | 25.5% |
| neither direction | 42.4% |
| optimal library size (knee) | k = 5 |

> 44% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 2/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress - CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | TRAF4 | 0.085 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 2 | RASGRP1 | 0.081 | CRISPRi or SM inhibitor (ligandable: High-Quality Ligand) |
| 3 | VCF2 | 0.080 | CRISPRi/shRNA (no SM tractability) |
| 4 | CNKSR2 | 0.079 | CRISPRi or antibody blockade (surface protein) |
| 5 | RMI1 | 0.078 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 6 | KYAT1 | 0.078 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 7 | ALDH3A2 | 0.075 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 8 | CENPS | 0.074 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |

## Activation recipe (induce - CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | TBC1D10A | 0.075 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 2 | NFATC2IP | 0.074 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 3 | COX7A2L | 0.073 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 4 | CTNNA1 | 0.070 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 5 | FZD3 | 0.065 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 6 | ZNF224 | 0.064 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 7 | CSTPP1 | 0.062 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 8 | ZNF583 | 0.062 | CRISPRa or ORF overexpression (activation, not SM-druggable) |

## Optimal next-screen library (knee at k=5)

`ZAP70 > SMG1 > CROT > WDFY1 > CD5`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*