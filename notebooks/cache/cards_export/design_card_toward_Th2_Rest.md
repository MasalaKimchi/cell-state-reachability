# Design card - toward_Th2 (Rest)

**Verdict:** partially reachable  |  **Confidence:** high

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample, full panel) | 0.643 |
| gene-panel bootstrap reach (85%, N=12) | mean 0.654, 95% CI 0.645–0.664 |
| reach cosine (held-out, honest) | 0.458 |
| in-sample overstatement | 0.184 |
| held-out significance (z) | 20.9 |
| knockdown-reachable (LOF) | 41.3% |
| activation-only (GOF, CRISPRa hypothesis) | 23.5% |
| neither direction | 35.2% |
| optimal library size (knee) | k = 4 |

> 36% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 2/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress - CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | MAPK14 | 0.207 | CRISPRi (SM inhibitor available: 28 drugs) |
| 2 | IRF1 | 0.180 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 3 | WWP2 | 0.173 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 4 | ILVBL | 0.166 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 5 | SIN3A | 0.157 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 6 | DENND10 | 0.146 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 7 | GALK1 | 0.139 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 8 | PPP5C | 0.135 | CRISPRi (SM inhibitor available: 1 drugs) |

## Activation recipe (induce - CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | ADCK2 | 0.141 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 2 | ICOS | 0.140 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 3 | CCND3 | 0.126 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 4 | ZBTB22 | 0.124 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 5 | C2CD5 | 0.118 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 6 | XAF1 | 0.117 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 7 | PLIN3 | 0.112 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 8 | ZNF627 | 0.110 | CRISPRa or ORF overexpression (activation, not SM-druggable) |

## Optimal next-screen library (knee at k=4)

`IRF1 > RSBN1L > DENND4B > SFXN1`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*