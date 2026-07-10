# Design card - toward_Th1 (Stim48hr)

**Verdict:** weakly reachable  |  **Confidence:** medium

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample, full panel) | 0.533 |
| reach cosine (held-out, honest) | 0.299 |
| in-sample overstatement | 0.234 |
| held-out significance (z) | 21.8 |
| knockdown-reachable (LOF) | 28.4% |
| activation-only (GOF, CRISPRa hypothesis) | 26.3% |
| neither direction | 45.3% |
| optimal library size (knee) | k = 4 |

> 48% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 3/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress - CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | FBXO32 | 0.230 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 2 | ICOS | 0.171 | CRISPRi or antibody blockade (surface; 4 drugs) |
| 3 | KLHL6 | 0.169 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 4 | RARA | 0.166 | CRISPRi (SM inhibitor available: 11 drugs) |
| 5 | ATG5 | 0.153 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 6 | KIF5B | 0.146 | CRISPRi (SM inhibitor available: 2 drugs) |
| 7 | APPBP2 | 0.130 | CRISPRi or antibody blockade (surface protein) |
| 8 | SNX4 | 0.130 | CRISPRi or antibody blockade (surface protein) |

## Activation recipe (induce - CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | GBP5 | 0.149 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 2 | ADD1 | 0.136 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 3 | SAMD8 | 0.123 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 4 | MYH9 | 0.122 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 5 | ZFP36 | 0.121 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 6 | ZNF213 | 0.120 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 7 | RAP1GDS1 | 0.117 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 8 | ALOX5 | 0.112 | CRISPRa or SM agonist (ligandable: Approved Drug) |

## Optimal next-screen library (knee at k=4)

`ATF7IP2 > FBXO32 > RARA > HPS1`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*