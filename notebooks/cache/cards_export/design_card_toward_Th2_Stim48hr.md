# Design card - toward_Th2 (Stim48hr)

**Verdict:** weakly reachable  |  **Confidence:** medium

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample, full panel) | 0.530 |
| reach cosine (held-out, honest) | 0.259 |
| in-sample overstatement | 0.271 |
| held-out significance (z) | 12.6 |
| knockdown-reachable (LOF) | 28.1% |
| activation-only (GOF, CRISPRa hypothesis) | 26.1% |
| neither direction | 45.8% |
| optimal library size (knee) | k = 4 |

> 48% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 2/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress - CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | IRF1 | 0.234 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 2 | TMSB15B | 0.181 | CRISPRi/shRNA (no SM tractability) |
| 3 | GFI1 | 0.165 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 4 | GBP5 | 0.162 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 5 | DNTTIP1 | 0.161 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 6 | TAFAZZIN | 0.158 | CRISPRi or antibody blockade (surface protein) |
| 7 | SERPINB1 | 0.151 | CRISPRi or antibody blockade (surface protein) |
| 8 | MAP3K10 | 0.141 | CRISPRi (SM inhibitor available: 1 drugs) |

## Activation recipe (induce - CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | FBXO32 | 0.183 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 2 | AMPD3 | 0.119 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 3 | ATG5 | 0.119 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 4 | ICOS | 0.119 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 5 | YWHAQ | 0.119 | CRISPRa or SM agonist (ligandable: High-Quality Pocket) |
| 6 | RNF139 | 0.113 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 7 | KIF5B | 0.112 | CRISPRa or SM agonist (ligandable: Approved Drug) |
| 8 | LIMA1 | 0.111 | CRISPRa or ORF overexpression (activation, not SM-druggable) |

## Optimal next-screen library (knee at k=4)

`RSBN1L > TMSB15B > IRF1 > TAFAZZIN`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*