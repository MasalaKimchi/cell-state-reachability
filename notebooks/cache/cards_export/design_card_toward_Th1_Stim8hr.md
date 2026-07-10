# Design card - toward_Th1 (Stim8hr)

**Verdict:** weakly reachable  |  **Confidence:** medium

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample, full panel) | 0.524 |
| reach cosine (held-out, honest) | 0.293 |
| in-sample overstatement | 0.231 |
| held-out significance (z) | 12.7 |
| knockdown-reachable (LOF) | 27.4% |
| activation-only (GOF, CRISPRa hypothesis) | 26.1% |
| neither direction | 46.5% |
| optimal library size (knee) | k = 6 |

> 49% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 2/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress - CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | AGTPBP1 | 0.166 | CRISPRi or antibody blockade (surface protein) |
| 2 | RARA | 0.151 | CRISPRi (SM inhibitor available: 11 drugs) |
| 3 | HERC1 | 0.141 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 4 | DPH1 | 0.138 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 5 | MAMLD1 | 0.135 | CRISPRi/shRNA (no SM tractability) |
| 6 | HELLS | 0.134 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 7 | ITPKB | 0.134 | CRISPRi or SM inhibitor (ligandable: High-Quality Ligand) |
| 8 | PLCL2 | 0.130 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |

## Activation recipe (induce - CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | ZBTB49 | 0.175 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 2 | AP5S1 | 0.164 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 3 | PEX16 | 0.129 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 4 | APOBEC3F | 0.128 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 5 | TXNDC5 | 0.122 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 6 | IFNGR1 | 0.118 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 7 | WDR54 | 0.118 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 8 | RTL8C | 0.113 | CRISPRa or ORF overexpression (activation, not SM-druggable) |

## Optimal next-screen library (knee at k=6)

`RARA > ICOS > PCNX3 > MAMLD1 > LAT2 > AGTPBP1`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*