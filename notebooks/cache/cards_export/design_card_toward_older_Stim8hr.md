# Design card - toward_older (Stim8hr)

**Verdict:** weakly reachable  |  **Confidence:** medium

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample, full panel) | 0.590 |
| reach cosine (held-out, honest) | 0.362 |
| in-sample overstatement | 0.228 |
| held-out significance (z) | 11.5 |
| knockdown-reachable (LOF) | 34.8% |
| activation-only (GOF, CRISPRa hypothesis) | 22.6% |
| neither direction | 42.6% |
| optimal library size (knee) | k = 4 |

> 39% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 1/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress - CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | DCP2 | 0.103 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 2 | ZNF579 | 0.089 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 3 | SNRNP40 | 0.084 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 4 | GALC | 0.083 | CRISPRi/shRNA (no SM tractability) |
| 5 | C3orf38 | 0.072 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 6 | NBN | 0.067 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 7 | PLEKHO1 | 0.067 | CRISPRi or antibody blockade (surface protein) |
| 8 | CALM1 | 0.067 | CRISPRi (SM inhibitor available: 2 drugs) |

## Activation recipe (induce - CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | KLF3 | 0.105 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 2 | CAMK1 | 0.086 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 3 | VANGL1 | 0.081 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 4 | ZNF260 | 0.079 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 5 | JMJD4 | 0.072 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 6 | CCDC43 | 0.071 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 7 | CD58 | 0.067 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 8 | IFT88 | 0.064 | CRISPRa or ORF overexpression (activation, not SM-druggable) |

## Optimal next-screen library (knee at k=4)

`CD3D > SMG1 > CROT > DCP2`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*