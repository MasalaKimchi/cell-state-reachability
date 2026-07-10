# Design card - toward_younger (Stim8hr)

**Verdict:** weakly reachable  |  **Confidence:** low

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample, full panel) | 0.588 |
| reach cosine (held-out, honest) | 0.349 |
| in-sample overstatement | 0.239 |
| held-out significance (z) | 4.7 |
| knockdown-reachable (LOF) | 34.6% |
| activation-only (GOF, CRISPRa hypothesis) | 22.8% |
| neither direction | 42.6% |
| optimal library size (knee) | k = 4 |

> 40% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 0/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress - CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | KLF3 | 0.112 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 2 | VANGL1 | 0.086 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 3 | CAMK1 | 0.077 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 4 | JMJD4 | 0.074 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 5 | VCPKMT | 0.073 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 6 | PAN3 | 0.073 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 7 | FBXO8 | 0.070 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 8 | CWC27 | 0.070 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |

## Activation recipe (induce - CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | C3orf38 | 0.084 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 2 | SNRNP40 | 0.064 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 3 | ABRAXAS1 | 0.062 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 4 | AMER1 | 0.062 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 5 | CD320 | 0.061 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 6 | SBF1 | 0.061 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 7 | CALM1 | 0.060 | CRISPRa or SM agonist (ligandable: Approved Drug) |
| 8 | RNF144A | 0.060 | CRISPRa or ORF overexpression (activation, not SM-druggable) |

## Optimal next-screen library (knee at k=4)

`TRIT1 > HVCN1 > PAN3 > ALKBH1`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*