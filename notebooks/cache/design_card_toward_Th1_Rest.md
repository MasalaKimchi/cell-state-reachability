# Design card — toward_Th1 (Rest)

**Verdict:** partially reachable  |  **Confidence:** high

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample) | 0.627 (95% CI 0.633–0.649) |
| reach cosine (held-out, honest) | 0.446 |
| in-sample overstatement | 0.180 |
| held-out significance (z) | 45.0 |
| knockdown-reachable (LOF) | 39.3% |
| activation-only (GOF, CRISPRa hypothesis) | 25.3% |
| neither direction | 35.4% |
| optimal library size (knee) | k = 5 |

> 39% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 1/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress — CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | ICOS | 0.254 | CRISPRi or antibody blockade (surface; 4 drugs) |
| 2 | C2CD5 | 0.187 | CRISPRi or antibody blockade (surface protein) |
| 3 | ADCK2 | 0.178 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 4 | VTI1A | 0.151 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 5 | FBXW7 | 0.138 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 6 | FAM91A1 | 0.131 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 7 | AKR1C3 | 0.122 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 8 | BCL9L | 0.121 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |

## Activation recipe (induce — CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | DENND10 | 0.135 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 2 | WWP2 | 0.132 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 3 | SP100 | 0.132 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 4 | ITPR3 | 0.128 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 5 | CREBZF | 0.120 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 6 | MAPK14 | 0.113 | CRISPRa or SM agonist (ligandable: Advanced Clinical) |
| 7 | IRF1 | 0.109 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 8 | GALK1 | 0.108 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |

## Optimal next-screen library (knee at k=5)

`LAT2 > APPBP2 > RARA > ICOS > SNAP23`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*
