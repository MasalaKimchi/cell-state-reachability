# Design card - toward_younger (Rest)

**Verdict:** partially reachable  |  **Confidence:** low

## Reachability

| metric | value |
|---|---|
| reach cosine (in-sample, full panel) | 0.626 |
| gene-panel bootstrap reach (85%, N=12) | mean 0.636, 95% CI 0.629–0.642 |
| reach cosine (held-out, honest) | 0.401 |
| in-sample overstatement | 0.225 |
| held-out significance (z) | 9.0 |
| knockdown-reachable (LOF) | 39.2% |
| activation-only (GOF, CRISPRa hypothesis) | 21.8% |
| neither direction | 39.0% |
| optimal library size (knee) | k = 4 |

> 36% of the addressable shift needs INDUCTION (CRISPRa/ORF); a knockdown-only screen structurally cannot reach it. 0/15 knockdown targets have an approved/clinical drug for repurposing.

## Knockdown recipe (suppress - CRISPRi/degrader/inhibitor)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | COPS7B | 0.095 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 2 | CCND3 | 0.094 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 3 | SEPTIN7 | 0.088 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 4 | CLINT1 | 0.086 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |
| 5 | HSPA1L | 0.086 | CRISPRi or SM inhibitor (ligandable: Structure with Ligand) |
| 6 | PBK | 0.084 | CRISPRi or SM inhibitor (ligandable: High-Quality Ligand) |
| 7 | AGL | 0.083 | CRISPRi or antibody blockade (surface protein) |
| 8 | UHRF2 | 0.082 | CRISPRi/degrader (no SM tractability; degrader PREDICTED only) |

## Activation recipe (induce - CRISPRa/ORF hypothesis; not testable in a knockdown screen)

| rank | gene | weight | delivery call |
|---|---|---|---|
| 1 | TYW5 | 0.083 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 2 | AGAP3 | 0.081 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 3 | ZNF544 | 0.081 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 4 | NAA80 | 0.077 | CRISPRa or SM agonist (ligandable: Structure with Ligand) |
| 5 | PDE4B | 0.069 | CRISPRa or SM agonist (ligandable: Approved Drug) |
| 6 | PRPF38B | 0.067 | CRISPRa or ORF overexpression (activation, not SM-druggable) |
| 7 | ACACA | 0.066 | CRISPRa or SM agonist (ligandable: Advanced Clinical) |
| 8 | SEC23IP | 0.065 | CRISPRa or ORF overexpression (activation, not SM-druggable) |

## Optimal next-screen library (knee at k=4)

`FAM98A > DEPDC5 > NFAT5 > COIL`

---
*in-sample reach cosine is optimistic; the held-out value is the honest estimate. The activation recipe is a sign-flip (CRISPRa) hypothesis, not a knockdown-screen result.*