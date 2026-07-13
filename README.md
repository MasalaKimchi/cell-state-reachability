# Cell-State Reachability

**Before designing a multi-gene perturbation, ask whether the intervention class can
reach the target direction at all.**

Cell-State Reachability turns measured Perturb-seq effects into a directional verdict, a
compact candidate panel, and a mathematical certificate of what the measured intervention
dictionary misses.

**[Watch the 3-minute demo](https://www.youtube.com/watch?v=GJbxLxYUBMo)**
· **[Explore the live result](https://raw.githack.com/MasalaKimchi/cell-state-reachability/5b96d535ba8e33e63c43f69b4a89b744ce7699e7/app/index.html)**
· **[Read the 2-minute pitch](./SUBMISSION.md)**
· [Open the paper](./manuscript/main.pdf)
· [Bring your own target](./notebooks/bring_your_own_target.ipynb)

![Cell-State Reachability at a glance: the Th2-to-Th1 result, signed modality decomposition, method outputs, and scope.](docs/figures/fig_at_a_glance.png)

## The result

In resting primary human CD4⁺ T cells, the flagship **Th2 → Th1** direction is **partly
directionally reachable by CRISPRi**.

| Result | What it means |
|---|---|
| **0.448 held-out cosine**; above all **60** shuffled targets | The signal generalizes to held-out genes and clears the target-shuffle control (plus-one empirical p = 1/61). |
| **39% measured LOF** | Knockdown can supply part of the target direction. |
| **25% sign-flipped GOF proxy** | The remaining direction motivates a CRISPRa or de-repression arm; it is not a measured activation effect. |
| **35% neither** | This share is unmatched by either the LOF cone or its sign-flipped proxy. |
| **1.1 × 10⁻¹¹ KKT optimality violation** | The cone projection is numerically certified; biological efficacy is not. |

Across all 12 transition–condition cases, knockdown is never the majority modality (mean
LOF fraction 0.34). The same operator runs without retuning on a held-out CEBPA state in a
K562 CRISPRa screen (cosine 0.878).

The practical answer is not “this target works.” It is: **test a focused knockdown panel,
add another modality, or change the perturbation dictionary before spending the next
experiment.**

## Why this is different

Most perturbation tools predict effects or rank candidate genes. This method asks a prior
question: can non-negative combinations of effects already measured in the relevant screen
point toward the desired transcriptomic state?

Under a first-order additivity approximation, the measured effects form a convex cone. The
target is projected onto that cone:

1. **Fit** the closest non-negative mixture of measured perturbation effects.
2. **Challenge** its alignment on held-out genes and shuffled targets.
3. **Compress** the mixture into a greedy sparse candidate panel.
4. **Certify** any whole-direction mismatch with the separating residual.

The full residual is the certificate. Its largest positive coordinates rank unmet readouts
for follow-up; no individual coordinate is, by itself, a proven activation target.

The contribution is the combination of measured-effect grounding, a target-specific
reachability verdict, an outside-the-cone certificate, and held-out/null calibration—not a
new list of Th1/Th2 regulators. A structured survey of 91 prior methods found no entry with
that full combination; this is a survey-bounded novelty claim.

## Try it on your target

The fastest hands-on route uses the committed Norman K562 effect bundle; no large download
or GPU is required.

```bash
git clone https://github.com/MasalaKimchi/cell-state-reachability.git
cd cell-state-reachability
python -m pip install numpy==2.4.6 scipy==1.17.1 matplotlib==3.11.0 jupyter
jupyter notebook notebooks/bring_your_own_target.ipynb
```

Replace the example target with signed gene weights or `up` / `down` gene lists. The
notebook returns the verdict, LOF/GOF-proxy/neither split, candidate panels, reachability
spectrum, and unmet-readout hypotheses.

For the full CD4⁺ T-cell analysis, obtain the source data, then follow the ordered workflow
in [notebooks/README.md](./notebooks/README.md). The 16.8 GB Tier-2 matrix is intentionally
not committed.

## What to open

| Goal | Start here |
|---|---|
| Watch the project story | [3-minute demo on YouTube](https://www.youtube.com/watch?v=GJbxLxYUBMo) |
| See the complete argument interactively | [Live walkthrough](https://raw.githack.com/MasalaKimchi/cell-state-reachability/5b96d535ba8e33e63c43f69b4a89b744ce7699e7/app/index.html) |
| Read the submission narrative | [SUBMISSION.md](./SUBMISSION.md) |
| Inspect the scientific result | [Paper](./manuscript/main.pdf) |
| Audit claims and limitations | [Validation report](./docs/VALIDATION_REPORT.md) |
| Inspect every result and literature anchor | [Technical dossier](./docs/Technical_Dossier.pdf) |
| Use the core API | [reachability.py](./reachability.py) |
| Verify the software invariants | `bash reproduce.sh` |

## Data and scope

The primary source is Zhu et al. 2025, *Genome-scale Perturb-seq in primary human CD4⁺ T
cells* ([dataset card](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq),
[preprint](https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1)). The derived
matrix contains 33,983 perturbation–condition profiles across 10,282 readout genes and
11,526 distinct targeted genes. Acquisition details and data grain are in
[data/README.md](./data/README.md).

- This is directional feasibility relative to one measured screen, not functional rescue.
- Multi-gene panels extrapolate from single perturbations under approximate additivity.
- The GOF share uses `-E` as a proxy; real CRISPRa effects need not mirror CRISPRi.
- Donor-collapsed effects do not permit true leave-one-donor-out validation.
- Every nomination is a wet-lab hypothesis, never a validated target.

## Repository guide

| Path | Purpose |
|---|---|
| [`reachability.py`](./reachability.py) | Cone fit, certificate, nulls, sparse spectrum, and design API |
| [`app/`](./app/) | Guided walkthrough and seven interactive explorers |
| [`notebooks/`](./notebooks/) | Ordered analyses plus the bring-your-own-target notebook |
| [`results/`](./results/) | Committed tables, fact ledgers, and design cards |
| [`manuscript/`](./manuscript/) | Paper, source, figures, and canonical fact file |
| [`docs/`](./docs/) | Validation report, technical dossier, and figures |

## License

MIT. Source data are MIT-licensed through the CZI Virtual Cells Platform. Please cite Zhu
et al. 2025.
