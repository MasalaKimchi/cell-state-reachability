# Cell-State Reachability

**Does a measured perturbation dictionary contain directional support for a target
transcriptomic change under a declared linear-combination model?**

This repository contains a small numerical core and a frozen retrospective case study.
It does not predict state conversion, prescribe an intervention, or validate a target.

![Measured target and effect profiles, model-relative projection, frozen retrospective challenge, and the boundary between a reasonable follow-up and untested biological claims.](docs/figures/fig_at_a_glance.png)

## Main finding

In donor-collapsed resting primary human CD4 data, screen-derived CRISPRi
differential-expression z-score profiles show split-stable alignment with a constructed
external Th1-vs-Th2 direction.

| Frozen result | Interpretation |
|---|---|
| **0.446 ± 0.010** | Mean ± SD across 12 fixed random-gene held-out splits. |
| **0.448** | Historical fixed-split held-out cosine. |
| **60/60 diagnostic shuffles below 0.448** | Plus-one empirical p = 1/61; diagnostic, not a precise tail estimate. |

The updated stress tests materially narrow the claim:

- the cone beats mean, PCA-1/PCA-5, and a matched random-cone mean, but PCA-20,
  unconstrained least squares, and ridge score higher;
- Rest is the strongest observed context; stimulated scores are weaker and Stim48hr is
  run-confounded;
- cell-cycle ablation, magnitude matching, and positive generator rescaling do not
  explain the registered score, but do not exclude all confounding;
- every admitted Rest generator was already source-flagged significant;
- only 9,831/25,672 target genes intersect the screen and only 38/50 strongest target DE
  genes survive; source-only held-out scores are 0.438 (Ota) and 0.342 (Höllbacher);
- on 126 Norman K562 CRISPRa doubles, one retired legacy threshold label and 100 staged
  `E`-to-`-E` proxy labels flip; these are sensitivity diagnostics, not modality claims.

The utility is therefore constrained, inspectable geometry—not predictive superiority.
See [updated findings](docs/FINDINGS.md), [method](docs/METHODS.md), and
[validation](docs/VALIDATION_REPORT.md).

## Run the maintained surface

```bash
python -m pip install -r requirements.txt
./reproduce.sh
```

The reproduction command runs the strict geometry tests, a synthetic demo, and checks
every frozen finding and artifact hash. Real-data inputs are intentionally not committed;
see [data acquisition](data/README.md).

```python
import numpy as np
from reachability import project_cone

effects = np.eye(4)                       # perturbations × genes
target = np.array([1.0, 0.0, -1.0, 0.0])
result = project_cone(effects, target)

print(result.geometry_status)             # outside_model_cone
print(result.cosine, result.kkt_violation)
```

The public API intentionally provides projection, KKT diagnostics, a model-relative
separator, held-out scoring, and conservative empirical p-values. Legacy biological
verdicts, signed activation proxies, “minimal” recipes, and experiment recommendations
were removed.

## Repository map

| Path | Role |
|---|---|
| [`reachability.py`](reachability.py) | Strict projection-only numerical core |
| [`tests/`](tests/) | Exact geometry, invariance, failure, and held-out tests |
| [`scripts/validate_findings.py`](scripts/validate_findings.py) | Frozen evidence and lineage validation |
| [`results/findings.json`](results/findings.json) | Canonical machine-readable findings |
| [`results/evidence/`](results/evidence/) | Selected supporting tables only |
| [`docs/FINDINGS.md`](docs/FINDINGS.md) | Human-readable updated results and limitations |
| [`docs/METHODS.md`](docs/METHODS.md) | Model and challenge semantics |
| [`ROADMAP.md`](ROADMAP.md) | Next strict-validation milestones |

## Hard boundaries

- The current case study uses donor-collapsed effects and random-gene splits.
- Its effect space is DE z-scores, not calibrated intervention dose.
- The diagnostic shuffle is not a structured biological null.
- The perturbations were not measured in established polarized Th2 cells.
- The target is restricted to the measured, sign-concordant cross-study subset.
- Direct combinations, matched CRISPRa, donor-held-out effects, chromatin, cytokines,
  viability, durability, and phenotype rescue remain untested.

## License

MIT. Source-data acquisition and licensing notes are in [data/README.md](data/README.md).
