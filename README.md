# Cell-State Reachability

**Does a measured perturbation dictionary contain directional support for a target
transcriptomic change under a declared non-negative linear-combination model?**

This repository provides a strict numerical core, a systemic synthetic stress harness,
a hash-bound retrospective reconstruction, and one independent cross-modality benchmark.
It does not predict state conversion, prescribe an intervention, or validate a target.

![Measured target and effect profiles, model-relative projection, source-bound challenges, and the boundary between a follow-up measurement and an untested biological claim.](docs/figures/fig_at_a_glance.png)

## Current evidence

The primary case study uses donor-collapsed `Rest` CRISPRi effects from primary human CD4
cells. Its target is a source-study-reused, selectively constructed cross-sectional
population contrast oriented toward the reported Th1 centroid; it is neither independent
nor a trajectory.

| Source-bound result | Interpretation |
|---|---|
| **0.444 ± 0.018** (range 0.417–0.473) | Mean ± SD over 12 lexicographically frozen random-gene splits; split variation, not donor uncertainty. |
| **25,672 → 11,616 → 7,960 → 6,188 genes** | Target union → shared sources → concordant signs → screen-measured registered target. |
| **+0.087 / +0.090 logFC cosine** | Ota→Höllbacher / reverse mean improvement over the better of a mean ray and best single atom across six correlated splits. Directional only: normalized RMSE does not improve. |
| **Arce Spearman 0.148 / 0.084 / 0.088** | Cross-study CRISPRi-transcript to CRISPR-KO CD25 ranking alignment in resting Teff / stimulated Teff / resting Treg; modest and context dependent. |

The systemic harness independently checks NNLS solutions, provenance and axis faults,
group leakage, exact uncertainty around familywise error, common-response confounding,
random-gene optimism, and sign-selection inflation.

The formerly displayed **0.446 ± 0.010** and **p = 1/61** were retired. Their deleted
pipeline depended on an unhashed `inputs.npz` whose gene order was not preserved, so the
multisplit table and shuffles cannot be source-reconstructed. The separately archived
fixed split (0.448154) does reproduce within `3e-10`; it is retained only as provenance in
the source report, not as the headline.

The utility is constrained, inspectable geometry and falsifiable transfer tests—not
predictive superiority. See [findings](docs/FINDINGS.md), [methods](docs/METHODS.md),
[technical validation](docs/VALIDATION_REPORT.md), and the
[expert-reviewed execution plan](docs/SCIENTIFIC_VALIDATION_PLAN.md).

## Run the maintained surface

```bash
python -m pip install -r requirements.txt
./reproduce.sh
```

The small reproduction path runs the numerical tests, demo, systemic harness check, and
artifact-lineage validation. External scientific data are gitignored. With the registered
inputs available:

```bash
python -m pip install -r requirements-external.txt
python scripts/run_source_reconstruction.py --check results/source_reconstruction.json
python scripts/run_arce_external_validation.py --check
```

Acquisition, hashes, licenses, and claim ceilings are in [data/README.md](data/README.md).

```python
import numpy as np
from reachability import project_cone

effects = np.eye(4)  # perturbations × genes
target = np.array([1.0, 0.0, -1.0, 0.0])
result = project_cone(effects, target)

print(result.geometry_status)  # outside_model_cone
print(result.cosine, result.kkt_violation)
```

The public API emits projection geometry, KKT diagnostics, a model-relative separator,
and held-out scores. It intentionally emits no biological verdict, recipe, dose, or
candidate ranking.

## Repository map

| Path | Role |
|---|---|
| [`reachability.py`](reachability.py) | Projection-only numerical core |
| [`validation.py`](validation.py) | Oracle, label/provenance, grouped-split, and multiplicity contracts |
| [`scripts/run_validation_harness.py`](scripts/run_validation_harness.py) | Deterministic systemic synthetic stress harness |
| [`scripts/run_source_reconstruction.py`](scripts/run_source_reconstruction.py) | Full-file-hash-bound target and cross-source reconstruction |
| [`scripts/run_arce_external_validation.py`](scripts/run_arce_external_validation.py) | Independent CRISPRi-transcript to CRISPR-KO CD25 benchmark |
| [`results/findings.json`](results/findings.json) | Canonical machine-readable findings |
| [`docs/SCIENTIFIC_VALIDATION_PLAN.md`](docs/SCIENTIFIC_VALIDATION_PLAN.md) | Ordered statistical, ML, and biological execution program |

## Hard boundaries

- The primary effects are donor-collapsed and random-gene splits are correlated.
- The target was constructed from sources reused by the source study.
- Source-transfer baselines are limited; directional gains do not imply magnitude accuracy.
- Arce S1 is aggregate CRISPR-KO screening, not donor-level validation.
- Established polarized Th2 cells, donor/guide holdout, direct combinations, matched
  CRISPRa, chromatin, durability, fitness, and functional state conversion remain untested.

## License

MIT. Source-data licenses and citation requirements are recorded separately in
[data/README.md](data/README.md).
