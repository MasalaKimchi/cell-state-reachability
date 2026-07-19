# Prospective active-acquisition recommender

`acquisition.py` turns CombiCone's per-combination triage score into a **batch
recommendation** for the next round of a combinatorial screen, and supports the
closed design loop (recommend → run → observe → refit → recommend).

## What it does

- **Relevance**: the training-free triage score `−agg_cos` (the single feature
  validated on Norman), or the ridge `TriageModel` once a labeled pilot exists.
- **Diversity**: greedy max-marginal-relevance (MMR). Each pick maximizes
  `relevance − diversity_weight · max_similarity_to_batch_so_far`. Similarity is
  the cosine of the candidates' predicted additive-effect directions (default)
  or gene-set Jaccard — so a batch spreads across distinct predicted biology,
  not B copies of the same hub gene.

CLI: `combicone recommend screen.h5ad --batch-size 10 --strategy diversified`.

## Retrospective validation (Norman, 131 doubles)

`simulate_campaign` replays the loop over the real screen, revealing a combo's
outcome only after it is "run". Enrichment = emergent hits discovered by the
recommender ÷ hits expected under random acquisition, at 20 combinations
acquired (≈15% of the screen):

| label regime | exploit (greedy) | diversified (w=0.5) |
|---|---|---|
| raw unreachable-residual (top tercile) | **1.96×** | 1.81× |
| noise-robust two-bar | 1.31× | 0.94× |

The honest reading:

1. On the **raw label**, the recommender front-loads discoveries at roughly
   **2×** random — consistent with the validated triage enrichment.
2. On the **strict two-bar label** the lift is modest (1.3× greedy), and the
   diversified strategy is *below* random in the first ~20 picks: diversity
   trades early exploitation for batch spread. That is the correct behavior for
   a real campaign that cannot afford a redundant batch, but it is not a
   free lunch, and we do not claim front-loading under every setting.
3. `diversity_weight` is the exploit↔explore knob. Set it to 0 for pure
   exploitation (reproduces top-B); raise it when redundancy in the batch is the
   binding constraint.

The recommender chooses what to **measure**; whether a discovered combination is
genuinely emergent is decided by `certify_emergence`, not by the acquisition
score. No probability or uncertainty is attached to the relevance score.

See `fig_acquisition.png` for the discovery curves and `tests/test_acquisition.py`
for the unit + planted-emergence tests.
