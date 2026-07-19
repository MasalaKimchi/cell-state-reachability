# Tutorial

[`tutorial.ipynb`](tutorial.ipynb) is a deterministic, reader-facing walkthrough of the
project's fail-closed reachability method. It runs on bounded synthetic examples and reads
only the frozen, data-free validation-harness summary. It requires no `slice/` directory,
download, credential, or external scientific dataset.

The tutorial covers:

- the portable effect-dictionary contract and safe NPZ round-trip;
- strict inside-cone status versus a declared cosine-threshold coverage rule;
- the model-relative infeasibility separator;
- why random-gene holdout is diagnostic only, with module holdout and frozen harness
  context;
- certificate acquisition order versus a separately computed realized order on supplied
  measured effect atoms; and
- fail-closed labeled target alignment, perturbation-labeled coefficients, and signed gap
  coordinates.

The notebook emits geometry and diagnostics, not a biological verdict, dose, recipe, gene
importance ranking, or wet-lab recommendation.

## Run

From the repository root:

```bash
python -m pip install -r requirements.txt
python -m pip install -e ".[tutorial]"
mkdir -p /tmp/cell-state-reachability-notebook
python -m jupyter nbconvert \
  --execute \
  --to notebook \
  --output-dir /tmp/cell-state-reachability-notebook \
  tutorial/tutorial.ipynb
python -m pytest -q tests/test_effect_dictionary.py
```

This is the clean-checkout CI path: it executes into `/tmp` without modifying the committed
notebook. The notebook is committed with bounded outputs from a successful top-to-bottom
execution. Its setup tolerates notebook runners that use `tutorial/` as the kernel
directory while keeping machine-specific absolute paths out of saved outputs.

## Effect-dictionary boundary

The reusable adapter lives at [`../effect_dictionary.py`](../effect_dictionary.py). Its
matrix is always `(perturbations, genes)`, with unique string labels on both axes. It reads
NPZ files with `allow_pickle=False` and supports dense or SciPy sparse **cell matrices**
when constructing a dense effect dictionary.

The builder computes pooled condition mean minus pooled control mean. It does not normalize
raw counts, construct replicate-aware pseudobulks, correct donor/batch effects, select
features, or estimate uncertainty. Cell counts are provenance fields, not biological
replicate counts. Those scientific choices must be made upstream.
