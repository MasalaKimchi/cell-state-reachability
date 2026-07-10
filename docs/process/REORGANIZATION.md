# Repository reorganization — structure & rationale

*One-time tidy of the repo from ~70 loose root files into grouped directories,
in preparation for public GitHub upload. Method behaviour is unchanged: the
test suite (11 tests) and `reachability._selftest()` pass identically before and
after, and every internal doc link resolves.*

## Final layout

| Path | Contents |
|---|---|
| `reachability.py` | The method (stays at root — import contract for `reproduce.sh`, tests, notebooks). |
| `reproduce.sh`, `environment.yml`, `requirements.txt` | One-command reproduction + pinned env. |
| `docs/` | Narrative writeups: `RESULTS.md` (start here), `NOVELTY.md`, `RELATED_WORK.md`, `CAUSAL.md`, `ROADMAP.md`. |
| `docs/figures/` | Doc figures (roadmap timeline). |
| `docs/process/` | Process notes: consolidation log, 5-day summary, reviewer-2 response, this file. |
| `scripts/` | Analysis drivers: `run_atlas`, `run_nulls`, `run_bootstrap`, `run_iv_compliance`, `run_a1_sensitivity`, `run_deg_weighted_eval`, `build_effect_matrices`, `build_nbB`, `a2_conditional_reachability_scaffold`, `_split_stability_worker`. |
| `notebooks/` | 01–09 analysis notebooks + `figures/` + READMEs. |
| `app/` | Interactive explorer: 8 self-contained HTML views (`index.html` is the hub) + inline data + `DEPLOY.md`; `previews/` holds PNG previews. |
| `results/` | Canonical output tables: atlas, modality triage, K562 demo, a-series (A1–A6), `references.csv`, dashboards. |
| `manuscript/` | LaTeX manuscript (`sections/`, `figures/`). |
| `data/` | Tier-1 supplementary CSVs + Tier-2 `GWCD4i.DE_stats.h5ad` (16.8 GB, gitignored). |
| `analysis_cache/` | Cached intermediates: `atlas_work/`, `nb_out/`, `czi_data/`, `czi_fig/`. |

## What moved

- **5 narrative docs** (`RESULTS`, `NOVELTY`, `RELATED_WORK`, `CAUSAL`, `ROADMAP`) root → `docs/`; **3 process docs** → `docs/process/`.
- **10 driver/​builder scripts** root → `scripts/`; each gained a repo-root anchor (`os.chdir` to the repo root computed from `__file__`) so their relative cache paths keep resolving regardless of the working directory.
- **~21 loose output files** root → `results/`.
- **8-file web app** + data + previews → `app/` (co-located so the `index.html` hub's relative links keep working).
- **`atlas_work/`, `nb_out/`, `czi_data/`, `czi_fig/`** → `analysis_cache/` (wholesale; internal structure preserved).

## What was removed (byte-identical duplicates / regenerable caches)

Verified identical to a surviving copy before removal: `ms_save2/` (== `manuscript/`), `deploy/index.html` (== `app/reachability_explorer.html`), root `50_discussion.tex` and several root `reviewer2_*` files (== `results/`/`manuscript/` copies), plus empty/scratch dirs.

## Path-coupling updates (so nothing breaks)

- **Scripts** — repo-root anchor added; cache-dir literals repointed to `analysis_cache/…`; redundant root-copy writes (`run_a1_sensitivity`, `a2_scaffold`) redirected into `results/`.
- **Notebooks 06/07/09** — `os.path.join(REPO, "atlas_work"|"nb_out")` → `…, "analysis_cache", …`; nb07's `../czi_data/` and `../czi_fig/` embeds repointed; nb06's A3 file search extended to look in `results/` (where `a3_kway_additivity_bound.csv` now lives).
- **Docs** — narrative-doc links gained the `docs/` prefix; command references gained `scripts/`; backticked cache-file provenance annotations gained `analysis_cache/`; the two truly broken README links and the `RESPONSE_TO_REVIEWER2 → ../CAUSAL.md` link fixed; README/CLAUDE repo-map trees rewritten to match the new layout.

## What is gitignored (regenerable / too large for GitHub)

Only large binaries are excluded — every small result table, per-cell JSON, figure and `.md` writeup under `analysis_cache/` and `notebooks/cache/` stays tracked:

- `data/*.h5ad|*.h5mu|*.csv` (Tier-1/Tier-2 raw data; refetch per `data/README.md`)
- `analysis_cache/**/inputs.npz` (772 MB), `…/cross_celltype_effects.npz` (17 MB), `…/first_stage_compliance.csv` (2 MB)
- `notebooks/cache/*.npz` — the three per-condition effect matrices `E_{Rest,Stim8hr,Stim48hr}.npz` (150–160 MB each, over GitHub's 100 MB limit; rebuilt by `scripts/run_deg_weighted_eval.py` and nb03/nb04)
- `notebooks/figures/`, `__pycache__/`, `.pytest_cache/`, `.venv/`, `*.log`, `.DS_Store`

After the reorg the staged tree is ~25 MB across ~261 files, with no single file above ~4 MB.
