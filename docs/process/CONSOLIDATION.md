# CONSOLIDATION — strategic-doc merge (record of what merged into what)

*A one-time record of the documentation consolidation that reduced the root-level strategic
markdown from 19 files to 7 canonical ones. Content was deduplicated, not summarized: every
unique claim, number, citation (7 DOIs, 1 PMID), and figure-embed marker (32 total) was
preserved verbatim and verified programmatically against the pre-merge sources.*

## The 7 canonical root docs (what survives)

| Doc | Role |
|---|---|
| `README.md` | Front door — framing, method, data, hackathon fit, where-to-start table. |
| `RESULTS.md` | Primary writeup — headline verdict, 12-cell atlas, modality triage, generalizability, design toolkit, method-improvement results. |
| `NOVELTY.md` | Novelty delta vs prior art + real-world/disease impact + field ("tide vs wave") positioning. |
| `RELATED_WORK.md` | The citation-grounded 91-prior-method survey (unchanged; already standalone). |
| `CAUSAL.md` | Design-based causal reframe + IV/compliance trust layer + A1–A6/B1–B4 agenda + validation ledger + adversarial dataset critique. |
| `ROADMAP.md` | The 3-day hackathon build plan (promoted from the former `ROADMAP_HACKATHON.md`). |
| `CLAUDE.md` | Operating manual — verified facts, guardrails, repo map, literature anchors (slimmed to current state). |

## What merged into what (13 files removed)

**→ NOVELTY.md** (novelty + impact + field positioning)
- `IMPACT.md` — the industry-scientist case (attrition economics, the two odds-moving levers, concrete users/decisions, quantified value proposition).
- `TIDE_VS_WAVE.md` — the overarching-problem framing, the tide-vs-wave call with its 2025–26 evidence, the five moves, the priority discipline.

**→ RESULTS.md** (spine + four absorbed sections)
- `MODALITY.md` — the signed LOF/GOF/neither modality triage against the druggable genome.
- `GENERALIZABILITY.md` — cross-dataset (K562 CRISPRa) and cross-cell-type (K562/RPE1) transfer, the input-contract framing, the application map.
- `DESIGN.md` — the `design_experiment()` toolkit and its calibration caveat.
- `METHOD_IMPROVEMENTS.md` — the nine method-improvement in-silico results and the methodological + experimental agendas.

**→ CAUSAL.md** (new unified trust-and-causal doc)
- `CAUSAL_INFERENCE.md` — the design-based reframe + the IV/compliance layer.
- `CAUSAL_RESEARCH_AGENDA.md` — the A1–A6 / B1–B4 agenda and the identifying-assumption stack.
- `VALIDATION_AND_EXPERIMENTS.md` — the assumption-by-assumption validation ledger and prioritized experiment plan.
- `COUNTERFACTUAL_EXPLANATION.md` — the verdict as a Wachter-2017-style counterfactual explanation.
- `REVIEWER2_DATASET_CRITIQUE.md` — the adversarial dataset appraisal (generators, donors, targets, provenance).

**→ ROADMAP.md** (promoted, replacing the superseded 7-day plan)
- `ROADMAP_HACKATHON.md` — promoted verbatim (footer cross-refs rewired to the 7-doc scheme); the old 7-day planning `ROADMAP.md` content it superseded was retired at the same time.

**Retired (superseded planning-era docs, content already in the current docs)**
- `ASSESSMENT.md` — pre-build technical review; its live facts (schema counts, literature verification, definition-of-done) now live in `CLAUDE.md` and `RESULTS.md`.

## Fidelity guarantee

The merge was executed as three parallel tracks over delimited source bundles, then verified:
- **Numbers:** every distinct numeric claim in each source cluster is present in its merged doc (verified by set-difference; the only deltas were the bundle's own delimiter char/line counts, not content).
- **Figures:** all 32 `{{artifact:…}}` / `art_…` embed markers preserved (4 in NOVELTY, 26 in RESULTS, 2 in CAUSAL).
- **Citations:** all 7 DOIs and 1 PMID preserved.
- **Links:** no surviving doc (or the manuscript `.tex`) contains a dangling markdown link to a removed file; provenance mentions of the former filenames are retained deliberately as prose so this decision trail stays visible.

The pre-merge source text is preserved in the saved artifacts `bundle_novelty.md`,
`bundle_results.md`, and `bundle_causal.md` should any original phrasing need to be recovered.


## Post-merge cross-reference sweep

A repo-wide scan (all `.md`, `.ipynb`, `.tex`, `.py`) found and repointed **19 stale
references** to merged-away files that the initial link-check (root docs only) had missed:

- `RELATED_WORK.md` (×3): companion intro, attrition-citation pointer, and see-also footer → `NOVELTY.md`.
- `data/README.md` (×1): "computed on it" pointer → `RESULTS.md`.
- `notebooks/README.md` (×2): design-toolkit → `RESULTS.md`; causal companions → `CAUSAL.md`.
- `FIVE_DAY_SUMMARY.md` (×1) → `NOVELTY.md`; `RESPONSE_TO_REVIEWER2.md` (×1) → `CAUSAL.md`.
- Notebooks 02 (×1), 04 (×2), 05 (×1), 07 (×1), 09 (×2): in-cell prose pointers → `RESULTS.md` / `CAUSAL.md`.
- Scripts `build_nbB.py` (×1), `run_iv_compliance.py` (×1),
  `a2_conditional_reachability_scaffold.py` (×1), `notebooks/_nb04_script.py` (×1):
  docstring/comment pointers → `RESULTS.md` / `CAUSAL.md`.

After the sweep, every surviving reference to a former filename is an intentional
"absorbs the former X.md" provenance note inside the merged docs and this file; all edited
notebooks remain valid JSON and all edited scripts still parse.
