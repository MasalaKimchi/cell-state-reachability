# Roadmap

The repository is intentionally small after consolidation. New work should earn its way
into the maintained surface through explicit gates.

## 1. Numerical certification

- Add an independent active-set oracle for small NNLS cases.
- Register well-conditioned, rank-deficient, duplicated, rescaled, and ill-conditioned
  fixtures.
- Compare objectives and fitted points rather than non-identifiable coefficients.
- Preserve atom-scale-invariant KKT and separator diagnostics.

Exit: exact fixtures and randomized oracle grids pass without reinstating legacy verdicts.

## 2. Labelled inputs

- Replace bare arrays at the public boundary with labelled gene axes and provenance.
- Require identifier namespace, modality, orientation, context, time, units, and hashes.
- Reject duplicate/missing genes and mixed effect scales.

Exit: deliberate reorder, namespace, unit, and orientation corruptions fail closed.

## 3. Systemic harness

- Add a single `load → align → fit → validate → report` runner.
- Emit immutable JSON, manifest, metrics, and Markdown report artifacts.
- Add fault injection for every stage and atomic failure semantics.
- Keep v0.1 single-process and data-free in pull-request CI.

Exit: a fresh clone produces a deterministic synthetic report in under five minutes.

## 4. Statistical contract

- Add module/pathway, source-study, donor, context, guide, and combination splits.
- Register structured nulls and maxT familywise correction.
- Gate false-positive rate, power, abstention, interval coverage, and approximation regret.
- Report equivalence sets unless independent units identify one candidate.

Exit: calibration cannot pass by always abstaining, and no dataset can choose its null or
baseline after outcomes are visible.

## 5. Frozen real-data ladder

1. Norman K562 CRISPRa singles-to-doubles with selection frozen before measured doubles.
2. Zhu primary-CD4 CRISPRi with log-fold-change-primary analysis and source-block target
   transfer.
3. Replogle K562/RPE1 with frozen cross-context support and coefficients.
4. Paired primary-T-cell CRISPRi/CRISPRa when available.

Exit: every suite has a dataset card, content hashes, explicit units, matched baselines,
and claim-specific approval.

## 6. Prospective biology

Before any conversion or efficacy language, preregister an experiment with an established
polarized starting state, power-derived independent donors, matched modalities, singles and
selected combinations, RNA plus protein/cytokine/fitness/durability readouts, and frozen
analysis.

This prospective program is not a software-release dependency; the corresponding
biological claims are blocked until it succeeds.
