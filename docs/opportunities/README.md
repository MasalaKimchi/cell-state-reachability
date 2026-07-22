# Opportunity execution — field-grounded extensions

This folder holds the deliverables from executing the nine opportunity initiatives in
`FIELD_OPPORTUNITIES.md` (the external strategy note). Each initiative was run as an
independent track; every result is labelled **COMPUTED** (from real repo data),
**SYNTHETIC** (from the planted `synth_triple_screen.npz` substrate), or **DESIGN-ONLY**
(a protocol/harness, not benchtop results). Full per-initiative status, numbers, and
caveats are in `EXECUTION_REPORT.md` at the repo root.

## Notes

| Note | Initiative | Status | One line |
|------|-----------|--------|----------|
| [drug_combination_triage.md](drug_combination_triage.md) | Drug combinations | done (reconciled) | NCI-ALMANAC single-agent geometry ports; training-free triage does **not** transfer to drug synergy (an honest negative). Not certified emergence on measured combos. |
| [conformal_certificate.md](conformal_certificate.md) | Conformal / selective prediction | done | Split-conformal layer giving the certificate a distribution-free false-certification rate ≤ α; coverage verified. |
| [causal_formalization.md](causal_formalization.md) | Design-based causal / uplift | done | Atoms as ATEs; emergence as a cone-membership counterfactual with explicit assumptions; uplift/HTE bridge + worked Norman example. |
| [order_aware_certification.md](order_aware_certification.md) | Order-aware theory + order-3 harness | done (theory) / synthetic (data) | Certified set is provably monotone as the cone is enriched (80→36 on synthetic triples); reusable order-3 harness; measured-triple wet spec. |
| [screen_design_decision_layer.md](screen_design_decision_layer.md) | Screen-design decision layer | done | Held-out-gene recovery reproduced from real Norman data (median rank 1); end-to-end CLI walkthrough. |
| [tcell_leave_one_donor_out.md](tcell_leave_one_donor_out.md) | CD4+ T-cell LODO | done (real data) | Public genome-scale CD4+ screen streamed; genuine leakage-free leave-one-donor-out (leakage contrast 0.62→0.00). |
| [prospective_validation_protocol.md](prospective_validation_protocol.md) | Prospective forward validation | design-only | Frozen pre-registered ranked list + power analysis (30 pairs/arm = 80% power for the 2.4× effect). |

## Code

- `../../conformal_certificate.py` — library module (repo root; imports `combicone`), the split-conformal / selective-prediction layer.
- [`../../scripts/opportunities/`](../../scripts/opportunities/) — runnable harnesses and reproduction scripts:
  `order3_harness.py`, `leave_one_donor_out.py`, `reproduce_screen_design_recovery.py`,
  `causal_formalization_example.py`, `extract_tcell_lodo_substrate.py`.

## Results

- [`../../results/opportunities/`](../../results/opportunities/) — machine-readable numbers for every track, plus
  `screen_design_walkthrough/` (CLI triage/certify outputs). The T-cell LODO substrate
  (`tcell_lodo_substrate.npz`, ~18 MB) is **not committed** — regenerate it with
  `scripts/opportunities/extract_tcell_lodo_substrate.py`.

Figures live in [`figures/`](figures/).
