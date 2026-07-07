"""Streamlit hypothesis explorer (CPU-only).

Run:  streamlit run app/explorer.py

Pick a target-state axis and stimulation condition -> view ranked minimal
perturbation sets -> expand a hypothesis card showing the confidence breakdown,
pathway enrichment, literature citations, and a permanent limitations banner.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import streamlit as st

sys.path.append(str(Path(__file__).resolve().parent.parent))
from src import counterfactual as cf  # noqa: E402
from src import data_loader, target_states  # noqa: E402

LIMITATIONS = (
    "⚠️ These are **falsifiable hypotheses**, not conclusions. CRISPRi effects are "
    "loss-of-function only; multi-gene sets assume additivity (no epistasis); a matched "
    "transcriptomic shift does not prove functional rescue; the target state is a "
    "transcriptomic proxy. Validate experimentally before acting."
)


@st.cache_data(show_spinner=True)
def _load(condition: str):
    d = data_loader.load_de_dictionary(condition=condition)
    return d.effects, d.perturbation_ids, d.genes


def main() -> None:
    st.set_page_config(page_title="Counterfactual Biology Explorer", layout="wide")
    st.title("Counterfactual Biology Explorer")
    st.caption("Minimal perturbation sets that shift CD4+ T-cell state — with confidence.")
    st.warning(LIMITATIONS)

    with st.sidebar:
        axis = st.selectbox("Target state axis", ["Th1/Th2 polarization", "Aging"])
        condition = st.selectbox("Condition", ["Stim8hr", "Stim48hr", "Rest"])
        k_max = st.slider("Max set size (k)", 1, 12, 6)
        solver = st.selectbox("Solver", ["greedy", "omp", "lasso"])

    try:
        E, pert_ids, genes = _load(condition)
    except Exception as exc:
        st.info(f"Load the DE_stats artifact first (see data/README.md). Details: {exc}")
        return

    if axis.startswith("Th1"):
        d = target_states.polarization_target(genes, "toward_Th1")
    else:
        d = target_states.aging_target(genes, "toward_young")

    solver_fn = {
        "greedy": cf.greedy_minimal_set,
        "omp": cf.omp_minimal_set,
        "lasso": cf.lasso_minimal_set,
    }[solver]
    result = solver_fn(E, d) if solver != "greedy" else cf.greedy_minimal_set(E, d, k_max)

    null = cf.random_null(E, d, k=max(1, len(result.indices)))
    p_emp = float(np.mean(null >= result.cosine))

    st.subheader(f"Nominated minimal set (k = {len(result.indices)})")
    st.metric("Cosine alignment with target", f"{result.cosine:.3f}",
              help=f"Empirical p vs random-perturbation null: {p_emp:.3g}")
    for local_i, w in zip(result.indices, result.weights or []):
        gene = str(pert_ids[local_i])
        with st.expander(f"{gene}  (weight {w:+.2f})"):
            st.write("Confidence and evidence panels populate here "
                     "(src/confidence.py, src/evidence.py — pathways now live in evidence.py).")
            st.caption("Knock-down nomination. Validate with an arrayed CRISPRi assay.")


if __name__ == "__main__":
    main()
