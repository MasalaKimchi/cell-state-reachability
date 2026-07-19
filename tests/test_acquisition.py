"""Tests for the prospective active-acquisition recommender."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import acquisition as aq  # noqa: E402


def _singles(seed=0, n=6, n_genes=40, clusters=True):
    rng = np.random.default_rng(seed)
    S = rng.normal(0, 1, (n, n_genes))
    if clusters and n >= 4:
        S[1] = S[0] + 0.05 * rng.normal(0, 1, n_genes)  # B ~ A
        S[3] = S[2] + 0.05 * rng.normal(0, 1, n_genes)  # D ~ C
    names = [chr(ord("A") + i) for i in range(n)]
    return S, names


def test_batch_size_and_order():
    S, names = _singles()
    b = aq.recommend_batch(S, names, batch_size=4)
    assert len(b.combos) == 4
    assert all(len(c) == 2 for c in b.combos)
    assert b.relevance.shape == (4,)
    assert b.novelty.shape == (4,)


def test_greedy_is_top_by_relevance():
    S, names = _singles()
    b = aq.recommend_batch(S, names, batch_size=5, strategy="greedy")
    # relevance must be non-increasing for a greedy top-k
    assert np.all(np.diff(b.relevance) <= 1e-9)


def test_diversified_differs_from_greedy_with_clusters():
    S, names = _singles(clusters=True)
    g = aq.recommend_batch(S, names, batch_size=4, strategy="greedy")
    d = aq.recommend_batch(S, names, batch_size=4, strategy="diversified", diversity_weight=0.8)
    assert set(g.combos) != set(d.combos)  # diversity changed the batch


def test_diversity_weight_zero_matches_greedy():
    S, names = _singles()
    g = aq.recommend_batch(S, names, batch_size=4, strategy="greedy")
    d = aq.recommend_batch(S, names, batch_size=4, strategy="diversified", diversity_weight=0.0)
    assert list(g.combos) == list(d.combos)


def test_excludes_measured():
    S, names = _singles()
    measured = [("A", "B"), ("C", "D")]
    b = aq.recommend_batch(S, names, batch_size=5, measured=measured)
    assert ("A", "B") not in b.combos
    assert ("C", "D") not in b.combos


def test_explicit_candidates():
    S, names = _singles()
    cands = [("A", "C"), ("B", "D"), ("E", "F")]
    b = aq.recommend_batch(S, names, batch_size=2, candidates=cands)
    assert all(c in {tuple(sorted(x)) for x in cands} for c in b.combos)


def test_labeled_switches_to_ridge():
    S, names = _singles(n=6)
    labeled = {("A", "B"): 5.0, ("C", "D"): 1.0, ("E", "F"): 3.0, ("A", "C"): 2.0}
    b = aq.recommend_batch(S, names, batch_size=3, labeled=labeled)
    assert b.model == "ridge"


def test_labeled_too_few_stays_training_free():
    S, names = _singles()
    b = aq.recommend_batch(S, names, batch_size=3, labeled={("A", "B"): 1.0})
    assert b.model == "training-free"


def test_novelty_bounds():
    S, names = _singles()
    b = aq.recommend_batch(S, names, batch_size=4, strategy="diversified")
    assert np.all(b.novelty <= 1.0 + 1e-9)
    assert np.all(b.novelty >= -1e-9)


def test_gene_jaccard_metric_runs():
    S, names = _singles()
    b = aq.recommend_batch(S, names, batch_size=4, diversity_metric="gene_jaccard")
    assert len(b.combos) == 4


def test_invalid_args():
    S, names = _singles()
    with pytest.raises(ValueError):
        aq.recommend_batch(S, names, batch_size=0)
    with pytest.raises(ValueError):
        aq.recommend_batch(S, names, batch_size=2, strategy="bogus")
    with pytest.raises(ValueError):
        aq.recommend_batch(S, names, batch_size=2, diversity_weight=2.0)


def test_as_rows():
    S, names = _singles()
    b = aq.recommend_batch(S, names, batch_size=3)
    rows = b.as_rows()
    assert len(rows) == 3
    assert rows[0]["run_order"] == 1
    assert "combination" in rows[0] and "relevance" in rows[0]


# --------------------------------------------------------------------------- #
# Campaign simulation
# --------------------------------------------------------------------------- #
def test_simulate_campaign_recovers_all():
    S, names = _singles(n=8)
    # every pair gets a random outcome; campaign must eventually acquire all
    import itertools
    combos = [tuple(sorted(c)) for c in itertools.combinations(names, 2)]
    rng = np.random.default_rng(1)
    outcomes = {c: float(rng.normal()) for c in combos}
    res = aq.simulate_campaign(S, names, outcomes, batch_size=5, refit=False)
    assert len(res["order"]) == len(combos)
    assert res["discovered"][-1] == res["n_emergent"]
    # discovery curve is non-decreasing
    assert np.all(np.diff(res["discovered"]) >= 0)


def test_simulate_campaign_planted_emergence_beats_random():
    """With emergence planted in dissimilar pairs, the recommender should
    front-load discoveries versus random acquisition (raw-label regime)."""
    rng = np.random.default_rng(2)
    n_genes = 60
    n = 12
    S = rng.normal(0, 1, (n, n_genes))
    names = [f"g{i}" for i in range(n)]
    import itertools
    combos = [tuple(sorted(c)) for c in itertools.combinations(names, 2)]
    # emergence label = low cosine (dissimilar) pairs, matching the real mechanism
    def cos(a, b):
        va, vb = S[names.index(a)], S[names.index(b)]
        return float(va @ vb / (np.linalg.norm(va) * np.linalg.norm(vb)))
    outcomes = {c: -cos(*c) for c in combos}
    thr = np.quantile(list(outcomes.values()), 2 / 3)
    emergent = {c: outcomes[c] >= thr for c in combos}
    res = aq.simulate_campaign(S, names, outcomes, is_emergent=emergent,
                               batch_size=8, strategy="greedy", refit=False)
    disc = res["discovered"]
    k = 16
    # random expectation at k acquisitions
    frac = res["n_emergent"] / res["n_total"]
    rand_at_k = frac * k
    assert disc[k - 1] > rand_at_k  # front-loaded
