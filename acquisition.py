"""Prospective active acquisition: which combinations to run next.

CombiCone's triage score ranks *every* unmeasured combination by predicted
emergence. A real screen, though, buys experiments a **batch** at a time, and the
top-B by score is usually a poor batch: the highest-scoring candidates cluster
(they share a hub gene, or point the same way in effect space), so running all B
of them re-measures the same biology. This module turns the per-combination
triage score into a **batch recommendation** that trades predicted emergence
against within-batch diversity, and supports the closed design loop
(recommend -> run -> observe -> refit -> recommend) that a prospective campaign
actually runs.

Two ingredients, both defensible and interpretable:

  * **Relevance** — the training-free triage score ``-agg_cos`` (the single
    validated feature on Norman), or the refined :class:`combicone.TriageModel`
    ridge score once a pilot has been labeled.
  * **Diversity** — a greedy max-marginal-relevance (MMR) selection. Each pick
    maximizes ``relevance - diversity_weight * max similarity to the batch so
    far``. Similarity is the cosine of the candidates' *predicted additive-effect
    directions* (default) or gene-set Jaccard — so the batch spreads across
    distinct predicted biology, not just distinct gene labels.

Honesty
-------
The relevance score is rank-only and training-free by default; no probability or
uncertainty is claimed. Diversity is a heuristic for batch efficiency, not a
guarantee. When a labeled pilot is supplied the ridge model is used, carrying the
same LOO-CV Spearman ~0.43 caveat as :func:`combicone.fit_triage_model`. This
module never invents outcomes: it recommends what to *measure*; certification of
what you find is :func:`combicone.certify_emergence`.

Public API
----------
``recommend_batch``    one batch of combinations to run next -> AcquisitionBatch
``AcquisitionBatch``   the recommendation + per-item relevance/novelty + rationale
``simulate_campaign``  retrospective design-loop simulation over labeled outcomes
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

import combicone as cc


__all__ = ["AcquisitionBatch", "recommend_batch", "simulate_campaign"]


# --------------------------------------------------------------------------- #
# Result container
# --------------------------------------------------------------------------- #
@dataclass
class AcquisitionBatch:
    """A recommended batch of combinations to run next.

    Attributes
    ----------
    combos : list[tuple[str, ...]]
        The recommended combinations, in run order (first = highest priority).
    relevance : np.ndarray
        Predicted-emergence score per recommended combo (higher = more promising).
    novelty : np.ndarray
        1 - (max similarity to the rest of the batch) per combo; higher = more
        distinct from its batch-mates. The first pick has novelty 1.0 by definition.
    strategy : str
        "diversified" (MMR) or "greedy" (top-B by relevance).
    diversity_weight : float
        The MMR trade-off used (0 = pure relevance, 1 = pure diversity).
    model : str
        "training-free" or "ridge".
    n_candidates : int
        How many unmeasured candidates were considered.
    scope : str
        Model-relative scope disclaimer, inherited from the triage engine.
    """

    combos: list[tuple[str, ...]]
    relevance: np.ndarray
    novelty: np.ndarray
    strategy: str
    diversity_weight: float
    model: str
    n_candidates: int
    scope: str = cc._SCOPE

    def as_rows(self) -> list[dict]:
        """One dict per recommended combo (for CSV/JSON emission)."""
        return [
            {
                "run_order": i + 1,
                "combination": "+".join(c),
                "relevance": float(self.relevance[i]),
                "novelty": float(self.novelty[i]),
            }
            for i, c in enumerate(self.combos)
        ]


# --------------------------------------------------------------------------- #
# Core recommender
# --------------------------------------------------------------------------- #
def _pairwise_similarity(
    combos: Sequence[tuple[str, ...]],
    singles: np.ndarray,
    index: Mapping[str, int],
    *,
    metric: str,
) -> np.ndarray:
    """Symmetric (n, n) similarity in [0, 1] between candidate combinations."""
    n = len(combos)
    if metric == "gene_jaccard":
        sets = [set(c) for c in combos]
        S = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                u = len(sets[i] | sets[j])
                val = len(sets[i] & sets[j]) / u if u else 0.0
                S[i, j] = S[j, i] = val
        return S
    # default: cosine of the predicted additive-effect direction (a + b + ...)
    preds = np.array([np.sum([singles[index[g]] for g in c], axis=0) for c in combos])
    norms = np.linalg.norm(preds, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    unit = preds / norms
    S = unit @ unit.T
    return np.clip(S, 0.0, 1.0)  # anti-correlation is not redundancy


def recommend_batch(
    singles: np.ndarray,
    single_names: Sequence[str],
    batch_size: int,
    *,
    candidates: Iterable[tuple[str, ...]] | None = None,
    measured: Iterable[tuple[str, ...]] | None = None,
    labeled: Mapping[tuple[str, ...], float] | None = None,
    order: int = 2,
    strategy: str = "diversified",
    diversity_weight: float = 0.5,
    diversity_metric: str = "effect_cosine",
    pairwise: str = "mean",
    gene_weights: np.ndarray | None = None,
    use_gap: bool = False,
) -> AcquisitionBatch:
    """Recommend the next ``batch_size`` combinations to run.

    Parameters
    ----------
    singles, single_names : as in :func:`combicone.triage_combinations`.
    batch_size : int
        Number of combinations to recommend (the experimental budget this round).
    candidates : iterable of name-tuples, optional
        Combinations eligible to be recommended. Defaults to all order-``order``
        combinations of the singles.
    measured : iterable of name-tuples, optional
        Combinations already run (excluded from the recommendation).
    labeled : mapping combo -> emergence label, optional
        Measured combinations with a known emergence score (e.g. the noise-robust
        ``z`` from certification). When given (>= 3), a ridge
        :class:`combicone.TriageModel` is fit and used as the relevance score;
        otherwise the training-free score is used.
    strategy : {"diversified", "greedy"}
        "diversified" = MMR (relevance traded against batch diversity);
        "greedy" = top-``batch_size`` by relevance alone.
    diversity_weight : float in [0, 1]
        MMR trade-off. 0 reproduces "greedy"; 1 ignores relevance after the first
        pick. Default 0.5.
    diversity_metric : {"effect_cosine", "gene_jaccard"}
        How batch redundancy is measured (see module docstring).

    Returns
    -------
    AcquisitionBatch
    """
    if batch_size < 1:
        raise ValueError("batch_size must be >= 1")
    if strategy not in {"diversified", "greedy"}:
        raise ValueError("strategy must be 'diversified' or 'greedy'")
    if not 0.0 <= diversity_weight <= 1.0:
        raise ValueError("diversity_weight must be in [0, 1]")

    singles, index = cc._resolve_singles(singles, single_names)

    # ---- candidate pool ---------------------------------------------------- #
    if candidates is None:
        candidates = [tuple(sorted(c)) for c in itertools.combinations(sorted(index), order)]
    else:
        candidates = [tuple(sorted(c)) for c in candidates]
    measured_set = {tuple(sorted(c)) for c in (measured or [])}
    pool = [c for c in candidates if c not in measured_set]
    if not pool:
        raise ValueError("no unmeasured candidates left to recommend")

    # ---- relevance score --------------------------------------------------- #
    model = None
    model_name = "training-free"
    if labeled:
        labeled = {tuple(sorted(k)): v for k, v in labeled.items()}
        if len(labeled) >= 3:
            model = cc.fit_triage_model(
                singles, single_names, labeled, pairwise=pairwise,
                gene_weights=gene_weights,
            )
            model_name = "ridge"
    tr = cc.triage_combinations(
        singles, single_names, pool, order=order, pairwise=pairwise,
        gene_weights=gene_weights, model=model, use_gap=use_gap or (model is not None),
    )
    relevance = np.asarray(tr.score, dtype=float)

    # min-max normalize relevance to [0, 1] for a stable MMR trade-off
    lo, hi = relevance.min(), relevance.max()
    rel_norm = (relevance - lo) / (hi - lo) if hi > lo else np.ones_like(relevance)

    k = min(batch_size, len(pool))

    # ---- selection --------------------------------------------------------- #
    if strategy == "greedy" or diversity_weight == 0.0:
        order_idx = list(np.argsort(-relevance, kind="stable")[:k])
        sim = _pairwise_similarity(pool, singles, index, metric=diversity_metric)
    else:
        sim = _pairwise_similarity(pool, singles, index, metric=diversity_metric)
        selected: list[int] = []
        remaining = set(range(len(pool)))
        # first pick = pure relevance
        first = int(np.argmax(rel_norm))
        selected.append(first)
        remaining.discard(first)
        while len(selected) < k and remaining:
            best_i, best_val = None, -np.inf
            for i in remaining:
                redundancy = max(sim[i, j] for j in selected)
                val = rel_norm[i] - diversity_weight * redundancy
                if val > best_val:
                    best_val, best_i = val, i
            selected.append(best_i)
            remaining.discard(best_i)
        order_idx = selected

    combos = [pool[i] for i in order_idx]
    rel_out = relevance[order_idx]
    # novelty = 1 - max similarity to the OTHER picks in the batch
    nov = []
    for pos, i in enumerate(order_idx):
        others = [order_idx[q] for q in range(len(order_idx)) if q != pos]
        max_sim = max((sim[i, j] for j in others), default=0.0)
        nov.append(1.0 - max_sim)
    return AcquisitionBatch(
        combos=combos,
        relevance=rel_out,
        novelty=np.asarray(nov, dtype=float),
        strategy=strategy,
        diversity_weight=diversity_weight,
        model=model_name,
        n_candidates=len(pool),
        scope=tr.scope,
    )


# --------------------------------------------------------------------------- #
# Retrospective design-loop simulation (for validation)
# --------------------------------------------------------------------------- #
def simulate_campaign(
    singles: np.ndarray,
    single_names: Sequence[str],
    outcomes: Mapping[tuple[str, ...], float],
    *,
    is_emergent: Mapping[tuple[str, ...], bool] | None = None,
    batch_size: int = 10,
    n_rounds: int | None = None,
    strategy: str = "diversified",
    diversity_weight: float = 0.5,
    refit: bool = True,
    seed: int = 0,
) -> dict:
    """Replay the design loop over a screen whose combo outcomes are known.

    Only singles and *already-revealed* outcomes drive each recommendation; a
    combo's outcome is revealed only after it is recommended and "run". Returns
    the cumulative discovery curve of emergent combinations, which is the honest
    retrospective test of whether the recommender front-loads emergent hits.

    Parameters
    ----------
    outcomes : mapping combo -> emergence label (e.g. noise-robust z)
        The measured combinations available to acquire, with their labels.
    is_emergent : mapping combo -> bool, optional
        Ground-truth "certified emergent" flag per combo. Defaults to top-tercile
        of ``outcomes`` if not given.
    batch_size, n_rounds : campaign budget.
    refit : bool
        If True, refit the ridge model on all revealed labels before each batch
        (closes the loop); if False, always use the training-free score.

    Returns
    -------
    dict with keys ``order`` (acquired combos in sequence), ``discovered``
    (cumulative emergent count after each acquisition), ``n_emergent``,
    ``n_total``, ``strategy``, ``batch_size``.
    """
    outcomes = {tuple(sorted(k)): float(v) for k, v in outcomes.items()}
    pool_all = list(outcomes)
    if is_emergent is None:
        thr = np.quantile(list(outcomes.values()), 2 / 3)
        is_emergent = {c: outcomes[c] >= thr for c in pool_all}
    else:
        is_emergent = {tuple(sorted(k)): bool(v) for k, v in is_emergent.items()}
    n_emergent = int(sum(is_emergent[c] for c in pool_all))

    measured: list[tuple[str, ...]] = []
    revealed: dict[tuple[str, ...], float] = {}
    acquired_order: list[tuple[str, ...]] = []
    discovered: list[int] = []
    hits = 0

    max_rounds = n_rounds if n_rounds is not None else (len(pool_all) + batch_size - 1) // batch_size
    for _ in range(max_rounds):
        remaining = [c for c in pool_all if c not in set(measured)]
        if not remaining:
            break
        batch = recommend_batch(
            singles, single_names, batch_size,
            candidates=remaining, measured=measured,
            labeled=(revealed if refit and len(revealed) >= 3 else None),
            strategy=strategy, diversity_weight=diversity_weight,
        )
        for combo in batch.combos:
            measured.append(combo)
            revealed[combo] = outcomes[combo]
            acquired_order.append(combo)
            hits += int(is_emergent[combo])
            discovered.append(hits)
    return {
        "order": acquired_order,
        "discovered": discovered,
        "n_emergent": n_emergent,
        "n_total": len(pool_all),
        "strategy": strategy,
        "batch_size": batch_size,
    }
