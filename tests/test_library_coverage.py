"""Exact-geometry and contract tests for the library coverage layer.

Mirrors the style of ``tests/test_reachability.py``: small hand-checkable
scenarios with known answers, plus the falsifiable contract that the cheap
certificate ranking recovers the expensive realized ranking.
"""

import hashlib
import json

import numpy as np
import pytest

import library_coverage as lc
from scripts import run_library_coverage_crossdataset as crossdataset


def test_coverage_fully_inside_cone():
    # Library is the standard basis; every non-negative target is reachable.
    effects = np.eye(4)
    targets = np.array([[1.0, 2.0, 3.0, 4.0], [0.0, 1.0, 0.0, 2.0]])
    cov = lc.coverage_report(effects, targets)
    assert cov.strict_inside_cone_fraction == 1.0
    assert cov.soft_coverage_fraction is None
    assert np.all(cov.strict_inside_cone)
    assert cov.soft_covered is None
    assert all(s == "inside_tolerance" for s in cov.statuses)
    np.testing.assert_allclose(cov.mean_cosine, 1.0, atol=1e-9)
    np.testing.assert_allclose(cov.mean_residual_fraction, 0.0, atol=1e-8)


def test_coverage_detects_unreachable_direction():
    # Library cannot produce anything negative in coord 2; a target that needs
    # it is outside the cone with a positive residual there.
    effects = np.eye(3)
    targets = np.array([[1.0, 0.0, -1.0]])   # needs -1 in a coord the cone can't reach
    cov = lc.coverage_report(effects, targets)
    assert cov.strict_inside_cone_fraction == 0.0
    assert cov.statuses[0] == "outside_model_cone"
    assert cov.residual_fractions[0] > 0.0


def test_reach_cosine_threshold_is_soft_bar():
    effects = np.eye(3)
    targets = np.array([[1.0, 0.0, -1.0]])   # best cosine ~0.707
    strict = lc.coverage_report(effects, targets)                     # geometric
    soft = lc.coverage_report(
        effects, targets, soft_cosine_threshold=0.5
    )
    assert strict.strict_inside_cone_fraction == 0.0
    assert soft.strict_inside_cone_fraction == 0.0
    assert soft.soft_coverage_fraction == 1.0  # 0.707 clears the 0.5 bar
    assert not soft.strict_inside_cone[0]
    assert soft.soft_covered[0]
    # Compatibility aliases preserve the former active-rule behavior without
    # erasing the explicit strict fields.
    assert soft.reachable_fraction == 1.0
    assert soft.reachable[0]


def test_reach_cosine_alias_and_threshold_validation():
    effects = np.eye(2)
    targets = np.array([[1.0, -1.0]])
    aliased = lc.coverage_report(effects, targets, reach_cosine=0.5)
    assert aliased.soft_cosine_threshold == 0.5
    assert aliased.soft_coverage_fraction == 1.0
    with pytest.raises(ValueError, match="not both"):
        lc.coverage_report(
            effects,
            targets,
            reach_cosine=0.5,
            soft_cosine_threshold=0.5,
        )
    for invalid in (-1.01, 1.01, np.nan, np.inf):
        with pytest.raises(ValueError, match=r"\[-1, 1\]"):
            lc.coverage_report(
                effects, targets, soft_cosine_threshold=invalid
            )


def test_redundancy_flags_duplicate_atom():
    # Atom 2 is an exact duplicate of atom 1; removing it must cost nothing.
    effects = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],   # duplicate of row 1
    ])
    targets = np.array([[1.0, 1.0, 0.0], [0.0, 2.0, 0.0]])
    red = lc.atom_redundancy(effects, targets)
    assert red.marginal_mean_cosine_loss[2] == pytest.approx(0.0, abs=1e-9)
    assert red.marginal_strict_inside_cone_fraction_loss[2] == pytest.approx(
        0.0, abs=1e-9
    )
    assert red.marginal_soft_coverage_fraction_loss is None
    np.testing.assert_array_equal(red.marginal_cosine_loss, red.marginal_mean_cosine_loss)
    np.testing.assert_array_equal(
        red.marginal_reach_loss, red.marginal_strict_inside_cone_fraction_loss
    )


def test_gap_directions_point_where_library_cannot_reach():
    effects = np.eye(3)                       # cannot make coord-2-negative
    targets = np.array([[1.0, 0.0, -1.0]])
    dirs, wts = lc.gap_directions(effects, targets)
    assert dirs.shape[0] == 1                 # one outside target
    # residual lives in coord 2 (the unreachable direction), sign negative
    assert abs(dirs[0, 2]) > 0.9
    assert wts[0] > 0.0


def test_gap_directions_use_weighted_dual_separator():
    # In a weighted metric the raw residual is not the separating normal.  The
    # acquisition layer must use the exact dual separator certified by the core.
    effects = np.array([[1.0, 1.0]])
    targets = np.array([[1.0, -1.0]])
    weights = np.array([1.0, 0.01])
    projection = lc.coverage_report(effects, targets, gene_weights=weights).results[0]
    dirs, _ = lc.gap_directions(effects, targets, gene_weights=weights)
    expected = projection.dual_separator / np.linalg.norm(projection.dual_separator)
    np.testing.assert_allclose(dirs[0], expected, atol=1e-12)


def test_certificate_prefers_gap_filling_candidate():
    # Library spans coords 0,1; a catalog target demands coord 2. The candidate
    # that supplies coord 2 must outscore one that only re-spans coords 0,1.
    effects = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    targets = np.array([[1.0, 1.0, 2.0]])                 # needs coord 2
    filler = np.array([0.0, 0.0, 1.0])                    # supplies coord 2
    redundant = np.array([1.0, 1.0, 0.0])                 # adds no new direction
    candidates = np.vstack([redundant, filler])
    rank = lc.rank_acquisitions(effects, targets, candidates)
    assert rank.order[0] == 1                              # filler wins
    assert rank.certificate_score[1] > rank.certificate_score[0]


def test_certificate_recovers_realized_ranking():
    # The falsifiable contract: on a random gap-filling scenario the cheap
    # certificate ranking agrees with the expensive realized ranking, and picks
    # the same best candidate.
    from scipy.stats import spearmanr
    rng = np.random.default_rng(0)
    D = 30
    span = np.arange(20)
    gap = np.arange(20, D)

    def in_span(n):
        values = np.zeros((n, D))
        values[:, span] = np.abs(rng.standard_normal((n, len(span))))
        return values

    effects = in_span(8) + 1e-9
    targets = in_span(40)
    for i in range(0, 40, 2):
        targets[i, rng.choice(gap)] += 2.0          # half demand a gap dim
    candidates = in_span(15) * 0.3
    for j in range(0, 15, 4):
        candidates[j] = 0.0
        candidates[j, rng.choice(gap)] = 1.0
    rank = lc.rank_acquisitions(effects, targets, candidates, realized=True)
    rho, _ = spearmanr(rank.certificate_score, rank.realized_mean_cosine_gain)
    assert rho > 0.6
    assert int(rank.certificate_order[0]) == int(
        np.argmax(rank.realized_mean_cosine_gain)
    )
    assert int(rank.realized_mean_cosine_order[0]) == int(
        np.argmax(rank.realized_mean_cosine_gain)
    )
    np.testing.assert_array_equal(rank.order, rank.certificate_order)
    np.testing.assert_array_equal(rank.realized_gain, rank.realized_mean_cosine_gain)
    np.testing.assert_array_equal(rank.realized_order, rank.realized_mean_cosine_order)


def test_certificate_order_stays_distinct_from_realized_order():
    # Regression for a former evaluation tautology: ``order`` used to switch to
    # realized-gain order whenever the comparator was requested.  Seed 7 gives
    # a stable case where the cheap and expensive top picks genuinely differ.
    rng = np.random.default_rng(7)
    effects = rng.normal(size=(3, 5))
    targets = rng.normal(size=(5, 5))
    candidates = rng.normal(size=(4, 5))
    rank = lc.rank_acquisitions(effects, targets, candidates, realized=True)
    assert int(np.argmax(rank.certificate_score)) == 1
    assert int(np.argmax(rank.realized_mean_cosine_gain)) == 0
    assert int(rank.certificate_order[0]) == 1
    assert int(rank.realized_mean_cosine_order[0]) == 0


def test_realized_mean_cosine_and_soft_coverage_gains_are_distinct():
    # Candidate 0 improves both targets a little but crosses no soft bar;
    # candidate 1 pushes one target across the bar. The API must not call the
    # primary mean-cosine objective a threshold-coverage gain.
    effects = np.array([[1.0, 0.0]])
    targets = np.array([[1.0, 2.0], [1.0, -2.0]])
    candidates = np.array([[0.0, 0.4], [0.0, 2.0]])
    rank = lc.rank_acquisitions(
        effects,
        targets,
        candidates,
        realized=True,
        soft_cosine_threshold=0.9,
    )
    assert rank.realized_mean_cosine_gain is not None
    assert rank.realized_strict_inside_cone_fraction_gain is not None
    assert rank.realized_soft_coverage_fraction_gain is not None
    assert rank.soft_cosine_threshold == 0.9
    # The arrays are distinct named objectives, even if this tiny example has
    # ties under one of them.
    assert not np.array_equal(
        rank.realized_mean_cosine_gain,
        rank.realized_soft_coverage_fraction_gain,
    )


def test_input_validation():
    effects = np.eye(3)
    with pytest.raises(ValueError):
        lc.coverage_report(effects, np.array([1.0, 0.0, 0.0]))   # 1D targets
    with pytest.raises(ValueError):
        lc.rank_acquisitions(effects, np.eye(3), np.array([1.0, 0.0, 0.0]))  # 1D candidates
    with pytest.raises(ValueError):
        lc.coverage_report(effects, np.empty((0, 3)))
    with pytest.raises(ValueError):
        lc.atom_redundancy(effects[:1], np.eye(3))
    with pytest.raises(ValueError):
        lc.rank_acquisitions(effects, np.eye(3), np.empty((0, 3)))


def test_crossdataset_cache_identity_fails_closed(tmp_path):
    cache = tmp_path / "tiny.npz"
    np.savez(cache, E=np.eye(2), perts=np.array(["A", "B"]))
    payload = cache.read_bytes()
    spec = {
        "path": "slice/tiny.npz",
        "bytes": len(payload),
        "sha256": "0" * 64,
    }
    with pytest.raises(ValueError, match="SHA-256 differs"):
        crossdataset.verify_input("tiny", spec, cache_dir=tmp_path)
    spec["sha256"] = hashlib.sha256(payload).hexdigest()
    path, identity = crossdataset.verify_input("tiny", spec, cache_dir=tmp_path)
    assert path == cache
    assert identity["hash_verified"] is True
    assert identity["sha256_actual"] == spec["sha256"]


def test_crossdataset_feature_selection_cannot_see_nonlibrary_rows():
    rng = np.random.default_rng(91)
    effects = rng.normal(size=(18, 9))
    source = np.arange(len(effects))
    library, catalog, candidates, unused = crossdataset._split_design(
        source, seed=7, library_size=5, catalog_size=5, candidate_size=4
    )
    before = crossdataset._library_only_features(effects, library, 4)

    mutated = effects.copy()
    hidden = np.concatenate([catalog, candidates, unused])
    mutated[hidden] = rng.normal(loc=10_000.0, scale=1_000.0, size=mutated[hidden].shape)
    after = crossdataset._library_only_features(mutated, library, 4)

    np.testing.assert_array_equal(before, after)
    assert crossdataset._index_sha256(before) == crossdataset._index_sha256(after)


def test_crossdataset_library_rows_do_control_feature_selection():
    effects = np.zeros((6, 4), dtype=float)
    library = np.array([0, 1, 2])
    effects[library, 0] = [0.0, 1.0, 2.0]
    effects[library, 1] = [0.0, 0.1, 0.2]
    assert crossdataset._library_only_features(effects, library, 1).tolist() == [0]
    effects[library, 1] = [0.0, 10.0, 20.0]
    assert crossdataset._library_only_features(effects, library, 1).tolist() == [1]


def test_crossdataset_signed_span_is_only_a_capacity_ceiling():
    library = np.array([[1.0, 0.0], [0.0, 1.0]])
    target = np.array([[-1.0, 1.0]])
    signed = crossdataset._signed_span_cosines(target, library)
    cone = lc.coverage_report(library, target)
    assert signed[0] == pytest.approx(1.0)
    assert cone.mean_cosine < 1.0
    assert cone.strict_inside_cone_fraction == 0.0


def test_crossdataset_capacity_dominance_fails_closed():
    crossdataset._assert_capacity_dominance(
        np.array([0.8]), np.array([0.7]), np.array([0.6]), np.array([0.9])
    )
    with pytest.raises(lc.InputError, match="signed-span"):
        crossdataset._assert_capacity_dominance(
            np.array([0.8]), np.array([0.7]), np.array([0.6]), np.array([0.5])
        )


def test_crossdataset_soft_threshold_sensitivity_is_monotone():
    fractions = crossdataset._soft_threshold_fractions(
        np.array([0.2, 0.5, 0.8]), [0.3, 0.5, 0.7, 0.9]
    )
    observed = list(fractions.values())
    assert observed == sorted(observed, reverse=True)
    assert observed == pytest.approx([2 / 3, 2 / 3, 1 / 3, 0.0])


def test_norman_acquisition_freezes_one_measured_row_per_gene_before_split():
    effects = np.arange(20, dtype=float).reshape(4, 5)
    labels = np.array(["ctrl+A", "A+ctrl", "ctrl+B", "C+ctrl"])
    sources = np.array([10, 11, 12, 13])
    representatives, chosen_labels, chosen_sources, meta = (
        crossdataset._norman_canonical_representatives(effects, labels, sources)
    )
    assert chosen_labels.tolist() == ["A+ctrl", "ctrl+B", "C+ctrl"]
    assert chosen_sources.tolist() == [11, 12, 13]
    np.testing.assert_array_equal(representatives[0], effects[1])
    assert meta["canonical_genes"] == 3
    assert meta["genes_with_both_cassette_positions"] == 1
    assert meta["opposite_position_fallbacks"] == 1

    permutation = np.random.default_rng(3).permutation(len(chosen_labels))
    library = {
        crossdataset._canonical_candidate_label(label)
        for label in chosen_labels[permutation[:2]]
    }
    candidates = {
        crossdataset._canonical_candidate_label(label)
        for label in chosen_labels[permutation[2:]]
    }
    assert library.isdisjoint(candidates)


def test_norman_sensitivity_reuses_representatives_and_one_fixed_catalog():
    labels = np.array(
        [
            "A+ctrl",
            "ctrl+A",
            "B+ctrl",
            "C+ctrl",
            "D+ctrl",
            "A+B",
            "C+D",
            "A+C",
        ]
    )
    effects = np.array(
        [
            [1.0, 0.0, 0.0],
            [100.0, 100.0, 100.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
        ]
    )
    report = crossdataset.norman_split_sensitivity(
        {"E": effects, "labels": labels},
        {
            "acquisition_library_size": 2,
            "sensitivity_catalog_doubles": 2,
            "sensitivity_catalog_seed": 91,
            "n_features": 3,
        },
        seeds=[3, 5],
        thresholds=[0.5],
    )
    protocol = report["protocol"]
    assert protocol["library_canonical_genes"] == 2
    assert protocol["canonical_representative_selection"][
        "genes_with_both_cassette_positions"
    ] == 1
    assert protocol["catalog_role"] == "fixed across every library-partition seed"
    assert len(
        {row["catalog_double_source_row_index_sha256"] for row in report["splits"]}
    ) == 1


def test_crossdataset_config_freezes_semantics_and_provenance():
    config = json.loads(crossdataset.DEFAULT_CONFIG.read_text(encoding="utf-8"))
    crossdataset._validate_config(config)
    assert config["analysis"]["soft_cosine_threshold"] == 0.5
    assert len(config["analysis"]["split_sensitivity"]["seeds"]) == 12
    assert "strict_inside_cone_definition" in config["analysis"]
    assert config["analysis"]["primary_realized_acquisition_objective"].startswith(
        "gain in mean catalog cosine"
    )
    assert "already-measured" in config["provenance"]["candidate_status"]
    discrepancy = config["provenance"]["norman"]["upstream_metadata_discrepancy"]
    assert config["provenance"]["norman"]["canonical_cell_system"] == "K562"
    assert discrepancy == {
        "field": "obs.cell_type",
        "redistributed_value": "A549",
        "canonical_value": "K562",
        "status": "FLAGGED",
        "handling": discrepancy["handling"],
    }
    for spec in config["inputs"].values():
        assert spec["bytes"] > 0
        assert len(spec["sha256"]) == 64


def test_frozen_crossdataset_report_exposes_fail_closed_contract():
    report = json.loads(crossdataset.DEFAULT_REPORT.read_text(encoding="utf-8"))
    config = json.loads(crossdataset.DEFAULT_CONFIG.read_text(encoding="utf-8"))
    assert report["status"] == "PASS"
    assert report["config_sha256"] == crossdataset.sha256_file(
        crossdataset.DEFAULT_CONFIG
    )
    assert "not coverage of an independent biological target catalog" in report[
        "claim_ceiling"
    ]
    assert report["analysis_contract"]["status_definition"].startswith(
        "PASS means every configured cache identity"
    )
    for name, dataset in report["datasets"].items():
        coverage = dataset["coverage"]
        assert "strict_inside_cone_fraction" in coverage
        assert "soft_coverage_fraction" in coverage
        assert coverage["soft_cosine_threshold"] == 0.5
        acquisition = dataset["acquisition"]
        assert acquisition["candidate_status"] == (
            "supplied_already_measured_effect_atoms"
        )
        assert "certificate_matches_realized_mean_cosine_top1" in acquisition
        assert acquisition["candidate_pool_effect_atom_labels_unique"] is True
        assert acquisition["candidate_pool_canonical_labels_unique"] is True
        assert (
            "certificate_matches_realized_mean_cosine_top1_effect_atom"
            not in acquisition
        )
        assert (
            "certificate_matches_realized_mean_cosine_top1_canonical_perturbation"
            not in acquisition
        )
        assert "top1_exact" not in acquisition
    norman = report["datasets"]["norman_k562_crispra"]
    assert norman["classification"] == (
        "retrospective_single_to_double_additivity_and_coverage"
    )
    assert {
        "paired_constituent_sum",
        "paired_constituent_two_atom_cone",
        "common_single_response_nonnegative_ray",
    } <= set(norman["comparators"])
    assert report["provenance"]["norman"]["canonical_cell_system"] == "K562"
    assert norman["acquisition"][
        "certificate_matches_realized_mean_cosine_top1"
    ] == (
        norman["acquisition"]["certificate_top1"][
            "source_effect_row_index_zero_based"
        ]
        == norman["acquisition"]["realized_mean_cosine_top1"][
            "source_effect_row_index_zero_based"
        ]
    )
    assert report["datasets"]["zhu_crispri_tcell"]["design"]["seed"] == 11
    assert norman["acquisition_design"]["seed"] == 3
    assert report["datasets"]["replogle_k562_essential_crispri"]["design"][
        "seed"
    ] == 5
    for name, identity in report["input_verification"].items():
        assert identity["hash_verified"] is True
        assert identity["sha256_actual"] == config["inputs"][name]["sha256"]
