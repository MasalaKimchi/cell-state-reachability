import numpy as np
import pytest

import reachability as rx


def test_exact_orthant_geometry():
    effects = np.eye(4)

    inside = rx.project_cone(effects, np.array([1.0, 2.0, 3.0, 4.0]))
    np.testing.assert_allclose(inside.fitted, [1, 2, 3, 4], atol=1e-12)
    assert inside.geometry_status == "inside_tolerance"
    assert inside.dual_separator is None
    assert inside.kkt_violation < 1e-10

    boundary = rx.project_cone(effects, np.array([1.0, 0.0, 0.0, 0.0]))
    np.testing.assert_allclose(boundary.fitted, [1, 0, 0, 0], atol=1e-12)
    assert boundary.geometry_status == "inside_tolerance"

    outside = rx.project_cone(effects, np.array([1.0, 0.0, -1.0, 0.0]))
    np.testing.assert_allclose(outside.fitted, [1, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(outside.residual, [0, 0, -1, 0], atol=1e-12)
    assert outside.geometry_status == "outside_model_cone"
    assert outside.separation_margin > 0


def test_weighted_separator_uses_metric_residual():
    effects = np.array([[1.0, 1.0]])
    target = np.array([2.0, 1.0])
    q = np.array([1.0, 4.0])
    result = rx.project_cone(effects, target, gene_weights=q)

    np.testing.assert_allclose(result.coefficients, [1.2], atol=1e-12)
    np.testing.assert_allclose(result.residual, [0.8, -0.2], atol=1e-12)
    np.testing.assert_allclose(result.dual_separator, [0.8, -0.8], atol=1e-12)
    assert effects @ result.dual_separator <= 1e-12
    assert abs(result.fitted @ result.dual_separator) < 1e-12
    assert target @ result.dual_separator > 0
    assert result.polarity_violation < 1e-12
    assert result.orthogonality_error < 1e-12


def test_zero_weight_coordinates_are_explicitly_excluded():
    effects = np.array([[1.0, 99.0], [0.0, -50.0]])
    target = np.array([1.0, 1_000.0])
    result = rx.project_cone(effects, target, gene_weights=np.array([1.0, 0.0]))
    np.testing.assert_allclose(result.fitted, [1.0, 0.0], atol=1e-12)
    np.testing.assert_allclose(result.residual, [0.0, 0.0], atol=1e-12)


@pytest.mark.parametrize("scale", [1e-6, 1.0, 1e6])
def test_target_scale_invariance(scale):
    effects = np.eye(3)
    target = np.array([1.0, -1.0, 0.5])
    base = rx.project_cone(effects, target)
    scaled = rx.project_cone(effects, scale * target)
    np.testing.assert_allclose(scaled.fitted / scale, base.fitted, rtol=1e-9, atol=1e-10)
    assert scaled.cosine == pytest.approx(base.cosine, abs=1e-10)
    assert scaled.residual_fraction == pytest.approx(base.residual_fraction, abs=1e-10)


def test_atom_scale_duplicate_and_zero_atom_do_not_change_fit():
    effects = np.array([[1.0, 0.0], [1.0, 1.0]])
    target = np.array([0.0, 1.0])
    base = rx.project_cone(effects, target)

    scaled = effects.copy()
    scaled[1] *= 7.0
    np.testing.assert_allclose(rx.project_cone(scaled, target).fitted, base.fitted, atol=1e-10)

    augmented = np.vstack([effects, effects[1], np.zeros(2)])
    np.testing.assert_allclose(rx.project_cone(augmented, target).fitted, base.fitted, atol=1e-10)


@pytest.mark.parametrize("maximum_scale", [1e4, 1e8, 1e12, 1e16])
def test_positive_atom_rescaling_stress_grid(maximum_scale):
    rng = np.random.default_rng(7)
    effects = rng.normal(size=(5, 8))
    target = rng.normal(size=8)
    base = rx.project_cone(effects, target)
    scales = np.geomspace(1.0, maximum_scale, effects.shape[0])
    stressed = rx.project_cone(effects * scales[:, None], target)

    np.testing.assert_allclose(stressed.fitted, base.fitted, rtol=1e-8, atol=1e-9)
    assert np.isfinite(stressed.kkt_violation)
    assert stressed.kkt_violation < 1e-8


def test_held_out_coefficients_are_frozen_before_scoring():
    effects = np.array([[1.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0]])
    target = np.array([2.0, 2.0, 0.0, 0.0])
    result = rx.held_out_alignment(effects, target, [0, 2], [1, 3])
    np.testing.assert_allclose(result.coefficients, [2.0, 0.0], atol=1e-12)
    assert result.fit_cosine == pytest.approx(1.0)
    assert result.held_out_cosine == pytest.approx(1.0)


@pytest.mark.parametrize(
    "effects,target,fit_idx,score_idx,weights",
    [
        (np.ones(4), np.ones(4), [0], [1], None),
        (np.eye(4), np.ones(3), [0], [1], None),
        (np.eye(4), np.ones(4), [0, 0], [1], None),
        (np.eye(4), np.ones(4), [0], [0], None),
        (np.eye(4), np.ones(4), [0], [4], None),
        (np.eye(4), np.ones(4), [0.9], [1], None),
        (np.eye(4), np.ones(4), [True], [1], None),
        (np.eye(4), np.ones(4), [0], [1], np.ones(3)),
        (np.eye(4), np.ones(4), [0], [1], [1.0, 0.0, 1.0, 1.0]),
    ],
)
def test_invalid_held_out_problems_fail_closed(effects, target, fit_idx, score_idx, weights):
    with pytest.raises(rx.InputError):
        rx.held_out_alignment(
            effects,
            target,
            fit_idx,
            score_idx,
            gene_weights=weights,
        )


def test_empirical_p_is_plus_one_and_counts_ties():
    assert rx.empirical_p(1.0, [0.1, 1.0, 2.0]) == pytest.approx(3 / 4)
    assert rx.empirical_p(3.0, [0.1, 1.0, 2.0]) == pytest.approx(1 / 4)


@pytest.mark.parametrize(
    "effects,target,weights",
    [
        (np.eye(2), np.zeros(2), None),
        (np.eye(2), np.array([np.nan, 1.0]), None),
        (np.eye(2), np.ones(3), None),
        (np.eye(2), np.ones(2), np.array([1.0, -1.0])),
        (np.eye(2), np.ones(2), np.zeros(2)),
    ],
)
def test_invalid_problems_fail_closed(effects, target, weights):
    with pytest.raises(rx.InputError):
        rx.project_cone(effects, target, gene_weights=weights)


@pytest.mark.parametrize("tolerance", [-1.0, np.nan, np.inf])
def test_invalid_separator_tolerance_fails_closed(tolerance):
    with pytest.raises(rx.InputError):
        rx.project_cone(np.eye(2), np.ones(2), separator_tolerance=tolerance)


@pytest.mark.parametrize("mask", [[1, 0], [np.nan, 0]])
def test_non_boolean_gene_masks_fail_closed(mask):
    with pytest.raises(rx.InputError):
        rx.project_cone(np.eye(2), np.ones(2), gene_mask=mask)


def test_legacy_decision_apis_are_removed():
    for name in (
        "signed_reachability",
        "design_experiment",
        "reachability_spectrum",
        "activation_certificate",
        "analytic_anisotropy_null",
    ):
        assert not hasattr(rx, name)
