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
    np.testing.assert_allclose(result.dual_separator, [0.2, -0.2], atol=1e-12)
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


@pytest.mark.parametrize("scale", [1e-250, 1e-154, 1e-6, 1.0, 1e6, 1e154])
def test_target_scale_invariance(scale):
    effects = np.eye(3)
    target = np.array([1.0, -1.0, 0.5])
    base = rx.project_cone(effects, target)
    scaled = rx.project_cone(effects, scale * target)
    np.testing.assert_allclose(scaled.fitted / scale, base.fitted, rtol=1e-9, atol=1e-10)
    assert scaled.cosine == pytest.approx(base.cosine, abs=1e-10)
    assert scaled.residual_fraction == pytest.approx(base.residual_fraction, abs=1e-10)
    assert np.isfinite(scaled.relative_objective)
    assert scaled.kkt_violation < 1e-8


def test_relative_objective_remains_finite_at_large_scale():
    result = rx.project_cone(np.eye(2), 1e200 * np.array([1.0, -1.0]))
    assert result.relative_objective == pytest.approx(0.25)
    assert np.isfinite(result.relative_objective)


def test_extreme_small_scale_separator_is_certified():
    result = rx.project_cone(np.eye(2), 1e-153 * np.array([1.0, -1.0]))
    assert result.geometry_status == "outside_model_cone"
    assert result.dual_separator is not None
    assert result.kkt_violation < 1e-8
    assert result.polarity_violation < 1e-8
    assert result.orthogonality_error < 1e-8
    assert result.separation_margin > 0


def test_representably_small_atom_is_not_discarded():
    result = rx.project_cone(np.array([[1e-154, 0.0]]), np.array([1.0, 1.0]))
    np.testing.assert_allclose(result.fitted, [1.0, 0.0], atol=1e-12)
    assert result.coefficients[0] == pytest.approx(1e154)


def test_large_common_gene_weight_scale_is_normalized_safely():
    result = rx.project_cone(
        np.eye(2),
        np.array([1e200, -1e200]),
        gene_weights=np.array([1e308, 1e308]),
    )
    np.testing.assert_allclose(result.fitted / 1e200, [1.0, 0.0], atol=1e-12)
    assert result.kkt_violation < 1e-8


def test_held_out_prediction_overflow_fails_closed():
    effects = np.array([[1e-200, 1e200]])
    target = np.ones(2)
    with pytest.raises(rx.InputError, match="held-out prediction is not representable"):
        rx.held_out_alignment(effects, target, [0], [1])


def test_nnls_iteration_failure_is_wrapped_as_input_error(monkeypatch):
    # Contract under test: when the NNLS solver raises RuntimeError (its
    # iteration-limit failure mode), project_cone wraps it as InputError.
    # We force that RuntimeError deterministically instead of relying on an
    # ill-conditioned matrix to exhaust the solver's internal iteration limit.
    # Whether a given matrix actually trips that limit depends on BLAS/LAPACK
    # rounding, which varies across CPU microarchitecture, Python build, and
    # scipy version, so an input-driven trigger is not portable across the CI
    # matrix (it passed on 3.10/3.12 but converged on 3.11).
    def _raise_iteration_failure(*args, **kwargs):
        raise RuntimeError("Maximum number of iterations reached.")

    monkeypatch.setattr(rx, "nnls", _raise_iteration_failure)

    rng = np.random.default_rng(22)
    effects = rng.normal(size=(7, 12))
    target = rng.normal(size=12)
    weights = np.exp(rng.uniform(-5, 5, size=12))
    with pytest.raises(rx.InputError, match="NNLS solver failed"):
        rx.project_cone(effects, target, gene_weights=weights)


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
    ):
        assert not hasattr(rx, name)


# --------------------------------------------------------------------------- #
# Closed-form anisotropy-corrected noise null
# --------------------------------------------------------------------------- #
def _mc_null_mean(effects, target, noise_sd, *, n=4000, seed=0):
    """Monte-Carlo reference for the emergence noise null mean."""
    rng = np.random.default_rng(seed)
    f0 = rx.project_cone(effects, target).fitted
    se = np.asarray(noise_sd, dtype=float)
    if se.ndim == 0:
        se = np.full(target.shape[0], float(se))
    vals = [
        rx.project_cone(effects, f0 + rng.normal(0.0, se)).residual_fraction
        for _ in range(n)
    ]
    return float(np.mean(vals)), float(np.std(vals))


def test_analytic_null_exists_and_returns_fields():
    an = rx.analytic_anisotropy_null(np.eye(4), np.array([1.0, 0.0, -1.0, 0.0]), 0.1)
    assert isinstance(an, rx.AnalyticNull)
    for field in ("null_mean", "null_sd", "p_value", "n_active", "gamma_shape"):
        assert np.isfinite(getattr(an, field))
    assert 0.0 < an.p_value <= 1.0


def test_analytic_null_matches_monte_carlo_on_stable_facet():
    """On a non-degenerate facet the closed form matches MC to a few percent."""
    rng = np.random.default_rng(0)
    effects = rng.normal(size=(6, 150))
    w = np.abs(rng.normal(size=6)) + 0.5
    f0 = w @ effects
    target = f0 + 3.0 * rng.normal(size=150)
    se = 0.15 * (1.0 + np.abs(rng.normal(size=150)))  # anisotropic
    an = rx.analytic_anisotropy_null(effects, target, se)
    mc_mean, mc_sd = _mc_null_mean(effects, target, se, n=6000)
    assert abs(an.null_mean / mc_mean - 1.0) < 0.05
    assert abs(an.null_sd / mc_sd - 1.0) < 0.20


def test_analytic_null_is_conservative_direction():
    """Analytic null >= MC null up to Monte-Carlo error (never materially inflates).

    The closed form measures distance to the active-atom *subspace*, which is
    locally inside the cone, so the analytic null is >= the true cone null. In
    general position the facet can be full-dimensional and the two coincide to
    within MC estimation error; on the biologically relevant cone-structured
    facets the analytic null is strictly and substantially larger. We assert the
    defensible property: it is never more than a fraction of a percent below MC.
    """
    rng = np.random.default_rng(1)
    ratios = []
    for _ in range(12):
        effects = rng.normal(size=(5, 80))
        w = np.abs(rng.normal(size=5)) + 0.3
        target = w @ effects + 2.0 * rng.normal(size=80)
        se = 0.1 * (1.0 + np.abs(rng.normal(size=80)))
        an = rx.analytic_anisotropy_null(effects, target, se)
        mc_mean, _ = _mc_null_mean(effects, target, se, n=2000)
        ratios.append(an.null_mean / mc_mean)
    ratios = np.array(ratios)
    # never materially below MC (deviations near equality are MC noise), and the
    # central tendency is conservative (>= MC).
    assert ratios.min() >= 0.98
    assert np.median(ratios) >= 1.0 - 1e-9


def test_analytic_null_captures_noise_anisotropy():
    """A scalar-se approximation misreads the null when leverage aligns with noise."""
    rng = np.random.default_rng(7)
    effects = 0.05 * rng.normal(size=(8, 240))
    effects[:, :40] += rng.normal(size=(8, 40)) * 2.0  # atoms load on genes 0:40
    se = 0.05 * np.ones(240)
    se[:40] = 0.8  # ...which are also the noisy genes -> high noise-leverage
    target = (np.abs(rng.normal(size=8)) + 0.5) @ effects + 1.5 * rng.normal(size=240)
    aniso = rx.analytic_anisotropy_null(effects, target, se)
    iso = rx.analytic_anisotropy_null(effects, target, np.full(240, np.sqrt((se**2).mean())))
    mc_mean, _ = _mc_null_mean(effects, target, se, n=5000)
    # the true (anisotropic) null is matched; the isotropic surrogate is not
    assert abs(aniso.null_mean / mc_mean - 1.0) < 0.05
    assert abs(iso.null_mean / mc_mean - 1.0) > 0.05


def test_analytic_null_is_deterministic():
    effects = np.eye(5)
    target = np.array([1.0, 0.5, -1.0, 0.2, 0.0])
    a = rx.analytic_anisotropy_null(effects, target, 0.1)
    b = rx.analytic_anisotropy_null(effects, target, 0.1)
    assert a.null_mean == b.null_mean and a.p_value == b.p_value


@pytest.mark.parametrize("bad", [-0.1, np.nan])
def test_analytic_null_rejects_bad_noise(bad):
    with pytest.raises(rx.InputError):
        rx.analytic_anisotropy_null(np.eye(3), np.ones(3), np.array([0.1, bad, 0.1]))


def test_analytic_null_respects_gene_weights():
    """Passing weights changes the metric the null is computed under."""
    rng = np.random.default_rng(3)
    effects = rng.normal(size=(4, 60))
    target = (np.abs(rng.normal(size=4))) @ effects + rng.normal(size=60)
    se = 0.2 * np.ones(60)
    w = np.abs(rng.normal(size=60)) + 0.1
    plain = rx.analytic_anisotropy_null(effects, target, se)
    weighted = rx.analytic_anisotropy_null(effects, target, se, gene_weights=w)
    assert plain.null_mean != weighted.null_mean


# --- analytic_certification_power (closed-form sample-size / power) -----------------

# A clear out-of-cone direction: the -1 coordinate is unreachable from the non-negative
# orthant spanned by eye(4), so there is real signal above a small measurement floor.
_CP_EFFECTS = np.eye(4)
_CP_TARGET = np.array([1.0, 0.0, -1.0, 0.0])


def test_certification_power_certifies_real_signal_and_curve_is_monotone():
    cp = rx.analytic_certification_power(_CP_EFFECTS, _CP_TARGET, 0.35)
    assert isinstance(cp, rx.CertificationPower)
    assert cp.certifiable and cp.r2_true_denoised > 0.0
    assert np.isfinite(cp.cells_multiplier_for_target_power)
    # power is monotone non-decreasing in the depth multiplier
    assert np.all(np.diff(np.array(cp.power_curve)) >= -1e-9)


def test_certification_power_mean_crossing_delivers_below_half_power():
    """Fix for the mislabeled 'median': the mean-crossing depth yields power < 0.5,
    and the calibrated target-power depth is deeper than it."""
    cp = rx.analytic_certification_power(_CP_EFFECTS, _CP_TARGET, 0.35)
    at_crossing = rx.analytic_certification_power(
        _CP_EFFECTS, _CP_TARGET, 0.35,
        curve_multipliers=np.array([cp.cells_multiplier_for_mean_crossing]),
    )
    assert at_crossing.power_curve[0] < 0.5
    assert cp.cells_multiplier_for_target_power >= cp.cells_multiplier_for_mean_crossing


def test_certification_power_below_floor_is_never_certifiable():
    """Residual at/below the noise floor -> not certifiable at any depth."""
    reachable_target = np.array([1.0, 0.002, 0.0, 0.0])  # inside the orthant
    cp = rx.analytic_certification_power(_CP_EFFECTS, reachable_target, 0.5)
    assert not cp.certifiable
    assert cp.cells_multiplier_for_target_power == float("inf")
    assert cp.cells_multiplier_for_mean_crossing == float("inf")


@pytest.mark.parametrize(
    "noise_sd, kwargs",
    [
        (0.1, {"bar": 1.0}),               # bar must exceed 1.0
        (0.1, {"target_power": 0.0}),      # power strictly in (0, 1)
        (0.1, {"target_power": 1.0}),
        (0.0, {}),                          # zero noise floor -> undefined, fail closed
    ],
)
def test_certification_power_fails_closed(noise_sd, kwargs):
    with pytest.raises(rx.InputError):
        rx.analytic_certification_power(_CP_EFFECTS, _CP_TARGET, noise_sd, **kwargs)
