"""
Test suite for the cell-state-reachability method (reachability.py).

Two layers:
  1. test_packaged_selftest — runs reachability._selftest(), the module's own
     end-to-end invariant battery (38 asserts). One green check here means the
     whole method reproduces: verdicts, KKT/Farkas certification, the signed
     LOF/GOF/neither decomposition, the greedy spectrum, held-out validation,
     the analytic anisotropy null, and DEG-weighting.
  2. Independent property tests — re-derive the load-bearing mathematical
     invariants from scratch on small fixtures, so the suite still has teeth
     if the internals of _selftest ever change. These encode the claims the
     manuscript actually makes.

Requires only numpy + scipy (the core method's only dependencies) + pytest.

    pytest -q                      # run everything
    pytest -q -k selftest          # just the packaged battery
    pytest -q -k property          # just the independent property tests
"""
import os
import sys

import numpy as np
import pytest

# make reachability.py importable regardless of where pytest is invoked
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import reachability as rx  # noqa: E402


# ---------------------------------------------------------------------------
# Layer 1 — the packaged self-test (the canonical reproduction check)
# ---------------------------------------------------------------------------
def test_packaged_selftest(capsys):
    """The module's own 38-assert invariant battery must pass end-to-end."""
    rx._selftest()
    out = capsys.readouterr().out
    assert "ALL SELF-TESTS PASSED" in out


# ---------------------------------------------------------------------------
# Fixtures for the independent property tests
# ---------------------------------------------------------------------------
@pytest.fixture
def cone_fixture():
    """A small non-negative generator dictionary E (P x G) and two targets:
    one built INSIDE the non-negative cone, one deliberately outside it."""
    rng = np.random.default_rng(0)
    P, G = 12, 40
    E = rng.standard_normal((P, G))
    # an in-cone target: a strictly non-negative mix of a few generators
    w_true = np.zeros(P)
    w_true[[1, 4, 7]] = [0.8, 1.3, 0.5]
    d_in = w_true @ E
    # an out-of-cone target: point AGAINST a generator (needs negative weight)
    d_out = -1.0 * E[2] + 0.2 * rng.standard_normal(G)
    return E, d_in, d_out, w_true


# ---------------------------------------------------------------------------
# Layer 2 — independent property tests
# ---------------------------------------------------------------------------
def test_property_incone_is_fully_reachable(cone_fixture):
    """A target that is literally a non-negative mix of rows must be recovered
    with cosine ~1 and ~0 residual."""
    E, d_in, _, _ = cone_fixture
    r = rx.reachability(E, d_in)
    assert r.reachable_cosine > 0.999
    assert r.residual_norm < 1e-3


def test_property_kkt_certified_at_optimum(cone_fixture):
    """The NNLS optimum must satisfy the KKT/Farkas conditions to ~machine
    precision — this is the certificate the method's guarantee rests on."""
    E, d_in, d_out, _ = cone_fixture
    for d in (d_in, d_out):
        r = rx.reachability(E, d)
        assert r.cert_max_violation < 1e-5


def test_property_right_triangle_identity(cone_fixture):
    """The headline geometric identity: at the NNLS optimum the residual is
    orthogonal to the fit, so reachable_cosine**2 + residual_norm**2 == 1
    for ANY target. This is what makes the 2D cone picture exact."""
    E, d_in, d_out, _ = cone_fixture
    for d in (d_in, d_out):
        r = rx.reachability(E, d)
        assert abs(r.reachable_cosine**2 + r.residual_norm**2 - 1.0) < 1e-6


def test_property_outcone_fits_worse(cone_fixture):
    """A target pointing against a generator cannot be reached as well as an
    in-cone target."""
    E, d_in, d_out, _ = cone_fixture
    r_in = rx.reachability(E, d_in)
    r_out = rx.reachability(E, d_out)
    assert r_out.reachable_cosine < r_in.reachable_cosine


def test_property_nonnegative_weights(cone_fixture):
    """The whole point of the cone: fitted weights are non-negative (no
    unrealizable 'anti-perturbations')."""
    E, d_in, _, _ = cone_fixture
    r = rx.reachability(E, d_in)
    assert np.all(r.weights >= -1e-9)


def test_property_signed_decomposition_sums_to_one(cone_fixture):
    """The signed LOF/GOF/neither decomposition is a 3-way orthogonal split of
    the target and must sum to 1."""
    E, _, d_out, _ = cone_fixture
    s = rx.signed_reachability(E, d_out)
    total = s.lof_fraction + s.gof_fraction + s.neither_fraction
    assert abs(total - 1.0) < 0.02


def test_property_pure_lof_target_is_knockdown_reachable(cone_fixture):
    """A target built as a non-negative (LOF) mix must be ~fully explained by
    the knockdown direction, needing ~no activation."""
    E, d_in, _, _ = cone_fixture
    s = rx.signed_reachability(E, d_in)
    assert s.lof_fraction > 0.99
    assert s.gof_fraction < 0.01


def test_property_negated_target_needs_activation(cone_fixture):
    """Negating an in-cone target should flip it into the activation (GOF)
    half of the joint cone."""
    E, d_in, _, _ = cone_fixture
    s = rx.signed_reachability(E, -d_in)
    assert s.gof_fraction > s.lof_fraction


def test_property_design_experiment_end_to_end(cone_fixture):
    """The one-call researcher API must return a coherent design card: verdict,
    a knee within k_max, and a library that mirrors the spectrum."""
    E, d_in, _, _ = cone_fixture
    names = [f"P{i}" for i in range(E.shape[0])]
    genes = [f"g{j}" for j in range(E.shape[1])]
    dz = rx.design_experiment(E, d_in, perturbation_names=names, readout_names=genes,
                              k_max=8, top=10, n_shuffles=10, seed=0)
    assert 1 <= dz.optimal_k <= 8
    assert len(dz.library) == dz.spectrum["k"].size
    assert dz.reachable_cosine > 0.9  # in-cone target designs a strong recipe


def test_property_spectrum_is_monotone(cone_fixture):
    """Cumulative reachable cosine must be non-decreasing as generators are
    added to the greedy recipe."""
    E, d_in, _, _ = cone_fixture
    spec = rx.reachability_spectrum(E, d_in, k_max=8)
    cos = np.asarray(spec["cosine"] if isinstance(spec, dict) else spec.cosine)
    assert np.all(np.diff(cos) >= -1e-9)
