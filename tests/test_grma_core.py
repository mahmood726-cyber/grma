"""Unit tests for the GRMA scientific core (grey_meta_v8.py).

These tests pin the *behaviour and invariants* of the estimator (input
validation, weight normalisation, reproducibility, fallback paths). They do NOT
assert any published manuscript number: the reference values used here were
read off the current implementation and are used only as regression anchors so
that future refactors cannot silently change the estimator.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from grey_meta_v8 import (  # noqa: E402
    GRMA,
    compare_methods,
    get_bcg_data,
    get_morris_data,
    valley_diagnostic,
)


# --------------------------------------------------------------------------- #
# Construction / parameter validation
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "kwargs",
    [
        {"zeta": 0.0},
        {"zeta": 1.5},
        {"zeta": -0.1},
        {"prec_cap": 0},
        {"prec_cap": -5},
        {"trim": 0.5},
        {"trim": -0.1},
        {"tukey_c": 0},
        {"tukey_c": -1},
    ],
)
def test_init_rejects_bad_params(kwargs):
    with pytest.raises(ValueError):
        GRMA(**kwargs)


def test_init_accepts_boundary_zeta():
    # zeta == 1 is allowed (0, 1]; zeta just above 0 allowed.
    GRMA(zeta=1.0)
    GRMA(zeta=1e-6)


# --------------------------------------------------------------------------- #
# Input validation on fit / bootstrap
# --------------------------------------------------------------------------- #
def test_fit_rejects_length_mismatch():
    g = GRMA()
    with pytest.raises(ValueError):
        g.fit([0.1, 0.2, 0.3], [0.05, 0.04])


def test_fit_rejects_empty():
    g = GRMA()
    with pytest.raises(ValueError):
        g.fit([], [])


def test_fit_rejects_nonpositive_variance():
    g = GRMA()
    with pytest.raises(ValueError):
        g.fit([0.1, 0.2, 0.3], [0.05, 0.0, 0.04])
    with pytest.raises(ValueError):
        g.fit([0.1, 0.2, 0.3], [0.05, -0.01, 0.04])


def test_fit_rejects_nonfinite():
    g = GRMA()
    with pytest.raises(ValueError):
        g.fit([0.1, np.nan, 0.3], [0.05, 0.04, 0.03])
    with pytest.raises(ValueError):
        g.fit([0.1, 0.2, 0.3], [0.05, np.inf, 0.03])


# --------------------------------------------------------------------------- #
# Core estimator invariants
# --------------------------------------------------------------------------- #
def test_weights_sum_to_one_and_estimate_in_range():
    yi, vi, _ = get_bcg_data()
    g = GRMA()
    fit = g.fit(yi, vi)
    w = np.asarray(fit["weights"])
    assert np.isclose(np.sum(w), 1.0, atol=1e-12)
    assert np.all(w >= 0)
    # Pooled estimate must lie within the observed effect range.
    assert yi.min() <= fit["estimate"] <= yi.max()
    # n_eff is bounded by k.
    assert 1.0 <= fit["n_eff"] <= fit["k"]
    assert fit["method"] == "GRMA"


def test_bcg_estimate_is_reproducible():
    # Regression anchor: value read from the current implementation.
    yi, vi, _ = get_bcg_data()
    fit = GRMA().fit(yi, vi)
    assert fit["estimate"] == pytest.approx(-0.5989525560283318, abs=1e-12)
    assert fit["w_max"] == pytest.approx(0.13598187171032014, abs=1e-12)


def test_noguard_differs_from_guard_on_bcg():
    yi, vi, _ = get_bcg_data()
    guarded = GRMA(effect_guard=True).fit(yi, vi)
    unguarded = GRMA(effect_guard=False).fit(yi, vi)
    assert guarded["method"] == "GRMA"
    assert unguarded["method"] == "GRMA_noguard"
    assert guarded["estimate"] != unguarded["estimate"]


def test_all_identical_studies_gives_equal_weights():
    # delta_max == 0 path: guard against 0/0 NaN.
    y = [0.5, 0.5, 0.5, 0.5]
    v = [0.1, 0.1, 0.1, 0.1]
    fit = GRMA().fit(y, v)
    assert fit["estimate"] == pytest.approx(0.5, abs=1e-12)
    assert np.allclose(fit["weights"], 0.25)


def test_k_less_than_3_uses_hk_fallback():
    fit = GRMA().fit([0.1, 0.2], [0.05, 0.04])
    assert fit["method"] == "HK_REML_fallback"
    assert np.isfinite(fit["estimate"])
    assert np.isfinite(fit["se"])
    assert fit["ci_lo"] < fit["ci_hi"]


def test_single_study_fallback_does_not_raise():
    # k == 1 routes through the HK/REML fallback. The point estimate is NaN
    # because SE/tau^2 are undefined for a single study; we only assert that
    # the public API returns the fallback method rather than raising. (A single
    # study is not a meaningful meta-analysis; this documents the boundary.)
    fit = GRMA().fit([0.42], [0.05])
    assert fit["method"] == "HK_REML_fallback"
    assert fit["k"] == 1


# --------------------------------------------------------------------------- #
# Bootstrap inference
# --------------------------------------------------------------------------- #
def test_bootstrap_is_seed_reproducible():
    yi, vi, _ = get_morris_data()
    g = GRMA()
    a = g.bootstrap_ci(yi, vi, B=200, seed=123)
    b = g.bootstrap_ci(yi, vi, B=200, seed=123)
    assert a["ci_lo_pct"] == b["ci_lo_pct"]
    assert a["ci_hi_pct"] == b["ci_hi_pct"]
    assert a["se"] == b["se"]


def test_bootstrap_ci_brackets_estimate():
    yi, vi, _ = get_bcg_data()
    ci = GRMA().bootstrap_ci(yi, vi, B=499, seed=7)
    assert ci["ci_lo_pct"] <= ci["estimate"] <= ci["ci_hi_pct"]
    assert 0.0 <= ci["pvalue"] <= 1.0
    assert ci["n_boot_ok"] > 0


def test_bootstrap_validates_inputs():
    with pytest.raises(ValueError):
        GRMA().bootstrap_ci([0.1, 0.2, 0.3], [0.05, 0.04])


# --------------------------------------------------------------------------- #
# Leave-one-out
# --------------------------------------------------------------------------- #
def test_leave_one_out_requires_k4():
    yi = [0.1, 0.2, 0.3]
    vi = [0.05, 0.04, 0.03]
    with pytest.raises(ValueError):
        GRMA().leave_one_out(yi, vi)


def test_leave_one_out_shapes_and_full_matches_fit():
    yi, vi, _ = get_bcg_data()
    g = GRMA()
    loo = g.leave_one_out(yi, vi)
    assert loo["est_loo"].shape == (len(yi),)
    assert np.all(np.isfinite(loo["est_loo"]))
    assert loo["est_full"] == pytest.approx(g.fit(yi, vi)["estimate"], abs=1e-12)
    assert 0 <= loo["idx_max_shift"] < len(yi)


# --------------------------------------------------------------------------- #
# Comparators + valley diagnostic
# --------------------------------------------------------------------------- #
def test_compare_methods_returns_all_methods():
    yi, vi, _ = get_bcg_data()
    res = compare_methods(yi, vi, seed=1)
    names = {r["method"] for r in res}
    assert {"IV_FE", "DL_RE", "PM_RE", "HK_REML", "HuberIV", "tRE_df4",
            "GRMA", "GRMA_noguard"} <= names


def test_valley_diagnostic_small_k_is_null():
    out = valley_diagnostic([0.1, 0.2], estimate=0.15)
    assert out["valley_flag"] is False
    assert out["valley_p"] == 1.0


def test_valley_diagnostic_reproducible():
    yi, _, _ = get_bcg_data()
    est = GRMA().fit(*get_bcg_data()[:2])["estimate"]
    a = valley_diagnostic(yi, est, n_perm=200, seed=99)
    b = valley_diagnostic(yi, est, n_perm=200, seed=99)
    assert a == b
