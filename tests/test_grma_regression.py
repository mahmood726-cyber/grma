"""Regression tests for grey_meta_v8.py.

Covers the previously-untested high-risk numerical paths and locks two fixed
correctness bugs:

- F2: the Hartung-Knapp SE floor mixed units (compared the dimensionless KNHA
  scaling factor q_hk against the variance se_wald**2), diverging from
  metafor::rma(test="knha") by up to ~50x for imprecise studies.
- F4: k == 1 divided by k-1 == 0 (and DL's C == 0), silently returning NaN
  SE/CI instead of the single-study Wald interval.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import grey_meta_v8  # noqa: E402
from grey_meta_v8 import (  # noqa: E402
    GRMA,
    compare_methods,
    get_bcg_data,
    simulate_scenario,
    valley_diagnostic,
    _hk_reml,
    _reml_random_effects,
)


# --------------------------------------------------------------------------
# F2: Hartung-Knapp SE must match metafor KNHA (no unit-mixed floor)
# --------------------------------------------------------------------------
def test_hk_se_matches_metafor_knha_no_unit_mixed_floor():
    # Imprecise, near-homogeneous studies: this is exactly where the old
    # `max(q_hk, se_wald**2)` floor kicked in and inflated the SE.
    y = np.array([0.30, 0.32])
    v = np.array([1.0, 1.0])
    hk = _hk_reml(y, v)

    # Independent metafor-KNHA closed form using the module's own tau2:
    reml = _reml_random_effects(y, v)
    tau2 = reml["tau2"]
    w = 1.0 / (v + tau2)
    est = np.sum(w * y) / np.sum(w)
    q_hk = np.sum(w * (y - est) ** 2) / (len(y) - 1)
    expected_se = float(np.sqrt(q_hk / np.sum(w)))

    # Before the fix this returned ~0.5 (the variance-scaled se_wald floor),
    # a ~50x divergence from the KNHA value (~0.01).
    assert hk["se"] == pytest.approx(expected_se, rel=1e-9)
    se_wald = float(np.sqrt(1.0 / np.sum(w)))
    assert hk["se"] < se_wald  # KNHA narrows here; the old floor forced >= se_wald


def test_hk_se_finite_and_brackets_estimate():
    hk = _hk_reml([0.2, 0.4], [0.02, 0.03])
    assert np.isfinite(hk["se"]) and hk["se"] > 0
    assert hk["ci_lo"] <= hk["estimate"] <= hk["ci_hi"]


# --------------------------------------------------------------------------
# F4: k == 1 must return a finite single-study Wald interval, not NaN
# --------------------------------------------------------------------------
def test_k1_fit_returns_finite_single_study_interval():
    g = GRMA()
    r = g.fit([0.5], [0.02])
    assert r["k"] == 1
    assert np.isfinite(r["se"])
    assert np.isfinite(r["ci_lo"]) and np.isfinite(r["ci_hi"])
    assert r["ci_lo"] <= r["estimate"] <= r["ci_hi"]
    # single-study Wald: se == sqrt(v)
    assert r["se"] == pytest.approx(np.sqrt(0.02))
    assert r["estimate"] == pytest.approx(0.5)


def test_k1_bootstrap_ci_returns_finite_se():
    g = GRMA()
    ci = g.bootstrap_ci([0.5], [0.02])
    assert np.isfinite(ci["se"])
    assert np.isfinite(ci["ci_lo_pct"]) and np.isfinite(ci["ci_hi_pct"])
    assert ci["ci_lo_pct"] <= ci["estimate"] <= ci["ci_hi_pct"]


# --------------------------------------------------------------------------
# F3: coverage for the untested high-risk paths
# --------------------------------------------------------------------------
def test_bootstrap_ci_percentile_and_bca_finite():
    yi, vi, _ = get_bcg_data()
    g = GRMA()
    ci = g.bootstrap_ci(yi, vi, B=299, bca=True, seed=12345)
    assert ci["n_boot_ok"] >= 50  # BCa branch actually exercised
    for key in ("estimate", "se", "ci_lo_pct", "ci_hi_pct", "ci_lo_bca", "ci_hi_bca", "pvalue"):
        assert np.isfinite(ci[key]), key
    assert ci["ci_lo_pct"] <= ci["estimate"] <= ci["ci_hi_pct"]
    assert ci["ci_lo_bca"] <= ci["ci_hi_bca"]
    assert 0.0 <= ci["pvalue"] <= 1.0


def test_leave_one_out_shapes_and_finite():
    yi, vi, _ = get_bcg_data()
    g = GRMA()
    loo = g.leave_one_out(yi, vi)
    k = len(yi)
    assert loo["est_loo"].shape == (k,)
    assert loo["w_max_loo"].shape == (k,)
    assert np.all(np.isfinite(loo["est_loo"]))
    assert 0 <= loo["idx_max_shift"] < k
    assert np.isfinite(loo["max_abs_delta_est"])


def test_leave_one_out_requires_k_at_least_4():
    g = GRMA()
    with pytest.raises(ValueError):
        g.leave_one_out([0.1, 0.2, 0.3], [0.01, 0.01, 0.01])


def test_valley_diagnostic_returns_probability():
    yi, vi, _ = get_bcg_data()
    g = GRMA()
    est = g.fit(yi, vi)["estimate"]
    vd = valley_diagnostic(yi, est, n_perm=199, seed=7)
    assert 0.0 < vd["valley_p"] <= 1.0
    assert isinstance(bool(vd["valley_flag"]), bool)


def test_compare_methods_returns_all_methods_finite():
    yi, vi, _ = get_bcg_data()
    results = compare_methods(yi, vi, seed=99)
    names = {r["method"] for r in results}
    assert {"IV_FE", "DL_RE", "PM_RE", "HK_REML", "HuberIV", "tRE_df4",
            "GRMA", "GRMA_noguard"} <= names
    for r in results:
        if "error" not in r:
            assert np.isfinite(r["estimate"]), r["method"]


def test_simulate_scenario_smoke():
    res = simulate_scenario(k=8, true_effect=0.5, tau2=0.05,
                            scenario_type="ideal", n_rep=20, seed=1)
    assert len(res) > 0
    for r in res:
        assert np.isfinite(r["bias"])
        assert np.isfinite(r["rmse"]) and r["rmse"] >= 0
        assert r["n_ok"] > 0
