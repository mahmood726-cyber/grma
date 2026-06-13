"""Minimal smoke tests for the GRMA Python engine (grey_meta_v8.py).

These guard the import path and the core estimator contracts:
- the estimator fits the embedded BCG dataset and returns a finite estimate,
- input validation rejects mismatched / empty / non-finite / non-positive inputs,
- the k < 3 HK-REML fallback path is taken and returns a finite estimate.

They are deliberately light (no bootstrap B sweeps) so the default suite stays fast.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from grey_meta_v8 import GRMA, get_bcg_data, _dl_random_effects  # noqa: E402


def test_fit_bcg_returns_finite_estimate():
    yi, vi, _ = get_bcg_data()
    g = GRMA()
    fit = g.fit(yi, vi)
    assert fit["k"] == len(yi)
    assert np.isfinite(fit["estimate"])
    # weights are a valid probability vector
    assert pytest.approx(1.0, abs=1e-9) == float(np.sum(fit["weights"]))
    assert 0.0 < fit["w_max"] <= 1.0


def test_validate_inputs_rejects_bad_data():
    g = GRMA()
    with pytest.raises(ValueError):
        g.fit([0.1, 0.2], [0.01])  # length mismatch
    with pytest.raises(ValueError):
        g.fit([], [])  # empty
    with pytest.raises(ValueError):
        g.fit([0.1, np.nan, 0.3], [0.01, 0.01, 0.01])  # non-finite effect
    with pytest.raises(ValueError):
        g.fit([0.1, 0.2, 0.3], [0.01, 0.0, 0.01])  # non-positive variance


def test_small_k_uses_hk_reml_fallback():
    g = GRMA()
    fit = g.fit([0.2, 0.4], [0.02, 0.03])
    assert fit["method"] == "HK_REML_fallback"
    assert np.isfinite(fit["estimate"])
    assert fit["ci_lo"] <= fit["estimate"] <= fit["ci_hi"]


def test_dl_random_effects_basic():
    yi, vi, _ = get_bcg_data()
    res = _dl_random_effects(yi, vi)
    assert np.isfinite(res["estimate"])
    assert res["tau2"] >= 0.0
    assert res["ci_lo"] <= res["estimate"] <= res["ci_hi"]
