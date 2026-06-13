"""
harness.py -- Truth-recovery yardstick for GRMA (grma).

GRMA's selling point is robustness: a Tukey-bisquare "effect guard" that
down-weights aberrant studies. The honest test is a known-truth simulation: under
a known true mu, does the guard let GRMA RECOVER the truth better than DL/REML
when a few studies are OUTLIERS, WITHOUT paying a meaningful price on clean data?

Methods compared (all from the app's own grey_meta_v8, run unchanged):
  DL_RE, PM_RE, HK_REML, GRMA (guard on), GRMA_noguard (guard off).

DGP: yi ~ N(mu, vi + tau2); in the 'outlier' scenario a fraction of studies are
shifted by a large aberrant offset. We measure coverage of the true mu, RMSE and
CI width. BCa bootstrap uses B=199 here for speed (the app ships B=999); coverage
is insensitive to B at this resolution.

Truth-first: every number is produced from seeded simulation. Run:
  python truth-recovery/harness.py --reps 200
"""
import sys, os, argparse
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from grey_meta_v8 import GRMA, _dl_random_effects, _reml_random_effects, _hk_reml

BASE_SEED = 20260613
B_BOOT = 199


def gen(rng, mu, k, tau2, scenario, frac=0.2, shift=1.5, se_lo=0.1, se_hi=0.4):
    se = np.exp(rng.uniform(np.log(se_lo), np.log(se_hi), k))
    v = se ** 2
    theta = mu + np.sqrt(tau2) * rng.standard_normal(k)
    y = theta + se * rng.standard_normal(k)
    if scenario == "outlier":
        n_out = max(1, int(round(frac * k)))
        idx = rng.choice(k, n_out, replace=False)
        y[idx] += shift   # aberrant positive contamination
    return y, v


def methods(y, v, seed):
    out = {}
    for name, fn in [("DL_RE", _dl_random_effects), ("PM_RE", _reml_random_effects), ("HK_REML", _hk_reml)]:
        try:
            r = fn(y, v); out[name] = (r["estimate"], r["ci_lo"], r["ci_hi"])
        except Exception:
            out[name] = None
    for name, guard in [("GRMA", True), ("GRMA_noguard", False)]:
        try:
            ci = GRMA(effect_guard=guard).bootstrap_ci(y, v, B=B_BOOT, bca=True, seed=seed)
            lo = ci.get("ci_lo_bca"); hi = ci.get("ci_hi_bca")
            if lo is None or not np.isfinite(lo):
                lo, hi = ci["ci_lo_pct"], ci["ci_hi_pct"]
            out[name] = (ci["estimate"], lo, hi)
        except Exception:
            out[name] = None
    return out


METHODS = ["DL_RE", "PM_RE", "HK_REML", "GRMA", "GRMA_noguard"]


def run_cell(mu, k, tau2, scenario, reps, seed0):
    acc = {m: {"cov": 0, "sq": 0.0, "bias": 0.0, "w": 0.0, "n": 0} for m in METHODS}
    rng = np.random.default_rng(seed0)
    for r in range(reps):
        y, v = gen(rng, mu, k, tau2, scenario)
        res = methods(y, v, seed=seed0 + r)
        for m in METHODS:
            o = res[m]
            if o is None or not np.isfinite(o[0]) or not np.isfinite(o[1]) or not np.isfinite(o[2]):
                continue
            est, lo, hi = o
            a = acc[m]; a["n"] += 1
            a["sq"] += (est - mu) ** 2; a["bias"] += est - mu; a["w"] += hi - lo
            if lo <= mu <= hi:
                a["cov"] += 1
    res = {}
    for m in METHODS:
        a = acc[m]
        res[m] = None if a["n"] == 0 else {
            "coverage": round(a["cov"] / a["n"], 3),
            "rmse": round((a["sq"] / a["n"]) ** 0.5, 4),
            "bias": round(a["bias"] / a["n"], 4),
            "width": round(a["w"] / a["n"], 3),
        }
    return res


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--reps", type=int, default=200)
    reps = ap.parse_args().reps
    import time; t0 = time.time()
    print(f"\n# Truth-recovery yardstick -- grma (GRMA robust pooling)")
    print(f"reps={reps}/cell  mu=0.3 tau2=0.05  B_boot={B_BOOT}  seed={BASE_SEED}\n")
    cells = [("clean", 8), ("clean", 15), ("outlier", 8), ("outlier", 15)]
    print(f"{'scenario':10s}{'k':>3s}  " + "".join(f"{m+' cov/rmse':>22s}" for m in ["DL_RE", "HK_REML", "GRMA", "GRMA_noguard"]))
    for scen, k in cells:
        r = run_cell(0.3, k, 0.05, scen, reps, BASE_SEED + k + (0 if scen == "clean" else 1000))
        row = f"{scen:10s}{k:>3d}  "
        for m in ["DL_RE", "HK_REML", "GRMA", "GRMA_noguard"]:
            d = r[m]; row += f"{(str(d['coverage'])+'/'+str(d['rmse'])):>22s}" if d else f"{'NA':>22s}"
        print(row)
    print(f"\n(coverage of true mu / RMSE of the point estimate. {time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
