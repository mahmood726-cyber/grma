"""
Cross-validate Python GRMA against R output.
Tests: BCG (with guard), BCG (no guard), Morris, all-identical, LOO.
Exit code 0 = all pass, 1 = failures.
"""
import numpy as np
import sys
sys.path.insert(0, ".")
from grey_meta_v8 import GRMA

pass_count = 0
fail_count = 0

def check(name, actual, expected, tol):
    global pass_count, fail_count
    diff = abs(actual - expected)
    if diff < tol:
        pass_count += 1
        print(f"  PASS: {name} (diff={diff:.2e})")
    else:
        fail_count += 1
        print(f"  FAIL: {name} (expected={expected:.10g}, got={actual:.10g}, diff={diff:.2e}, tol={tol:.0e})")

def check_array(name, actual, expected, tol):
    global pass_count, fail_count
    diff = float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))
    if diff < tol:
        pass_count += 1
        print(f"  PASS: {name} (max_diff={diff:.2e})")
    else:
        fail_count += 1
        print(f"  FAIL: {name} (max_diff={diff:.2e}, tol={tol:.0e})")

# ============================================================================
# BCG data (shared across tests)
# ============================================================================
yi_bcg = np.array([
    -0.606135803570315, -0.42128298435055, -1.01819072135754,
    -0.551986956723663, -1.62089822359839, 1.3959407381717,
    -0.786115585818864, -1.39541182955542, 0.0120206014527326,
    -0.471746035838696, -1.40121013928348, -0.340849646366211,
    0.446634682274447
])
vi_bcg = np.array([
    0.250909090909091, 0.0494711164157345, 0.116517390655322,
    0.00604363757840157, 0.223017247572315, 0.222884491114701,
    0.00690561845590876, 0.103426577894663, 0.00396145005951926,
    0.0564328376980611, 0.0729945632778398, 0.0124119515149481,
    0.53250449419195
])

# Morris data
yi_morris = np.array([1.2, 0.99, 0.78, 0.54, 1.10, 0.85, 1.45, 0.62])
vi_morris = np.array([0.10, 0.08, 0.12, 0.15, 0.09, 0.11, 0.07, 0.14])

# ============================================================================
# Test 1: BCG with guard (R reference values)
# ============================================================================
print("Test 1: BCG with guard")
r_estimate = -0.598952556028332
r_w_max = 0.13598187171032
r_n_eff = 10.5618584094101
r_anchor_y = -0.551986956723663
r_anchor_p = 252.432817522722
r_weights = np.array([
    0.092050102584, 0.0959567634236, 0.0712584239572,
    0.13598187171, 0.0370588247625, 0.00239786009813,
    0.117484806589, 0.0501889897206, 0.103147929726,
    0.0980158253963, 0.0513028355602, 0.106608426328,
    0.0385473401446
])

g = GRMA(effect_guard=True)
fit = g.fit(yi_bcg, vi_bcg)
check("estimate", fit['estimate'], r_estimate, 1e-10)
check("w_max", fit['w_max'], r_w_max, 1e-10)
check("n_eff", fit['n_eff'], r_n_eff, 1e-6)
check("anchor_y", fit['anchor_y'], r_anchor_y, 1e-12)
check("anchor_p", fit['anchor_p'], r_anchor_p, 1e-6)
check_array("weights", fit['weights'], r_weights, 1e-8)
print()

# ============================================================================
# Test 2: BCG without guard
# ============================================================================
print("Test 2: BCG without guard (R reference values)")
# R reference: grma_meta(yi, vi, effect_guard = FALSE, n_boot=50, bca=FALSE, seed=1)
r_ng_estimate = -0.529553878551645
r_ng_w_max = 0.116019690475892
r_ng_n_eff = 12.1097572748081
g_ng = GRMA(effect_guard=False)
fit_ng = g_ng.fit(yi_bcg, vi_bcg)
check("no-guard estimate", fit_ng['estimate'], r_ng_estimate, 1e-10)
check("no-guard w_max", fit_ng['w_max'], r_ng_w_max, 1e-10)
check("no-guard n_eff", fit_ng['n_eff'], r_ng_n_eff, 1e-6)
# Without guard, study 6 (positive outlier) pulls estimate toward zero
check("no-guard pulled toward zero", 1.0 if fit_ng['estimate'] > fit['estimate'] else 0.0, 1.0, 0.5)
print()

# ============================================================================
# Test 3: Morris dataset
# ============================================================================
print("Test 3: Morris dataset (R reference values)")
# R reference: grma_meta(yi_m, vi_m, n_boot=50, bca=FALSE, seed=1)
r_morris_estimate = 0.957955145630988
r_morris_w_max = 0.191195052260139
r_morris_n_eff = 7.42304986760445
r_morris_anchor_y = 0.92
r_morris_weights = np.array([
    0.11306110913794, 0.19119505226014, 0.13023213915622,
    0.077609892687721, 0.14522539883156, 0.15191740705615,
    0.099228900099626, 0.091530100770642
])
g_morris = GRMA(effect_guard=True)
fit_m = g_morris.fit(yi_morris, vi_morris)
check("Morris estimate", fit_m['estimate'], r_morris_estimate, 1e-10)
check("Morris w_max", fit_m['w_max'], r_morris_w_max, 1e-10)
check("Morris n_eff", fit_m['n_eff'], r_morris_n_eff, 1e-6)
check("Morris anchor_y", fit_m['anchor_y'], r_morris_anchor_y, 1e-12)
check_array("Morris weights", fit_m['weights'], r_morris_weights, 1e-8)
check("Morris weights sum to 1", float(np.sum(fit_m['weights'])), 1.0, 1e-10)
print()

# ============================================================================
# Test 4: All-identical studies (regression test for NaN fix)
# ============================================================================
print("Test 4: All-identical studies (NaN regression test)")
g_id = GRMA(effect_guard=True)
fit_id = g_id.fit([1.0, 1.0, 1.0], [0.1, 0.1, 0.1])
check("identical: estimate", fit_id['estimate'], 1.0, 1e-10)
check_array("identical: equal weights", fit_id['weights'], [1/3, 1/3, 1/3], 1e-10)
check("identical: w_max", fit_id['w_max'], 1/3, 1e-10)
check("identical: n_eff", fit_id['n_eff'], 3.0, 1e-10)

# All-identical with varying variances
fit_id2 = g_id.fit([0.5, 0.5, 0.5, 0.5, 0.5], [0.01, 0.02, 0.03, 0.04, 0.05])
check("identical-yi: estimate", fit_id2['estimate'], 0.5, 1e-10)
print()

# ============================================================================
# Test 5: Leave-one-out (BCG)
# ============================================================================
print("Test 5: Leave-one-out (BCG)")
g_loo = GRMA(effect_guard=True)
loo = g_loo.leave_one_out(yi_bcg, vi_bcg)
check("LOO est_full matches fit", loo['est_full'], fit['estimate'], 1e-10)
check("LOO 13 entries", float(len(loo['est_loo'])), 13.0, 0.5)
check("LOO all finite", 1.0 if np.all(np.isfinite(loo['est_loo'])) else 0.0, 1.0, 0.5)
check("LOO max_shift > 0", 1.0 if loo['max_abs_delta_est'] > 0 else 0.0, 1.0, 0.5)
# Max-shift study: with guard, the outlier (idx 5) has near-zero weight,
# so removing a high-precision study (idx 8, v=0.004) shifts more
check("LOO max-shift study valid idx", 1.0 if 0 <= loo['idx_max_shift'] < 13 else 0.0, 1.0, 0.5)
print()

# ============================================================================
# Test 6: BCa operator consistency (strict < not <=)
# ============================================================================
print("Test 6: BCa z0 uses strict < (canonical Efron)")
import inspect
source = inspect.getsource(g.bootstrap_ci)
check("BCa uses strict <",
      1.0 if "boot_ok < est0" in source and "boot_ok <= est0" not in source else 0.0,
      1.0, 0.5)
print()

# ============================================================================
# Summary
# ============================================================================
print(f"=== {pass_count} PASSED, {fail_count} FAILED ===")
if fail_count > 0:
    print("CROSS-VALIDATION: SOME FAILURES")
    sys.exit(1)
else:
    print("CROSS-VALIDATION: ALL PASS")
