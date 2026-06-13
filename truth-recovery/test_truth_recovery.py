"""
test_truth_recovery.py -- measured invariants for the GRMA truth-recovery
yardstick. Seeded; no hand-entered numbers. Reduced reps so it runs in ~90s
(the full grid is in harness.py / REPORT.md).

Run:  python truth-recovery/test_truth_recovery.py    (exit 0 = all pass)
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from harness import run_cell

REPS = 120
ok = True


def check(name, cond, detail):
    global ok
    print(f"{'PASS' if cond else 'FAIL'}  {name}  {detail}")
    if not cond:
        ok = False


# Outlier scenario (the central robustness claim)
r_out = run_cell(0.3, 12, 0.05, "outlier", REPS, 20260613 + 5)
# Clean scenario (the no-harm check)
r_cln = run_cell(0.3, 12, 0.05, "clean", REPS, 20260613 + 6)

check(
    "guard recovers truth under outliers (GRMA RMSE << DL RMSE)",
    r_out["GRMA"]["rmse"] < 0.6 * r_out["DL_RE"]["rmse"],
    f"(GRMA {r_out['GRMA']['rmse']} vs DL {r_out['DL_RE']['rmse']})",
)
check(
    "the GUARD is what provides robustness (guard coverage >> no-guard under outliers)",
    r_out["GRMA"]["coverage"] > r_out["GRMA_noguard"]["coverage"] + 0.15,
    f"(GRMA {r_out['GRMA']['coverage']} vs GRMA_noguard {r_out['GRMA_noguard']['coverage']})",
)
check(
    "GRMA keeps near-nominal coverage under outliers where DL/REML cannot",
    r_out["GRMA"]["coverage"] > 0.85,
    f"(GRMA coverage under outliers {r_out['GRMA']['coverage']})",
)
check(
    "clean-data cost is modest (GRMA still covers, RMSE within ~40% of DL)",
    r_cln["GRMA"]["coverage"] > 0.85 and r_cln["GRMA"]["rmse"] < 1.4 * r_cln["DL_RE"]["rmse"],
    f"(GRMA clean cov {r_cln['GRMA']['coverage']}, rmse {r_cln['GRMA']['rmse']} vs DL {r_cln['DL_RE']['rmse']})",
)

print("\nAll measured invariants hold." if ok else "\nSOME INVARIANTS FAILED.")
sys.exit(0 if ok else 1)
