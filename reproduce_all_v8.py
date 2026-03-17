#!/usr/bin/env python3
"""
GRMA v8 — Reproduce All Tables
================================

Single command: python reproduce_all_v8.py

Produces:
  Table_bias_rmse_v8_py.csv       — Simulation bias/RMSE (Python capsule)
  Table_applied_v8_py.csv         — Applied results (BCG + Morris)
  Table_diagnostics_v8_py.csv     — GRMA diagnostics
  Table_loo_v8_py.csv             — Leave-one-out influence
  Table_zeta_sensitivity_bcg_v8_py.csv — Zeta sensitivity (BCG)
  Table_defaults_v8.csv           — Default parameter table
  Table_artifact_map_v8.csv       — Manuscript-to-script-to-CSV mapping

Note: The definitive simulation results (2000 reps x 25 scenarios) come from
the R implementation in Pairwise70. This Python capsule reproduces a lightweight
simulation (250 reps, 4 scenarios) for verification and accessibility.

Dependencies: numpy, scipy (see requirements.txt)
"""

import csv
import os
import sys
import time

import numpy as np

# Ensure we can import from the same directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grey_meta_v8 import (
    GRMA,
    compare_methods,
    get_bcg_data,
    get_morris_data,
    simulate_scenario,
    valley_diagnostic,
)

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
SEED = 20260213


def write_csv(filename, rows, fieldnames):
    path = os.path.join(OUT_DIR, filename)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {filename} ({len(rows)} rows)")


def main():
    t0 = time.time()
    print("=" * 60)
    print("GRMA v8 — Reproduce All Tables (Python capsule)")
    print("=" * 60)
    print()

    # ==================================================================
    # 1. Default parameters table
    # ==================================================================
    print("[1/7] Default parameters...")
    defaults = [
        {"parameter": "zeta", "default": "0.5",
         "meaning": "Grey resolution coefficient in GRC denominator"},
        {"parameter": "norm_method", "default": "robust_minmax",
         "meaning": "Fitted robust scaling to [0,1]"},
        {"parameter": "anchor_mode", "default": "median",
         "meaning": "Effect anchor uses median(y)"},
        {"parameter": "trim", "default": "0.1",
         "meaning": "Trim fraction if anchor_mode=trimmed_mean"},
        {"parameter": "prec_cap", "default": "1e6",
         "meaning": "Upper cap on precision before log transform"},
        {"parameter": "effect_guard", "default": "True",
         "meaning": "Apply Tukey bisquare guard"},
        {"parameter": "tukey_c", "default": "4.685",
         "meaning": "Bisquare tuning constant"},
        {"parameter": "guard_power", "default": "1",
         "meaning": "Exponent on guard term"},
    ]
    write_csv("Table_defaults_v8.csv", defaults,
              ["parameter", "default", "meaning"])

    # ==================================================================
    # 2. Simulation: Bias & RMSE
    # ==================================================================
    print("[2/7] Simulation (250 reps x 4 scenarios)...")
    scenarios = [
        ("ideal", {"k": 15, "true_effect": 0.5, "tau2": 0.05,
                   "scenario_type": "ideal"}),
        ("outlier", {"k": 15, "true_effect": 0.5, "tau2": 0.05,
                     "scenario_type": "outlier",
                     "outlier_spec": {"n": 2, "shift": 3.0}}),
        ("high_prec", {"k": 15, "true_effect": 0.5, "tau2": 0.05,
                       "scenario_type": "high_prec"}),
        ("pub_bias", {"k": 15, "true_effect": 0.5, "tau2": 0.05,
                      "scenario_type": "pub_bias"}),
    ]

    sim_rows = []
    for sc_name, sc_params in scenarios:
        results = simulate_scenario(
            n_rep=250, seed=SEED, **sc_params
        )
        for r in results:
            sim_rows.append({
                "scenario": sc_name,
                "method": r["method"],
                "bias": f"{r['bias']:.7g}",
                "rmse": f"{r['rmse']:.7g}",
                "n_ok": r["n_ok"],
            })

    write_csv("Table_bias_rmse_v8_py.csv", sim_rows,
              ["scenario", "method", "bias", "rmse", "n_ok"])

    # ==================================================================
    # 3. Applied results (BCG + Morris)
    # ==================================================================
    print("[3/7] Applied results (BCG + Morris)...")
    datasets = [
        get_bcg_data(),
        get_morris_data(),
    ]

    applied_rows = []
    for yi, vi, label in datasets:
        results = compare_methods(yi, vi, seed=SEED)
        for r in results:
            applied_rows.append({
                "dataset": label,
                "method": r["method"],
                "estimate": f"{r['estimate']:.7g}",
                "se": f"{r.get('se', float('nan')):.7g}",
                "ci_lo": f"{r['ci_lo']:.7g}",
                "ci_hi": f"{r['ci_hi']:.7g}",
                "ci_lo_bca": f"{r.get('ci_lo_bca', float('nan')):.7g}",
                "ci_hi_bca": f"{r.get('ci_hi_bca', float('nan')):.7g}",
            })

    write_csv("Table_applied_v8_py.csv", applied_rows,
              ["dataset", "method", "estimate", "se", "ci_lo", "ci_hi",
               "ci_lo_bca", "ci_hi_bca"])

    # ==================================================================
    # 4. Diagnostics
    # ==================================================================
    print("[4/7] Diagnostics...")
    diag_rows = []
    for yi, vi, label in datasets:
        g = GRMA(effect_guard=True)
        fit = g.fit(yi, vi)
        vd = valley_diagnostic(yi, fit["estimate"], seed=SEED)
        diag_rows.append({
            "dataset": label,
            "k": fit["k"],
            "anchor_y": f"{fit['anchor_y']:.7g}",
            "anchor_p": f"{fit['anchor_p']:.7g}",
            "zeta": "0.5",
            "w_max": f"{fit['w_max']:.7g}",
            "n_eff": f"{fit['n_eff']:.5g}",
            "valley_flag": str(vd["valley_flag"]),
            "valley_p": f"{vd['valley_p']:.6g}",
        })

    write_csv("Table_diagnostics_v8_py.csv", diag_rows,
              ["dataset", "k", "anchor_y", "anchor_p", "zeta",
               "w_max", "n_eff", "valley_flag", "valley_p"])

    # ==================================================================
    # 5. Leave-one-out
    # ==================================================================
    print("[5/7] Leave-one-out influence...")
    loo_rows = []
    for yi, vi, label in datasets:
        g = GRMA(effect_guard=True)
        loo = g.leave_one_out(yi, vi)
        loo_rows.append({
            "dataset": label,
            "k": len(yi),
            "est_full": f"{loo['est_full']:.7g}",
            "max_abs_delta_est": f"{loo['max_abs_delta_est']:.7g}",
            "idx_max_shift": loo["idx_max_shift"],
            "y_idx": f"{loo['y_at_max']:.7g}",
            "v_idx": f"{loo['v_at_max']:.7g}",
            "maxw_full": f"{loo['w_max_full']:.7g}",
            "max_abs_delta_maxw": f"{loo['max_abs_delta_maxw']:.7g}",
            "idx_maxw_shift": loo["idx_maxw_shift"],
        })

    write_csv("Table_loo_v8_py.csv", loo_rows,
              ["dataset", "k", "est_full", "max_abs_delta_est",
               "idx_max_shift", "y_idx", "v_idx", "maxw_full",
               "max_abs_delta_maxw", "idx_maxw_shift"])

    # ==================================================================
    # 6. Zeta sensitivity (BCG)
    # ==================================================================
    print("[6/7] Zeta sensitivity (BCG)...")
    yi_bcg, vi_bcg, _ = get_bcg_data()
    zeta_rows = []
    for z in np.arange(0.1, 1.0, 0.1):
        g = GRMA(zeta=z, effect_guard=True)
        ci = g.bootstrap_ci(yi_bcg, vi_bcg, B=999, bca=True, seed=SEED)
        fit = g.fit(yi_bcg, vi_bcg)
        zeta_rows.append({
            "zeta": f"{z:.1f}",
            "estimate": f"{ci['estimate']:.7g}",
            "ci_lo_pct": f"{ci['ci_lo_pct']:.7g}",
            "ci_hi_pct": f"{ci['ci_hi_pct']:.7g}",
            "ci_lo_bca": f"{ci.get('ci_lo_bca', float('nan')):.7g}",
            "ci_hi_bca": f"{ci.get('ci_hi_bca', float('nan')):.7g}",
            "w_max": f"{fit['w_max']:.7g}",
            "n_eff": f"{fit['n_eff']:.5g}",
        })

    write_csv("Table_zeta_sensitivity_bcg_v8_py.csv", zeta_rows,
              ["zeta", "estimate", "ci_lo_pct", "ci_hi_pct",
               "ci_lo_bca", "ci_hi_bca", "w_max", "n_eff"])

    # ==================================================================
    # 7. Artifact mapping table
    # ==================================================================
    print("[7/7] Artifact mapping...")
    artifact_rows = [
        {"manuscript_item": "Table 1 (Bias/RMSE, R sim)",
         "script": "run_grma_simulation.R",
         "output_csv": "Table_bias_rmse_v8.csv",
         "note": "Definitive: 2000 reps x 25 scenarios"},
        {"manuscript_item": "Table 1 (Bias/RMSE, Python verification)",
         "script": "reproduce_all_v8.py",
         "output_csv": "Table_bias_rmse_v8_py.csv",
         "note": "Verification: 250 reps x 4 scenarios"},
        {"manuscript_item": "Table 2 (Coverage, R sim)",
         "script": "run_grma_simulation.R",
         "output_csv": "Table_coverage_v8.csv",
         "note": "Definitive: 2000 reps x 25 scenarios"},
        {"manuscript_item": "Table 3 (Applied results)",
         "script": "applied_examples.R",
         "output_csv": "Table_applied_v8.csv",
         "note": "R (authoritative); Table_applied_v8_py.csv (Python)"},
        {"manuscript_item": "Table 4 (Diagnostics)",
         "script": "applied_examples.R",
         "output_csv": "Table_diagnostics_v8.csv",
         "note": "R (authoritative)"},
        {"manuscript_item": "Table 5 (LOO summary)",
         "script": "applied_examples.R",
         "output_csv": "Table_loo_v8.csv",
         "note": "R (authoritative)"},
        {"manuscript_item": "Table S1 (Zeta sensitivity)",
         "script": "applied_examples.R",
         "output_csv": "Table_zeta_sensitivity_bcg_v8.csv",
         "note": "R (authoritative)"},
        {"manuscript_item": "Table S2 (Guard ablation)",
         "script": "run_grma_simulation.R",
         "output_csv": "Table_ablation_v8.csv",
         "note": "Definitive: derived from simulation"},
        {"manuscript_item": "Table S0 (Defaults)",
         "script": "reproduce_all_v8.py",
         "output_csv": "Table_defaults_v8.csv",
         "note": "Static table"},
    ]
    write_csv("Table_artifact_map_v8.csv", artifact_rows,
              ["manuscript_item", "script", "output_csv", "note"])

    # ==================================================================
    # Summary
    # ==================================================================
    elapsed = time.time() - t0
    print()
    print("=" * 60)
    print(f"All tables produced in {elapsed:.1f} seconds.")
    print("=" * 60)

    # Quick cross-check: print BCG estimates for verification
    print()
    print("--- BCG Quick Check (Python) ---")
    yi_bcg, vi_bcg, _ = get_bcg_data()
    g = GRMA()
    fit = g.fit(yi_bcg, vi_bcg)
    ci = g.bootstrap_ci(yi_bcg, vi_bcg, B=999, bca=True, seed=SEED)
    print(f"  GRMA estimate: {fit['estimate']:.8f}")
    print(f"  w_max: {fit['w_max']:.8f}")
    print(f"  n_eff: {fit['n_eff']:.4f}")
    print(f"  Percentile CI: [{ci['ci_lo_pct']:.6f}, {ci['ci_hi_pct']:.6f}]")
    print(f"  BCa CI:        [{ci['ci_lo_bca']:.6f}, {ci['ci_hi_bca']:.6f}]")


if __name__ == "__main__":
    main()
