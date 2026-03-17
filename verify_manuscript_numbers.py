"""
Verify specific numbers cited in the GRMA manuscript text.
Reads Table_bias_rmse_v8.csv, Table_coverage_v8.csv, Table_ablation_v8.csv.
"""
import csv
import sys
import io

# Force UTF-8 stdout on Windows
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import os
BASE = os.path.dirname(os.path.abspath(__file__))

def read_csv(filename):
    rows = []
    with open(os.path.join(BASE, filename), newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

# -------------------------------------------------------------------
# Load data
# -------------------------------------------------------------------
bias_rmse = read_csv("Table_bias_rmse_v8.csv")
coverage = read_csv("Table_coverage_v8.csv")
ablation = read_csv("Table_ablation_v8.csv")

# Build lookup: (scenario, method) -> row for bias_rmse
br = {}
for row in bias_rmse:
    key = (row['scenario'], row['method'])
    br[key] = row

cov = {}
for row in coverage:
    key = (row['scenario'], row['method'])
    cov[key] = row

abl = {}
for row in ablation:
    abl[row['scenario']] = row

# -------------------------------------------------------------------
print("=" * 70)
print("A. OUTLIER BIAS REDUCTION CLAIMS (O2, O4, O5)")
print("=" * 70)

for sc in ['O2', 'O4', 'O5']:
    grma_bias = abs(float(br[(sc, 'GRMA')]['bias']))
    hksj_bias = abs(float(br[(sc, 'HKSJ')]['bias']))
    reduction_pct = (1 - grma_bias / hksj_bias) * 100
    print(f"  {sc}: GRMA |bias| = {grma_bias:.4f}, HKSJ |bias| = {hksj_bias:.4f}, "
          f"reduction = {reduction_pct:.1f}%")

# -------------------------------------------------------------------
print()
print("=" * 70)
print("B. OUTLIER RMSE RATIO GRMA/HKSJ (mean across O1-O5)")
print("=" * 70)

o_scenarios = ['O1', 'O2', 'O3', 'O4', 'O5']
ratios = []
for sc in o_scenarios:
    grma_rmse = float(br[(sc, 'GRMA')]['rmse'])
    hksj_rmse = float(br[(sc, 'HKSJ')]['rmse'])
    ratio = grma_rmse / hksj_rmse
    ratios.append(ratio)
    print(f"  {sc}: GRMA RMSE = {grma_rmse:.4f}, HKSJ RMSE = {hksj_rmse:.4f}, ratio = {ratio:.3f}")
print(f"  Mean ratio across O1-O5: {sum(ratios)/len(ratios):.3f}")

# -------------------------------------------------------------------
print()
print("=" * 70)
print("C. COVERAGE EXCLUDING PB SCENARIOS (mean GRMA, HKSJ, REML)")
print("=" * 70)

all_scenarios = ['B1','B2','B3','B4','B5',
                 'H1','H2','H3','H4','H5',
                 'S1','S2','S3',
                 'D1','D2',
                 'O1','O2','O3','O4','O5']  # excludes PB1-PB5

for method in ['GRMA', 'HKSJ', 'REML']:
    coverages = []
    for sc in all_scenarios:
        key = (sc, method)
        if key in cov:
            coverages.append(float(cov[key]['coverage']))
    mean_cov = sum(coverages) / len(coverages)
    print(f"  {method}: mean coverage = {mean_cov:.4f} ({mean_cov*100:.2f}%) over {len(coverages)} scenarios")

# -------------------------------------------------------------------
print()
print("=" * 70)
print("D. ABLATION TABLE FOR O1-O5")
print("=" * 70)

print(f"  {'Scenario':<10} {'|bias| guard':>13} {'|bias| noguard':>15} "
      f"{'rmse guard':>11} {'rmse noguard':>13} {'guard_helps_bias':>17} {'guard_helps_rmse':>17}")
print("  " + "-" * 95)

for sc in o_scenarios:
    r = abl[sc]
    bg = abs(float(r['bias_guard']))
    bng = abs(float(r['bias_noguard']))
    rg = float(r['rmse_guard'])
    rng = float(r['rmse_noguard'])
    helps_bias = bg < bng
    helps_rmse = rg < rng
    print(f"  {sc:<10} {bg:>13.6f} {bng:>15.6f} {rg:>11.6f} {rng:>13.6f} "
          f"{'YES' if helps_bias else 'no':>17} {'YES' if helps_rmse else 'no':>17}")

# -------------------------------------------------------------------
print()
print("=" * 70)
print("E. HOW MANY OF 25 SCENARIOS DOES GUARD HELP BIAS? RMSE?")
print("=" * 70)

all_25 = sorted(abl.keys())
help_bias_count = 0
help_rmse_count = 0
help_bias_list = []
help_rmse_list = []

for sc in all_25:
    r = abl[sc]
    bg = abs(float(r['bias_guard']))
    bng = abs(float(r['bias_noguard']))
    rg = float(r['rmse_guard'])
    rng = float(r['rmse_noguard'])
    if bg < bng:
        help_bias_count += 1
        help_bias_list.append(sc)
    if rg < rng:
        help_rmse_count += 1
        help_rmse_list.append(sc)

print(f"  Guard helps bias (|bias| guard < |bias| noguard): {help_bias_count}/25")
print(f"    Scenarios: {', '.join(help_bias_list)}")
print(f"  Guard helps RMSE (rmse guard < rmse noguard): {help_rmse_count}/25")
print(f"    Scenarios: {', '.join(help_rmse_list)}")

# -------------------------------------------------------------------
print()
print("=" * 70)
print("F. ABLATION O2 AND O5 BIAS REDUCTION %")
print("=" * 70)

for sc in ['O2', 'O5']:
    r = abl[sc]
    bg = abs(float(r['bias_guard']))
    bng = abs(float(r['bias_noguard']))
    reduction = (1 - bg / bng) * 100
    print(f"  {sc}: |bias| guard = {bg:.6f}, |bias| noguard = {bng:.6f}, "
          f"reduction = {reduction:.1f}%")

# -------------------------------------------------------------------
print()
print("=" * 70)
print("G. RMSE PREMIUM FOR UNCONTAMINATED CATEGORIES (B, H, S, D)")
print("=" * 70)

categories = {
    'B': ['B1','B2','B3','B4','B5'],
    'H': ['H1','H2','H3','H4','H5'],
    'S': ['S1','S2','S3'],
    'D': ['D1','D2'],
}

all_ratios = []
for cat, scenarios in categories.items():
    cat_ratios = []
    for sc in scenarios:
        grma_rmse = float(br[(sc, 'GRMA')]['rmse'])
        hksj_rmse = float(br[(sc, 'HKSJ')]['rmse'])
        ratio = grma_rmse / hksj_rmse
        cat_ratios.append(ratio)
        all_ratios.append(ratio)
    mean_r = sum(cat_ratios) / len(cat_ratios)
    print(f"  {cat}: mean GRMA/HKSJ RMSE ratio = {mean_r:.3f} "
          f"(range {min(cat_ratios):.3f} - {max(cat_ratios):.3f})")

overall_min = min(all_ratios)
overall_max = max(all_ratios)
overall_mean = sum(all_ratios) / len(all_ratios)
print(f"  Overall (B+H+S+D): mean = {overall_mean:.3f}, range = {overall_min:.3f} - {overall_max:.3f}")

# -------------------------------------------------------------------
print()
print("=" * 70)
print("EXTRA: Full bias/RMSE for all outlier scenarios, all methods")
print("=" * 70)
for sc in o_scenarios:
    print(f"\n  {sc} ({br[(sc,'GRMA')].get('name','')}):")
    for method in ['GRMA', 'HKSJ', 'REML', 'RBM', 'WRD', 'GRMA_noguard']:
        key = (sc, method)
        if key in br:
            b = float(br[key]['bias'])
            r = float(br[key]['rmse'])
            print(f"    {method:<15} bias={b:+.6f}  |bias|={abs(b):.6f}  rmse={r:.6f}")
