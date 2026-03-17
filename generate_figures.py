"""
Generate publication-quality figures for the GRMA v8 manuscript.
Requires: matplotlib, numpy, pandas (or csv fallback).

Produces:
  Fig1_bias_rmse_by_scenario.png   — Bias and RMSE panel plot (all 25 scenarios)
  Fig2_coverage_comparison.png     — Coverage by scenario category
  Fig3_guard_ablation.png          — Guard vs no-guard for O1-O5
  Fig4_benchmark_scatter.png       — GRMA vs RE scatter (Pairwise70, 4572 analyses)
  Fig5_benchmark_kband.png         — Agreement metrics by k-band
"""

import csv
import os
import sys
import json
from collections import defaultdict

import numpy as np

# Force UTF-8 stdout on Windows
sys.stdout = __import__('io').TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

# Try matplotlib
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
except ImportError:
    print("ERROR: matplotlib required. Install with: pip install matplotlib")
    sys.exit(1)

BASE = os.path.dirname(os.path.abspath(__file__))

# ── Style ──
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 9,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

# Color palette (colorblind-safe, adapted from Okabe-Ito)
COLORS = {
    'GRMA': '#0072B2',       # blue
    'HKSJ': '#D55E00',       # vermillion
    'REML': '#CC79A7',       # reddish purple
    'WRD': '#009E73',        # bluish green
    'RBM': '#F0E442',        # yellow
    'GRMA_noguard': '#56B4E9',  # sky blue
}

# ── Load CSVs ──
def load_csv(filename):
    rows = []
    with open(os.path.join(BASE, filename), newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

bias_rmse = load_csv("Table_bias_rmse_v8.csv")
coverage = load_csv("Table_coverage_v8.csv")
ablation = load_csv("Table_ablation_v8.csv")

# Build lookup dicts
br = {}
for row in bias_rmse:
    br[(row['scenario'], row['method'])] = row

cov = {}
for row in coverage:
    cov[(row['scenario'], row['method'])] = row

abl = {}
for row in ablation:
    abl[row['scenario']] = row

# Scenario ordering
scenario_order = [
    'B1','B2','B3','B4','B5',
    'H1','H2','H3','H4','H5',
    'O1','O2','O3','O4','O5',
    'PB1','PB2','PB3','PB4','PB5',
    'S1','S2','S3',
    'D1','D2',
]

category_spans = {
    'Baseline': (0, 5),
    'Heterogeneity': (5, 5),
    'Outlier': (10, 5),
    'Pub. bias': (15, 5),
    'Small study': (20, 3),
    'Distribution': (23, 2),
}

methods_main = ['GRMA', 'HKSJ', 'REML', 'WRD', 'RBM']


# =====================================================================
# FIGURE 1: Bias and RMSE across all 25 scenarios (2-panel)
# =====================================================================
def fig1_bias_rmse():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    x = np.arange(len(scenario_order))
    width = 0.15
    offsets = np.array([-2, -1, 0, 1, 2]) * width

    for ax, metric, title in [(ax1, 'bias', 'A. Absolute Bias'),
                               (ax2, 'rmse', 'B. RMSE')]:
        for i, method in enumerate(methods_main):
            vals = []
            for sc in scenario_order:
                key = (sc, method)
                if key in br:
                    v = float(br[key][metric])
                    vals.append(abs(v) if metric == 'bias' else v)
                else:
                    vals.append(0)
            ax.bar(x + offsets[i], vals, width, label=method,
                   color=COLORS[method], edgecolor='white', linewidth=0.3,
                   zorder=3)

        ax.set_ylabel('|Bias|' if metric == 'bias' else 'RMSE')
        ax.set_title(title, loc='left', fontweight='bold')
        ax.grid(axis='y', alpha=0.3, zorder=0)

        # Category shading
        for cat, (start, count) in category_spans.items():
            if cat == 'Outlier':
                ax.axvspan(start - 0.4, start + count - 0.6, alpha=0.08,
                          color='#0072B2', zorder=0)

    ax2.set_xticks(x)
    ax2.set_xticklabels(scenario_order, rotation=45, ha='right')
    ax1.legend(ncol=5, loc='upper right', framealpha=0.9)

    # Category labels at bottom
    for cat, (start, count) in category_spans.items():
        mid = start + count / 2 - 0.5
        ax2.annotate(cat, xy=(mid, 0), xytext=(mid, -0.06),
                    textcoords=('data', 'axes fraction'),
                    ha='center', fontsize=7, fontstyle='italic',
                    annotation_clip=False)

    fig.suptitle('Figure 1. Simulation performance across 25 scenarios (2000 replicates)',
                 fontsize=11, fontweight='bold', y=1.01)
    plt.tight_layout()
    out = os.path.join(BASE, 'Fig1_bias_rmse_by_scenario.png')
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# =====================================================================
# FIGURE 2: Coverage comparison by scenario category
# =====================================================================
def fig2_coverage():
    fig, ax = plt.subplots(figsize=(8, 4.5))

    methods_cov = ['GRMA', 'HKSJ', 'REML', 'RBM']
    colors_cov = [COLORS[m] for m in methods_cov]

    cats = list(category_spans.keys())
    x = np.arange(len(cats))
    width = 0.18
    offsets = np.array([-1.5, -0.5, 0.5, 1.5]) * width

    for i, method in enumerate(methods_cov):
        means = []
        for cat in cats:
            start, count = category_spans[cat]
            scenarios = scenario_order[start:start+count]
            vals = []
            for sc in scenarios:
                key = (sc, method)
                if key in cov:
                    vals.append(float(cov[key]['coverage']))
            means.append(np.mean(vals) if vals else 0)
        ax.bar(x + offsets[i], means, width, label=method,
               color=colors_cov[i], edgecolor='white', linewidth=0.3, zorder=3)

    ax.axhline(y=0.95, color='black', linestyle='--', linewidth=0.8,
               label='Nominal 95%', zorder=2)
    ax.set_ylabel('Mean Coverage')
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylim(0.0, 1.05)
    ax.legend(ncol=5, loc='lower left', framealpha=0.9)
    ax.grid(axis='y', alpha=0.3, zorder=0)
    ax.set_title('Figure 2. Mean coverage by scenario category (2000 replicates)',
                 fontsize=11, fontweight='bold')

    plt.tight_layout()
    out = os.path.join(BASE, 'Fig2_coverage_comparison.png')
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# =====================================================================
# FIGURE 3: Guard ablation (O1-O5)
# =====================================================================
def fig3_guard_ablation():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))

    o_scenarios = ['O1', 'O2', 'O3', 'O4', 'O5']
    o_labels = ['O1\n(mild)', 'O2\n(extreme)', 'O3\n(opposing)',
                'O4\n(multiple)', 'O5\n(clustered)']
    x = np.arange(len(o_scenarios))
    width = 0.35

    # Panel A: Bias
    bias_guard = [abs(float(abl[sc]['bias_guard'])) for sc in o_scenarios]
    bias_noguard = [abs(float(abl[sc]['bias_noguard'])) for sc in o_scenarios]

    ax1.bar(x - width/2, bias_guard, width, label='GRMA (guard)',
            color=COLORS['GRMA'], edgecolor='white', linewidth=0.3)
    ax1.bar(x + width/2, bias_noguard, width, label='GRMA (no guard)',
            color=COLORS['GRMA_noguard'], edgecolor='white', linewidth=0.3)
    ax1.set_ylabel('|Bias|')
    ax1.set_title('A. Absolute Bias', fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(o_labels)
    ax1.legend(framealpha=0.9)
    ax1.grid(axis='y', alpha=0.3)

    # Mark O3 as within MC noise
    ax1.annotate('\u2020', xy=(2, max(bias_guard[2], bias_noguard[2])),
                 ha='center', va='bottom', fontsize=12, color='gray')

    # Panel B: RMSE
    rmse_guard = [float(abl[sc]['rmse_guard']) for sc in o_scenarios]
    rmse_noguard = [float(abl[sc]['rmse_noguard']) for sc in o_scenarios]

    ax2.bar(x - width/2, rmse_guard, width, label='GRMA (guard)',
            color=COLORS['GRMA'], edgecolor='white', linewidth=0.3)
    ax2.bar(x + width/2, rmse_noguard, width, label='GRMA (no guard)',
            color=COLORS['GRMA_noguard'], edgecolor='white', linewidth=0.3)
    ax2.set_ylabel('RMSE')
    ax2.set_title('B. RMSE', fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(o_labels)
    ax2.legend(framealpha=0.9)
    ax2.grid(axis='y', alpha=0.3)

    fig.suptitle('Figure 3. Guard ablation: outlier scenarios (2000 replicates)',
                 fontsize=11, fontweight='bold', y=1.02)
    plt.tight_layout()
    out = os.path.join(BASE, 'Fig3_guard_ablation.png')
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# =====================================================================
# FIGURE 4: Pairwise70 benchmark scatter (GRMA vs RE)
# =====================================================================
def fig4_benchmark_scatter():
    # Load benchmark results
    bench_file = os.path.join(BASE, 'pairwise70_benchmark_grma', 'analysis_results.csv')
    mu_re_vals = []
    mu_grma_vals = []
    k_vals = []

    with open(bench_file, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['mu_grma'] and row['mu_re']:
                try:
                    mu_re = float(row['mu_re'])
                    mu_grma = float(row['mu_grma'])
                    k = int(row['k'])
                    if k >= 3 and abs(mu_re) < 5 and abs(mu_grma) < 5:
                        mu_re_vals.append(mu_re)
                        mu_grma_vals.append(mu_grma)
                        k_vals.append(k)
                except (ValueError, KeyError):
                    pass

    mu_re_arr = np.array(mu_re_vals)
    mu_grma_arr = np.array(mu_grma_vals)
    k_arr = np.array(k_vals)

    fig, ax = plt.subplots(figsize=(6, 6))

    # Color by k-band
    k3_mask = k_arr == 3
    k49_mask = (k_arr >= 4) & (k_arr <= 9)
    k10_mask = k_arr >= 10

    ax.scatter(mu_re_arr[k3_mask], mu_grma_arr[k3_mask],
               s=4, alpha=0.15, color='#CC79A7', label=f'k=3 (n={k3_mask.sum()})', zorder=2)
    ax.scatter(mu_re_arr[k49_mask], mu_grma_arr[k49_mask],
               s=6, alpha=0.2, color='#0072B2', label=f'k=4-9 (n={k49_mask.sum()})', zorder=3)
    ax.scatter(mu_re_arr[k10_mask], mu_grma_arr[k10_mask],
               s=10, alpha=0.4, color='#D55E00', label=f'k>=10 (n={k10_mask.sum()})', zorder=4)

    # Identity line
    lims = [-4, 4]
    ax.plot(lims, lims, 'k--', linewidth=0.8, alpha=0.5, zorder=1)

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel('Random-Effects estimate (log OR)')
    ax.set_ylabel('GRMA estimate (log OR)')
    ax.set_aspect('equal')
    ax.legend(loc='upper left', framealpha=0.9)
    ax.grid(alpha=0.2)

    # Correlation annotation
    r = np.corrcoef(mu_re_arr, mu_grma_arr)[0, 1]
    ax.text(0.97, 0.03, f'r = {r:.3f} (n = {len(mu_re_arr)})',
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_title('Figure 4. GRMA vs. RE estimates\n(Pairwise70 Cochrane benchmark, k >= 3)',
                 fontsize=10, fontweight='bold')

    plt.tight_layout()
    out = os.path.join(BASE, 'Fig4_benchmark_scatter.png')
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# =====================================================================
# FIGURE 5: Benchmark k-band agreement metrics
# =====================================================================
def fig5_benchmark_kband():
    # Load from summary.json
    with open(os.path.join(BASE, 'pairwise70_benchmark_grma', 'summary.json'),
              encoding='utf-8') as f:
        summary = json.load(f)

    kband_data = summary['grma_comparison']['by_k_band']
    bands = ['k_3', 'k_4_to_9', 'k_ge_10']
    band_labels = ['k = 3', 'k = 4-9', 'k >= 10']
    n_vals = [kband_data[b]['n'] for b in bands]
    corr_vals = [kband_data[b]['correlation_grma_re'] for b in bands]
    shift_vals = [kband_data[b]['median_abs_shift'] for b in bands]
    wmax_vals = [kband_data[b]['median_w_max'] for b in bands]
    neff_vals = [kband_data[b]['median_n_eff'] for b in bands]

    fig, axes = plt.subplots(2, 2, figsize=(8, 6))

    x = np.arange(3)
    bar_color = '#0072B2'

    # A: Correlation
    ax = axes[0, 0]
    ax.bar(x, corr_vals, color=bar_color, edgecolor='white')
    ax.set_ylim(0.9, 1.0)
    ax.set_ylabel('Pearson r')
    ax.set_title('A. Correlation (GRMA vs RE)', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(band_labels)
    ax.axhline(1.0, color='gray', linestyle=':', linewidth=0.5)
    for i, v in enumerate(corr_vals):
        ax.text(i, v + 0.001, f'{v:.3f}', ha='center', va='bottom', fontsize=8)

    # B: Median absolute shift
    ax = axes[0, 1]
    ax.bar(x, shift_vals, color='#D55E00', edgecolor='white')
    ax.set_ylabel('Median |shift| (log OR)')
    ax.set_title('B. Median absolute shift', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(band_labels)
    for i, v in enumerate(shift_vals):
        ax.text(i, v + 0.002, f'{v:.3f}', ha='center', va='bottom', fontsize=8)

    # C: Median w_max
    ax = axes[1, 0]
    ax.bar(x, wmax_vals, color='#009E73', edgecolor='white')
    ax.set_ylabel('Median w_max')
    ax.set_title('C. Maximum weight (w_max)', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(band_labels)
    ax.axhline(0.5, color='gray', linestyle=':', linewidth=0.5)
    for i, v in enumerate(wmax_vals):
        ax.text(i, v + 0.01, f'{v:.3f}', ha='center', va='bottom', fontsize=8)

    # D: Median n_eff
    ax = axes[1, 1]
    ax.bar(x, neff_vals, color='#CC79A7', edgecolor='white')
    ax.set_ylabel('Median n_eff')
    ax.set_title('D. Effective sample size (n_eff)', fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(band_labels)
    for i, v in enumerate(neff_vals):
        ax.text(i, v + 0.2, f'{v:.1f}', ha='center', va='bottom', fontsize=8)

    fig.suptitle('Figure 5. GRMA weight diagnostics by study count\n'
                 '(Pairwise70 Cochrane benchmark)',
                 fontsize=11, fontweight='bold', y=1.03)
    plt.tight_layout()
    out = os.path.join(BASE, 'Fig5_benchmark_kband.png')
    fig.savefig(out)
    plt.close(fig)
    print(f"  Saved {out}")


# =====================================================================
# Main
# =====================================================================
if __name__ == '__main__':
    print("Generating GRMA v8 manuscript figures...")
    fig1_bias_rmse()
    fig2_coverage()
    fig3_guard_ablation()
    fig4_benchmark_scatter()
    fig5_benchmark_kband()
    print("Done. All figures saved to:", BASE)
