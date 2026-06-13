# Truth-recovery yardstick — grma (Grey Relational Meta-Analysis)

**Verdict: STRONG VALIDATION of the robustness claim — with the honest finding
that the robustness comes specifically from the Tukey-bisquare GUARD, not the
grey-relational weighting, at a modest clean-data efficiency cost.**

## Method
GRMA is sold as a robust pooling estimator (grey-relational similarity + a
Tukey-bisquare "effect guard" that down-weights aberrant studies). The honest test
is a known-truth simulation: under a known `mu`, does the guard let GRMA RECOVER
the truth better than DL/REML when a few studies are OUTLIERS, without a meaningful
price on clean data? The harness uses the app's OWN `grey_meta_v8` (GRMA +
`_dl_random_effects` + `_reml_random_effects` + `_hk_reml`), run unchanged. DGP:
`yi ~ N(mu, vi+tau2)`; the `outlier` scenario shifts 20% of studies by +1.5. BCa
bootstrap B=199 (app ships 999; coverage is insensitive at this resolution).
250 reps/cell, mu=0.3, tau2=0.05.

## Results (coverage of true mu / RMSE of the point estimate)

| scenario | k  | DL_RE | HK_REML | **GRMA (guard)** | GRMA_noguard |
|----------|----|-------|---------|------------------|--------------|
| clean    | 8  | 0.888 / 0.110 | 0.924 / 0.110 | 0.888 / 0.128 | 0.892 / 0.111 |
| clean    | 15 | 0.916 / 0.085 | 0.928 / 0.086 | 0.920 / 0.099 | 0.892 / 0.087 |
| outlier  | 8  | 0.880 / **0.379** | 0.964 / 0.383 | 0.932 / **0.168** | 0.660 / 0.311 |
| outlier  | 15 | 0.704 / **0.314** | 0.796 / 0.317 | 0.904 / **0.128** | 0.524 / 0.257 |

## Findings (all measured)
1. **VALIDATION — the guard delivers real outlier-robustness.** Under 20% outlier
   contamination, GRMA (guard) recovers the true mu with **2–2.5× lower RMSE** than
   DL/REML (0.168/0.128 vs 0.379/0.314) and keeps near-nominal coverage
   (0.93/0.90) exactly where DL collapses (RMSE explodes; coverage falls to 0.70 at
   k=15). The Tukey-bisquare guard does what it claims.
2. **The GUARD is the active ingredient, not the grey-relational weighting.**
   GRMA_noguard (grey-relational weights, guard off) is *not* robust — under
   outliers it has high RMSE (0.31/0.26) and badly under-covers (0.66/0.52, worse
   than plain DL). So the robustness is attributable to the bisquare guard
   specifically; the grey-relational layer alone does not provide it. (Honest, and
   useful for the manuscript's attribution of the effect.)
3. **Modest clean-data cost.** On clean data GRMA (guard) carries slightly higher
   RMSE (0.128 vs DL 0.110 at k=8, ~16%) for maintained coverage (~0.89–0.92) —
   the standard robust-estimator efficiency trade-off, honestly quantified. It does
   not break or over-shrink when there are no outliers.
4. HK_REML has the best clean-data coverage (the HKSJ widening) but no
   outlier-robustness — its point estimate is dragged just like DL's; only its
   wider CI keeps coverage up, at the cost of a wrong center.

## What did NOT transfer
NPE/conformal machinery is estimator-of-distribution specific; GRMA is a robust
point+bootstrap estimator, so the known-truth bias/RMSE/coverage harness (with the
outlier-contamination stress the repo's Cochrane benchmark does not isolate)
transferred. The shipped Python engine is run unchanged; no new dependency.

## Reproduce
```
python truth-recovery/harness.py --reps 250
python truth-recovery/test_truth_recovery.py
```
