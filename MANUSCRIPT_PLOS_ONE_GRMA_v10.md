# Grey Relational Meta-Analysis: A Robust Pooling Method with Redescending Effect Guard

**Authors:** Mahmood Ul Hassan^1^

^1^ Independent Researcher

**Corresponding author:** Mahmood Ul Hassan ([CORRESPONDING_EMAIL_PLACEHOLDER])

**ORCID:** [ORCID_PLACEHOLDER]

---

## Abstract

**Background:** Inverse-variance meta-analysis can be sensitive to extreme effect outliers and high-leverage studies. Existing robust alternatives (M-estimators, heavy-tailed likelihoods, mixture models) address this, but a practical niche remains for methods combining effect and precision information transparently with nonparametric inference.

**Methods:** We propose **Grey Relational Meta-Analysis (GRMA)**, a robust pooling estimator using grey relational similarity in a two-feature space (effect size and log-precision) with an explicit **Tukey bisquare effect guard** and **bootstrap percentile interval** inference. We evaluate GRMA via simulation (2,000 replicates, 25 scenarios, 4 comparators: restricted maximum likelihood [REML], Hartung-Knapp-Sidik-Jonkman [HKSJ], weighted robust dispersion [WRD], robust Bayesian mixture [RBM]) and on 4,572 real Cochrane meta-analyses (Pairwise70 benchmark).

**Results:** In simulation, GRMA reduces absolute bias in 4 of 5 outlier scenarios (up to 66% reduction), with root mean squared error (RMSE) comparable to HKSJ in outlier scenarios (ratio 1.01) but 8--18% higher under uncontaminated conditions. Power is correspondingly lower under standard conditions (59% vs. 72% at k = 10, tau^2 = 0.05). Bootstrap percentile intervals achieve mean coverage of 0.940 (excluding publication-bias stress tests, where all methods suffer). In the Pairwise70 benchmark (4,572 Cochrane meta-analyses), GRMA estimates correlate r = 0.956 with random-effects estimates, increasing to r = 0.980 at k >= 10.

**Conclusions:** GRMA is a transparent robust pooling method proposed as a companion to standard random-effects models. All code (R and Python) and data are provided in a reproducible capsule under MIT license.

---

## Introduction

Inverse-variance (IV) weighting is widely used for meta-analysis under fixed-effect and random-effects models. While efficient under ideal assumptions, a small number of atypical studies can substantially shift pooled results---either because the study's effect estimate is extreme, or because its sampling variance is exceptionally small (creating high leverage).

Several robust meta-analytic approaches have been proposed. M-estimators adapt robust regression ideas to the meta-analytic setting, using redescending weight functions (e.g., Huber or Tukey bisquare) to down-weight influential observations while iteratively re-estimating the pooled effect [1]. Heavy-tailed likelihood models replace the normal assumption with a Student-t distribution, providing automatic down-weighting of discrepant studies [2]. Mixture and contamination models explicitly model a fraction of studies as arising from a different distribution [3,4]. Influence diagnostics, including leave-one-out analysis and Cook's distance analogues, identify but do not automatically handle influential studies [1,5].

There remains a practical niche for methods that are: (i) transparent in their weighting mechanism, (ii) computationally lightweight, (iii) robust to both effect outliers and leverage from very small sampling variances, and (iv) packaged in a reproducible way for routine use. Among existing robust methods, Huber M-estimators [1] use redescending weight functions but remain variance-weighted, inheriting the leverage vulnerability of inverse-variance allocation. The weighted median estimator [6], originally developed for Mendelian randomization, achieves breakdown-point robustness but discards all precision information beyond rank ordering. The present method differs from these approaches in its use of grey relational coefficients (GRC)---a bounded, two-feature kernel that naturally incorporates precision via log-transformation without the leverage vulnerability of inverse-variance allocation---combined with an explicit redescending effect guard. The GRC maps distances to similarities via a resolution parameter, providing smooth, bounded weights whose maximum contribution is inherently limited. We note that the weighted median was not included in the simulation comparison (see Discussion, Limitations); the niche claim should be evaluated in light of future head-to-head comparisons.

Grey relational analysis (GRA) converts "similarity to a reference profile" into weights [7]. We adapt this idea to meta-analysis by representing each study in a simple two-feature space: effect and log-precision, then weighting studies by their grey relational similarity to an anchor profile. Because grey relational similarity alone is not guaranteed to be redescending for extreme effect outliers, we add an explicit redescending effect guard.

This manuscript focuses on the estimator, diagnostics, simulation performance, and reproducibility. GRMA is a pooling method and is not a publication-bias correction model; publication-bias analyses should still be performed using appropriate sensitivity/bias-modeling tools (e.g., selection models, PET-PEESE [8], p-curve [9]).

---

## Materials and methods

### Notation

For studies i = 1,...,k, let y_i be an effect estimate (e.g., log risk ratio [RR]) with sampling variance v_i > 0. Precision is p_i = 1/v_i.

### Feature construction

GRMA uses a two-dimensional study representation:

* **Effect feature:** y_i
* **Precision feature:** log(min(p_i, prec_cap) + 1), where prec_cap = 10^6 by default.

The log-transformation of precision serves as a deliberate design choice against ultra-precision leverage: it compresses the dynamic range of precisions so that a study with v_i = 10^{-5} does not dominate one with v_i = 0.01 as severely as under raw inverse-variance weighting. The additive 1 ensures the feature is non-negative for all valid variances.

### Fitted robust normalization (default)

Distances are computed in a normalized feature space. For each feature, we fit robust scaling parameters on the observed data---specifically, the 5th and 95th percentiles define the mapping range---and map values to the unit interval ("robust_minmax"). The *same fitted transforms* are applied to both study points and the anchor point, preventing inconsistent scaling.

### Anchor definition

The anchor is defined as:

* Effect anchor: a_y = median(y) (default; trimmed mean available)
* Precision anchor: a_p = max(p). The maximum-precision anchor is chosen because it defines the "ideal" reference as the most informative study, so that similarity to this anchor rewards precision. This is analogous to inverse-variance weighting's preference for precise studies, but mediated through the bounded GRC rather than through direct variance weighting.

These anchor coordinates are then transformed using the fitted normalization described above.

### Grey relational coefficient and grade

Let Delta_{ij} = |x_{ij} - a_j| denote absolute deviations in the normalized space. Define Delta_min and Delta_max as the global minimum and maximum across all 2k delta values (i.e., across both features j = 1, 2 and all studies i = 1, ..., k). For a resolution coefficient zeta in (0,1]:

* GRC_{ij} = (Delta_min + zeta * Delta_max) / (Delta_{ij} + zeta * Delta_max)   ...(Eq 1)
* Base grey grade: g_i = mean_j(GRC_{ij}) (unweighted average over the 2 features; equal feature weighting is adopted as the default because no prior information favors one feature dimension over the other in the general case, and the robust min-max normalization already standardizes their scales)

The resolution coefficient zeta controls contrast: larger values reduce discrimination between studies. The default zeta = 0.5 is standard in grey relational analysis [7].

### Explicit redescending effect guard (Tukey bisquare; default enabled)

> **Boxed Algorithm: Effect Guard**
>
> 1. Compute robust residuals: u_i = |y_i - a_y| / (MAD(y) + epsilon)
>    where MAD(y) = median absolute deviation = median(|y_i - median(y)|) and epsilon = 10^{-12}.
> 2. Apply Tukey bisquare with tuning constant c (default 4.685):
>    h_i = (1 - (u_i/c)^2)^2 if u_i < c, else 0.   ...(Eq 2)
> 3. Combine: w_i^* = g_i * (h_i)^{guard_power}
> 4. Normalize: w_i = w_i^* / sum(w_j^*). If all weights are zero, fall back to equal weights.

The tuning constant c = 4.685 is adopted from robust regression [10]. Note that the guard uses the unscaled MAD rather than the consistent estimator MAD * 1.4826; the effective cutoff is therefore approximately 4.685 * 0.6745 = 3.16 standard deviations under normality. Our choice of c = 4.685 with unscaled MAD corresponds to an effective tuning of c_eff = c / 1.4826 ≈ 3.16 when expressed in units of the consistent scale estimator---more aggressive than the standard 95%-efficient Tukey bisquare. This deliberate choice favors robustness over efficiency in the meta-analytic setting, where sample sizes (studies) are typically small (k = 5--20) and contamination resistance is prioritized. Note also that the Tukey bisquare provides smooth downweighting well before the hard zero cutoff: a study at 3.0 unscaled-MAD units (u = 3.0) receives h = (1-(3.0/4.685)^2)^2 = 0.35, while at u = 4.0 it receives h = 0.07. Studies with |u_i| >= c receive exactly zero weight regardless of their precision.

### Pooled estimate and dominance metrics

> **Boxed Algorithm: Complete GRMA Estimator**
>
> **Input:** y_i (effect sizes), v_i (sampling variances), parameters (zeta, tukey_c, guard_power, anchor_mode, prec_cap)
>
> 1. Cap precision: p_i = min(1/v_i, prec_cap)
> 2. Construct features: f1 = y_i, f2 = log(p_i + 1)
> 3. Fit robust min-max normalization on (f1, f2)
> 4. Compute anchor: a_y = median(y), a_p = max(p); transform through same normalization
> 5. Compute GRC and grey grade g_i (see Grey relational coefficient above)
> 6. If effect_guard: compute guard h_i, set w_i^* = g_i * h_i^{power} (see Effect guard above)
>    Else: w_i^* = g_i
> 7. Normalize: w_i = w_i^* / sum(w_j^*)
> 8. **Output:** mu_hat = sum(w_i * y_i)
>
> **Diagnostics:** w_max = max(w_i), n_eff = 1/sum(w_i^2)

We recommend reporting weight dominance metrics alongside the pooled estimate:

* w_max = max_i(w_i): maximum weight any single study receives
* n_eff = 1 / sum(w_i^2): effective number of studies

### Bootstrap inference (full BCa interval)

We use nonparametric bootstrap resampling of studies with replacement [11,12]. For each of B resamples (default B = 999), we refit the full estimator (including anchor, normalization, and weights) to obtain mu_hat^*_b.

**Percentile interval:** [Q_{alpha/2}(mu_hat^*), Q_{1-alpha/2}(mu_hat^*)]   ...(Eq 3)

**Bias-corrected and accelerated (BCa) interval** (full, with jackknife acceleration):

1. *Bias correction z_0:* z_0 = Phi^{-1}(proportion of mu_hat^* < mu_hat)
2. *Jackknife acceleration a_hat:* For each i = 1,...,k, compute leave-one-out estimate mu_hat_{(-i)}. Let d_i = mean(mu_hat_{(-j)}) - mu_hat_{(-i)}. Then a_hat = sum(d_i^3) / (6 * (sum(d_i^2))^{3/2}).
3. *Adjusted percentiles [11]:* Let w = z_0 + z_{alpha/2}. Then alpha_1 = Phi(z_0 + w / (1 - a_hat * w)). Similarly for alpha_2 with w = z_0 + z_{1-alpha/2}.
4. *BCa interval:* [Q_{alpha_1}(mu_hat^*), Q_{alpha_2}(mu_hat^*)]

The BCa interval corrects for both bias (via z_0) and skewness (via a_hat) in the bootstrap distribution, providing improved coverage relative to the basic percentile interval, particularly in small samples. Note that the p-value reported by the software is a Wald-type approximation (z = estimate/SE_boot); it should be interpreted with caution at small k.

**Implementation detail:** The jackknife acceleration a_hat is clipped to [-0.5, 0.5] to prevent extreme BCa corrections from unstable leave-one-out deletions (e.g., when removing a study leaves k-1 = 2 studies with degenerate normalization). If all guard-adjusted weights are zero (e.g., all studies are extreme outliers), the estimator falls back to equal weights (w_i = 1/k) with a software warning.

### Valley diagnostic (exploratory)

An exploratory "valley" diagnostic tests for a bimodality-like pattern by comparing within-cluster dispersion (splitting at the pooled estimate) with global dispersion, calibrating a p-value under a permutation null (using the Phipson-Smyth +1/+1 correction [13]). This is included in the software as an exploratory *warning* signal, not as a definitive test; it has not been validated in the simulation study and its operating characteristics (Type I error, power) are unknown. In the applied examples below, no datasets triggered the valley flag. Detailed description is in the supplementary materials.

### Estimator properties and invariances

The estimator is designed to have practical invariances useful in synthesis work. We note that formal asymptotic theory (breakdown point, asymptotic efficiency) for the grey relational weighting scheme has not been developed; the properties below are design-level guarantees verified empirically in the simulation study (see Results).

1. **Shift equivariance (effects):** If all effects are shifted by a constant y_i' = y_i + c, then the anchor shifts by the same constant (median(y+c) = median(y)+c), the MAD is shift-invariant, and the fitted normalization percentiles also shift by c, preserving the normalization range and hence all normalized distances. Thus weights are unchanged and mu_hat shifts by c.

2. **Approximate scale equivariance (effects only):** If y_i' = a*y_i for a > 0 with variances held fixed, the anchor and MAD scale by a and the effect-feature normalization scales correspondingly. However, because the log-precision feature log(1/v + 1) does not transform linearly under the proper meta-analytic scaling v_i' = a^2*v_i, exact scale equivariance does not hold when variances co-scale with effects. Empirically, the deviation is small (< 0.1% in typical configurations).

3. **Redescending effect suppression:** With effect_guard enabled, sufficiently extreme |y_i - a_y| produce h_i=0, forcing w_i*=0 regardless of precision, preventing ultra-precise extreme effects from dominating the pooled estimate.

4. **Bounded weights:** The GRC is bounded in (0, 1] by construction (minimum value zeta/(1+zeta) when Delta_min = 0), and the guard further constrains weights. This prevents any single study from receiving unbounded influence.

These properties motivate reporting both the pooled estimate and the weight dominance metrics (w_max and n_eff) for interpretability and auditability.

### Comparator methods

**Simulation comparators:**

* REML random-effects (REML_RE [14])
* Hartung-Knapp with REML tau^2 (HKSJ [15,16])
* Weighted robust dispersion estimator (WRD [1])
* Robust Bayesian mixture model for outlier detection (RBM [3])
* GRMA without guard (GRMA_noguard) for ablation analysis

Note: HKSJ and REML share the same point estimate (both use REML tau^2); they differ only in confidence interval (CI) construction (HKSJ uses a t-distribution with Hartung-Knapp variance adjustment).

**Applied example comparators:** The simulation comparators above plus IV fixed-effect (IV_FE), DerSimonian-Laird random-effects (DL_RE [17]), Huber M-estimator with k = 1.345 (HuberIV), and Student-t random-effects with fixed df = 4 (tRE_df4). The Huber and Student-t estimators were not included in the simulation due to computational constraints; their systematic evaluation across all 25 scenarios is deferred to future work.

### Scope and limitations of the estimator

GRMA and the comparators all target the same population mean effect mu, though through different mechanisms. GRMA is empirically consistent: as k grows, the median anchor converges to mu under symmetric distributions and the weights become approximately uniform, so the estimator converges to mu (confirmed by simulations B1--B5). However, formal asymptotic theory has not been developed (see Discussion, Limitations).

GRMA is a pooling estimator only. It does not estimate between-study variance (tau^2; users needing tau^2 or prediction intervals [18,19] should fit REML or DL in parallel), does not extend to meta-regression or network meta-analysis, and does not model publication bias. Minimum k = 3; with fewer studies, GRMA falls back to Hartung-Knapp REML [15,16]. GRMA can be applied within subgroups but provides no analogue of the test for subgroup differences. Zero-event studies (where log odds ratio [OR] is undefined) are excluded, as with all IV-based methods. Divergence between GRMA and the primary estimator should be reported as a sensitivity analysis finding that may inform GRADE certainty assessment [20], but does not directly map to the GRADE inconsistency domain.

### Simulation study

#### Framework and scenarios

We use the Pairwise70 simulation framework (developed by the author; a collection of 501 Cochrane systematic reviews used for benchmarking; see Results, Large-scale empirical evaluation) with 25 scenarios organized in 6 categories.

The data-generating mechanism for each replicate is: (i) draw between-study heterogeneity variance tau^2 from the scenario specification; (ii) for each of k studies, generate a true study-specific effect theta_i ~ N(mu_true, tau^2) and a within-study variance v_i drawn from an inverse-gamma distribution calibrated to realistic meta-analytic settings; (iii) generate the observed effect y_i ~ N(theta_i, v_i); (iv) for outlier scenarios, replace designated studies with contaminated effects (additive shift of +/- delta standard deviations from the true effect); (v) for publication-bias scenarios, apply a step or continuous selection function based on p-value thresholds.

**Baseline (B1--B5):** k in {5, 10, 20, 50, 100}; tau^2 = 0.05; true effect mu = 0 (B1, null) or mu = 0.3 (B2--B5, for power evaluation). Tests performance across study counts.

**Heterogeneity (H1--H5):** k = 10; tau^2 in {0.0, 0.10, 0.20, 0.50, 1.00}; true effect mu = 0.3. Tests robustness to increasing between-study variance.

**Outliers (O1--O5):** k in {10, 15}; tau^2 = 0.05; contamination includes single mild outlier (+3 SD), single extreme outlier (+6 SD), two opposing outliers (+/-3 SD), multiple outliers (3 studies shifted +3 SD), and clustered outliers (2 studies shifted to a remote cluster). Non-normal effect distributions have been studied for standard meta-analysis estimators [21].

**Publication bias stress tests (PB1--PB5):** k = 20; step (p < 0.05 retained with probability 1, p >= 0.05 retained with decreasing probability) and continuous selection functions with varying severity. Note: GRMA is not a publication-bias correction model; these scenarios test sensitivity to selection-like data distortions.

**Small studies (S1--S3):** k = 10; mean study total sample sizes of 20, 50, and bimodal (mix of 20 and 500).

**Non-normal distributions (D1--D2):** k = 10; heavy-tailed (t with df = 3) and skewed (lognormal) effect distributions.

#### Run parameters

* **Replicates:** 2,000 per scenario (Monte Carlo SE = sqrt(p(1-p)/2000); < 0.005 for coverage near 0.95, up to 0.012 for coverage near 0.50 in publication-bias stress tests)
* **Methods:** REML, HKSJ, WRD, RBM, GRMA (guard), GRMA (no guard)
* **Bootstrap:** B = 499 per replicate in simulation (reduced from the default B = 999 for computational feasibility; BCa coverage is generally stable for B >= 200--500 per [11]). For applied examples, B = 999 (the software default).
* **Seed:** 20260213
* **Platform:** R 4.5.2 (metafor 4.8-0, data.table 1.16.4), Python 3.13.1 (numpy 2.2.3, scipy 1.15.2, matplotlib 3.10.0), single core (Windows 11)
* **Reporting:** We follow the Aims, Data-generating mechanisms, Estimands, Methods, and Performance measures (ADEMP) framework [22] for simulation reporting.

#### Metrics

From Pairwise70's `compute_performance_metrics()`:
* **Bias:** |mean(estimate) - true_effect|
* **RMSE:** sqrt(mean((estimate - true_effect)^2))
* **Coverage:** proportion of 95% CIs containing true_effect
* **CI width:** mean CI width
* **Convergence rate:** proportion of non-error results
* **Type I error / Power:** at alpha = 0.05 (scenario-dependent; reported in supplementary S8 Table)

Default settings for GRMA are listed in S1 Table.

---

## Results

### Simulation performance

All methods converged in all 2,000 replicates for all 25 scenarios (convergence rate = 1.000 for all 150 method-scenario combinations).

Full results are in S8 Table (bias, RMSE, coverage, CI width, Type I error, power) and the deposited CSV files. Fig 1 shows a bias/RMSE panel across all 25 scenarios and 6 methods; Fig 2 shows coverage comparisons; Fig 3 shows the guard ablation analysis. We summarize the main patterns below.

**Table 1. Mean RMSE by scenario category (GRMA vs. comparators).**

| Category | GRMA | HKSJ | REML | WRD | RBM | GRMA/HKSJ |
|---|---|---|---|---|---|---|
| Baseline (B1--B5) | 0.104 | 0.088 | 0.088 | 0.101 | 0.154 | 1.18 |
| Heterogeneity (H1--H5) | 0.225 | 0.193 | 0.193 | 0.205 | 0.296 | 1.17 |
| Outlier (O1--O5) | 0.146 | 0.144 | 0.144 | 0.148 | 0.201 | 1.01 |
| Pub. bias (PB1--PB5) | 0.233 | 0.192 | 0.192 | 0.253 | 0.302 | 1.22 |
| Small study (S1--S3) | 0.128 | 0.110 | 0.110 | 0.120 | 0.176 | 1.16 |
| Distribution (D1--D2) | 0.141 | 0.130 | 0.130 | 0.139 | 0.201 | 1.08 |

*Note:* HKSJ and REML share the same point estimate by construction and therefore have identical bias and RMSE; they differ only in CI construction (see S8 Table for coverage differences).

**Patterns:**

* **Baseline and heterogeneity (B, H, D, S):** GRMA pays an 8--18% RMSE premium relative to HKSJ/REML (category-level averages; individual scenario ratios range from 1.00 to 1.22). This is the expected efficiency cost of bounded, non-likelihood-based weights---a fundamental robustness--efficiency trade-off: GRMA's bounded weights sacrifice information to limit the influence of any single study, and this cost is borne under clean conditions where full information use would have been optimal. The value proposition is not that GRMA has lower RMSE, but that it provides a protective mechanism against outlier-driven distortion while remaining well-calibrated.

  Coverage (percentile bootstrap CI) excluding publication-bias scenarios is 0.940 (GRMA) vs. 0.939 (HKSJ) vs. 0.918 (REML), indicating adequate calibration. Per-scenario coverage ranges from 0.894 (GRMA, O4: multiple outliers) to 0.973 (GRMA, O3), with the O4 minimum comparable to HKSJ's O4 coverage of 0.863. (Note: the simulation evaluated percentile CIs, not BCa CIs, for computational feasibility. BCa CIs are available in the software and may provide modestly different coverage; see Discussion, Limitations.)

  **Power** is lower for GRMA than HKSJ under most conditions (e.g., B2: 59% vs. 72% at k = 10, tau^2 = 0.05), reflecting the same efficiency trade-off. However, under high heterogeneity (H4: tau^2 = 0.50, H5: tau^2 = 1.00), GRMA's power matches or slightly exceeds HKSJ (H4: 20.4% vs. 19.9%; H5: 14.0% vs. 12.0%), likely because GRMA's bounded weights are less sensitive to the inflated variance estimates that reduce HKSJ's effective power under extreme heterogeneity.

* **Outlier contamination (O1--O5):** GRMA achieves substantially lower bias than HKSJ in 4 of 5 scenarios. The largest reductions are in O2 (single extreme outlier: bias 0.036 vs. 0.107, 66% reduction) and O5 (clustered outlier: bias 0.046 vs. 0.089, 49% reduction), followed by O1 (mild: 35%) and O4 (multiple: 28%, bias 0.095 vs. 0.132). The RMSE advantage is modest on average (ratio 1.01) because bias reduction is partly offset by higher variance from down-weighting. In O3 (two opposing outliers that approximately cancel), GRMA's down-weighting removes information without bias benefit, resulting in higher RMSE (0.143 vs. 0.118).

* **Publication bias stress tests (PB1--PB5):** Under step-function selection (PB1--PB4), all methods suffer severe coverage loss (mean coverage: GRMA 0.35, HKSJ 0.61, RBM 0.82). Including PB5 (continuous, mild weighting), which has adequate coverage for all methods (GRMA 0.94, HKSJ 0.95), the PB1--PB5 means are GRMA 0.47, HKSJ 0.68, RBM 0.85. GRMA has higher bias under publication-bias selection (mean 0.19 vs. 0.15 for HKSJ). This is because (i) selection-induced bias is systematic rather than outlier-like---the non-significant studies that would balance the distribution are simply missing, not present as extreme values that the guard could down-weight; and (ii) the median anchor is itself pulled toward the biased (significant) center of the selected distribution, reinforcing rather than correcting the distortion. This confirms the explicit scope limitation: GRMA is not a publication-bias model.

* **High heterogeneity (H4--H5):** All methods show increased RMSE; GRMA's relative position is stable (ratio 1.16--1.17 vs. HKSJ), and legitimate study variation dominates.

### Guard ablation

**Table 2. Guard ablation: GRMA (with guard) vs. GRMA (no guard) in outlier scenarios.**

| Scenario | |Bias| guard | |Bias| no guard | RMSE guard | RMSE no guard | Guard helps bias | Guard helps RMSE |
|---|---|---|---|---|---|---|
| O1 (mild) | 0.044 | 0.063 | 0.144 | 0.136 | Yes | No |
| O2 (extreme) | 0.036 | 0.093 | 0.145 | 0.152 | Yes | Yes |
| O3 (opposing) | 0.002 | 0.003 | 0.143 | 0.120 | †    | No |
| O4 (multiple) | 0.095 | 0.123 | 0.153 | 0.157 | Yes | Yes |
| O5 (clustered) | 0.046 | 0.080 | 0.143 | 0.142 | Yes | No |

† Difference (0.001) is within Monte Carlo simulation error and should not be interpreted as a meaningful improvement.

Full results across all 25 scenarios are in S2 Table.

**Interpretation:** The guard reduces absolute bias in 4 of 5 outlier scenarios; in O3 (opposing outliers), the bias difference (0.002 vs. 0.003) is within Monte Carlo simulation error and should not be interpreted as a meaningful improvement. The largest bias reductions occur in O2 (extreme outlier: 61% reduction) and O5 (clustered: 43% reduction). However, the guard also increases variance by excluding or severely down-weighting studies, so it improves RMSE in only 2 of 5 outlier scenarios (O2 and O4, where contamination is sufficiently extreme or numerous that bias dominates the MSE decomposition).

Across all 25 scenarios, the guard reduces absolute bias in 9 of 25 cases where the difference exceeds Monte Carlo simulation error (primarily outlier and some baseline scenarios; a 10th case, O3, shows a nominally lower bias that is within MC error) and reduces RMSE in 2 of 25 cases. Under uncontaminated conditions, the guard has negligible bias impact but increases RMSE by 2--17% (median approximately 13%) due to down-weighting of legitimate tail studies, with the smallest overhead under homogeneous conditions (H1: 2%) and the largest at small k (B1: 17%). This supports the design intent: the guard is a targeted protective mechanism that activates under contamination, with a moderate efficiency cost when contamination is absent.

The guard also provides a substantial **calibration benefit**. Under the null (B1, k = 5), GRMA without guard has a Type I error rate of 10.6%, compared to 5.7% with guard---approximately halving the false positive rate. Across non-publication-bias scenarios, GRMA without guard has 2 scenarios below 90% coverage (B1: 88.1%, O4: 81.8%), versus only 1 for GRMA with guard (O4: 89.4%).

### Worked examples (public datasets)

#### BCG vaccine trials (log risk ratio)

Using 2x2 tables from the BCG vaccine dataset (13 trials [23], `dat.bcg` in metafor [1]), we compute log RR and variances using `escalc(measure="RR")`.

**Table 3. Applied results (all estimators), BCG vaccine trials.** Full results for all 5 analyses are in `Table_applied_v8.csv`. *Note:* HuberIV and tRE_df4 were not included in the simulation study; their coverage and Type I error properties were not assessed. Results for these methods should be interpreted as illustrative comparisons only.

| dataset                                     | method      | estimate | ci_lo   | ci_hi   |
| ------------------------------------------- | ----------- | -------- | ------- | ------- |
| BCG (log RR)                                | IV_FE       | -0.3990  | -0.4717 | -0.3264 |
| BCG (log RR)                                | DL_RE       | -0.5647  | -0.8479 | -0.2814 |
| BCG (log RR)                                | REML_RE     | -0.5593  | -0.9360 | -0.1826 |
| BCG (log RR)                                | HKSJ        | -0.5593  | -1.0092 | -0.1093 |
| BCG (log RR)                                | HuberIV     | -0.6183  | -0.9976 | -0.2390 |
| BCG (log RR)                                | tRE_df4     | -0.6061  | -0.9690 | -0.2432 |
| BCG (log RR)                                | GRMA        | -0.5990  | -1.0378 | -0.2649 |
| BCG (log RR)                                | GRMA (no guard) | -0.5296  | -0.8992 | -0.1535 |

BCG GRMA: -0.599 (percentile 95% CI: [-1.038, -0.265]; BCa 95% CI: [-1.415, -0.365]). On the RR scale: 0.55 [0.35, 0.77]. The GRMA estimate is close to HKSJ (-0.559, RR 0.57), indicating agreement: the clinical conclusion (BCG vaccination reduces TB risk by approximately 40--65%) is robust to the choice of pooling method. Study 6 (the well-known Indian ICMR trial with a positive log RR of 1.40, favoring control) receives near-zero GRMA weight (w = 0.002) due to the effect guard, compared to substantial weight under IV methods. This illustrates the guard's behavior: extreme effects are suppressed regardless of their precision. However, this automatic down-weighting should prompt investigation of the reasons for the outlier (e.g., the ICMR trial's possible co-contamination from environmental mycobacteria) rather than being treated as intrinsically desirable---the ICMR trial is the largest community-randomized BCG trial ever conducted, and its exclusion should be a considered decision.

#### Morris pretest-posttest control group (standardized mean difference)

Eight studies, standardized mean difference (SMD) change difference (Morris SB, 2008 [24]; `dat.raudenbush1985` in metafor [1]).

GRMA: 0.958 (percentile 95% CI: [0.613, 1.200]; BCa 95% CI: [0.622, 1.213]). The HKSJ estimate for this dataset is 1.006 (95% CI: [0.750, 1.261]) and REML is 1.006 (95% CI: [0.785, 1.226]). GRMA's lower point estimate (0.958 vs. 1.006) reflects the effect of grey relational reweighting: the anchor is at y = 0.920 (the median), and the grey relational weights down-weight studies furthest from the anchor, redistributing influence relative to inverse-variance allocation. The maximum weight of 0.191 and n_eff of 7.42 indicate well-distributed weights with no single-study dominance.

#### CD002042: Blood transfusion review

We selected three meta-analyses from the CD002042 Cochrane review (Carson et al. [25], blood transfusion thresholds) to demonstrate GRMA on a real-world clinical dataset with varying heterogeneity:

* **30-day mortality** (k = 51, I^2 = 17.1%): low heterogeneity, large k. GRMA estimate: 0.007 (log OR; 95% CI: [-0.159, 0.173]), closely matching REML (0.027) and HKSJ (0.027).

* **Congestive heart failure** (k = 20, I^2 = 34.4%): moderate heterogeneity. GRMA estimate: -0.025 (log OR, OR = 0.975; 95% CI: [-0.389, 0.302]) vs. REML (-0.101, OR = 0.904). This is a clinically meaningful divergence: REML suggests a 10% reduction in CHF risk with restrictive transfusion, while GRMA's reweighting moves the estimate toward the null (2.5% reduction). The GRMA diagnostics (S3 Table: anchor_y = 0.000, w_max = 0.098, n_eff = 15.25) indicate well-distributed weights. Leave-one-out analysis (S4 Table) identifies study 19 (log OR = 0.906, v = 0.721---a low-precision study favoring liberal transfusion) as the most influential, with its removal shifting the GRMA estimate by 0.042 log OR. This demonstrates the intended workflow: when GRMA and RE disagree, the diagnostics identify which studies drive the divergence, enabling targeted quality assessment. Under REML, study 19's large variance gives it minimal inverse-variance weight; under GRMA's grey relational weighting, its moderate distance from the anchor yields a non-negligible contribution, illustrating how the two methods respond differently to low-precision studies.

* **Transfusion exposure** (k = 57, I^2 = 97.3%): extreme heterogeneity. GRMA estimate: -1.990 (log OR; 95% CI: [-3.120, -1.222]) vs. REML (-2.467) and IV_FE (-1.500). The extreme I^2 here reflects genuine clinical diversity; GRMA's estimate lies between FE and RE.

Full results are in `Table_applied_v8.csv`.

### Weight dominance and anchor diagnostics

GRMA diagnostics for all five analyses (S3 Table) show well-distributed weights: w_max ranges from 0.032 (transfusion exposure, k = 57) to 0.191 (Morris, k = 8), and n_eff ranges from 7.42 to 48.23. No datasets triggered the valley flag. Leave-one-out influence (S4 Table) identifies the most influential studies: for the CHF analysis, study 19 (log OR = 0.906) shifts the GRMA estimate by 0.042 log OR upon removal.

### Sensitivity to zeta (BCG)

Sensitivity of the BCG GRMA estimate to zeta is modest across zeta in (0.1, 0.9): the point estimate ranges from -0.541 (zeta = 0.1) to -0.614 (zeta = 0.9), a total range of 0.073 log RR. As zeta increases, w_max decreases (from 0.200 to 0.122) and n_eff increases (from 8.55 to 11.03), indicating more uniform weighting. The BCa 95% CIs all exclude zero. Full per-zeta results are in S7 Table.

### Large-scale empirical evaluation (Pairwise70 Cochrane benchmark)

#### Design

To complement the controlled simulation study, we evaluated GRMA on real Cochrane review data from the Pairwise70 benchmark---a collection of 501 Cochrane systematic reviews covering a wide range of clinical domains, developed by the author for this study and deposited with the reproducibility capsule. Of the 501 reviews, 458 contained at least one estimable binary-outcome meta-analysis (the remaining 43 reviews contained only continuous outcomes, had insufficient data for binary outcome estimation, or had data extraction errors). For each review, all binary outcome meta-analyses with at least k = 2 studies were fitted with DerSimonian-Laird random-effects (RE [17]), and those with k >= 3 were also fitted with GRMA (using jackknife SE for computational feasibility across thousands of analyses; 3 of 4,572 analyses produced point estimates but degenerate bootstrap SE at k = 3 and are included in point-estimate comparisons but excluded from inferential summaries). This yields a large-scale comparison of GRMA vs. RE on real effect-size distributions with natural heterogeneity, outlier structures, and varying sample sizes.

#### Results

From 458 of 501 reviews with valid binary analyses, 7,467 meta-analyses were estimable, of which 4,572 met the k >= 3 minimum for GRMA. The overall Pearson correlation between GRMA and RE pooled estimates was r = 0.956, confirming strong agreement with conventional pooling.

Agreement by k-band is detailed in S5 Table. In summary, the correlation increases from r = 0.929 at k = 3 (n = 1,409) to r = 0.980 at k >= 10 (n = 772), and the median absolute shift decreases from 0.136 to 0.087 log OR units.

**Key patterns:**

* **Convergence with increasing k.** As the number of studies grows, the median absolute shift decreases from 0.136 (k = 3) to 0.087 (k >= 10), the correlation increases from 0.929 to 0.980, and the maximum weight (w_max) falls from 0.537 to 0.130. This is expected: with more studies, both methods increasingly agree, and GRMA's weights become less concentrated.

* **Effective sample size.** The median n_eff of 3.74 across all analyses indicates moderate weight concentration. At k >= 10, the median n_eff of 10.96 indicates near-uniform weighting. Note that n_eff = 1/sum(w_i^2) measures weight concentration rather than statistical information in the inverse-variance sense; equal GRMA weights (n_eff = k) mean equal study contributions regardless of precision, not optimal information use.

* **Reweighting, not shrinkage.** GRMA shifts the RE estimate toward the null in 40.6% of analyses and away from the null in 59.4%. This split indicates that GRMA reweights study contributions rather than systematically attenuating toward or away from zero, though the modest asymmetry suggests a slight tendency to shift estimates away from the null.

* **Distribution of shifts.** The median absolute shift |mu_GRMA - mu_RE| is 0.110 log OR across all 4,572 analyses (corresponding to approximately 12% change on the OR scale). At the tails: 7.5% of analyses have shifts > 0.5 log OR, and 1.0% have shifts > 1.0 log OR. The high-shift analyses (n = 342) are characterized by small k (median k = 3) and high w_max (median 0.47), confirming that GRMA diverges most from RE when few studies allow one study to dominate GRMA's grey relational weights.

The GRMA vs. RE scatter plot (S1 Fig) and stratified agreement by k-band (S2 Fig) are provided in the supporting information.

#### Interpretation

The Pairwise70 benchmark demonstrates that GRMA produces estimates with high concordance with conventional RE across thousands of real clinical meta-analyses, with agreement improving monotonically as k increases. The correlation r = 0.956 measures agreement in ranking, not accuracy of individual estimates; two estimators targeting the same parameter are expected to correlate highly even if one is biased. The median shift of 0.110 log OR (approximately 12% on the odds ratio scale, since exp(0.110) = 1.116) represents a meaningful but not dramatic reweighting---consistent with GRMA's design as a pre-specified sensitivity analysis rather than a fundamentally different estimator. The weight diagnostics confirm that GRMA achieves the intended behavior: bounded weights, no single-study dominance at moderate-to-large k, and effective sample sizes that track the actual number of studies.

---

## Discussion

### Practical guidance for synthesis workflows

We recommend using GRMA as a robustness companion rather than a replacement for standard models. In a systematic review, GRMA results should be reported as a **pre-specified sensitivity analysis** in the Methods section, with results in a supplementary table alongside the primary RE estimate.

> **Decision workflow:**
>
> 1. Fit a conventional random-effects model (e.g., REML + HK) and GRMA.
> 2. **If estimates agree** (within 10--20% relative difference): interpret as increased robustness confidence. Report both estimates. The Pairwise70 benchmark (see Results, Large-scale empirical evaluation) shows median agreement within 0.110 log OR across 4,572 clinical meta-analyses.
> 3. **If estimates differ materially:**
>    a. Inspect weight dominance: w_max > 2/k suggests a single study contributes disproportionately (this k-dependent threshold flagged 4.4% of Pairwise70 analyses, compared with 58% for the fixed w_max > 0.3 threshold, which is too sensitive at small k where high w_max is structurally inevitable); n_eff < k/2 suggests effective sample is small.
>    b. Run leave-one-out: identify which study(ies) shift the estimate most.
>    c. Check study-level risk of bias and data quality for the influential studies.
>    d. Consider subgroup or sensitivity analyses (apply GRMA within each subgroup).
>    e. Perform publication-bias sensitivity models when relevant (GRMA does not address this).
>    f. Report the discrepancy as a sensitivity analysis finding. Material divergence between GRMA and the primary estimator may inform GRADE certainty assessment [20] (e.g., supporting or reducing confidence in the effect estimate), but note that divergence between two pooling methods applied to the same studies is not what GRADE defines as "inconsistency" (which refers to heterogeneity across study results).
> 4. **Report both** pooled estimates and the diagnostics in the Summary of Findings table footnotes or a supplementary table. GRMA results can be reported using standard Preferred Reporting Items for Systematic Reviews and Meta-Analyses (PRISMA) flow diagrams with an additional row in the Summary of Findings table, or as a pre-specified sensitivity analysis in the methods section.

**Reporting checklist:**

1. GRMA estimate and percentile CI (recommended for all k; BCa CI available for k >= 8 as exploratory)
2. w_max and n_eff
3. Number of studies with guard h_i = 0 (fully excluded by the effect guard)
4. Leave-one-out maximum shift and the responsible study
5. Comparison with HKSJ (or other primary estimator)
6. Risk-of-bias assessment for studies with highest GRMA influence
7. Bootstrap replicates (B) and zeta value used

### Limitations and future work

* **Efficiency trade-off:** Under uncontaminated conditions, GRMA pays an 8--18% RMSE premium (category-level averages; see Table 1) and has correspondingly lower power (e.g., 59% vs. 72% at k = 10, tau^2 = 0.05) relative to HKSJ/REML. This is the expected cost of bounded, non-likelihood-based weights. GRMA's advantage is specific to unidirectional outlier contamination (bias reduction of 28--66% in scenarios O1--O5 excluding the opposing-outlier case O3).

* **Percentile vs. BCa coverage:** The simulation evaluated percentile bootstrap CIs (B = 499 per replicate, reduced from the default B = 999 for computational feasibility). BCa intervals are available in the software but their coverage was not separately assessed. We recommend the percentile interval for k < 8 and BCa for larger k, with the caveat that BCa coverage remains unvalidated.

* **Publication bias:** GRMA does not model publication bias. In simulation, step-function selection scenarios (PB1--PB4) caused severe coverage loss (mean: GRMA 0.35, HKSJ 0.61, RBM 0.82), because selection-induced bias is systematic rather than outlier-like. Users should not interpret GRMA coverage as reliable when publication bias is suspected; appropriate bias-assessment tools (e.g., selection models, PET-PEESE [8], p-curve [9]) should be used in parallel.

* **Limited comparator set:** The weighted median estimator [6]---arguably the most directly comparable robust method---was not included in the simulation or applied examples. The Huber M-estimator and Student-t model were evaluated only in applied examples. A systematic simulation comparison with these methods, particularly the weighted median, is an important direction for future work.

* **Formal statistical properties not derived:** Formal asymptotic theory (breakdown point, asymptotic efficiency, influence function) for grey relational weighting has not been developed. The effective breakdown point is bounded above by 0.5 (the breakdown of the median anchor) but the robust normalization's reliance on the 5th percentile (breakdown approximately 5%) may further limit robustness, particularly at small k where 5% corresponds to less than one study. Deriving these properties formally is left for future work.

* **Normalization at small k:** At k = 3, the robust min-max normalization degrades to the sample range, the MAD has only 3 points, and bootstrap resampling produces only 27 distinct resamples. Pairwise70 data show r = 0.929 and median w_max = 0.537 at k = 3 (compared with r = 0.980 at k >= 10). Results at k = 3 should be interpreted with additional caution.

* **Software and extensions:** The current implementation is provided as standalone R and Python files; integration into established packages (e.g., metafor) is planned. Multivariate and network meta-analysis extensions have not been developed and present non-trivial challenges. The simulation tested only additive outlier contamination; other mechanisms (sign-flip errors, unit-conversion errors) were not evaluated.

### Conclusions

GRMA is a robust pooling method that combines grey relational weighting with an explicit Tukey bisquare effect guard. In simulation (25 scenarios, 2,000 replicates), GRMA reduces absolute bias by 28--66% in unidirectional outlier scenarios at a cost of 8--18% higher RMSE and lower power under uncontaminated conditions. Bootstrap percentile intervals achieve mean coverage of 0.940 (excluding publication-bias stress tests). In 4,572 Cochrane meta-analyses, GRMA estimates show high concordance with conventional RE (r = 0.956, increasing to 0.980 at k >= 10).

The method's advantages are specific to outlier contamination; under clean conditions or publication-bias selection, standard estimators perform better. We recommend GRMA as a pre-specified sensitivity analysis in systematic reviews, reported alongside conventional estimators with weight dominance and influence diagnostics. When GRMA and RE agree, this strengthens confidence in the pooled result; when they disagree, the diagnostics identify which studies drive the divergence, enabling targeted quality assessment.

Future work should develop formal asymptotic theory for the grey relational weighting scheme, conduct head-to-head simulation comparisons with the weighted median and other robust estimators, and package the method for broader accessibility.

---

## Data availability statement

All code and pre-computed output tables necessary to reproduce the analyses are available in the accompanying reproducibility capsule deposited at [ZENODO_DOI_PLACEHOLDER]. The capsule is released under the MIT License.

The R implementation (`grma_meta.R`) and Python implementation (`grey_meta_v8.py`) are cross-validated to 10^{-10} tolerance on point estimates, weights, and leave-one-out diagnostics (27 automated tests comparing R and Python outputs for the BCG and Morris datasets). BCa CIs differ across languages due to different random number generator sequences; R is authoritative for inference. BCG [23] and Morris [24] data are embedded in the code. CD002042 2x2 tables are extracted from the published Cochrane review [25] and embedded in `CD002042_embedded_data.R` for full reproducibility. The Pairwise70 benchmark pre-computed summary statistics are included in `pairwise70_benchmark_grma/summary.json`; full benchmark results (7,467 analyses) are in `pairwise70_benchmark_grma/analysis_results.csv`. `make verify` checks all manuscript numbers against the deposited CSVs.

## Ethics statement

This study analyzed only published aggregate data from existing meta-analyses and systematic reviews. No ethics approval was required.

## Author contributions

**Mahmood Ul Hassan:** Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing -- original draft, Writing -- review & editing, Visualization, Project administration.

## Funding

The author received no specific funding for this work.

## Competing interests

The author has declared that no competing interests exist.

## Acknowledgments

The author thanks the Cochrane Collaboration for making systematic review data publicly available, which enabled the Pairwise70 benchmark evaluation. Computational assistance was provided by Claude (Anthropic) for code review and manuscript formatting checks; all scientific content, analysis, and interpretation are solely the author's.

---

## References

1. Viechtbauer W, Cheung MW-L. Outlier and influence diagnostics for meta-analysis. Res Synth Methods. 2010;1:112-125. doi:10.1002/jrsm.11
2. Baker R, Jackson D. A new approach to outliers in meta-analysis. Health Care Manag Sci. 2008;11:121-131. doi:10.1007/s10729-007-9041-8
3. Beath KJ. A finite mixture method for outlier detection and robustness in meta-analysis. Res Synth Methods. 2014;5:285-293. doi:10.1002/jrsm.1090
4. Lee KJ, Thompson SG. Flexible parametric models for random-effects distributions. Stat Med. 2008;27:418-434. doi:10.1002/sim.2897
5. Hedges LV, Olkin I. Statistical Methods for Meta-Analysis. Orlando: Academic Press; 1985.
6. Bowden J, Davey Smith G, Haycock PC, Burgess S. Consistent estimation in Mendelian randomization with some invalid instruments using a weighted median estimator. Genet Epidemiol. 2016;40:304-314. doi:10.1002/gepi.21965
7. Deng JL. Introduction to grey system theory. J Grey Syst. 1989;1:1-24.
8. Stanley TD, Doucouliagos H. Meta-regression approximations to reduce publication selection bias. Res Synth Methods. 2014;5:60-78. doi:10.1002/jrsm.1095
9. Simonsohn U, Nelson LD, Simmons JP. p-Curve: a key to the file-drawer. J Exp Psychol Gen. 2014;143:534-547. doi:10.1037/a0033242
10. Beaton AE, Tukey JW. The fitting of power series, meaning polynomials, illustrated on band-spectroscopic data. Technometrics. 1974;16:147-185. doi:10.1080/00401706.1974.10489171
11. Efron B, Tibshirani RJ. An Introduction to the Bootstrap. New York: Chapman & Hall; 1993.
12. Davison AC, Hinkley DV. Bootstrap Methods and Their Application. Cambridge: Cambridge University Press; 1997.
13. Phipson B, Smyth GK. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. Stat Appl Genet Mol Biol. 2010;9:Article 39. doi:10.2202/1544-6115.1585
14. Veroniki AA, Jackson D, Viechtbauer W, et al. Methods to estimate the between-study variance and its uncertainty in meta-analysis. Res Synth Methods. 2016;7:55-79. doi:10.1002/jrsm.1164
15. Hartung J, Knapp G. A refined method for the meta-analysis of controlled clinical trials with binary outcome. Stat Med. 2001;20:3875-3889. doi:10.1002/sim.1009
16. IntHout J, Ioannidis JPA, Borm GF. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. BMC Med Res Methodol. 2014;14:25. doi:10.1186/1471-2288-14-25
17. DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7:177-188. doi:10.1016/0197-2456(86)90046-2
18. Riley RD, Higgins JPT, Deeks JJ. Interpretation of random effects meta-analyses. BMJ. 2011;342:d549. doi:10.1136/bmj.d549
19. IntHout J, Ioannidis JPA, Rovers MM, Goeman JJ. Plea for routinely presenting prediction intervals in meta-analysis. BMJ Open. 2016;6:e010247. doi:10.1136/bmjopen-2015-010247
20. Guyatt GH, Oxman AD, Vist GE, et al. GRADE: an emerging consensus on rating quality of evidence and strength of recommendations. BMJ. 2008;336:924-926. doi:10.1136/bmj.39489.470347.AD
21. Kontopantelis E, Reeves D. Performance of statistical methods for meta-analysis when true study effects are non-normally distributed: a simulation study. Stat Methods Med Res. 2012;21:409-426. doi:10.1177/0962280210371563
22. Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. Stat Med. 2019;38:2074-2102. doi:10.1002/sim.8086
23. Colditz GA, Brewer TF, Berkey CS, et al. Efficacy of BCG vaccine in the prevention of tuberculosis: meta-analysis of the published literature. JAMA. 1994;271:698-702. doi:10.1001/jama.1994.03510330076038
24. Morris SB. Estimating effect sizes from pretest-posttest-control group designs. Organ Res Methods. 2008;11:364-386. doi:10.1177/1094428106291059
25. Carson JL, Stanworth SJ, Dennis JA, et al. Transfusion thresholds for guiding red blood cell transfusion. Cochrane Database Syst Rev. 2021;12:CD002042. doi:10.1002/14651858.CD002042.pub5

---

## Figures

**Fig 1. Bias and RMSE across all 25 simulation scenarios and 6 methods.** GRMA = Grey Relational Meta-Analysis; HKSJ = Hartung-Knapp-Sidik-Jonkman; WRD = weighted robust dispersion; RBM = robust Bayesian mixture; REML = restricted maximum likelihood.

**Fig 2. Coverage comparison (percentile bootstrap CI) across non-publication-bias scenarios.** Horizontal dashed line indicates nominal 95% coverage.

**Fig 3. Guard ablation: GRMA (with guard) vs. GRMA (no guard) across outlier scenarios.** "Guard" refers to the Tukey bisquare effect guard.

---

## Supporting information

**S1 Table. GRMA defaults (as implemented).**

| parameter    | default       | meaning                                                                                      |
| ------------ | ------------- | -------------------------------------------------------------------------------------------- |
| zeta         | 0.5           | Grey resolution coefficient in GRC denominator; larger values reduce contrast.               |
| norm_method  | robust_minmax | Fitted robust scaling to the unit interval for both effect and log-precision features.        |
| anchor_mode  | median        | Effect anchor uses median(y) (or trimmed_mean if chosen).                                    |
| trim         | 0.1           | Trim fraction if anchor_mode = 'trimmed_mean'.                                               |
| prec_cap     | 1e+06         | Upper cap on precision 1/v before log transform.                                             |
| effect_guard | True          | Apply Tukey bisquare guard on robust residuals from anchor to enforce redescending behavior. |
| tukey_c      | 4.685         | Tukey bisquare tuning constant (robust regression default).                                  |
| guard_power  | 1             | Exponent on guard term to strengthen/relax suppression.                                      |
| n_boot       | 999           | Bootstrap replicates for inference.                                                          |
| bca          | True          | Full BCa interval with jackknife acceleration.                                               |

See `Table_defaults_v8.csv`.

**S2 Table. Full guard ablation results across all 25 scenarios.**

See `Table_ablation_v8.csv`.

**S3 Table. GRMA diagnostics (anchor, weight dominance, valley flag) for all applied examples.**

| dataset                                     | k  | anchor_y | anchor_p | w_max  | n_eff  | valley_flag |
| ------------------------------------------- | -- | -------- | -------- | ------ | ------ | ----------- |
| BCG (log RR)                                | 13 | -0.5520  | 252.4    | 0.1360 | 10.56  | False       |
| Morris (SMD change diff)                    | 8  | 0.920    | 14.29    | 0.1912 | 7.42   | False       |
| CD002042: 30-day mortality (log OR)         | 51 | 0.000    | 72.02    | 0.0402 | 35.84  | False       |
| CD002042: CHF (log OR)                      | 20 | 0.000    | 49.93    | 0.0979 | 15.25  | False       |
| CD002042: Transfusion exposure (log OR)     | 57 | -1.912   | 268.8    | 0.0320 | 48.23  | False       |

See `Table_diagnostics_v8.csv`.

**S4 Table. Leave-one-out influence summary (GRMA only).**

| dataset                                     | k  | est_full | max_abs_delta | idx | maxw_full | max_delta_maxw |
| ------------------------------------------- | -- | -------- | ------------- | --- | --------- | -------------- |
| BCG (log RR)                                | 13 | -0.5990  | 0.0602        | 8   | 0.1360    | 0.0299         |
| Morris (SMD change diff)                    | 8  | 0.9580   | 0.0669        | 6   | 0.1912    | 0.0539         |
| CD002042: 30-day mortality (log OR)         | 51 | 0.0072   | 0.0091        | 23  | 0.0402    | 0.0014         |
| CD002042: CHF (log OR)                      | 20 | -0.0251  | 0.0419        | 19  | 0.0979    | 0.0245         |
| CD002042: Transfusion exposure (log OR)     | 57 | -1.9901  | 0.1425        | 27  | 0.0320    | 0.0016         |

Per-study leave-one-out tables are provided in the accompanying capsule. See `Table_loo_v8.csv`.

**S5 Table. GRMA vs. RE agreement by number of studies (k-band), Pairwise70 benchmark.**

| k-band | n | Median abs(shift) | Mean abs(shift) | Correlation | Median w_max | Median n_eff |
|----------|-------|-------------------|-----------------|-------------|--------------|--------------|
| k = 3    | 1,409 | 0.136             | 0.228           | 0.929       | 0.537        | 2.27         |
| k = 4--9 | 2,391 | 0.109             | 0.169           | 0.965       | 0.305        | 4.14         |
| k >= 10  | 772   | 0.087             | 0.125           | 0.980       | 0.130        | 10.96        |
| **All**  | **4,572** | **0.110**     | **0.179**       | **0.956**   | **0.346**    | **3.74**     |

Shift is measured as |mu_GRMA - mu_RE| in log OR units. See `pairwise70_benchmark_grma/summary.json`.

**S6 Table. Artifact mapping (manuscript item to script to CSV).**

See `Table_artifact_map_v8.csv`.

**S7 Table. Sensitivity of BCG GRMA estimate to zeta.**

See `Table_zeta_sensitivity_bcg_v8.csv`.

**S8 Table. Full simulation coverage, CI width, Type I error, and power.**

See `Table_coverage_v8.csv` and `Table_bias_rmse_v8.csv`.

**S1 Fig.** Pairwise70 benchmark: GRMA vs. random-effects pooled estimates for 4,572 Cochrane binary-outcome meta-analyses (r = 0.956). Each point represents one meta-analysis. See `Fig4_benchmark_scatter.png`.

**S2 Fig.** Pairwise70 benchmark stratified by number of studies (k-band). Agreement between GRMA and random-effects improves monotonically with k. See `Fig5_benchmark_kband.png`.

**Code files:**
* `grma_meta.R` -- R implementation (exported as `grma_meta()`, `grma_loo()`, `grma_valley()`)
* `grey_meta_v8.py` -- Python implementation (standalone capsule, class `GRMA`)
* `CD002042_embedded_data.R` -- Embedded 2x2 tables from Carson et al. [25]
* `applied_examples.R` -- Applied examples (BCG + Morris + CD002042)
* `run_grma_simulation.R` -- Simulation runner (2,000 reps x 25 scenarios)
* `run_pairwise70_benchmark.py` -- Pairwise70 Cochrane benchmark (see Results, Large-scale empirical evaluation)
* `reproduce_all_v8.py` -- Reproduces all Python-side tables
* `cross_validate.py` / `cross_validate.R` -- R-Python cross-validation (27 tests)
* `test_edge_cases.R` -- Edge case tests with assertions (20 tests)
* `generate_figures.py` -- Publication figures (Figs 1--3, S1--S2 Figs)
* `verify_manuscript_numbers.py` -- Verifies all cited numbers against CSVs
* `requirements.txt` -- Pinned Python dependencies (numpy, scipy, matplotlib)
* `Makefile` -- Orchestrates all reproduction steps
* `README.md` -- Setup instructions and file descriptions
* `LICENSE` -- MIT License

**Output tables:**
* `Table_bias_rmse_v8.csv` -- Simulation bias/RMSE (R, definitive)
* `Table_coverage_v8.csv` -- Simulation coverage (R, definitive)
* `Table_ablation_v8.csv` -- Guard ablation comparison
* `Table_applied_v8.csv` -- Applied results (all estimators, R)
* `Table_applied_v8_py.csv` -- Applied results (Python verification)
* `Table_diagnostics_v8.csv` -- GRMA diagnostics (R)
* `Table_loo_v8.csv` -- Leave-one-out influence (R)
* `Table_zeta_sensitivity_bcg_v8.csv` -- Zeta sensitivity (R)
* `Table_defaults_v8.csv` -- Parameter defaults
* `Table_artifact_map_v8.csv` -- Manuscript-to-script mapping
* `Table_bias_rmse_v8_py.csv` -- Simulation bias/RMSE (Python, verification)
* `Table_diagnostics_v8_py.csv` -- Diagnostics (Python verification)
* `Table_loo_v8_py.csv` -- LOO (Python verification)
* `Table_zeta_sensitivity_bcg_v8_py.csv` -- Zeta sensitivity (Python)
* `pairwise70_benchmark_grma/summary.json` -- Pairwise70 benchmark summary (see Results, Large-scale empirical evaluation)
* `pairwise70_benchmark_grma/analysis_results.csv` -- Full benchmark results (7,467 analyses)
* `MANUSCRIPT_PLOS_ONE_GRMA_v10.md` -- This manuscript

