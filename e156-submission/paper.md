Mahmood Ahmad
Tahir Heart Institute
author@example.com

Grey Relational Meta-Analysis: A Robust Pooling Method with Redescending Effect Guard

Can grey relational analysis provide a transparent and robust alternative to inverse-variance meta-analytic pooling? GRMA was evaluated via simulation with 2,000 replicates across 25 scenarios against four comparators and benchmarked on 4,572 real Cochrane meta-analyses from the Pairwise70 dataset. The method represents each study in a two-feature space of effect size and log-precision, computing grey relational coefficients with an explicit Tukey bisquare redescending effect guard and bootstrap percentile interval inference. GRMA reduced median absolute bias of the pooled RR by up to 66 percent in four of five outlier scenarios, with bootstrap 95% CI coverage of 0.940. In the Pairwise70 benchmark, GRMA estimates correlated at r equals 0.956 with standard random-effects estimates, increasing to 0.980 at k of at least 10. GRMA offers a transparent companion to standard random-effects models for detecting and down-weighting influential outliers in routine synthesis. However, a limitation is that statistical power under uncontaminated conditions was 13 percentage points lower than HKSJ.

Outside Notes

Type: methods
Primary estimand: Bias reduction (percent)
App: GRMA v10.0 (R + Python)
Data: 4,572 Cochrane meta-analyses (Pairwise70), 2,000-rep simulation
Code: https://github.com/mahmood726-cyber/grma
Version: 10.0
Validation: DRAFT

References

1. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.
2. Higgins JPT, Thompson SG, Deeks JJ, Altman DG. Measuring inconsistency in meta-analyses. BMJ. 2003;327(7414):557-560.
3. Cochrane Handbook for Systematic Reviews of Interventions. Version 6.4. Cochrane; 2023.

AI Disclosure

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI (Claude, Anthropic) was used as a constrained synthesis engine operating on structured inputs and predefined rules for infrastructure generation, not as an autonomous author. The 156-word body was written and verified by the author, who takes full responsibility for the content. This disclosure follows ICMJE recommendations (2023) that AI tools do not meet authorship criteria, COPE guidance on transparency in AI-assisted research, and WAME recommendations requiring disclosure of AI use. All analysis code, data, and versioned evidence capsules (TruthCert) are archived for independent verification.
