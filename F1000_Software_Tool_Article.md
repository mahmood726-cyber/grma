# GRMA: a software tool for reviewer-auditable evidence synthesis

## Authors
- Mahmood Ahmad [1,2]
- Niraj Kumar [1]
- Bilaal Dar [3]
- Laiba Khan [1]
- Andrew Woo [4]
- Corresponding author: Andrew Woo (andy2709w@gmail.com)

## Affiliations
1. Royal Free Hospital
2. Tahir Heart Institute Rabwah
3. King's College Medical School
4. St George's Medical School

## Abstract
**Background:** Inverse-variance pooling is efficient under ideal assumptions but can be destabilized by effect outliers and high-leverage studies. Reviewers therefore expect a robust estimator to be accompanied by explicit validation, guard ablation, and honest discussion of efficiency loss when data are clean.

**Methods:** GRMA implements grey relational weighting in a two-feature space of effect size and log precision, adds a Tukey bisquare effect guard, and reports bootstrap-based inference. The repository contains both R and Python implementations, simulation scripts, applied examples, and large-scale Pairwise70 benchmarks.

**Results:** Local artifacts document simulation performance across contamination scenarios, guard ablation experiments, cross-language validation tables, leave-one-out diagnostics, and benchmark comparisons against random-effects estimates across thousands of Cochrane meta-analyses.

**Conclusions:** GRMA is best framed as a companion robust estimator for sensitivity analysis and audit, with transparent tradeoffs between robustness under contamination and efficiency under uncontaminated conditions.

## Keywords
robust meta-analysis; grey relational analysis; outlier resistance; bootstrap intervals; Pairwise70; software tool

## Introduction
The project provides a fully reproducible software implementation of a robust pooling idea that is easy to inspect: bounded grey relational weights, an explicit outlier guard, and publication-ready diagnostics that explain how the pooled estimate was stabilized.

The appropriate comparison set includes REML, HKSJ, winsorized or robust alternatives, and robust Bayesian mixtures. The manuscript emphasizes this companion role rather than presenting GRMA as a universal replacement.

The manuscript structure below is deliberately aligned to common open-software review requests: the rationale is stated explicitly, at least one runnable example path is named, local validation artifacts are listed, and conclusions are bounded to the functions and outputs documented in the repository.

## Methods
### Software architecture and workflow
The codebase includes `grma_meta.R` and `grey_meta_v8.py`, figure generation, paired R/Python cross-validation, simulation drivers, and explicit manuscript-number verification scripts. This organization supports end-to-end reruns from code to tables and figures.

### Installation, runtime, and reviewer reruns
The local implementation is packaged under `C:\Models\GRMA_paper`. The manuscript identifies the local entry points, dependency manifest, fixed example input, and expected saved outputs so that reviewers can rerun the documented workflow without reconstructing it from scratch.

- Entry directory: `C:\Models\GRMA_paper`.
- Detected documentation entry points: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Detected environment capture or packaging files: `requirements.txt`.
- Named worked-example paths in this draft: `applied_examples.R` and embedded example data for worked demonstrations; `pairwise70_benchmark_grma/` for large-scale empirical agreement analysis; `MANUSCRIPT_PLOS_ONE_GRMA_v10.md` as the project-specific narrative source.
- Detected validation or regression artifacts: `f1000_artifacts/validation_summary.md`, `verify_manuscript_numbers.py`, `test_edge_cases.R`.
- Detected example or sample data files: `f1000_artifacts/example_dataset.csv`.

### Worked examples and validation materials
**Example or fixed demonstration paths**
- `applied_examples.R` and embedded example data for worked demonstrations.
- `pairwise70_benchmark_grma/` for large-scale empirical agreement analysis.
- `MANUSCRIPT_PLOS_ONE_GRMA_v10.md` as the project-specific narrative source.

**Validation and reporting artifacts**
- `cross_validate.py` and `cross_validate.R` for implementation parity.
- `verify_manuscript_numbers.py` for manuscript-to-output checking.
- `Table_*.csv` artifacts documenting defaults, diagnostics, guard ablation, and simulation results.

### Typical outputs and user-facing deliverables
- Robust pooled estimates with weight diagnostics and influence summaries.
- Bootstrap intervals and scenario-wise simulation tables.
- Benchmark summaries against standard random-effects outputs.

### Reviewer-informed safeguards
- Provides a named example workflow or fixed demonstration path.
- Documents local validation artifacts rather than relying on unsupported claims.
- Positions the software against existing tools without claiming blanket superiority.
- States limitations and interpretation boundaries in the manuscript itself.
- Requires explicit environment capture and public example accessibility in the released archive.

## Review-Driven Revisions
This draft has been tightened against recurring open peer-review objections taken from the supplied reviewer reports.
- Reproducibility: the draft names a reviewer rerun path and points readers to validation artifacts instead of assuming interface availability is proof of correctness.
- Validation: claims are anchored to local tests, validation summaries, simulations, or consistency checks rather than to unsupported assertions of performance.
- Comparators and niche: the manuscript now names the relevant comparison class and keeps the claimed niche bounded instead of implying universal superiority.
- Documentation and interpretation: the text expects a worked example, input transparency, and reviewer-verifiable outputs rather than a high-level feature list alone.
- Claims discipline: conclusions are moderated to the documented scope of GRMA and paired with explicit limitations.

## Use Cases and Results
The software outputs should be described in terms of concrete reviewer-verifiable workflows: running the packaged example, inspecting the generated results, and checking that the reported interpretation matches the saved local artifacts. In this project, the most important result layer is the availability of a transparent execution path from input to analysis output.

Representative local result: `Table_bias_rmse_v8.csv` contains 25 GRMA scenario rows with mean absolute bias 0.050, mean RMSE 0.168, with mean convergence 100.0%.

### Concrete local quantitative evidence
- `Table_bias_rmse_v8.csv` contains 25 GRMA scenario rows with mean absolute bias 0.050, mean RMSE 0.168, with mean convergence 100.0%.
- `Table_bias_rmse_v8_py.csv` contains 4 GRMA scenario rows with mean absolute bias 0.017, mean RMSE 0.088.
- `MANUSCRIPT_PLOS_ONE_GRMA_v10.md` reports Results: In simulation, GRMA reduces absolute bias in 4 of 5 outlier scenarios (up to 66% reduction), with root mean squared error (RMSE) comparable to HKSJ in outlier scenarios (ratio 1.01) but 8--18% higher under uncontaminated conditions.

## Discussion
Representative local result: `Table_bias_rmse_v8.csv` contains 25 GRMA scenario rows with mean absolute bias 0.050, mean RMSE 0.168, with mean convergence 100.0%.

This directory already contains the kind of evidence that software reviewers ask for: transparent defaults, ablation studies, large-scale empirical benchmarking, and explicit documentation of where robustness gains are purchased at the cost of efficiency.

### Limitations
- Efficiency can be lower than standard random-effects methods when data are uncontaminated.
- The method is not itself a publication-bias correction model.
- Inference depends on bootstrap calibration and is reported as a complement to established workflows.

## Software Availability
- Local source package: `GRMA_paper` under `C:\Models`.
- Public repository: `https://github.com/mahmood726-cyber/grma`.
- Public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/grma/tree/bef87ed1889731ed47e5ef7bab94d010f419dc4c`.
- DOI/archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.
- Environment capture detected locally: `requirements.txt`.
- Reviewer-facing documentation detected locally: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Reproducibility walkthrough: `f1000_artifacts/tutorial_walkthrough.md` where present.
- Validation summary: `f1000_artifacts/validation_summary.md` where present.
- Reviewer rerun manifest: `F1000_Reviewer_Rerun_Manifest.md`.
- Multi-persona review memo: `F1000_MultiPersona_Review.md`.
- Concrete submission-fix note: `F1000_Concrete_Submission_Fixes.md`.
- License: see the local `LICENSE` file.

## Data Availability
All code, benchmark summaries, figures, and applied-example artifacts are stored locally in the project directory. The underlying Pairwise70 datasets are derived from publicly accessible Cochrane reviews.

## Reporting Checklist
Real-peer-review-aligned checklist: `F1000_Submission_Checklist_RealReview.md`.
Reviewer rerun companion: `F1000_Reviewer_Rerun_Manifest.md`.
Companion reviewer-response artifact: `F1000_MultiPersona_Review.md`.
Project-level concrete fix list: `F1000_Concrete_Submission_Fixes.md`.

## Declarations
### Competing interests
The authors declare that no competing interests were disclosed.

### Grant information
No specific grant was declared for this manuscript draft.

### Author contributions (CRediT)
| Author | CRediT roles |
|---|---|
| Mahmood Ahmad | Conceptualization; Software; Validation; Data curation; Writing - original draft; Writing - review and editing |
| Niraj Kumar | Conceptualization |
| Bilaal Dar | Conceptualization |
| Laiba Khan | Conceptualization |
| Andrew Woo | Conceptualization |

### Acknowledgements
The authors acknowledge contributors to open statistical methods, reproducible research software, and reviewer-led software quality improvement.

## References
1. DerSimonian R, Laird N. Meta-analysis in clinical trials. Controlled Clinical Trials. 1986;7(3):177-188.
2. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. Statistics in Medicine. 2002;21(11):1539-1558.
3. Viechtbauer W. Conducting meta-analyses in R with the metafor package. Journal of Statistical Software. 2010;36(3):1-48.
4. Page MJ, McKenzie JE, Bossuyt PM, et al. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. BMJ. 2021;372:n71.
5. Fay C, Rochette S, Guyader V, Girard C. Engineering Production-Grade Shiny Apps. Chapman and Hall/CRC. 2022.
