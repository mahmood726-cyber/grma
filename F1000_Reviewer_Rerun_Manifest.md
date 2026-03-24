# GRMA: reviewer rerun manifest

This manifest is the shortest reviewer-facing rerun path for the local software package. It lists the files that should be sufficient to recreate one worked example, inspect saved outputs, and verify that the manuscript claims remain bounded to what the repository actually demonstrates.

## Reviewer Entry Points
- Project directory: `C:\Models\GRMA_paper`.
- Preferred documentation start points: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Detected public repository root: `https://github.com/mahmood726-cyber/grma`.
- Detected public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/grma/tree/bef87ed1889731ed47e5ef7bab94d010f419dc4c`.
- Detected public archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.
- Environment capture files: `requirements.txt`.
- Validation/test artifacts: `f1000_artifacts/validation_summary.md`, `verify_manuscript_numbers.py`, `test_edge_cases.R`.

## Worked Example Inputs
- Manuscript-named example paths: `applied_examples.R` and embedded example data for worked demonstrations; `pairwise70_benchmark_grma/` for large-scale empirical agreement analysis; `MANUSCRIPT_PLOS_ONE_GRMA_v10.md` as the project-specific narrative source; f1000_artifacts/example_dataset.csv.
- Auto-detected sample/example files: `f1000_artifacts/example_dataset.csv`.

## Expected Outputs To Inspect
- Robust pooled estimates with weight diagnostics and influence summaries.
- Bootstrap intervals and scenario-wise simulation tables.
- Benchmark summaries against standard random-effects outputs.

## Minimal Reviewer Rerun Sequence
- Start with the README/tutorial files listed below and keep the manuscript paths synchronized with the public archive.
- Create the local runtime from the detected environment capture files if available: `requirements.txt`.
- Run at least one named example path from the manuscript and confirm that the generated outputs match the saved validation materials.
- Quote one concrete numeric result from the local validation snippets below when preparing the final software paper.

## Local Numeric Evidence Available
- `Table_bias_rmse_v8.csv` contains 25 GRMA scenario rows with mean absolute bias 0.050, mean RMSE 0.168, with mean convergence 100.0%.
- `Table_bias_rmse_v8_py.csv` contains 4 GRMA scenario rows with mean absolute bias 0.017, mean RMSE 0.088.
- `MANUSCRIPT_PLOS_ONE_GRMA_v10.md` reports Results: In simulation, GRMA reduces absolute bias in 4 of 5 outlier scenarios (up to 66% reduction), with root mean squared error (RMSE) comparable to HKSJ in outlier scenarios (ratio 1.01) but 8--18% higher under uncontaminated conditions.
