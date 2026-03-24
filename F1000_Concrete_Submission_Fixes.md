# GRMA: concrete submission fixes

This file converts the multi-persona review into repository-side actions that should be checked before external submission of the F1000 software paper for `GRMA_paper`.

## Detectable Local State
- Documentation files detected: `README.md`, `f1000_artifacts/tutorial_walkthrough.md`.
- Environment lock or container files detected: `requirements.txt`.
- Package manifests detected: none detected.
- Example data files detected: `f1000_artifacts/example_dataset.csv`.
- Validation artifacts detected: `f1000_artifacts/validation_summary.md`, `verify_manuscript_numbers.py`, `test_edge_cases.R`.
- Detected public repository root: `https://github.com/mahmood726-cyber/grma`.
- Detected public source snapshot: Fixed public commit snapshot available at `https://github.com/mahmood726-cyber/grma/tree/bef87ed1889731ed47e5ef7bab94d010f419dc4c`.
- Detected public archive record: No project-specific DOI or Zenodo record URL was detected locally; archive registration pending.

## High-Priority Fixes
- Check that the manuscript's named example paths exist in the public archive and can be run without repository archaeology.
- Confirm that the cited repository root (`https://github.com/mahmood726-cyber/grma`) resolves to the same fixed public source snapshot used for submission.
- Archive the tagged release and insert the Zenodo DOI or record URL once it has been minted; no project-specific archive DOI was detected locally.
- Reconfirm the quoted benchmark or validation sentence after the final rerun so the narrative text matches the shipped artifacts.

## Numeric Evidence Available To Quote
- `Table_bias_rmse_v8.csv` contains 25 GRMA scenario rows with mean absolute bias 0.050, mean RMSE 0.168, with mean convergence 100.0%.
- `Table_bias_rmse_v8_py.csv` contains 4 GRMA scenario rows with mean absolute bias 0.017, mean RMSE 0.088.
- `MANUSCRIPT_PLOS_ONE_GRMA_v10.md` reports Results: In simulation, GRMA reduces absolute bias in 4 of 5 outlier scenarios (up to 66% reduction), with root mean squared error (RMSE) comparable to HKSJ in outlier scenarios (ratio 1.01) but 8--18% higher under uncontaminated conditions.

## Manuscript Files To Keep In Sync
- `F1000_Software_Tool_Article.md`
- `F1000_Reviewer_Rerun_Manifest.md`
- `F1000_MultiPersona_Review.md`
- `F1000_Submission_Checklist_RealReview.md` where present
- README/tutorial files and the public repository release metadata
