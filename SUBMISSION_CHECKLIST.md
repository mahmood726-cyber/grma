# PLOS ONE Submission Checklist -- GRMA v10.0

> Generated 2026-03-24. This checklist documents exactly what the author must do before submission.

## Status: MANUSCRIPT READY (3 placeholders remain)

---

## Required Author Actions (must complete before submission)

### 1. Fill Placeholders in `MANUSCRIPT_PLOS_ONE_GRMA_v10.md`

| Placeholder | Location (line) | Action |
|---|---|---|
| `[CORRESPONDING_EMAIL_PLACEHOLDER]` | Line 7 | Replace with your real email address |
| `[ORCID_PLACEHOLDER]` | Line 9 | Replace with your ORCID iD (format: 0000-0000-0000-0000) |
| `[ZENODO_DOI_PLACEHOLDER]` | Line 401 | Replace with Zenodo DOI after depositing the reproducibility capsule |

### 2. Zenodo Deposit

1. Create a new Zenodo deposit (or link via GitHub release)
2. Upload the full reproducibility capsule (all files listed in the manuscript's Supporting Information)
3. Set license to MIT
4. Set author to "Mahmood Ul Hassan" (sole author, matching manuscript)
5. Copy the minted DOI and paste into the manuscript at `[ZENODO_DOI_PLACEHOLDER]`

**Note:** The existing `.zenodo.json` was created for the F1000 submission and lists 5 creators with a different name format ("Ahmad, Mahmood"). For the PLOS ONE submission, the Zenodo deposit should list the sole author as "Mahmood Ul Hassan" to match the manuscript. Either update `.zenodo.json` or configure the deposit manually.

### 3. Update README Citation

In `README.md` line 90, the citation says `[DOI to be added upon publication]`. After obtaining the Zenodo DOI, update this line as well.

---

## Verification Completed (no action needed)

### Manuscript Quality
- [x] Abstract: 234 words (within PLOS ONE 300-word limit)
- [x] References: 25 references, sequential [1]-[25] in first-appearance order
- [x] All 25 references cited in body text
- [x] Reference list complete (entries 1-25 all present)
- [x] No TODOs, FIXMEs, or placeholder text (besides the 3 listed above)
- [x] No formatting issues (double spaces are intentional equation labels only)
- [x] Total lines: 577

### Artifact Completeness (38/38 files present)
- [x] R implementation: `grma_meta.R`
- [x] Python implementation: `grey_meta_v8.py`
- [x] Cross-validation: `cross_validate.py` + `cross_validate.R` (27 tests, R-Python parity)
- [x] Edge case tests: `test_edge_cases.R` (20 tests)
- [x] Simulation tables: `Table_bias_rmse_v8.csv`, `Table_coverage_v8.csv`
- [x] Applied tables: `Table_applied_v8.csv` (R) + `Table_applied_v8_py.csv` (Python)
- [x] Diagnostics: `Table_diagnostics_v8.csv` + `Table_diagnostics_v8_py.csv`
- [x] Leave-one-out: `Table_loo_v8.csv` + `Table_loo_v8_py.csv`
- [x] Guard ablation: `Table_ablation_v8.csv`
- [x] Zeta sensitivity: `Table_zeta_sensitivity_bcg_v8.csv` + `Table_zeta_sensitivity_bcg_v8_py.csv`
- [x] Defaults: `Table_defaults_v8.csv`
- [x] Artifact map: `Table_artifact_map_v8.csv`
- [x] Benchmark: `pairwise70_benchmark_grma/summary.json` + `analysis_results.csv`
- [x] Figures: Fig1-Fig5 PNG files (all 5 present)
- [x] Manuscript number verifier: `verify_manuscript_numbers.py`
- [x] Makefile: orchestrates all reproduction steps
- [x] Requirements: `requirements.txt` (pinned versions)
- [x] R session info: `R_session_info.txt`
- [x] License: MIT

### Sections Present (PLOS ONE requirements)
- [x] Title
- [x] Authors + affiliations
- [x] Abstract (Background / Methods / Results / Conclusions)
- [x] Introduction
- [x] Materials and methods (notation, algorithm, simulation design, ADEMP)
- [x] Results (simulation, ablation, worked examples, benchmark)
- [x] Discussion (practical guidance, limitations, conclusions)
- [x] Data availability statement
- [x] Ethics statement
- [x] Author contributions (CRediT format)
- [x] Funding statement
- [x] Competing interests
- [x] Acknowledgments
- [x] References
- [x] Figures + captions
- [x] Supporting information (S1-S8 Tables, S1-S2 Figs, code files, output tables)

---

## Minor Notes (non-blocking)

1. **data.table version discrepancy**: Manuscript line 193 says `data.table 1.16.4` but `R_session_info.txt` says `data.table 1.18.0`. Verify which version was actually used during the simulation run and reconcile if needed.

2. **`.zenodo.json` author mismatch**: The `.zenodo.json` file (untracked) was created for F1000 and lists 5 creators with "Ahmad, Mahmood" as first author. For the PLOS ONE GRMA paper (sole author "Mahmood Ul Hassan"), either update this file or configure Zenodo manually.

3. **F1000 files (untracked)**: Several F1000-related files are untracked in git (`F1000_Concrete_Submission_Fixes.md`, `F1000_MultiPersona_Review.md`, `F1000_Reviewer_Rerun_Manifest.md`, `f1000_artifacts/`). These are separate from the PLOS ONE submission and can remain untracked.
