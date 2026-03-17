# Grey Relational Meta-Analysis (GRMA)

## Installation
Use the dependency files in this directory (for example `requirements.txt`, `environment.yml`, `DESCRIPTION`, or equivalent project-specific files) to create a clean local environment before running analyses.
Document any package-version mismatch encountered during first run.

Reproducibility capsule for: *"Grey Relational Meta-Analysis with an Explicit Redescending Effect Guard and Bootstrap Inference"*

## Overview

GRMA is a robust pooling estimator for meta-analysis that uses grey-relational similarity in a two-feature space (effect size and log-precision) with an explicit Tukey bisquare redescending effect guard and full BCa bootstrap inference.

## System Requirements

- **R** >= 4.3.0 with packages: `metafor`, `data.table`
- **Python** >= 3.10 with packages: `numpy`, `scipy`, `matplotlib` (see `requirements.txt`)
- **OS**: Tested on Windows 11; should work on Linux/macOS
- **RAM**: 4 GB minimum (8 GB recommended for simulation)
- **Disk**: ~50 MB for full capsule with outputs

## Quick Start

```bash
# Install Python dependencies
pip install -r requirements.txt

# Run Python-side reproducibility (tables + verification): ~2 min
python reproduce_all_v8.py

# Run edge case tests (R): ~10 sec
Rscript test_edge_cases.R

# Run applied examples (R, requires metafor): ~1 min
Rscript applied_examples.R

# Generate figures (Python, requires matplotlib): ~30 sec
python generate_figures.py

# Verify manuscript numbers against tables: ~5 sec
python verify_manuscript_numbers.py
```

## File Structure

### Core Implementation
| File | Language | Description |
|------|----------|-------------|
| `grma_meta.R` | R | GRMA estimator (exported as `grma_meta()`, `grma_loo()`, `grma_valley()`) |
| `grey_meta_v8.py` | Python | GRMA estimator (class `GRMA`, standalone) |

### Scripts
| File | Runtime | Description |
|------|---------|-------------|
| `reproduce_all_v8.py` | ~2 min | Reproduces all Python-side tables |
| `applied_examples.R` | ~1 min | Five applied examples (BCG, Morris, CD002042) |
| `run_grma_simulation.R` | ~8-12 h | Full 2000-rep x 25-scenario simulation (R) |
| `run_pairwise70_benchmark.py` | ~30 min | Pairwise70 Cochrane benchmark (4,572 analyses) |
| `generate_figures.py` | ~30 sec | Publication figures (Figs 1-5) |
| `verify_manuscript_numbers.py` | ~5 sec | Checks manuscript claims against CSVs |
| `cross_validate.py` / `.R` | ~10 sec | R-Python cross-validation |
| `test_edge_cases.R` | ~10 sec | Edge case tests with assertions |

### Output Tables
See `Table_artifact_map_v8.csv` for the full mapping from manuscript items to scripts and output files.

### Datasets
- **BCG** (Colditz et al. 1994): Embedded in `grey_meta_v8.py` and `applied_examples.R`
- **Morris** (2008): Embedded in both implementations
- **CD002042** (Carson et al. 2021): 2x2 tables embedded in `CD002042_embedded_data.R` (extracted from published Cochrane review)

### Development Scripts (not part of the reproducibility capsule)
- `dev_gen_tables_from_checkpoint.R` — Generate tables from intermediate simulation checkpoint
- `dev_resume_simulation.R` — Resume simulation from checkpoint

These are development utilities used during the simulation phase and are not required for reproduction.

### External Dependencies
The simulation (`run_grma_simulation.R`) requires the **Pairwise70** framework, which is not included in this capsule. The applied examples are fully self-contained (all datasets embedded). The pre-computed results (all CSV/JSON tables) are provided for verification. The Pairwise70 summary statistics are deposited in `pairwise70_benchmark_grma/`. Contact the corresponding author for full Pairwise70 access.

## Cross-Validation

The R implementation is authoritative. The Python implementation is cross-validated against R to 10^-10 tolerance on point estimates and weights (BCG dataset). Bootstrap CIs differ across languages due to different RNG sequences; R is authoritative for inference.

## License

MIT License. See `LICENSE`.

## Citation

Hassan MU. Grey Relational Meta-Analysis: A Robust Pooling Method with Redescending Effect Guard. PLOS ONE. 2026. [DOI to be added upon publication]
