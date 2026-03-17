# GRMA Reproducibility Makefile
# Usage: make all  (runs everything reproducible without Pairwise70)
# Usage: make full (requires Pairwise70 framework)

PYTHON = python
RSCRIPT = Rscript

.PHONY: all full test tables figures verify crossval clean

# Default: reproduce everything that doesn't need external dependencies
all: test crossval tables figures verify
	@echo ""
	@echo "=== All reproducible steps completed ==="

# Full reproduction (requires Pairwise70 framework)
full: all applied simulation benchmark
	@echo ""
	@echo "=== Full reproduction completed (including Pairwise70) ==="

# Edge case tests (R)
test:
	@echo "--- Running edge case tests ---"
	cd "$(CURDIR)" && $(RSCRIPT) test_edge_cases.R

# Python-side tables (standalone, no external deps)
tables:
	@echo "--- Reproducing Python tables ---"
	$(PYTHON) reproduce_all_v8.py

# Generate publication figures
figures:
	@echo "--- Generating figures ---"
	$(PYTHON) generate_figures.py

# Verify manuscript numbers
verify:
	@echo "--- Verifying manuscript numbers ---"
	$(PYTHON) verify_manuscript_numbers.py

# Applied examples (requires metafor; CD002042 requires Pairwise70 data)
applied:
	@echo "--- Running applied examples ---"
	$(RSCRIPT) applied_examples.R

# Full simulation (requires Pairwise70 framework; ~8-12 hours)
simulation:
	@echo "--- Running full simulation (this takes 8-12 hours) ---"
	$(RSCRIPT) run_grma_simulation.R

# Pairwise70 benchmark (requires Pairwise70 data; ~30 min)
benchmark:
	@echo "--- Running Pairwise70 benchmark ---"
	$(PYTHON) run_pairwise70_benchmark.py

# Cross-validation (R then Python)
crossval:
	@echo "--- Cross-validating R vs Python ---"
	$(RSCRIPT) cross_validate.R
	$(PYTHON) cross_validate.py

clean:
	@echo "--- Cleaning generated files ---"
	rm -f Table_*_py.csv Table_defaults_v8.csv
	rm -f Fig*.png
	rm -rf __pycache__
