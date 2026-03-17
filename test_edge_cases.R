# Edge case tests for GRMA with programmatic assertions
# Usage: Rscript test_edge_cases.R
# Exit code 0 = all pass, 1 = failure

source("grma_meta.R")
library(metafor)

pass_count <- 0
fail_count <- 0

assert_equal <- function(actual, expected, tol = 1e-6, label = "") {
  if (is.na(actual) && is.na(expected)) {
    pass_count <<- pass_count + 1
    cat(sprintf("  PASS: %s\n", label))
    return(invisible(TRUE))
  }
  if (is.na(actual) || is.na(expected) || abs(actual - expected) > tol) {
    fail_count <<- fail_count + 1
    cat(sprintf("  FAIL: %s (expected %.8f, got %.8f)\n", label, expected, actual))
    return(invisible(FALSE))
  }
  pass_count <<- pass_count + 1
  cat(sprintf("  PASS: %s\n", label))
  invisible(TRUE)
}

assert_true <- function(cond, label = "") {
  if (!isTRUE(cond)) {
    fail_count <<- fail_count + 1
    cat(sprintf("  FAIL: %s\n", label))
    return(invisible(FALSE))
  }
  pass_count <<- pass_count + 1
  cat(sprintf("  PASS: %s\n", label))
  invisible(TRUE)
}

# ---- Test 1: k=3 (minimum) ----
cat("Test 1: k=3 minimum\n")
yi <- c(0.5, 0.6, 0.7)
vi <- c(0.01, 0.02, 0.03)
res <- grma_meta(yi, vi, n_boot = 50, bca = FALSE, seed = 1)
assert_true(res$method == "GRMA", "k=3 uses GRMA method")
assert_true(is.finite(res$estimate), "k=3 estimate is finite")
assert_true(res$w_max > 0 && res$w_max <= 1, "k=3 w_max in (0,1]")
assert_true(abs(sum(res$weights) - 1.0) < 1e-10, "k=3 weights sum to 1")

# ---- Test 2: k=2 (should fall back to HKSJ) ----
cat("Test 2: k=2 fallback to HKSJ\n")
yi2 <- c(0.5, 0.6)
vi2 <- c(0.01, 0.02)
res2 <- suppressWarnings(grma_meta(yi2, vi2, n_boot = 50, seed = 1))
assert_true(res2$method == "HKSJ_fallback", "k=2 falls back to HKSJ")
assert_true(is.finite(res2$estimate), "k=2 estimate is finite")

# ---- Test 3: All identical yi (was NaN before fix) ----
cat("Test 3: All-identical effects\n")
set.seed(42)
yi3 <- rep(0.5, 10)
vi3 <- runif(10, 0.01, 0.05)
res3 <- grma_meta(yi3, vi3, n_boot = 50, bca = FALSE, seed = 1)
assert_true(is.finite(res3$estimate), "identical yi: estimate is finite")
assert_equal(res3$estimate, 0.5, tol = 1e-10, label = "identical yi: estimate = 0.5")
assert_true(res3$method == "GRMA", "identical yi: uses GRMA method")

# ---- Test 3b: All identical yi AND vi ----
cat("Test 3b: All-identical effects AND variances\n")
yi3b <- rep(1.0, 5)
vi3b <- rep(0.1, 5)
res3b <- grma_meta(yi3b, vi3b, n_boot = 50, bca = FALSE, seed = 1)
assert_true(is.finite(res3b$estimate), "all-identical: estimate is finite")
assert_equal(res3b$estimate, 1.0, tol = 1e-10, label = "all-identical: estimate = 1.0")
assert_true(all(abs(res3b$weights - 0.2) < 1e-10), "all-identical: equal weights")

# ---- Test 4: Single extreme outlier ----
cat("Test 4: Single extreme outlier\n")
yi4 <- c(rep(0.5, 9), 10.0)
vi4 <- rep(0.02, 10)
res4 <- grma_meta(yi4, vi4, n_boot = 50, bca = FALSE, seed = 1)
assert_true(is.finite(res4$estimate), "outlier: estimate is finite")
assert_true(res4$weights[10] < 0.01, "outlier: study 10 weight near zero (guarded)")
assert_true(res4$guard_values[10] < 1e-10, "outlier: guard zeroed study 10")

# ---- Test 5: No guard ----
cat("Test 5: No guard\n")
res5 <- grma_meta(yi4, vi4, effect_guard = FALSE, n_boot = 50, bca = FALSE, seed = 1)
assert_true(res5$method == "GRMA_noguard", "no guard: correct method label")
assert_true(res5$estimate > res4$estimate, "no guard: estimate pulled toward outlier")
assert_true(is.null(res5$guard_values), "no guard: guard_values is NULL")

# ---- Test 6: tau2=0 homogeneous ----
cat("Test 6: Homogeneous studies\n")
set.seed(99)
yi6 <- rnorm(20, mean = 0.3, sd = 0.1)
vi6 <- rep(0.01, 20)
res6 <- grma_meta(yi6, vi6, n_boot = 50, bca = FALSE, seed = 1)
assert_true(is.finite(res6$estimate), "homogeneous: estimate is finite")
assert_true(res6$n_eff > 15, "homogeneous: n_eff high (near k=20)")

# ---- Summary ----
cat(sprintf("\n=== %d PASSED, %d FAILED ===\n", pass_count, fail_count))
if (fail_count > 0) {
  quit(status = 1)
} else {
  cat("All edge case tests passed.\n")
}
