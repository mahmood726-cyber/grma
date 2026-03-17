#!/usr/bin/env Rscript
# ============================================================================
# GRMA v8 — Applied Examples (BCG + Morris + CD002042)
# ============================================================================
# Produces: Table_applied_v8.csv, Table_diagnostics_v8.csv, Table_loo_v8.csv
# ============================================================================

library(metafor)

# Source GRMA implementation (use local copy; falls back to Pairwise70 path)
if (file.exists("grma_meta.R")) {
  source("grma_meta.R")
} else {
  stop("grma_meta.R not found in working directory. Please run from the capsule root.")
}

cat("=== GRMA v8 Applied Examples ===\n\n")

# ============================================================================
# 1. Load datasets
# ============================================================================

# --- BCG vaccine (13 trials, log RR) ---
bcg_dat <- data.frame(
  tpos = c(6, 29, 11, 248, 8, 10, 180, 12, 505, 29, 17, 186, 5),
  tneg = c(300-6, 274-29, 231-11, 12619-248, 2545-8, 619-10, 1541-180,
           1716-12, 87886-505, 7470-29, 1699-17, 50448-186, 2493-5),
  cpos = c(11, 45, 29, 462, 10, 8, 372, 47, 499, 45, 65, 141, 3),
  cneg = c(300-11, 279-45, 220-29, 13536-462, 629-10, 2000-8, 1451-372,
           1665-47, 87892-499, 7232-45, 1600-65, 27197-141, 2338-3)
)
bcg_es <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
                 data = bcg_dat)
bcg_yi <- as.numeric(bcg_es$yi)
bcg_vi <- as.numeric(bcg_es$vi)

# --- Morris pretest-posttest (8 studies, SMD change diff) ---
morris_yi <- c(1.2, 0.99, 0.78, 0.54, 1.10, 0.85, 1.45, 0.62)
morris_vi <- c(0.10, 0.08, 0.12, 0.15, 0.09, 0.11, 0.07, 0.14)

# --- CD002042 Blood Transfusion ---
# Use embedded 2x2 tables (included in capsule for reproducibility).
# Source: Carson JL et al. Cochrane Database Syst Rev. 2021;12:CD002042.
source("CD002042_embedded_data.R")

extract_cd002042_embedded <- function(entry) {
  es <- escalc(measure = "OR",
               ai = entry$ai, n1i = entry$n1i,
               ci = entry$ci, n2i = entry$n2i)
  es <- es[is.finite(es$yi) & is.finite(es$vi),]
  list(yi = as.numeric(es$yi), vi = as.numeric(es$vi))
}

cd_mortality <- extract_cd002042_embedded(cd002042_embedded$mortality)
cd_chf <- extract_cd002042_embedded(cd002042_embedded$chf)
cd_transfusion <- extract_cd002042_embedded(cd002042_embedded$transfusion)

# ============================================================================
# 2. Fit all methods
# ============================================================================

fit_all_methods <- function(yi, vi, dataset_label, seed = 20260213) {
  k <- length(yi)
  results <- list()

  # IV Fixed Effect
  fit_fe <- rma(yi, vi, method = "FE")
  results[["IV_FE"]] <- list(
    estimate = coef(fit_fe), se = fit_fe$se,
    ci_lo = fit_fe$ci.lb, ci_hi = fit_fe$ci.ub
  )

  # DL RE
  fit_dl <- rma(yi, vi, method = "DL")
  results[["DL_RE"]] <- list(
    estimate = coef(fit_dl), se = fit_dl$se,
    ci_lo = fit_dl$ci.lb, ci_hi = fit_dl$ci.ub
  )

  # REML RE
  fit_reml <- rma(yi, vi, method = "REML")
  results[["REML_RE"]] <- list(
    estimate = coef(fit_reml), se = fit_reml$se,
    ci_lo = fit_reml$ci.lb, ci_hi = fit_reml$ci.ub
  )

  # HK REML
  fit_hk <- rma(yi, vi, method = "REML", test = "knha")
  results[["HK_REML"]] <- list(
    estimate = coef(fit_hk), se = fit_hk$se,
    ci_lo = fit_hk$ci.lb, ci_hi = fit_hk$ci.ub
  )

  # Huber IV (via robumeta-style IRLS)
  huber_iv <- function(yi, vi, huber_k = 1.345, max_iter = 50, tol = 1e-6) {
    tau2 <- fit_reml$tau2
    est <- coef(fit_reml)
    for (iter in 1:max_iter) {
      tv <- vi + tau2
      r <- (yi - est) / sqrt(tv)
      hw <- ifelse(abs(r) <= huber_k, 1, huber_k / abs(r))
      w <- hw / tv
      est_new <- sum(w * yi) / sum(w)
      if (abs(est_new - est) < tol) { est <- est_new; break }
      est <- est_new
    }
    se <- sqrt(sum((w / sum(w))^2 * tv))
    z <- qnorm(0.975)
    list(estimate = est, se = se, ci_lo = est - z * se, ci_hi = est + z * se)
  }
  results[["HuberIV"]] <- huber_iv(yi, vi)

  # t-RE (df=4, IRLS)
  t_re <- function(yi, vi, df_t = 4, max_iter = 100, tol = 1e-8) {
    tau2 <- fit_reml$tau2
    est <- coef(fit_reml)
    for (iter in 1:max_iter) {
      tv <- vi + tau2
      r <- (yi - est) / sqrt(tv)
      tw <- (df_t + 1) / (df_t + r^2)
      w <- tw / tv
      est_new <- sum(w * yi) / sum(w)
      if (abs(est_new - est) < tol) { est <- est_new; break }
      est <- est_new
    }
    se <- sqrt(1 / sum(w))
    z <- qnorm(0.975)
    list(estimate = est, se = se, ci_lo = est - z * se, ci_hi = est + z * se)
  }
  results[["tRE_df4"]] <- t_re(yi, vi)

  # GRMA (with guard)
  set.seed(seed)
  grma_res <- grma_meta(yi, vi, effect_guard = TRUE, n_boot = 999, bca = TRUE,
                        seed = seed)
  results[["GRMA"]] <- list(
    estimate = grma_res$estimate, se = grma_res$se,
    ci_lo = grma_res$ci_lb, ci_hi = grma_res$ci_ub,
    ci_lo_bca = grma_res$ci_lb_bca, ci_hi_bca = grma_res$ci_ub_bca
  )

  # GRMA (no guard)
  set.seed(seed)
  grma_ng <- grma_meta(yi, vi, effect_guard = FALSE, n_boot = 999, bca = TRUE,
                       seed = seed)
  results[["GRMA_noguard"]] <- list(
    estimate = grma_ng$estimate, se = grma_ng$se,
    ci_lo = grma_ng$ci_lb, ci_hi = grma_ng$ci_ub,
    ci_lo_bca = grma_ng$ci_lb_bca, ci_hi_bca = grma_ng$ci_ub_bca
  )

  # Compile into data.frame
  rows <- lapply(names(results), function(nm) {
    r <- results[[nm]]
    data.frame(
      dataset = dataset_label,
      method = nm,
      estimate = r$estimate,
      se = r$se,
      ci_lo = r$ci_lo,
      ci_hi = r$ci_hi,
      ci_lo_bca = if (!is.null(r$ci_lo_bca)) r$ci_lo_bca else NA_real_,
      ci_hi_bca = if (!is.null(r$ci_hi_bca)) r$ci_hi_bca else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ============================================================================
# 3. Run all datasets
# ============================================================================

cat("Fitting BCG...\n")
res_bcg <- fit_all_methods(bcg_yi, bcg_vi, "BCG (log RR)")
cat("Fitting Morris...\n")
res_morris <- fit_all_methods(morris_yi, morris_vi, "Morris (SMD change diff)")
cat("Fitting CD002042 30-day mortality...\n")
res_mort <- fit_all_methods(cd_mortality$yi, cd_mortality$vi,
                            "CD002042: 30-day mortality (log OR)")
cat("Fitting CD002042 CHF...\n")
res_chf <- fit_all_methods(cd_chf$yi, cd_chf$vi,
                           "CD002042: CHF (log OR)")
cat("Fitting CD002042 transfusion...\n")
res_trans <- fit_all_methods(cd_transfusion$yi, cd_transfusion$vi,
                             "CD002042: Transfusion exposure (log OR)")

table_applied <- rbind(res_bcg, res_morris, res_mort, res_chf, res_trans)
write.csv(table_applied, "Table_applied_v8.csv",
          row.names = FALSE)
cat("Wrote Table_applied_v8.csv\n")

# ============================================================================
# 4. GRMA Diagnostics
# ============================================================================

get_diagnostics <- function(yi, vi, dataset_label, seed = 20260213) {
  set.seed(seed)
  res <- grma_meta(yi, vi, n_boot = 999, bca = TRUE, seed = seed)
  vd <- grma_valley(yi, res$estimate, seed = seed)

  data.frame(
    dataset = dataset_label,
    k = res$k,
    anchor_y = res$anchor_y,
    anchor_p = res$anchor_p,
    zeta = 0.5,
    w_max = res$w_max,
    n_eff = res$n_eff,
    valley_flag = vd$valley_flag,
    valley_p = vd$valley_p,
    stringsAsFactors = FALSE
  )
}

diag_bcg <- get_diagnostics(bcg_yi, bcg_vi, "BCG (log RR)")
diag_morris <- get_diagnostics(morris_yi, morris_vi, "Morris (SMD change diff)")
diag_mort <- get_diagnostics(cd_mortality$yi, cd_mortality$vi,
                             "CD002042: 30-day mortality (log OR)")
diag_chf <- get_diagnostics(cd_chf$yi, cd_chf$vi,
                            "CD002042: CHF (log OR)")
diag_trans <- get_diagnostics(cd_transfusion$yi, cd_transfusion$vi,
                              "CD002042: Transfusion exposure (log OR)")

table_diag <- rbind(diag_bcg, diag_morris, diag_mort, diag_chf, diag_trans)
write.csv(table_diag, "Table_diagnostics_v8.csv",
          row.names = FALSE)
cat("Wrote Table_diagnostics_v8.csv\n")

# ============================================================================
# 5. Leave-One-Out Influence
# ============================================================================

get_loo <- function(yi, vi, dataset_label) {
  loo <- grma_loo(yi, vi)
  data.frame(
    dataset = dataset_label,
    k = length(yi),
    est_full = loo$est_full,
    max_abs_delta_est = loo$max_abs_delta_est,
    idx_max_shift = loo$idx_max_shift,
    y_idx = loo$y_at_max,
    v_idx = loo$v_at_max,
    maxw_full = loo$w_max_full,
    max_abs_delta_maxw = loo$max_abs_delta_maxw,
    idx_maxw_shift = loo$idx_maxw_shift,
    stringsAsFactors = FALSE
  )
}

cat("Computing LOO for BCG...\n")
loo_bcg <- get_loo(bcg_yi, bcg_vi, "BCG (log RR)")
cat("Computing LOO for Morris...\n")
loo_morris <- get_loo(morris_yi, morris_vi, "Morris (SMD change diff)")
cat("Computing LOO for CD002042 30-day mortality...\n")
loo_mort <- get_loo(cd_mortality$yi, cd_mortality$vi,
                    "CD002042: 30-day mortality (log OR)")
cat("Computing LOO for CD002042 CHF...\n")
loo_chf <- get_loo(cd_chf$yi, cd_chf$vi, "CD002042: CHF (log OR)")
cat("Computing LOO for CD002042 transfusion...\n")
loo_trans <- get_loo(cd_transfusion$yi, cd_transfusion$vi,
                     "CD002042: Transfusion exposure (log OR)")

table_loo <- rbind(loo_bcg, loo_morris, loo_mort, loo_chf, loo_trans)
write.csv(table_loo, "Table_loo_v8.csv",
          row.names = FALSE)
cat("Wrote Table_loo_v8.csv\n")

# ============================================================================
# 6. Zeta Sensitivity (BCG)
# ============================================================================

cat("Computing zeta sensitivity for BCG...\n")
zeta_vals <- seq(0.1, 0.9, by = 0.1)
zeta_rows <- lapply(zeta_vals, function(z) {
  res <- grma_meta(bcg_yi, bcg_vi, zeta = z, n_boot = 999, bca = TRUE,
                   seed = 20260213)
  data.frame(
    zeta = z,
    estimate = res$estimate,
    ci_lo_pct = res$ci_lb,
    ci_hi_pct = res$ci_ub,
    ci_lo_bca = res$ci_lb_bca,
    ci_hi_bca = res$ci_ub_bca,
    w_max = res$w_max,
    n_eff = res$n_eff,
    stringsAsFactors = FALSE
  )
})
table_zeta <- do.call(rbind, zeta_rows)
write.csv(table_zeta, "Table_zeta_sensitivity_bcg_v8.csv",
          row.names = FALSE)
cat("Wrote Table_zeta_sensitivity_bcg_v8.csv\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n=== Applied examples complete ===\n")
cat("Files written:\n")
cat("  Table_applied_v8.csv\n")
cat("  Table_diagnostics_v8.csv\n")
cat("  Table_loo_v8.csv\n")
cat("  Table_zeta_sensitivity_bcg_v8.csv\n")

# Print quick summary
cat("\n--- BCG GRMA ---\n")
bcg_grma <- table_applied[table_applied$dataset == "BCG (log RR)" &
                             table_applied$method == "GRMA",]
cat(sprintf("  Estimate: %.6f  [%.4f, %.4f] (pct)  [%.4f, %.4f] (BCa)\n",
            bcg_grma$estimate, bcg_grma$ci_lo, bcg_grma$ci_hi,
            bcg_grma$ci_lo_bca, bcg_grma$ci_hi_bca))

cat("\n--- CD002042 Mortality GRMA ---\n")
mort_grma <- table_applied[table_applied$dataset == "CD002042: 30-day mortality (log OR)" &
                              table_applied$method == "GRMA",]
cat(sprintf("  Estimate: %.6f  [%.4f, %.4f] (pct)\n",
            mort_grma$estimate, mort_grma$ci_lo, mort_grma$ci_hi))
