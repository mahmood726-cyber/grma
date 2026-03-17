#!/usr/bin/env Rscript
# Cross-validation: output BCG yi/vi from R for Python comparison
library(metafor)
# Source local copy (capsule); fall back to Pairwise70 if needed
if (file.exists("grma_meta.R")) {
  source("grma_meta.R")
} else {
  stop("grma_meta.R not found in working directory. Please run from the capsule root.")
}

# BCG data (same 2x2 tables)
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
yi <- as.numeric(bcg_es$yi)
vi <- as.numeric(bcg_es$vi)

# Print yi/vi for Python to use
cat("yi_r = [", paste(sprintf("%.15g", yi), collapse = ", "), "]\n")
cat("vi_r = [", paste(sprintf("%.15g", vi), collapse = ", "), "]\n")

# Fit GRMA (point estimate only, no bootstrap randomness)
res <- grma_meta(yi, vi, n_boot = 50, bca = FALSE, seed = 99999)
cat(sprintf("\nR GRMA estimate: %.15g\n", res$estimate))
cat(sprintf("R w_max:         %.15g\n", res$w_max))
cat(sprintf("R n_eff:         %.15g\n", res$n_eff))
cat(sprintf("R anchor_y:      %.15g\n", res$anchor_y))
cat(sprintf("R anchor_p:      %.15g\n", res$anchor_p))
cat("R weights:", paste(sprintf("%.12g", res$weights), collapse = ", "), "\n")
