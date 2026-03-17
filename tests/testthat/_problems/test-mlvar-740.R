# Extracted from test-mlvar.R:740

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Nestimate", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
.compare_mlvar <- function(seed) {
  d <- simulate_data("mlvar", seed = seed)
  vars <- attr(d, "vars")
  nv <- length(vars)

  # Our implementation (no standardization to match mlVAR scale=FALSE)
  our <- mlvar(d, vars = vars, id = "id", day = "day", beep = "beep",
               standardize = FALSE)

  # mlVAR: lmer with fixed temporal + fixed contemporaneous, no scaling
  ref <- suppressWarnings(mlVAR::mlVAR(
    d, vars = vars, idvar = "id", dayvar = "day", beepvar = "beep",
    lags = 1, estimator = "lmer", temporal = "fixed",
    contemporaneous = "fixed", scale = FALSE, verbose = FALSE
  ))

  ref_B <- ref$results$Beta$mean[, , 1]

  # P-values
  our_pvals <- matrix(NA, nv, nv)
  for (k in seq_len(nv)) our_pvals[k, ] <- our$coefs[[k]]$p
  ref_pvals <- ref$results$Beta$P[, , 1]

  # Residual correlations
  prepared <- Nestimate:::.mlvar_prepare_data(d, vars, "id", "day", "beep")
  lag_result <- Nestimate:::.mlvar_build_lag_pairs(prepared, vars, "id",
                                                  "day", "beep", 1L)
  centered <- Nestimate:::.mlvar_within_center(lag_result$Y, lag_result$X,
                                              lag_result$id_vec, FALSE)
  n_subjects <- length(unique(lag_result$id_vec))
  temp <- Nestimate:::.mlvar_temporal_ols(centered$Y, centered$X,
                                         n_subjects, vars)
  our_rcor <- cor(temp$residuals)
  ref_rcor <- ref$results$Theta$cor$mean

  # Person mean correlations
  pm <- aggregate(. ~ id, data = d[, c("id", vars)], FUN = mean)
  our_pmcor <- cor(as.matrix(pm[, vars]))
  ref_pmcor <- ref$results$Omega_mu$cor$mean

  list(
    B_max_diff       = max(abs(our$temporal - ref_B)),
    B_cor            = cor(as.vector(our$temporal), as.vector(ref_B)),
    p_max_diff       = max(abs(our_pvals - ref_pvals)),
    rcor_max_diff    = max(abs(our_rcor - ref_rcor)),
    rcor_cor         = cor(our_rcor[upper.tri(our_rcor)],
                           ref_rcor[upper.tri(ref_rcor)]),
    pmcor_max_diff   = max(abs(our_pmcor - ref_pmcor)),
    pmcor_cor        = cor(our_pmcor[upper.tri(our_pmcor)],
                           ref_pmcor[upper.tri(ref_pmcor)]),
    seed             = seed,
    d                = nv
  )
}

# test -------------------------------------------------------------------------
skip_if_not_installed("mlVAR")
seeds <- seq(201, 220)
results <- lapply(seeds, .compare_mlvar)
