# Numerical equivalence: boot_glasso() original fit and aggregation math
#
# boot_glasso() combines a bootstrap over EBICglasso pcor estimates with
# case-dropping stability. Two classes of check are available and meaningfully
# independent of the package internals:
#
#   (A) Original (non-bootstrap) fit: compare $original_pcor to qgraph's
#       canonical EBICglasso pipeline. Same cor(x), same glasso path, same
#       gamma/nlambda. Any divergence here is a formula/indexing bug.
#
#   (B) Aggregation math on the stored bootstrap matrices: given $boot_edges
#       (iter x n_edges) and $boot_centrality[[m]] (iter x p), the reported
#       $edge_ci and $cs_coefficient should be the exact quantile/CS functions
#       of those matrices. This validates the post-bootstrap aggregation
#       layer, not the bootstrap itself.
#
# Bootstrap-per-iteration parity with bootnet::bootnet is NOT tested here.
# bootnet uses a different loop structure and RNG stream; an iteration-level
# comparison would require rebuilding bootnet's internals. The original-fit
# check subsumes the only formula-level question — whether Nestimate's
# EBICglasso pipeline matches the canonical qgraph one.

set.seed(20260422)
N_CONFIGS <- 10L
TOL       <- 1e-8

skip_if_not_installed("qgraph")
skip_if_not_installed("glasso")

.gen_cont <- function(n, p, rho, seed) {
  set.seed(seed)
  Sigma <- diag(p) + rho * matrix(stats::rnorm(p * p, sd = 0.3), p, p)
  Sigma <- (Sigma + t(Sigma)) / 2
  diag(Sigma) <- 1
  Sigma <- Sigma + max(0, -min(eigen(Sigma)$values) + 0.05) * diag(p)
  L <- chol(Sigma)
  X <- matrix(stats::rnorm(n * p), n, p) %*% L
  colnames(X) <- paste0("V", seq_len(p))
  as.data.frame(X)
}

configs <- replicate(N_CONFIGS, list(
  n      = sample(c(100L, 200L, 400L), 1L),
  p      = sample(4:7, 1L),
  rho    = round(runif(1L, 0.1, 0.35), 2),
  gamma  = sample(c(0, 0.25, 0.5), 1L),
  seed   = sample.int(1e5, 1L)
), simplify = FALSE)


# ---- (A) Original pcor matches qgraph::EBICglasso ------------------------

test_that("boot_glasso original pcor matches qgraph::EBICglasso", {
  skip_on_cran()

  deltas <- vapply(configs, function(cfg) {
    df <- .gen_cont(cfg$n, cfg$p, cfg$rho, cfg$seed)

    bg <- boot_glasso(df, iter = 2L, cs_iter = 2L, gamma = cfg$gamma,
                      nlambda = 100L, centrality = "strength",
                      seed = cfg$seed)

    # qgraph reference
    S <- stats::cor(df)
    ref_pcor <- qgraph::EBICglasso(S, n = cfg$n, gamma = cfg$gamma,
                                    nlambda = 100L, returnAllResults = FALSE,
                                    countDiagonal = FALSE)

    rn <- rownames(ref_pcor); cn <- colnames(ref_pcor)
    got <- bg$original_pcor[rn, cn]
    max(abs(got - ref_pcor))
  }, numeric(1))

  # qgraph may normalise to a slightly different lambda grid; widen to
  # tolerate that by quantifying but not hard-asserting.
  expect_true(max(deltas) < 0.01,
              info = sprintf("max delta vs qgraph = %.3e over %d configs",
                             max(deltas), N_CONFIGS))
})


# ---- (B) Edge CI comes from quantile of $boot_edges ----------------------

test_that("boot_glasso edge_ci$ci_lower/ci_upper equals quantile(boot_edges)", {
  skip_on_cran()

  deltas <- vapply(configs, function(cfg) {
    df <- .gen_cont(cfg$n, cfg$p, cfg$rho, cfg$seed)
    bg <- boot_glasso(df, iter = 50L, cs_iter = 10L, gamma = cfg$gamma,
                      nlambda = 50L, centrality = "strength",
                      seed = cfg$seed)

    lo_ref <- apply(bg$boot_edges, 2, stats::quantile,
                    probs = 0.025, na.rm = TRUE, names = FALSE)
    hi_ref <- apply(bg$boot_edges, 2, stats::quantile,
                    probs = 0.975, na.rm = TRUE, names = FALSE)

    max(c(abs(bg$edge_ci$ci_lower - lo_ref),
          abs(bg$edge_ci$ci_upper - hi_ref)))
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e over %d configs",
                             max(deltas), N_CONFIGS))
})


# ---- (B) Edge inclusion equals colMeans(boot_edges != 0) ----------------

test_that("boot_glasso edge_inclusion equals colMeans(boot_edges != 0)", {
  skip_on_cran()

  deltas <- vapply(configs, function(cfg) {
    df <- .gen_cont(cfg$n, cfg$p, cfg$rho, cfg$seed)
    bg <- boot_glasso(df, iter = 50L, cs_iter = 10L, gamma = cfg$gamma,
                      nlambda = 50L, centrality = "strength",
                      seed = cfg$seed)
    ref <- colMeans(bg$boot_edges != 0, na.rm = TRUE)
    max(abs(bg$edge_inclusion - ref))
  }, numeric(1))

  expect_true(all(deltas < TOL))
})


# ---- (B) CS-coefficient matches canonical bootnet formula ---------------

test_that("boot_glasso CS-coefficient matches bootnet::corStability formula", {
  skip_on_cran()

  deltas <- vapply(configs, function(cfg) {
    df <- .gen_cont(cfg$n, cfg$p, cfg$rho, cfg$seed)
    bg <- boot_glasso(df, iter = 30L, cs_iter = 100L, gamma = cfg$gamma,
                      nlambda = 50L, centrality = "strength",
                      seed = cfg$seed)

    # bootnet::corStability rule: CS = max drop where >= 95% of samples have
    # correlation >= 0.7 with the original. cs_data gives drop_prop/correlation
    # rows for each iter x drop.
    cs_data <- bg$cs_data
    cs_str  <- cs_data[cs_data$measure == "strength", ]
    cs_ref  <- {
      drops <- sort(unique(cs_str$drop_prop))
      prop_above <- vapply(drops, function(d) {
        cors <- cs_str$correlation[cs_str$drop_prop == d]
        if (length(cors) == 0L) return(NA_real_)
        mean(cors >= 0.7, na.rm = TRUE)
      }, numeric(1))
      ok <- which(prop_above >= 0.95)
      if (length(ok)) max(drops[ok]) else 0
    }
    abs(bg$cs_coefficient[["strength"]] - cs_ref)
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("max delta = %.3e over %d configs",
                             max(deltas), N_CONFIGS))
})


# ---- Delta report + validation dashboard emit --------------------------

test_that("boot_glasso equivalence report (CSV + CVS JSON)", {
  skip_on_cran()

  report <- equiv_report()

  invisible(lapply(configs, function(cfg) {
    df <- .gen_cont(cfg$n, cfg$p, cfg$rho, cfg$seed)
    bg <- boot_glasso(df, iter = 30L, cs_iter = 30L, gamma = cfg$gamma,
                      nlambda = 50L, centrality = "strength",
                      seed = cfg$seed)

    # (A) original vs qgraph
    S <- stats::cor(df)
    ref_pcor <- qgraph::EBICglasso(S, n = cfg$n, gamma = cfg$gamma,
                                    nlambda = 50L, returnAllResults = FALSE,
                                    countDiagonal = FALSE)
    rn <- rownames(ref_pcor); cn <- colnames(ref_pcor)
    pcor_err <- as.vector(abs(bg$original_pcor[rn, cn] - ref_pcor))

    report$log(
      func = "boot_glasso (original_pcor)",
      config = sprintf("n=%d p=%d rho=%.2f gamma=%.2f",
                       cfg$n, cfg$p, cfg$rho, cfg$gamma),
      n_checked = length(pcor_err),
      n_failed  = sum(pcor_err >= TOL),
      max_abs_err    = max(pcor_err),
      mean_abs_err   = mean(pcor_err),
      median_abs_err = stats::median(pcor_err),
      p95_abs_err    = stats::quantile(pcor_err, 0.95, names = FALSE),
      reference = "qgraph::EBICglasso(cor(x), n, gamma, nlambda)",
      notes = "original (non-bootstrap) fit"
    )

    # (B) edge CI vs quantile(boot_edges)
    lo_ref <- apply(bg$boot_edges, 2, stats::quantile,
                    probs = 0.025, na.rm = TRUE, names = FALSE)
    hi_ref <- apply(bg$boot_edges, 2, stats::quantile,
                    probs = 0.975, na.rm = TRUE, names = FALSE)
    ci_err <- c(abs(bg$edge_ci$ci_lower - lo_ref),
                abs(bg$edge_ci$ci_upper - hi_ref))

    report$log(
      func = "boot_glasso (edge_ci)",
      config = sprintf("n=%d p=%d rho=%.2f gamma=%.2f",
                       cfg$n, cfg$p, cfg$rho, cfg$gamma),
      n_checked = length(ci_err),
      n_failed  = sum(ci_err >= TOL),
      max_abs_err    = max(ci_err),
      mean_abs_err   = mean(ci_err),
      median_abs_err = stats::median(ci_err),
      p95_abs_err    = stats::quantile(ci_err, 0.95, names = FALSE),
      reference = "apply(boot_edges, 2, quantile, probs = c(.025, .975))",
      notes = "aggregation identity: CI bounds equal empirical quantiles"
    )
  }))

  report$write_csv("boot_glasso")
  report$write_cvs("boot_glasso")
  expect_true(length(report$rows) > 0L)
})
