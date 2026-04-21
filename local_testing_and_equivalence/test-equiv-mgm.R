# Numerical equivalence: Nestimate's "mgm" estimator vs mgm::mgm()
#
# Per-value comparison of the symmetric weighted adjacency matrix across a
# random grid of (seed, n, p_g, p_c) configurations. Defaults match
# mgm::mgm() with lambdaSel = "EBIC", ruleReg = "AND", threshold = "LW",
# overparameterize = FALSE, scale = TRUE.

set.seed(20260413)
N_CONFIGS <- 30L
TOL <- 1e-10

skip_if_not_installed("mgm")
skip_if_not_installed("glmnet")

.gen_mixed <- function(seed, n, p_g, p_c, k = 3) {
  set.seed(seed)
  p <- p_g + p_c
  Z <- matrix(stats::rnorm(n * p), n, p)
  M <- diag(p) + 0.4 * matrix(stats::rnorm(p * p, sd = 0.3), p, p)
  Z <- Z %*% M

  out <- as.data.frame(Z)
  for (j in seq_len(p_c)) {
    col <- p_g + j
    out[[col]] <- as.integer(cut(
      Z[, col],
      breaks = stats::quantile(Z[, col],
                               probs = seq(0, 1, length.out = k + 1)),
      include.lowest = TRUE, labels = FALSE
    ))
  }
  colnames(out) <- paste0("V", seq_len(p))
  out
}

cfgs <- expand.grid(
  seed = seq_len(5),
  n    = c(150, 300),
  p_g  = c(2, 4),
  p_c  = c(0, 2, 3)
)
cfgs <- cfgs[(cfgs$p_g + cfgs$p_c) >= 3, ]
cfgs <- cfgs[seq_len(min(N_CONFIGS, nrow(cfgs))), ]

test_that("mgm estimator matches mgm::mgm() at machine precision", {
  invisible(lapply(seq_len(nrow(cfgs)), function(i) {
    cfg <- cfgs[i, ]
    dat <- .gen_mixed(cfg$seed, cfg$n, cfg$p_g, cfg$p_c)
    type  <- c(rep("g", cfg$p_g), rep("c", cfg$p_c))
    level <- c(rep(1L, cfg$p_g), rep(3L, cfg$p_c))

    nest_out <- Nestimate:::.estimator_mgm(
      data      = dat,
      type      = type,
      level     = level,
      lambdaGam = 0.25,
      ruleReg   = "AND",
      threshold = "LW",
      scale     = TRUE
    )
    nest_wadj <- nest_out$matrix

    ref <- suppressWarnings(suppressMessages(mgm::mgm(
      data       = as.matrix(dat),
      type       = type,
      level      = level,
      lambdaSel  = "EBIC",
      lambdaGam  = 0.25,
      ruleReg    = "AND",
      threshold  = "LW",
      overparameterize = FALSE,
      scale      = TRUE,
      pbar       = FALSE,
      signInfo   = FALSE,
      warnings   = FALSE
    )))
    ref_wadj <- ref$pairwise$wadj

    expect_equal(
      unname(nest_wadj), unname(ref_wadj),
      tolerance = TOL,
      label = sprintf("cfg(seed=%d,n=%d,p_g=%d,p_c=%d)",
                       cfg$seed, cfg$n, cfg$p_g, cfg$p_c)
    )
  }))
})


test_that("build_network dispatches method = 'mgm' to .estimator_mgm", {
  set.seed(7)
  dat <- .gen_mixed(1, 200, 3, 2)
  # Make categorical columns factors so the auto-detector picks them up
  dat$V4 <- factor(dat$V4)
  dat$V5 <- factor(dat$V5)
  net <- build_network(dat, method = "mgm")
  expect_true(inherits(net, "netobject"))
  expect_equal(dim(net$weights), c(5, 5))
  expect_true(isSymmetric(net$weights))
})


test_that("aliases 'mixed' and 'mixed_graphical' resolve to mgm", {
  expect_equal(Nestimate:::.resolve_method_alias("mixed"),           "mgm")
  expect_equal(Nestimate:::.resolve_method_alias("mixed_graphical"), "mgm")
})
