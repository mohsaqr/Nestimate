# Equivalence tests: moderated MGM against mgm::mgm(..., moderators) + mgm::condition()
#
# Validates .mgm_estimate_moderated() + condition_moderated() against the
# mgm reference across random configurations spanning:
#   n in {150, 200, 300}, p_g in {3, 4, 5, 6}, p_c in {0, 1, 2}
#   Binary moderator always in the last column.
#
# Each configuration compares condition(mod=0) and condition(mod=1) element-by-element.
# All 30 configs must match at tolerance 1e-10.

skip_if_not_installed("mgm")
skip_if_not_installed("glmnet")

TOL <- 1e-10

.gen_mod_data <- function(seed, n, p_g, p_c = 0L) {
  set.seed(seed)
  Z <- matrix(stats::rnorm(n * p_g), n, p_g)
  Z <- Z + 0.3 * matrix(rowSums(Z), n, p_g)
  G <- sample(0:1, n, replace = TRUE)
  Z[, p_g] <- Z[, p_g] + ifelse(G == 1, 0.6 * Z[, 1], -0.1 * Z[, 1])
  dat <- as.data.frame(Z)
  if (p_c > 0L) {
    for (k in seq_len(p_c)) {
      n_lev <- sample(2:4, 1)
      dat[[paste0("C", k)]] <- sample(seq_len(n_lev), n, replace = TRUE)
    }
  }
  dat$G <- G
  type <- c(rep("g", p_g), rep("c", p_c), "c")
  level <- c(
    rep(1L, p_g),
    if (p_c > 0L)
      vapply(dat[, (p_g + 1):(p_g + p_c), drop = FALSE],
             function(x) length(unique(x)), integer(1))
    else integer(0),
    2L
  )
  colnames(dat) <- paste0("V", seq_len(ncol(dat)))
  list(dat = dat, type = type, level = level, mod_idx = ncol(dat))
}

.run_equiv <- function(seed, n, p_g, p_c) {
  d <- .gen_mod_data(seed, n, p_g, p_c)
  fit <- Nestimate:::.mgm_estimate_moderated(
    d$dat, d$type, d$level, moderator = d$mod_idx,
    lambdaGam = 0.25, ruleReg = "AND", threshold = "LW", scale = TRUE
  )
  ours0 <- Nestimate:::condition_moderated(fit, mod_value = 0)
  ours1 <- Nestimate:::condition_moderated(fit, mod_value = 1)

  ref <- suppressWarnings(suppressMessages(mgm::mgm(
    as.matrix(d$dat), type = d$type, level = d$level,
    moderators = d$mod_idx, lambdaSel = "EBIC", lambdaGam = 0.25,
    ruleReg = "AND", threshold = "LW", overparameterize = FALSE,
    scale = TRUE, pbar = FALSE, signInfo = FALSE, warnings = FALSE
  )))
  ref0 <- mgm::condition(ref, values = stats::setNames(list(0), as.character(d$mod_idx)))$pairwise$wadj
  ref1 <- mgm::condition(ref, values = stats::setNames(list(1), as.character(d$mod_idx)))$pairwise$wadj
  list(d0 = max(abs(ours0 - ref0)), d1 = max(abs(ours1 - ref1)))
}

# 30 baked configs: 10 continuous-only, 10 with 1 categorical, 10 with 2 categorical
set.seed(20260415)
cfgs <- rbind(
  data.frame(seed = sample.int(100000, 10), n = sample(c(150, 200, 300), 10, TRUE),
             p_g = sample(3:6, 10, TRUE), p_c = 0L),
  data.frame(seed = sample.int(100000, 10), n = sample(c(150, 200, 300), 10, TRUE),
             p_g = sample(3:5, 10, TRUE), p_c = 1L),
  data.frame(seed = sample.int(100000, 10), n = sample(c(200, 300), 10, TRUE),
             p_g = sample(3:4, 10, TRUE), p_c = 2L)
)

test_that("condition_moderated matches mgm::condition across 30 configs", {
  for (i in seq_len(nrow(cfgs))) {
    res <- .run_equiv(cfgs$seed[i], cfgs$n[i], cfgs$p_g[i], cfgs$p_c[i])
    expect_true(res$d0 < TOL,
      info = sprintf("cfg %d (seed=%d n=%d p_g=%d p_c=%d) cond(0) delta=%.2e",
                     i, cfgs$seed[i], cfgs$n[i], cfgs$p_g[i], cfgs$p_c[i], res$d0))
    expect_true(res$d1 < TOL,
      info = sprintf("cfg %d (seed=%d n=%d p_g=%d p_c=%d) cond(1) delta=%.2e",
                     i, cfgs$seed[i], cfgs$n[i], cfgs$p_g[i], cfgs$p_c[i], res$d1))
  }
})
