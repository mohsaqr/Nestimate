# Regression tests for the third batch of review fixes (2026-10-26).
testthat::skip_on_cran()

.mc_seqs <- function() {
  data.frame(T1 = c("A", "C", "B", "A"), T2 = c("A", "D", "A", "A"),
             T3 = c("B", "C", "D", "A"), T4 = c("D", "A", "C", "B"))
}
.mc_clusters <- function() list(G1 = c("A", "B"), G2 = c("C", "D"))

test_that("build_mcml() decodes a tna object's integer-coded data", {
  skip_if_not_installed("tna")
  from_tna <- build_mcml(tna::tna(.mc_seqs()), .mc_clusters(), type = "raw")
  from_seq <- build_mcml(.mc_seqs(), .mc_clusters(), type = "raw")
  expect_equal(from_tna$macro$weights, from_seq$macro$weights)
})

test_that("summary.mcml(): no within layers gives NA; else the macro diagonal", {
  with_layers <- summary(build_mcml(.mc_seqs(), .mc_clusters(), type = "raw"))
  m <- build_mcml(.mc_seqs(), .mc_clusters(), type = "raw")
  expect_equal(with_layers$within_total, unname(diag(m$macro$weights)))
  without <- summary(build_mcml(.mc_seqs(), .mc_clusters(), type = "raw",
                                compute_within = FALSE))
  expect_identical(without$cluster, c("G1", "G2"))
  expect_true(all(is.na(without$within_total)))
})

test_that("entropy_trajectory() windows run through time", {
  set.seed(1)
  d <- data.frame(Actor = rep(c("A", "B"), each = 40),
                  Time = c(101:140, 1:40),
                  Action = sample(c("x", "y", "z"), 80, TRUE))
  tr <- entropy_trajectory(d, action = "Action", actor = "Actor",
                           time = "Time", window = 10)$trajectory
  expect_false(is.unsorted(tr$time_start))
  expect_true(all(tr$time_end >= tr$time_start))
  expect_identical(tr$time_start[1], 2L)
})

test_that("Markov measures refuse correlation-type networks", {
  expect_error(transition_entropy(build_network(srl_strategies, "cor")),
               class = "nestimate_not_markov")
  expect_error(passage_time(build_network(srl_strategies, "glasso")),
               class = "nestimate_not_markov")
  expect_error(passage_time(matrix(c(0.5, -0.5, 0.5, 0.5), 2)),
               class = "nestimate_not_markov")
  # a co-occurrence network is a random walk: still accepted
  expect_s3_class(suppressWarnings(transition_entropy(
    build_network(trajectories, method = "co_occurrence"))),
    "net_transition_entropy")
})

test_that("path_dependence() keeps states whose names contain 'NA'", {
  set.seed(1)
  d <- as.data.frame(matrix(sample(c("FINAL", "ANALYSIS", "Plan", "Read"),
                                   1200, TRUE), 200, 6))
  expect_setequal(path_dependence(d, order = 2, min_count = 1)$states,
                  c("ANALYSIS", "FINAL", "Plan", "Read"))
  # windows that touch a real NA are dropped: only A->B and B->C remain...
  k <- .pd_kgram_counts(data.frame(a = c("A", "A"), b = c("B", NA),
                                   c = c("C", "C")), 2)
  expect_identical(as.vector(k), c(1L, 1L))
  expect_identical(names(k), c("A\x01B", "B\x01C"))
})

test_that("sequence_compare() reads numeric states of a netobject as states", {
  set.seed(1)
  d <- as.data.frame(matrix(sample(c(1L, 3L, 5L), 300, TRUE), 60, 5))
  g <- rep(c("X", "Y"), 30)
  from_net <- summary(sequence_compare(build_network(d, "relative"),
                                       group = g, sub = 1, test = "none",
                                       min_freq = 1))
  from_df <- summary(sequence_compare(d, group = g, sub = 1, test = "none",
                                      min_freq = 1))
  expect_setequal(from_net$pattern, c("1", "3", "5"))
  expect_equal(from_net[order(from_net$pattern), c("freq_X", "freq_Y")],
               from_df[order(from_df$pattern), c("freq_X", "freq_Y")],
               ignore_attr = TRUE)
})

test_that("build_mmm(): factor covariates counted by design column", {
  set.seed(11)
  seqs <- as.data.frame(matrix(sample(c("A", "B", "C"), 150 * 12, TRUE),
                               150, 12))
  cov <- data.frame(z = factor(sample(c("a", "b", "c", "d"), 150, TRUE)))
  fit <- suppressWarnings(build_mmm(seqs, k = 3, covariates = cov,
                                    n_starts = 2, seed = 1))
  # transitions 3*3*2 + initial 3*2 + mixing 2 * (intercept + 3 dummies)
  expect_identical(as.integer(fit$n_params), 32L)
})

test_that("build_mmm(): an unknown first state adds no initial-state term", {
  counts <- diag(4)
  init <- c(1L, 2L, NA, NA)
  post <- matrix(c(0.6, 0.4, 0.3, 0.7, 0.5, 0.5, 0.2, 0.8), 4, byrow = TRUE)
  fit <- .mmm_em(counts, init, 2L, 1L, 1e-9, 0.01, 2L, init_posterior = post)
  # hand-computed M-step then E-step for one iteration
  p_un <- crossprod(counts, post) + 0.01
  dim(p_un) <- c(2, 2, 2)
  from_sums <- colSums(p_un)
  dim(p_un) <- c(4, 2)
  p <- p_un / from_sums[rep(1:2, each = 2), ]
  ind <- matrix(0, 4, 2)
  ind[cbind(1:2, 1:2)] <- 1
  i_un <- crossprod(ind, post) + 0.01
  i_p <- i_un / rep(colSums(i_un), each = 2)
  rows <- counts %*% log(p) + rbind(log(i_p[1, ]), log(i_p[2, ]), 0, 0) +
    rep(log(colMeans(post)), each = 4)
  top <- apply(rows, 1, max)
  expect_equal(fit$ll, sum(top + log(rowSums(exp(rows - top)))))
})

test_that("set_state_colors() keeps a netobject_group a list of networks", {
  g <- build_network(group_regulation_long, method = "relative",
                     actor = "Actor", action = "Action", time = "Time",
                     group = "Achiever")
  g2 <- set_state_colors(g, c(plan = "#0072B2"))
  expect_identical(names(g2), names(g))
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  res <- sequence_plot(g2, type = "index")
  expect_identical(res$palette[res$levels == "plan"], "#0072B2")
  expect_s3_class(summary(g2), "data.frame")
})

test_that("plots with a one-sequence panel draw", {
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(sequence_plot(trajectories[1, , drop = FALSE]))
  expect_no_error(sequence_plot(trajectories[1:20, ], type = "index",
                                group = c(rep("A", 19), "B")))
  expect_no_error(sequence_plot(data.frame(s = "a", d = c(1, 1, 2),
                                           a = c("x", "y", "x")),
                                actor = "s", action = "a", periods = "d"))
})

test_that("mosaic_plot() draws networks and tidy data frames; plot_mosaic is gone", {
  g <- build_network(group_regulation_long, method = "relative",
                     actor = "Actor", action = "Action", time = "Time",
                     group = "Achiever")
  expect_s3_class(mosaic_plot(g, seed = 1), "ggplot")
  expect_error(mosaic_plot(g, x_var = "a"), class = "nestimate_bad_input")
  df <- data.frame(group = rep(c("A", "B"), each = 2),
                   state = rep(c("s1", "s2"), 2), count = c(3, 1, 2, 4))
  expect_s3_class(mosaic_plot(df, x_var = "group", y_var = "state",
                              weight = "count"), "ggplot")
  expect_error(mosaic_plot(df, x_var = "group", y_var = "state",
                           weight = "count", bogus = 1),
               "unsupported argument")
  # No longer exported: it masked tna::plot_mosaic() when both were attached
  expect_false("plot_mosaic" %in% getNamespaceExports("Nestimate"))
})

test_that("cluster_network(cluster_by = 'mmm') passes covariate_effect on", {
  set.seed(1)
  d <- as.data.frame(matrix(sample(c("A", "B", "C"), 60 * 8, TRUE), 60, 8))
  cv <- data.frame(g = rnorm(60))
  nets <- cluster_network(d, k = 2, cluster_by = "mmm", covariates = cv,
                          covariate_effect = "posthoc", n_starts = 2, seed = 1)
  fit <- build_mmm(d, k = 2, covariates = cv, covariate_effect = "posthoc",
                   n_starts = 2, seed = 1)
  expect_identical(unname(vapply(nets, \(n) nrow(n$data), 1L)),
                   as.integer(table(fit$assignments)))
})
