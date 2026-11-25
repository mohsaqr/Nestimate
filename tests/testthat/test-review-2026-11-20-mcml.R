# Regression tests for the 2026-11-20 review of MCML, clustering and the
# landscape / period plots. Every expected value is hand-computed from the
# definition, not read off the code.
testthat::skip_on_cran()

# ---- K1 / K3: a resampled mcml layer is recounted as it was fitted ----
# These rely on .precompute_per_sequence() honouring the `symmetrize` and
# `row_weight` parameters that as_tna() records on an mcml layer.

test_that("K1: an undirected mcml resamples symmetrised adjacent counts", {
  # four identical sequences A, B, C: every resample is the data, so the
  # bootstrap mean is the fit, (W + t(W)) / 2 of the adjacent pairs
  d <- data.frame(T1 = rep("A", 4), T2 = "B", T3 = "C")
  m <- build_mcml(d, list(A = "A", B = "B", C = "C"),
                  type = "cooccurrence", compute_within = FALSE)
  fit <- matrix(c(0, 2, 0, 2, 0, 2, 0, 2, 0), 3,
                dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  expect_equal(m$macro$weights, fit)
  net <- as_tna(m)$macro
  expect_identical(net$method, "frequency")
  expect_true(isTRUE(net$params$symmetrize))
  b <- bootstrap_network(m, iter = 5, seed = 1)$macro
  # A and C are never adjacent: no phantom A-C edge
  expect_equal(unname(b$mean), unname(fit))
})

test_that("K1: per-sequence counts of a directed = FALSE mcml sum to its fit", {
  d <- data.frame(T1 = c("A", "C", "B", "A"), T2 = c("B", "A", "D", "C"),
                  T3 = c("C", "D", "A", "B"))
  m <- build_mcml(d, list(G1 = c("A", "B"), G2 = c("C", "D")),
                  type = "raw", directed = FALSE)
  net <- as_tna(m)$macro
  per <- Nestimate:::.precompute_per_sequence(net$data, net$method,
                                              net$params, c("G1", "G2"))
  # by hand: G1->G1 1 (A->B), G1->G2 4 (B->C, A->D, B->D, A->C),
  # G2->G1 3 (C->A, D->A, C->B), G2->G2 0; symmetrised (4 + 3) / 2 = 3.5
  expect_equal(matrix(colSums(per), 2, byrow = TRUE),
               matrix(c(1, 3.5, 3.5, 0), 2))
  expect_equal(unname(m$macro$weights), matrix(c(1, 3.5, 3.5, 0), 2))
})

test_that("K3: a weighted edge list is resampled with its weights", {
  d <- data.frame(src = rep("A", 4), dst = "B", w = 10)
  m <- build_mcml(d, list(A = "A", B = "B"), from = "src", to = "dst",
                  weight = "w", type = "raw", compute_within = FALSE)
  expect_equal(m$macro$weights[["A", "B"]], 40)
  expect_warning(b <- bootstrap_network(m, iter = 5, seed = 1),
                 class = "nestimate_edgelist_resampling")
  # identical rows: every resample sums four weights of 10
  expect_equal(b$macro$mean[["A", "B"]], 40)
})

# ---- K2: a windowed network is recounted with its window ----

test_that("K2: build_mcml() of a windowed network counts block pairs", {
  d <- data.frame(T1 = rep("A", 3), T2 = "B", T3 = "C", T4 = "D")
  n <- build_network(d, method = "frequency", params = list(window_size = 2))
  m <- build_mcml(n, as.list(stats::setNames(LETTERS[1:4], LETTERS[1:4])),
                  type = "raw")
  # windows {T1, T2} -> {T3, T4}: A and B each to C and D, 3 sequences
  expect_equal(m$macro$weights[["A", "C"]], 3)
  expect_equal(m$macro$weights[["A", "B"]], 0)
  expect_equal(m$macro$weights, n$weights)
})

test_that("K2: windowed clusters aggregate, resample and re-partition by window", {
  set.seed(2)
  d <- as.data.frame(matrix(sample(c("a", "b", "c", "d"), 50, TRUE), 10))
  n <- build_network(d, method = "frequency", params = list(window_size = 2))
  members <- list(X = c("a", "b"), Y = c("c", "d"))
  m <- build_mcml(n, members, type = "raw")
  of <- c(a = "X", b = "X", c = "Y", d = "Y")
  W <- n$weights
  by_cluster <- t(rowsum(t(rowsum(W, of[rownames(W)])), of[colnames(W)]))
  expect_equal(unname(m$macro$weights), unname(by_cluster))
  expect_equal(m$clusters$X$weights, W[c("a", "b"), c("a", "b")])
  # the resampling units count the same window pairs
  net <- as_tna(m)$macro
  per <- Nestimate:::.precompute_per_sequence(net$data, net$method,
                                              net$params, c("X", "Y"))
  expect_equal(unname(matrix(colSums(per), 2, byrow = TRUE)),
               unname(by_cluster))
  # re-partitioning keeps the window
  opened <- build_mcml(m, expand = "X")
  expect_equal(opened$macro$weights[c("a", "b"), c("a", "b")],
               W[c("a", "b"), c("a", "b")])
})

# ---- K4: distance clustering of a one-hot network ----

test_that("K4: cluster_network() clusters a one-hot network's actors", {
  s <- rep(c("A", "A", "B", "B", "B", "B", "A", "A"), 2)
  d <- data.frame(id = rep(1:4, each = 4), A = as.integer(s == "A"),
                  B = as.integer(s == "B"))
  n <- build_network(d, method = "frequency", actor = "id",
                     codes = c("A", "B"), window_size = 2)
  cls <- build_clusters(n, k = 2, method = "complete")
  r <- cluster_network(n, k = 2, cluster_by = "complete")
  # each cluster's network is the sum of its actors' own networks
  per_actor <- lapply(split(d, d$id), function(g) {
    build_network(g, method = "frequency", codes = c("A", "B"),
                  window_size = 2)$weights
  })
  expected <- lapply(1:2, function(cl) {
    Reduce(`+`, per_actor[cls$assignments == cl])
  })
  expect_equal(unname(r[[1]]$weights), unname(expected[[1]]))
  expect_equal(unname(r[[2]]$weights), unname(expected[[2]]))
})

# ---- K5: build_mcml_pc() keeps the network's correlation ----

test_that("K5: build_mcml_pc() re-estimates with the network's cor_method", {
  d <- data.frame(a = 1:20, b = (1:20)^4, c = rep(0:1, 10))
  # a and b are perfectly rank-correlated, so c's neighbours are singular:
  # that warning (class nestimate_singular_neighbours) is expected here
  quiet_singular <- function(expr) {
    withCallingHandlers(expr, nestimate_singular_neighbours = function(w) {
      invokeRestart("muffleWarning")
    })
  }
  n <- quiet_singular(build_network(d, method = "cor", cor_method = "spearman"))
  f <- quiet_singular(build_mcml_pc(
    n, list(X = c("a", "b"), Y = "c"), signed = FALSE,
    weighting = function(W, data, nodes) c(W[1, 2], 1)))
  expect_identical(f$meta$cor_method, "spearman")
  # a and b are monotone: Spearman 1 in every resample, so the weight of a
  # is 1 / (1 + 1) = 0.5 in each, and so is the bootstrap mean
  ls <- quiet_singular(loading_stability(f, iter = 5, seed = 1))
  a_row <- subset(quiet_singular(summary(ls)), node == "a")
  expect_equal(a_row$weight, 0.5)
  expect_equal(a_row$boot_mean, 0.5)
})

# ---- K7: free parameters of the covariate mixing ----

test_that("K7: an unused factor level adds no mixing parameter", {
  d <- data.frame(T1 = rep(c("A", "B"), each = 6))
  d$T2 <- d$T1
  cv <- data.frame(g = factor(rep(c("x", "y"), 6), levels = c("x", "y", "z")))
  expect_warning(
    f <- build_mmm(d, covariates = cv, estimator = "multinom",
                   n_starts = 1, seed = 1),
    class = "nestimate_multinom_ml")
  # k = 2, 2 states: 2 * 2 * 1 transition + 2 * 1 initial + (2 - 1) * rank,
  # rank(model.matrix(~ g)) = 2 (the column of level z is all zero)
  expect_equal(qr(stats::model.matrix(~ g, cv))$rank, 2L)
  expect_equal(f$n_params, 4 + 2 + 1 * 2)
  expect_equal(f$BIC, -2 * f$log_likelihood + 8 * log(12))
  # the same count against nnet::multinom(): tests/equivalence/test-build_mmm.R
})

# ---- K8: an undirected edge counts once ----

test_that("K8: cluster_summary(directed = FALSE) counts each edge once", {
  w <- matrix(c(3, 2, 2, 0), 2, dimnames = list(c("A", "B"), c("A", "B")))
  m <- cluster_summary(w, list(G = c("A", "B")), directed = FALSE)
  # self-loop A 3 + edge A-B 2
  expect_equal(m$macro$weights[["G", "G"]], 5)
  expect_equal(summary(m)$within_total, 5)
  expect_equal(sum(as_tna(m)$G$edges$weight), 5)
  # density: 5 over the 3 possible undirected edges (A-A, A-B, B-B)
  dens <- cluster_summary(w, list(G = c("A", "B")), directed = FALSE,
                          method = "density")
  expect_equal(dens$macro$weights[["G", "G"]], 5 / 3)
  # directed, A->B and B->A are two edges
  md <- cluster_summary(w, list(G = c("A", "B")), directed = TRUE)
  expect_equal(md$macro$weights[["G", "G"]], 7)
  expect_equal(summary(md)$within_total, 7)
})

test_that("K8: undirected transitions keep within_total = macro diagonal", {
  # each transition counts once: A->B within G1 gives A-B 0.5 and B-A 0.5
  d <- data.frame(T1 = c("A", "C"), T2 = c("B", "C"), T3 = c("A", "A"))
  m <- build_mcml(d, list(G1 = c("A", "B"), G2 = "C"), type = "raw",
                  directed = FALSE)
  # within G1: A->B, B->A = 2 transitions; within G2: C->C = 1
  expect_equal(unname(diag(m$macro$weights)), c(2, 1))
  expect_equal(summary(m)$within_total, c(2, 1))
})

# ---- K9: the post-hoc multinomial model predicts on new data ----

test_that("K9: the covariate model predicts from the original columns", {
  y <- factor(c(rep(1, 30), rep(2, 10), rep(1, 12), rep(2, 28)))
  x <- c(rep(2015, 40), rep(2016, 40))
  fit <- .multinom_standardised(y, stats::model.matrix(~ x), y ~ x,
                                data.frame(y = y, x = x))
  p <- predict(fit$model, newdata = data.frame(x = c(2015, 2016)),
               type = "probs")
  # P(y = 2 | x) is the cell share: 10 / 40 and 28 / 40
  expect_equal(unname(p), c(10 / 40, 28 / 40), tolerance = 1e-5)
  expect_equal(unname(stats::coef(fit$model)),
               unname(fit$estimates$coefficients[1, ]), tolerance = 1e-6)

  d <- data.frame(T1 = rep(c("A", "B"), each = 6), x = rep(1:6, 2))
  d$T2 <- d$T1
  expect_warning(f <- build_clusters(d, k = 2, covariates = "x",
                                     estimator = "multinom"),
                 class = "nestimate_multinom_ml")
  new <- predict(f$covariates$model, newdata = data.frame(x = c(1, 3, 6)),
                 type = "probs")
  expect_length(new, 3L)
  # x takes each value once per cluster, so the fitted probability is 1/2
  expect_equal(unname(new), rep(0.5, 3), tolerance = 1e-4)
})

# ---- K10: the within-cluster networks keep counts under "density" ----

test_that("K10: density keeps the within-cluster counts", {
  d <- data.frame(T1 = rep("A", 3), T2 = rep("B", 3))
  m <- build_mcml(d, list(G = c("A", "B")), method = "density",
                  type = "raw")
  expect_equal(m$clusters$G$weights[["A", "B"]], 3)
})

# ---- B7: labels reach the source macro_network() expands ----

test_that("B7: macro_network(expand = ) works on a labelled mcml", {
  d <- data.frame(T1 = c("A", "B", "B"), T2 = c("B", "A", "B"))
  mc <- build_mcml(d, clusters = list(G1 = "A", G2 = "B"),
                   labels = c(A = "Alpha"))
  r <- macro_network(mc, expand = "G1")
  # transitions A->B, B->A, B->B: Alpha -> G2 1; G2 -> Alpha 1/2, G2 -> G2 1/2
  expect_equal(r$weights[["Alpha", "G2"]], 1)
  expect_equal(r$weights[["G2", "Alpha"]], 0.5)
  expect_equal(r$weights[["G2", "G2"]], 0.5)
})

# ---- A11: one definition of the classification error ----

test_that("A11: overall classification error agrees with the per-cluster one", {
  old <- structure(list(
    k = 2L, assignments = 1:2, mixing = c(0.5, 0.5),
    posterior = rbind(c(0.9, 0.1), c(0.2, 0.8)),
    quality = list(classification_error = 0)
  ), class = "net_mmm")
  z <- cluster_diagnostics(old)
  # mean(1 - max posterior) = mean(0.1, 0.2)
  expect_equal(z$overall$classification_error, 0.15)
  expect_equal(z$per_cluster$class_err_pct, c(10, 20))
})

# ---- P1: landscapes refuse windowed networks ----

test_that("P1: a windowed network has no one-step landscape", {
  d <- expand.grid(T1 = c("A", "B"), T2 = c("A", "B"), T3 = c("A", "B"),
                   T4 = c("A", "B"), stringsAsFactors = FALSE)
  n <- build_network(d, method = "relative", params = list(window_size = 2L))
  expect_error(stability_landscape(n, iter = 0), class = "nestimate_bad_input")
  expect_error(landscape_trend(n, window = 2, iter = 0),
               class = "nestimate_bad_input")
  # the same network unwindowed: 16 sequences x 3 adjacent pairs = 48
  n1 <- build_network(d, method = "relative")
  z <- stability_landscape(n1, iter = 0)
  expect_equal(sum(summary(z)$n_from), 48)
})

# ---- P4: a group_tna gives one landscape per group ----

test_that("P4: landscapes split a group_tna by group", {
  skip_if_not_installed("tna")
  set.seed(4)
  d <- as.data.frame(matrix(sample(c("A", "B", "C"), 240, TRUE), 30))
  grp <- rep(c("G1", "G2"), each = 15)
  m <- tna::group_tna(d, group = grp)
  expect_identical(unique(summary(stability_landscape(m, iter = 0))$group),
                   c("G1", "G2"))
  expect_identical(
    unique(summary(landscape_trend(m, window = 4, iter = 0))$group),
    c("G1", "G2"))
})

# ---- P5: one closed class has a unique stationary distribution ----

test_that("P5: a transient state has an infinite potential, not NA", {
  d <- data.frame(T1 = rep("A", 2), T2 = "B", T3 = "C", T4 = "B")
  z <- summary(landscape_trend(d, depth = "occupancy", window = 4, iter = 0))
  wells <- subset(z, feature == "well")
  # A transient, B <-> C closed: pi = (0, 1/2, 1/2)
  expect_equal(wells$potential, c(Inf, log(2), log(2)))
  ridges <- subset(z, feature == "ridge")
  # B | C: (1/2 * 1 + 1/2 * 1) / 2 = 1/2; A | B and C | A never crossed
  expect_equal(ridges$potential, c(NA, log(2), NA))
  circ <- summary(stability_landscape(d, iter = 0), table = "circulation")
  expect_equal(circ$current, 0)
  # two closed classes: no unique distribution
  two <- data.frame(T1 = c("A", "B"), T2 = c("A", "B"), T3 = c("A", "B"))
  expect_true(is.na(Nestimate:::.ls_current(matrix(c(1, 0, 0, 0, 1, 0,
                                                     0, 0, 1), 3))))
  expect_null(Nestimate:::.ls_stationary(diag(2)))
})

# ---- P7: one state gives a count matrix of the right shape ----

test_that("P7: one-state sequences count without recycling", {
  s <- matrix(c("A", "A", NA, "A", "A", "A"), 2, byrow = TRUE)
  counts <- Nestimate:::.ls_counts(s, "A")
  # one A->A pair in row 1, two in row 2
  expect_equal(counts, matrix(c(1, 2), 2, 1))
})

# ---- P8: a group that varies inside one period is refused ----

test_that("P8: periods refuse a group that changes within one period", {
  d <- data.frame(id = "s1", day = 1, act = c("A", "B", "A"),
                  condition = c("X", "X", "Y"))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(sequence_plot(d, type = "index", actor = "id", action = "act",
                             periods = "day", group = "condition",
                             legend = "none"),
               class = "nestimate_bad_input")
})

# ---- P10: the trend line breaks at a point that is not estimable ----

test_that("P10: no line joins across an unestimable window", {
  d <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"), T3 = c("B", "A"),
                  T4 = c("A", "B"))
  z <- landscape_trend(d, window = 2, iter = 0)
  p <- plot(z)
  is_line <- vapply(p$layers, function(l) inherits(l$geom, "GeomLine"),
                    logical(1))
  built <- ggplot2::ggplot_build(p)$data[[which(is_line)]]
  # windows at 1.5, 2.5 (both states never left: Inf), 3.5: no line group
  # holds points on both sides of 2.5
  spans <- tapply(built$x, built$group, function(v) min(v) < 2.5 && max(v) > 2.5)
  expect_false(any(spans))
})
