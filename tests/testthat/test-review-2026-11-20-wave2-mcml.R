testthat::skip_on_cran()

# Regression tests for the second wave of the 2026-11-20 MCML review.
# Expected values are hand-computed from the definitions.

# ---- C3: a non-sum aggregation is applied within each resample ----

test_that("C3: bootstrap of a mean-aggregated edge-list mcml applies the mean", {
  e <- data.frame(f = c("A", "A"), t = c("B", "B"), w = c(2, 2))
  m <- build_mcml(e, from = "f", to = "t", weight = "w",
                  clusters = list(X = "A", Y = "B"), type = "raw",
                  method = "mean", compute_within = FALSE)
  macro <- as_tna(m)$macro
  expect_equal(macro$weights[["X", "Y"]], 2)
  expect_identical(macro$params$aggregation, "mean")
  expect_warning(b <- bootstrap_network(macro, iter = 5, seed = 1),
                 class = "nestimate_edgelist_resampling")
  # every resample holds two rows of weight 2: their mean is 2, not 4
  expect_equal(b$mean[["X", "Y"]], 2)
  expect_equal(b$sd[["X", "Y"]], 0)
})

test_that("C3: max, min, median and geomean of equal weights resample to the weight", {
  e <- data.frame(f = rep("A", 3), t = rep("B", 3), w = c(2, 2, 2))
  boot_mean <- vapply(c("max", "min", "median", "geomean"), function(meth) {
    m <- build_mcml(e, from = "f", to = "t", weight = "w",
                    clusters = list(X = "A", Y = "B"), type = "raw",
                    method = meth, compute_within = FALSE)
    # the edge-list warning is tested above
    b <- withCallingHandlers(
      bootstrap_network(as_tna(m)$macro, iter = 4, seed = 2),
      nestimate_edgelist_resampling = function(w) {
        invokeRestart("muffleWarning")
      })
    b$mean[["X", "Y"]]
  }, numeric(1))
  expect_equal(unname(boot_mean), rep(2, 4))
})

test_that("C3: a mean-aggregated sequence mcml gives 1 per observed cluster pair", {
  # every sequence holds A -> B, so X -> Y is observed in every resample
  d <- data.frame(T1 = rep("A", 4), T2 = "B", T3 = c("B", "A", "B", "A"))
  m <- build_mcml(d, list(X = "A", Y = "B"), type = "raw", method = "mean",
                  compute_within = FALSE)
  expect_equal(m$macro$weights[["X", "Y"]], 1)
  b <- bootstrap_network(m, iter = 5, seed = 1)$macro
  expect_equal(b$mean[["X", "Y"]], 1)
  expect_equal(b$sd[["X", "Y"]], 0)
})

test_that("C3: each replicate is the aggregate of the rows it draws", {
  e <- data.frame(f = c("A", "C", "A", "C", "B", "D"),
                  t = c("B", "D", "D", "B", "A", "C"),
                  w = c(1, 3, 2, 6, 4, 5))
  clusters <- list(X = c("A", "B"), Y = c("C", "D"))
  of <- c(A = "X", B = "X", C = "Y", D = "Y")
  cells <- paste(of[e$f], of[e$t])
  replicates <- function(meth, type) {
    m <- build_mcml(e, from = "f", to = "t", weight = "w", clusters = clusters,
                    type = type, method = meth, compute_within = FALSE)
    net <- as_tna(m)$macro
    Nestimate:::.local_seed(11)
    reps <- Nestimate:::.bootstrap_aggregated(
      net$data, net$method, net$params, c("X", "Y"), scaling = NULL,
      threshold = 0, iter = 3, is_relative = net$method == "relative",
      include_zeros = TRUE)
    set.seed(11)
    draws <- lapply(1:3, function(i) sample.int(6, 6, replace = TRUE))
    list(reps = reps, draws = draws)
  }
  by_hand <- function(idx, meth, n_possible = 1) {
    vapply(c("X X", "Y X", "X Y", "Y Y"), function(cell) {
      net_aggregate_weights(e$w[idx][cells[idx] == cell], meth,
                            n_possible = n_possible)
    }, numeric(1))
  }
  # mean, raw: column-major cells XX, YX, XY, YY
  r <- replicates("mean", "raw")
  expected <- t(vapply(r$draws, by_hand, numeric(4), meth = "mean"))
  expect_equal(unname(r$reps), unname(expected))
  # median, tna: the medians, then each row normalised
  r <- replicates("median", "tna")
  expected <- t(vapply(r$draws, function(idx) {
    w <- matrix(by_hand(idx, "median"), 2)
    as.vector(w / ifelse(rowSums(w) == 0, 1, rowSums(w)))
  }, numeric(4)))
  expect_equal(unname(r$reps), unname(expected))
  # density, raw: the sums divided by n_i * n_j = 4
  r <- replicates("density", "raw")
  expected <- t(vapply(r$draws, by_hand, numeric(4), meth = "density",
                       n_possible = 4))
  expect_equal(unname(r$reps), unname(expected))
})

test_that("C3: .aggregate_repeated() is net_aggregate_weights() of the repeats", {
  values <- c(1, 3, 0, 2)
  times <- c(2, 1, 3, 0)
  # rep(values, times) = 1 1 3 0 0 0; zeros are left out: 1 1 3
  hand <- c(mean = 5 / 3, median = 1, max = 3, min = 1,
            geomean = 3^(1 / 3), density = 5 / 4)
  got <- vapply(names(hand), function(meth) {
    Nestimate:::.aggregate_repeated(values, times, meth, n_possible = 4)
  }, numeric(1))
  expect_equal(got, hand)
  # an even number of repeats: the median averages the middle two
  expect_equal(Nestimate:::.aggregate_repeated(c(1, 5), c(1, 1), "median"), 3)
})

test_that("C3: a sum mcml records no aggregation, so its counts are kept", {
  d <- data.frame(src = rep("A", 4), dst = "B", w = 10)
  m <- build_mcml(d, list(A = "A", B = "B"), from = "src", to = "dst",
                  weight = "w", type = "raw", compute_within = FALSE)
  expect_identical(as_tna(m)$macro$params,
                   list(cols = c("from", "to"), row_weight = "weight"))
})

test_that("C3: actor is refused for a non-sum aggregation", {
  d <- data.frame(T1 = c("A", "B", "A"), T2 = c("B", "A", "A"),
                  T3 = c("A", "B", "B"))
  m <- build_mcml(d, list(X = "A", Y = "B"), type = "raw", method = "mean",
                  compute_within = FALSE)
  expect_error(bootstrap_network(as_tna(m)$macro, iter = 3, actor = "T1"),
               class = "nestimate_actor_unsupported")
})

# ---- C7: a windowed network without between-window pairs ----

test_that("C7: build_mcml() of a one-window-per-sequence network is all zero", {
  d <- data.frame(T1 = c("A", "A"), T2 = c("B", "B"))
  n <- build_network(d, "frequency", params = list(window_size = 2L))
  m <- build_mcml(n, clusters = list(X = "A", Y = "B"), type = "raw")
  expect_true(inherits(m, "mcml"))
  expect_equal(m$macro$weights,
               matrix(0, 2, 2, dimnames = list(c("X", "Y"), c("X", "Y"))))
  expect_equal(nrow(m$edges), 0L)
  expect_equal(unname(m$clusters$X$weights), matrix(0, 1, 1))
})

# ---- C10: a Kendall network is read as it is when nothing is re-estimated ----

kendall_data <- function() {
  data.frame(a = 1:8, b = c(1, 3, 2, 4, 6, 5, 8, 7),
             c = c(4, 1, 6, 8, 2, 7, 5, 3), d = c(5, 2, 7, 8, 1, 6, 4, 3))
}

test_that("C10: average + subnetwork averages and slices a Kendall network", {
  n <- build_cor(kendall_data(), params = list(cor_method = "kendall"))
  clusters <- list(X = c("a", "b"), Y = c("c", "d"))
  f <- build_mcml_pc(n, clusters = clusters, aggregation = "average",
                     within = "subnetwork")
  W <- n$weights
  expect_identical(f$meta$cor_method, "kendall")
  expect_equal(f$macro$weights[["X", "Y"]], mean(W[c("a", "b"), c("c", "d")]))
  expect_equal(f$macro$weights[["X", "X"]], W[["a", "b"]])
  expect_equal(f$clusters$X$weights, W[c("a", "b"), c("a", "b")])
  expect_equal(f$clusters$Y$weights, W[c("c", "d"), c("c", "d")])
  # the same call on a Pearson network carrying these weights agrees
  p <- build_cor(kendall_data())
  p$weights <- W
  g <- build_mcml_pc(p, clusters = clusters, aggregation = "average",
                     within = "subnetwork")
  expect_equal(f$macro$weights, g$macro$weights)
  expect_equal(f$clusters$X$weights, g$clusters$X$weights)
})

test_that("C10: the other descriptive and score paths pass with subnetwork", {
  n <- build_cor(kendall_data(), params = list(cor_method = "kendall"))
  clusters <- list(X = c("a", "b"), Y = c("c", "d"))
  fits <- lapply(c("escoufier", "cancor", "scaled"), function(agg) {
    build_mcml_pc(n, clusters = clusters, aggregation = agg,
                  within = "subnetwork")
  })
  expect_true(all(vapply(fits, inherits, logical(1), what = "mcml_pc")))
})

test_that("C10: Kendall is still refused where a correlation is re-estimated", {
  n <- build_cor(kendall_data(), params = list(cor_method = "kendall"))
  clusters <- list(X = c("a", "b"), Y = c("c", "d"))
  expect_error(
    build_mcml_pc(n, clusters = clusters, aggregation = "average"),
    class = "nestimate_bad_input")
  expect_error(
    build_mcml_pc(n, clusters = clusters, aggregation = "scaled",
                  within = "subnetwork", weighting = "factor"),
    class = "nestimate_bad_input")
  # an explicit cor_method is honoured
  f <- build_mcml_pc(n, clusters = clusters, aggregation = "average",
                     cor_method = "spearman")
  expect_identical(f$meta$cor_method, "spearman")
})

test_that("only bootstrap_network() resamples a non-sum mcml layer", {
  set.seed(1)
  edges <- data.frame(f = sample(c("A", "B", "C", "D"), 60, TRUE),
                      t = sample(c("A", "B", "C", "D"), 60, TRUE),
                      w = round(runif(60, 1, 9)))
  clusters <- list(X = c("A", "B"), Y = c("C", "D"))
  macro <- as_tna(build_mcml(edges, from = "f", to = "t", weight = "w",
                             clusters = clusters, type = "raw",
                             method = "mean"))$macro
  expect_warning(boot <- bootstrap_network(macro, iter = 20, seed = 1),
                 class = "nestimate_edgelist_resampling")
  expect_s3_class(boot, "net_bootstrap")
  expect_error(permutation(macro, macro, iter = 20, seed = 1),
               class = "nestimate_bad_input")
  expect_error(casedrop_reliability(macro, iter = 10, seed = 1),
               class = "nestimate_bad_input")
  expect_error(centrality_stability(macro, iter = 10, seed = 1),
               class = "nestimate_bad_input")
  expect_error(entropy_bayes(macro, draws = 100, seed = 1),
               class = "nestimate_bad_input")
})
