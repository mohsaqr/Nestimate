# Regression tests for build_network() input handling and the start/end
# fast path (review of 2026-10-26).
testthat::skip_on_cran()

.toy_wide <- function() {
  data.frame(T1 = c("A", "B", "C"), T2 = c("B", "C", "A"),
             T3 = c("C", "A", "B"))
}

test_that("actor names the ID column in wide sequence data", {
  with_id <- build_network(cbind(id = c("p1", "p2", "p3"), .toy_wide()),
                           method = "relative", actor = "id")
  without <- build_network(.toy_wide(), method = "relative")
  expect_identical(rownames(with_id$weights), c("A", "B", "C"))
  expect_equal(with_id$weights, without$weights)
  expect_identical(with_id$metadata$id, c("p1", "p2", "p3"))
})

test_that("actor names the ID column in association data", {
  set.seed(1)
  d <- data.frame(id = 1:50, x = rnorm(50), y = rnorm(50), z = rnorm(50))
  net <- build_network(d, method = "cor", actor = "id")
  expect_identical(rownames(net$weights), c("x", "y", "z"))
  expect_equal(net$weights,
               build_network(d[c("x", "y", "z")], method = "cor")$weights)
  y <- as.data.frame(matrix(rbinom(600, 1, 0.5), 200))
  y$pid <- 1:200
  expect_identical(rownames(build_network(y, method = "ising",
                                          actor = "pid")$weights),
                   c("V1", "V2", "V3"))
})

test_that("a numeric column in character sequence data is metadata", {
  net <- build_network(cbind(.toy_wide(), age = c(20, 30, 40)),
                       method = "relative")
  expect_identical(rownames(net$weights), c("A", "B", "C"))
  expect_equal(net$weights, build_network(.toy_wide(), "relative")$weights)
  expect_identical(net$metadata$age, c(20, 30, 40))
})

test_that("a character column does not turn numeric data into sequences", {
  set.seed(2)
  d <- data.frame(x = rnorm(50), y = rnorm(50), z = rnorm(50),
                  gender = sample(c("F", "M"), 50, TRUE))
  expect_identical(rownames(build_network(d, method = "cor")$weights),
                   c("x", "y", "z"))
  expect_identical(rownames(build_network(d, method = "glasso",
                                          metadata_cols = "gender")$weights),
                   c("x", "y", "z"))
  # all-categorical data are still sequences: one node per state
  seqs <- as.data.frame(matrix(sample(c("A", "B", "C"), 500, TRUE), 50, 10))
  expect_identical(rownames(build_network(seqs, method = "cor")$weights),
                   c("A", "B", "C"))
})

test_that("labels relabel $data, $initial and keep resampling intact", {
  set.seed(2)
  s <- as.data.frame(matrix(sample(c("A", "B", "C"), 300, TRUE), 30, 10))
  plain <- build_network(s, method = "relative")
  named <- build_network(s, method = "relative",
                         labels = c(A = "Alpha", B = "Beta"))
  expect_identical(names(named$initial), c("Alpha", "Beta", "C"))
  expect_true(all(unlist(named$data) %in% c("Alpha", "Beta", "C", NA)))
  b_plain <- bootstrap_network(plain, iter = 20, seed = 1)
  b_named <- bootstrap_network(named, iter = 20, seed = 1)
  expect_equal(unname(b_named$mean), unname(b_plain$mean))

  g <- build_network(srl_strategies, method = "glasso")
  gl <- build_network(srl_strategies, method = "glasso",
                      labels = stats::setNames(c("one", "two"),
                                               names(srl_strategies)[1:2]))
  expect_equal(predictability(gl)$RMSE, predictability(g)$RMSE)
  expect_equal(predictability(gl, data = srl_strategies)$RMSE,
               predictability(g)$RMSE)
})

test_that("millisecond timestamps: unit forwarded, wrong unit warned", {
  d <- data.frame(id = rep(c("u1", "u2"), each = 4),
                  act = rep(c("a", "b", "c", "a"), 2),
                  t = 1.7e12 + rep(0:3, 2) * 60000)
  net <- build_network(d, "frequency", actor = "id", action = "act",
                       time = "t", unix_time_unit = "milliseconds")
  # a->b, b->c, c->a once per user
  expect_equal(unname(net$weights),
               matrix(c(0, 2, 0, 0, 0, 2, 2, 0, 0), 3, byrow = TRUE))
  expect_false("unix_time_unit" %in% names(net$params))
  expect_warning(build_network(d, "frequency", actor = "id", action = "act",
                               time = "t"),
                 class = "nestimate_time_unit")
})

test_that("fast-path counts reproduce networks with start/end markers", {
  rebuild <- function(net) {
    st <- rownames(net$weights)
    counts <- .precompute_per_sequence(net$data, net$method, net$params, st)
    m <- matrix(colSums(counts), length(st), byrow = TRUE,
                dimnames = list(st, st))
    if (net$method == "relative") {
      rs <- rowSums(m)
      m[rs > 0, ] <- m[rs > 0, ] / rs[rs > 0]
    }
    m
  }
  configs <- list(list("relative", start = TRUE),
                  list("relative", end = TRUE),
                  list("relative", start = TRUE, end = TRUE),
                  list("frequency", start = TRUE, end = TRUE),
                  list("co_occurrence", start = TRUE))
  invisible(lapply(configs, function(cfg) {
    net <- do.call(build_network,
                   c(list(trajectories, method = cfg[[1]]), cfg[-1]))
    expect_equal(rebuild(net), net$weights, tolerance = 1e-12)
  }))
  b <- bootstrap_network(build_network(trajectories, method = "relative",
                                       start = TRUE), iter = 30, seed = 1)
  expect_true(all(b$sd["Start", c("Active", "Average", "Disengaged")] > 0))
})

test_that("level = 'both' keeps actor; group errors on NA; group keeps labels", {
  set.seed(1)
  d <- data.frame(id = rep(1:20, each = 10), x = rnorm(200), y = rnorm(200),
                  z = rnorm(200))
  ml <- build_network(d, method = "cor", actor = "id", level = "both")
  expect_s3_class(ml, "netobject_ml")
  expect_identical(rownames(ml$between$weights), c("x", "y", "z"))

  s <- data.frame(T1 = c("A", "B", "C", "A"), T2 = c("B", "C", "A", "C"),
                  g = c("x", "y", NA, "x"))
  expect_error(build_network(s, "frequency", group = "g"),
               class = "nestimate_bad_input")
  s$g[3] <- "y"
  grouped <- build_network(s, "frequency", group = "g",
                           labels = c(A = "Alpha"))
  expect_identical(grouped[[1]]$nodes$label, c("Alpha", "B", "C"))
})
