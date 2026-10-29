# Regression cases verified independently in the October 26 review.
testthat::skip_on_cran()

test_that("undefined association refits refuse the whole permutation test", {
  set.seed(2)
  d <- as.data.frame(matrix(rnorm(56), 8L, 7L))
  net <- build_network(d, method = "pcor")
  failure <- tryCatch(permutation(net, net, iter = 199L, seed = 1L),
                      error = identity)
  expect_s3_class(failure, "nestimate_permutation_failed")
  expect_identical(failure$n_failed, 194L)
  expect_identical(failure$iter, 199L)
  expect_match(conditionMessage(failure), "No p-values were computed")

  # This ordinary control has a statistic defined on every split.
  set.seed(26)
  ordinary <- build_network(as.data.frame(matrix(rnorm(240), 60L, 4L)),
                             method = "pcor")
  result <- permutation(ordinary, ordinary, iter = 29L, seed = 1L)
  expect_equal(result$global$observed, c(0, 0))
  expect_equal(result$global$p_value, c(1, 1))
  expect_true(all(result$p_values == 1))
})

test_that("non-finite refits do not return apparently finite permutation p-values", {
  make <- function(seed) {
    set.seed(seed)
    data.frame(rare1 = as.numeric(seq_len(12L) == 1L),
               rare2 = as.numeric(seq_len(12L) == 2L),
               rare3 = as.numeric(seq_len(12L) == 3L),
               ordinary = rnorm(12L))
  }
  x <- build_network(make(26L), method = "pcor")
  y <- build_network(make(27L), method = "pcor")
  failure <- suppressWarnings(tryCatch(
    permutation(x, y, iter = 39L, seed = 1L), error = identity))
  expect_s3_class(failure, "nestimate_permutation_failed")
  expect_gt(failure$n_failed, 0L)
  expect_lte(failure$n_failed, 39L)
})

test_that("bootstrap weighting agrees with independent full-estimator draws", {
  d <- data.frame(T1 = c("A", "A", "B", "B"),
                   T2 = c("A", "B", "A", "B"),
                   T3 = c("B", NA, "B", "A"),
                   T4 = c("C", NA, NA, "C"))
  invisible(lapply(c("relative", "frequency", "co_occurrence"), function(method) {
    net <- build_network(d, method = method, params = list(weighted = TRUE))
    boot <- bootstrap_network(net, iter = 19L, seed = 26L)
    states <- rownames(net$weights)
    set.seed(26L)
    direct <- vapply(seq_len(19L), function(i) {
      take <- sample.int(nrow(d), nrow(d), replace = TRUE)
      fitted <- build_network(d[take, , drop = FALSE], method = method,
                               params = list(weighted = TRUE, alphabet = states))
      as.vector(fitted$weights[states, states])
    }, numeric(length(states)^2L))
    expect_equal(as.vector(boot$mean), rowMeans(direct), info = method)
    expect_equal(as.vector(boot$ci_lower),
                   apply(direct, 1L, quantile, probs = 0.025), info = method)
    expect_equal(as.vector(boot$ci_upper),
                   apply(direct, 1L, quantile, probs = 0.975), info = method)
  }))
})

test_that("concat resamples merged sequences with their boundary joins", {
  raw <- data.frame(T1 = c("A", "B", "A", "C", "C", "A"),
                     T2 = c("B", "C", "A", "A", "B", "C"))
  # Independently construct the three sequences produced by concat=2.
  merged <- data.frame(S1 = raw$T1[c(1L, 3L, 5L)],
                       S2 = raw$T2[c(1L, 3L, 5L)],
                       S3 = raw$T1[c(2L, 4L, 6L)],
                       S4 = raw$T2[c(2L, 4L, 6L)])
  invisible(lapply(c("relative", "frequency", "co_occurrence"), function(method) {
    a <- build_network(raw, method, start = TRUE, end = TRUE,
                        params = list(concat = 2L, weighted = TRUE))
    b <- build_network(merged, method, start = TRUE, end = TRUE,
                        params = list(weighted = TRUE))
    expect_equal(a$weights, b$weights)
    ba <- bootstrap_network(a, iter = 19L, seed = 26L)
    bb <- bootstrap_network(b, iter = 19L, seed = 26L)
    expect_equal(ba$mean, bb$mean, info = method)
    expect_equal(ba$ci_lower, bb$ci_lower, info = method)
    pa <- permutation(a, a, iter = 19L, seed = 26L)
    pb <- permutation(b, b, iter = 19L, seed = 26L)
    expect_equal(pa$p_values, pb$p_values, info = method)
    expect_equal(pa$null_sd, pb$null_sd, info = method)
  }))
})

test_that("concatenated actor units stay within their actor", {
  raw <- data.frame(id = rep(paste0("p", 1:4), each = 4L),
                     T1 = rep(c("A", "B", "C", "A"), 4L),
                     T2 = rep(c("B", "A", "A", "C"), 4L))
  odd <- seq.int(1L, nrow(raw), by = 2L)
  merged <- data.frame(id = raw$id[odd], S1 = raw$T1[odd], S2 = raw$T2[odd],
                       S3 = raw$T1[odd + 1L], S4 = raw$T2[odd + 1L])
  a <- build_network(raw, "frequency", actor = "id", params = list(concat = 2L))
  b <- build_network(merged, "frequency", actor = "id")
  ba <- bootstrap_network(a, iter = 19L, seed = 26L, actor = "id")
  bb <- bootstrap_network(b, iter = 19L, seed = 26L, actor = "id")
  expect_equal(ba$mean, bb$mean)
  expect_identical(ba$n_actors, 4L)
  expect_identical(ba$clustering$n_sequences, 8L)
  pa <- permutation(a, a, iter = 19L, seed = 26L, actor = "id")
  pb <- permutation(b, b, iter = 19L, seed = 26L, actor = "id")
  expect_equal(pa$p_values, pb$p_values)
  expect_equal(pa$null_sd, pb$null_sd)

  raw$id <- rep(paste0("p", 1:4), 4L)
  crossed <- build_network(raw, "frequency", actor = "id",
                            params = list(concat = 2L))
  expect_error(bootstrap_network(crossed, iter = 3L, actor = "id"),
                class = "nestimate_actor_misaligned")
  expect_error(permutation(crossed, crossed, iter = 3L, actor = "id"),
                class = "nestimate_actor_misaligned")
})

test_that("paired concat tests validate the actual resampling units", {
  d <- data.frame(T1 = rep(c("A", "B"), 3L), T2 = rep(c("B", "A"), 3L))
  x <- build_network(d, "frequency", params = list(concat = 2L))
  y <- build_network(d, "frequency", params = list(concat = 3L))
  expect_error(permutation(x, y, iter = 3L, paired = TRUE),
                "after concatenation", class = "nestimate_bad_input")
  weighted <- build_network(d, "frequency", params = list(weighted = TRUE))
  unweighted <- build_network(d, "frequency")
  expect_error(permutation(weighted, unweighted, iter = 3L),
                "same `weighted`", class = "nestimate_bad_input")
})
