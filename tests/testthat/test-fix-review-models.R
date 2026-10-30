# Regression oracles for the package review of 2026-10-30.
testthat::skip_on_cran()

test_that('post-hoc covariates do not alter sequence distances or clustering', {
  d <- data.frame(V1 = rep(c('A', 'B'), each = 8),
                  V2 = rep(c('B', 'A'), each = 8),
                  outcome = rep(c('High', 'Low'), 8))
  reference <- build_clusters(d[c('V1', 'V2')], k = 2)
  invisible(lapply(list('outcome', ~ outcome, d['outcome']), function(covariates) {
    fit <- build_clusters(d, k = 2, covariates = covariates, estimator = 'chisq')
    expect_equal(fit$distance, reference$distance)
    expect_equal(fit$assignments, reference$assignments)
    expect_named(fit$data, c('V1', 'V2'))
    expect_false(is.null(fit$covariates))
  }))
  # A separately supplied data frame can share a name with a state position.
  external <- data.frame(V1 = rep(c('High', 'Low'), 8))
  fit <- build_clusters(d[c('V1', 'V2')], k = 2,
                        covariates = external, estimator = 'chisq')
  expect_equal(fit$distance, reference$distance)
  expect_named(fit$data, c('V1', 'V2'))
  expect_error(build_clusters(d['outcome'], k = 2, covariates = 'outcome'),
               'sequence columns')
})

test_that('loading stability replays equal and custom composite weights', {
  set.seed(91)
  n <- 120
  z <- rnorm(n)
  q <- rnorm(n)
  d <- data.frame(a = z + rnorm(n, sd = .2), b = z + rnorm(n, sd = .4),
                  c = z + rnorm(n, sd = 2), d = q + rnorm(n, sd = .2),
                  e = q + rnorm(n, sd = .4), f = q + rnorm(n, sd = 2))
  cl <- list(one = c('a', 'b', 'c'), two = c('d', 'e', 'f'))
  w <- setNames(c(1, 2, 3, 3, 2, 1), names(d))
  schemes <- list('equal', w, function(W_block, data_block, nodes) {
    seq_along(nodes)
  })
  invisible(lapply(schemes, function(scheme) {
    fit <- suppressWarnings(build_mcml_pc(d, cl, aggregation = 'composite',
                            weighting = scheme, method = 'cor', signed = FALSE))
    out <- loading_stability(fit, iter = 3, seed = 8)
    expect_equal(out$summary$boot_mean, fit$loadings$weight)
    expect_equal(out$summary$boot_sd, rep(0, 6))
    expect_equal(out$summary$ci_lower, fit$loadings$weight)
    expect_equal(out$summary$ci_upper, fit$loadings$weight)
  }))
  # Data-dependent weights are checked against a complete refit of a known draw.
  fit <- suppressWarnings(build_mcml_pc(d, cl, aggregation = 'composite',
                          weighting = 'pca', method = 'cor', signed = FALSE))
  set.seed(19)
  rows <- sample.int(n, n, replace = TRUE)
  expected <- suppressWarnings(build_mcml_pc(d[rows, ], cl,
                            aggregation = 'composite', weighting = 'pca',
                            method = 'cor', signed = FALSE))
  out <- loading_stability(fit, iter = 2, seed = 19)
  expect_equal(unname(out$boot_weights[1, ]), expected$loadings$weight)
  # Old custom-weight objects cannot replay an unavailable function/vector.
  fit$meta$weighting <- 'custom (function)'
  fit$meta$weighting_spec <- NULL
  expect_error(loading_stability(fit, iter = 2), 'Rebuild')
})

test_that('compare_model aligns named nodes and rejects different node sets', {
  W <- matrix(c(0, 1, 2, 3, 0, 4, 5, 6, 0), 3, 3, byrow = TRUE,
               dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  expected <- compare_model(W, W, network = FALSE)
  perm <- c(3, 1, 2)
  fit <- compare_model(W, W[perm, perm], network = FALSE)
  expect_equal(fit$summary_metrics, expected$summary_metrics)
  expect_equal(fit$difference_matrix, expected$difference_matrix)
  expect_equal(compare_model(W, unname(W), network = FALSE)$difference_matrix,
                expected$difference_matrix)
  expect_equal(compare_model(unname(W), W, network = FALSE)$summary_metrics,
                expected$summary_metrics)
  different <- W
  dimnames(different) <- list(c('A', 'B', 'D'), c('A', 'B', 'D'))
  expect_error(compare_model(W, different, network = FALSE), 'same nodes')
  duplicate <- W
  dimnames(duplicate) <- list(c('A', 'A', 'C'), c('A', 'A', 'C'))
  expect_error(compare_model(W, duplicate, network = FALSE), 'unique')
  # Row and column orders are independently labelled.
  expect_equal(compare_model(W, W[perm, ], network = FALSE)$difference_matrix,
                expected$difference_matrix)
})

test_that('ordinary Gaussian outcome intervals agree with linear model t intervals', {
  d <- data.frame(x = 1:6)
  d$y <- .6 * d$x + c(1, -1, 0, 0, -1, 1)
  fit <- outcome_model(d, 'y', 'x', family = 'gaussian')
  oracle <- stats::confint(stats::lm(y ~ x, d))
  expect_equal(fit$effects$ci_lower, unname(oracle[, 1]))
  expect_equal(fit$effects$ci_upper, unname(oracle[, 2]))
  expect_lt(fit$effects$ci_lower[2], 0)
  expect_gt(fit$effects$p_value[2], .05)
  binary <- data.frame(x = rep(0:1, each = 6), y = rep(c(0, 0, 1, 1, 0, 1), 2))
  b <- outcome_model(binary, 'y', 'x', family = 'binomial')
  z <- qnorm(.975)
  expect_equal(b$effects$ci_lower,
                b$effects$estimate - z * b$effects$std_error)
})

test_that('range normalization is invariant to small common weight scales', {
  invisible(lapply(c(1, 1e-12, 1e-100), function(multiplier) {
    expect_equal(Nestimate:::.range01_tna(c(1, 2, 3) * multiplier), c(0, .5, 1))
    expect_true(all(is.nan(Nestimate:::.range01_tna(rep(multiplier, 3)))))
  }))
  expect_true(all(is.nan(Nestimate:::.range01_tna(rep(0, 3)))))
  W <- matrix(c(1e12, 1, 0, 0, 1e12, 2, 3, 0, 1e12), 3, 3,
               byrow = TRUE, dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  fit <- build_network(W, method = 'relative')
  out <- suppressMessages(net_centrality(fit, measures = 'InStrength', normalize = TRUE))
  expect_equal(out$InStrength, c(1, 0, .5), tolerance = 1e-10)
})
