# Clustered association tests behind mosaic_plot(residuals = ) and
# mosaic_analysis(test = ). Hand-computed values and invariants only; the
# agreement with survey, geepack, aod and stats::mantelhaen.test is checked
# in local_testing_and_equivalence/test-equiv-nested-association.R.

# A small clustered design: rows x columns x clusters.
.na_array <- function() {
  counts <- c(
    5, 1, 2, 4,   # cluster 1 (2 x 2, column-major)
    3, 3, 1, 5,
    6, 0, 2, 2,
    2, 2, 3, 4,
    4, 1, 1, 6
  )
  array(counts, c(2, 2, 5),
        dimnames = list(c("a", "b"), c("x", "y"), paste0("k", 1:5)))
}

# "Rao-Scott residuals equal chisq.test()'s with one row per cluster" moved to
# tests/equivalence/test-mosaic_test.R (a comparison with stats::chisq.test()).

test_that("Rao-Scott is symmetric in rows and columns and ignores cluster order", {
  arr <- .na_array()
  a <- Nestimate:::.nested_association(arr, "rao_scott")
  b <- Nestimate:::.nested_association(aperm(arr, c(2, 1, 3)), "rao_scott")
  d <- Nestimate:::.nested_association(arr[, , 5:1], "rao_scott")
  expect_equal(a$test$statistic, b$test$statistic, tolerance = 1e-12)
  expect_equal(a$test$df, b$test$df, tolerance = 1e-12)
  expect_equal(a$residuals, t(b$residuals), tolerance = 1e-12)
  expect_equal(a$test[c("statistic", "df", "df2", "p_value")],
               d$test[c("statistic", "df", "df2", "p_value")],
               tolerance = 1e-12)
})

test_that("Rao-Scott F is Pearson X2 over tr(Delta) for a 2 x 2 table", {
  arr <- .na_array()
  fit <- Nestimate:::.nested_association(arr, "rao_scott")
  tab <- rowSums(arr, dims = 2)
  N <- sum(tab)
  E <- outer(rowSums(tab), colSums(tab)) / N
  x2 <- sum((tab - E)^2 / E)
  # 2 x 2: one contrast, so d0 = 1 and the F is X2 / delta.
  expect_equal(fit$test$df, 1, tolerance = 1e-12)
  expect_equal(fit$test$statistic, x2 / fit$test$design_effect,
               tolerance = 1e-12)
  expect_equal(fit$test$df2, dim(arr)[3] - 1, tolerance = 1e-12)
})

test_that("Donner-Klar matches the hand computation on a balanced design", {
  # Two groups, two clusters each, four rows per cluster, binary outcome.
  y <- list(g1k1 = c(1, 1, 1, 0), g1k2 = c(1, 0, 0, 0),
            g2k1 = c(0, 0, 0, 1), g2k2 = c(1, 1, 0, 0))
  grp <- c("g1", "g1", "g2", "g2")
  arr <- array(0, c(2, 2, 4), dimnames = list(c("g1", "g2"), c("yes", "no"),
                                               names(y)))
  arr[cbind(match(grp, c("g1", "g2")), 1, 1:4)] <- vapply(y, sum, 1)
  arr[cbind(match(grp, c("g1", "g2")), 2, 1:4)] <- 4 - vapply(y, sum, 1)
  fit <- Nestimate:::.nested_association(arr, "donner_klar")

  # One-way ANOVA of the 0/1 outcome, clusters within groups.
  obs <- unlist(y)
  cl <- rep(names(y), each = 4)
  g <- rep(grp, each = 4)
  cl_mean <- ave(obs, cl)
  g_mean <- ave(obs, g)
  msw <- sum((obs - cl_mean)^2) / (16 - 4)
  msc <- sum((cl_mean - g_mean)^2) / (4 - 2)
  rho <- (msc - msw) / (msc + (4 - 1) * msw)
  tab <- rowSums(arr, dims = 2)
  E <- outer(rowSums(tab), colSums(tab)) / 16
  x2 <- sum((tab - E)^2 / E)
  expect_equal(fit$test$icc, rho, tolerance = 1e-12)
  expect_equal(fit$test$statistic, x2 / (1 + 3 * max(rho, 0)),
               tolerance = 1e-12)
  expect_equal(fit$test$df, 1)
})

test_that("Donner-Klar reduces to Pearson when clusters do not differ", {
  # Every cluster of a group has the same composition: rho <= 0, used as 0.
  arr <- array(0, c(2, 3, 4), dimnames = list(c("a", "b"), c("x", "y", "z"),
                                               NULL))
  arr[1, , 1] <- c(2, 1, 1); arr[1, , 2] <- c(2, 1, 1)
  arr[2, , 3] <- c(1, 1, 2); arr[2, , 4] <- c(1, 1, 2)
  fit <- Nestimate:::.nested_association(arr, "donner_klar")
  tab <- rowSums(arr, dims = 2)
  E <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  expect_lte(fit$test$icc, 0)
  expect_equal(fit$test$statistic, sum((tab - E)^2 / E), tolerance = 1e-12)
  expect_equal(fit$test$design_effect, 1, tolerance = 1e-12)
})

test_that("CMH residuals are summed O - E over summed hypergeometric variances", {
  arr <- .na_array()
  fit <- Nestimate:::.nested_association(arr, "cmh")
  per_k <- lapply(seq_len(dim(arr)[3]), function(k) {
    t_k <- arr[, , k]
    m <- sum(t_k); ri <- rowSums(t_k); cj <- colSums(t_k)
    list(oe = t_k[1, 1] - ri[1] * cj[1] / m,
         v = ri[1] * cj[1] * (m - ri[1]) * (m - cj[1]) / (m^2 * (m - 1)))
  })
  oe <- sum(vapply(per_k, `[[`, 1, "oe"))
  v <- sum(vapply(per_k, `[[`, 1, "v"))
  expect_equal(fit$residuals[1, 1], oe / sqrt(v), tolerance = 1e-12)
  # 2 x 2 x K without correction: the statistic is the residual squared.
  expect_equal(fit$test$statistic, (oe / sqrt(v))^2, tolerance = 1e-12)
})

test_that("GEE Wald matches the sandwich computed by hand for a 2 x 2 table", {
  arr <- .na_array()
  fit <- Nestimate:::.nested_association(arr, "gee")
  tab <- rowSums(arr, dims = 2)
  ni <- rowSums(tab)
  p <- tab[, 1] / ni                      # P(x | a), P(x | b)
  # Cluster score sums for the two proportions.
  u <- t(vapply(seq_len(dim(arr)[3]), function(k) {
    (arr[, 1, k] - rowSums(arr[, , k]) * p) / ni
  }, numeric(2)))
  V <- crossprod(u)
  wald <- (p[1] - p[2])^2 / (V[1, 1] + V[2, 2] - 2 * V[1, 2])
  expect_equal(fit$test$statistic, unname(wald), tolerance = 1e-12)
  expect_equal(fit$test$df, 1)
})

test_that("clustered tests refuse designs they cannot handle, by class", {
  arr <- .na_array()   # both variables vary within clusters
  expect_error(Nestimate:::.nested_association(arr, "donner_klar"),
               class = "nestimate_bad_input")
  grp <- array(0, c(2, 2, 4))
  grp[1, , 1:2] <- c(3, 1, 2, 2); grp[2, , 3:4] <- c(1, 3, 2, 2)
  expect_error(Nestimate:::.nested_association(grp, "cmh"),
               class = "nestimate_bad_input")
  expect_error(Nestimate:::.nested_association(arr[, , 1, drop = FALSE],
                                               "rao_scott"),
               class = "nestimate_bad_input")
  expect_error(Nestimate:::.nested_association(arr + 0.5, "gee"),
               class = "nestimate_bad_input")
})

test_that("nested cluster levels reduce to the outermost for design-based tests", {
  arr <- .na_array()
  lv <- data.frame(team = c("t1", "t1", "t2", "t2", "t3"),
                   student = paste0("s", 1:5))
  by_team <- array(c(rowSums(arr[, , 1:2], dims = 2),
                     rowSums(arr[, , 3:4], dims = 2), arr[, , 5]),
                   c(2, 2, 3), dimnames = dimnames(arr)[1:2])
  a <- Nestimate:::.nested_association(arr, "rao_scott", levels = lv)
  b <- Nestimate:::.nested_association(by_team, "rao_scott")
  expect_equal(a$test$statistic, b$test$statistic, tolerance = 1e-12)
  expect_identical(a$test$clusters, 3L)
})

test_that("mosaic_plot draws clustered residuals with the test as subtitle", {
  skip_on_cran()
  data(group_regulation_long, package = "Nestimate")
  net <- build_network(group_regulation_long, method = "frequency",
                       format = "long", actor = "Actor", action = "Action",
                       order = "Time")
  expect_output(p <- mosaic_plot(net, nested = TRUE, actor = "Group",
                               report = TRUE), "rao_scott")
  expect_true(inherits(p, "ggplot"))
  expect_match(p$labels$subtitle, "^Rao-Scott: F\\(")
  expect_match(p$labels$subtitle, "200 actors")
  expect_error(mosaic_plot(net, nested = TRUE, nested_method = "donner_klar"),
               class = "nestimate_bad_input")
  expect_error(mosaic_plot(net, nested = TRUE, nested_method = "gee", actor = "nope"),
               class = "nestimate_actor_missing")
  expect_warning(mosaic_plot(net, actor = "Group", seed = 1),
                 class = "nestimate_ignored_argument")

  grp <- build_network(group_regulation_long, method = "frequency",
                       format = "long", actor = "Actor", action = "Action",
                       order = "Time", group = "Achiever")
  expect_output(q <- mosaic_plot(grp, nested = TRUE,
                                 nested_method = "donner_klar",
                                 style = "flat", report = TRUE))
  expect_match(q$labels$subtitle, "^Donner-Klar: X2\\(8\\)")
  expect_error(mosaic_plot(grp, nested = TRUE, nested_method = "cmh"),
               class = "nestimate_bad_input")
})

test_that("mosaic_plot takes a three-way table and refuses a two-way one", {
  arr <- as.table(.na_array())
  expect_output(p <- mosaic_plot(arr, nested = TRUE, nested_method = "gee",
                               report = TRUE))
  expect_match(p$labels$subtitle, "^GEE Wald")
  expect_error(mosaic_plot(as.table(rowSums(.na_array(), dims = 2)),
                           nested = TRUE, nested_method = "gee"),
               class = "nestimate_bad_input")
})

test_that("mosaic_analysis reports a clustered test in its tidy tables", {
  d <- data.frame(
    cl = rep(paste0("k", 1:5), times = apply(.na_array(), 3, sum)),
    stringsAsFactors = FALSE)
  arr <- .na_array()
  rows <- lapply(seq_len(5), function(k) {
    t_k <- arr[, , k]
    idx <- which(t_k > 0, arr.ind = TRUE)
    reps <- rep(seq_len(nrow(idx)), t_k[idx])
    data.frame(v1 = rownames(t_k)[idx[reps, 1]], v2 = colnames(t_k)[idx[reps, 2]])
  })
  d <- cbind(d, do.call(rbind, rows))
  res <- mosaic_analysis(d, "v1", "v2", min_count = 1, nested = TRUE,
                         actor = "cl")
  fit <- Nestimate:::.nested_association(arr, "rao_scott")
  expect_identical(res$stats$test, "Rao-Scott")
  expect_identical(res$stats$clusters, 5L)
  expect_equal(res$stats$p_value, signif(fit$test$p_value, 4))
  expect_equal(summary(res)$residual,
               round(as.vector(fit$residuals), 3))
  expect_true(all(c("df2", "clusters", "design_effect", "icc") %in%
                    names(res$stats)))
  expect_error(mosaic_analysis(d, "v1", "v2", min_count = 1, nested = TRUE, nested_method = "gee"),
               class = "nestimate_bad_input")
  expect_warning(mosaic_analysis(d, "v1", "v2", min_count = 1, actor = "cl"),
                 class = "nestimate_ignored_argument")
})

test_that("GLMM gives a Wald test on the interaction and finite residuals", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  set.seed(3)
  arr <- array(stats::rpois(3 * 3 * 12, 4) + 1, c(3, 3, 12))
  fit <- Nestimate:::.nested_association(arr, "glmm")
  expect_identical(fit$test$df, 4L)
  expect_true(fit$test$statistic > 0)
  expect_true(all(is.finite(fit$residuals)))
  # A planted interaction is found: cell [1, 1] tripled in every cluster.
  planted <- arr
  planted[1, 1, ] <- planted[1, 1, ] * 3
  hit <- Nestimate:::.nested_association(planted, "glmm")
  expect_gt(hit$residuals[1, 1], 2)
  expect_gt(hit$test$statistic, fit$test$statistic)
})

test_that("mosaic_test returns one tidy row per test, matching the engine", {
  arr <- as.table(.na_array())
  out <- mosaic_test(arr)
  expect_true(is.data.frame(out))
  expect_identical(out$method,
                   c("chisq", "rao_scott", "donner_klar", "cmh", "gee", "glmm"))
  expect_identical(names(out),
                   c("method", "test", "statistic", "df", "df2", "p_value",
                     "n", "clusters", "design_effect", "icc", "note"))
  rs <- Nestimate:::.nested_association(unclass(arr), "rao_scott")$test
  expect_equal(out$statistic[out$method == "rao_scott"], rs$statistic,
               tolerance = 1e-12)
  tab <- rowSums(unclass(arr), dims = 2)
  E <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  expect_equal(out$statistic[out$method == "chisq"], sum((tab - E)^2 / E),
               tolerance = 1e-12)
  # Both variables vary within actors: Donner-Klar cannot run, says why.
  dk <- out[out$method == "donner_klar", ]
  expect_true(is.na(dk$statistic))
  expect_match(dk$note, "^Donner-Klar needs")
  expect_true(all(is.na(out$note[out$method %in% c("chisq", "rao_scott", "cmh")])))
  expect_match(out$note[out$method == "gee"], "liberal")
})

test_that("mosaic_test raises the error when one inapplicable test is asked", {
  expect_error(mosaic_test(as.table(.na_array()), method = "donner_klar"),
               class = "nestimate_bad_input")
  expect_error(mosaic_test(as.table(rowSums(.na_array(), dims = 2))),
               class = "nestimate_bad_input")
  expect_error(mosaic_test(1:3), class = "nestimate_bad_input")
  expect_error(mosaic_test(as.table(.na_array()), method = "bogus"),
               class = "nestimate_bad_input")
})

test_that("mosaic_test cells give the residuals the mosaic draws", {
  arr <- as.table(.na_array())
  cells <- mosaic_test(arr, method = c("chisq", "cmh"), level = "cells")
  expect_identical(names(cells),
                   c("method", "row", "column", "observed", "residual"))
  expect_identical(nrow(cells), 8L)
  cmh <- Nestimate:::.nested_association(unclass(arr), "cmh")$residuals
  expect_equal(cells$residual[cells$method == "cmh"], as.vector(cmh),
               tolerance = 1e-12)
})

test_that("mosaic_test reads networks and data frames with actor", {
  skip_on_cran()
  data(group_regulation_long, package = "Nestimate")
  net <- build_network(group_regulation_long, method = "frequency",
                       format = "long", actor = "Actor", action = "Action",
                       order = "Time")
  out <- mosaic_test(net, method = c("chisq", "rao_scott"), actor = "Group")
  expect_identical(out$clusters, c(NA, 200L))
  expect_equal(out$n, rep(sum(net$weights), 2))
  cells <- mosaic_test(net, method = "rao_scott", actor = "Group",
                       level = "cells")
  expect_identical(names(cells)[2:3], c("from", "to"))
  expect_equal(sum(cells$observed[cells$from == "plan" &
                                    cells$to == "consensus"]),
               net$weights["plan", "consensus"])
  df_out <- mosaic_test(group_regulation_long, var1 = "Course",
                        var2 = "Action", actor = "Actor",
                        method = c("rao_scott", "donner_klar"))
  expect_identical(df_out$clusters, c(2000L, 2000L))
  expect_true(all(df_out$p_value < 0.05))
})

test_that("nested = FALSE keeps the old default; ignored arguments warn", {
  arr <- as.table(.na_array())
  tab <- as.table(rowSums(.na_array(), dims = 2))
  # Without nested = TRUE a three-way table is refused, as in 1.0.12 (which
  # failed on it); the two-way table draws as before, with no subtitle.
  expect_error(mosaic_plot(arr, seed = 1), class = "nestimate_bad_input")
  expect_null(mosaic_plot(tab, seed = 1)$labels$subtitle)
  # nested = TRUE defaults to Rao-Scott.
  expect_output(r <- mosaic_plot(arr, nested = TRUE, report = TRUE),
                "rao_scott")
  expect_match(r$labels$subtitle, "^Rao-Scott")
  expect_null(mosaic_plot(arr, nested = TRUE)$labels$subtitle)
  expect_error(mosaic_plot(arr, nested = NA), class = "nestimate_bad_input")
  expect_error(mosaic_plot(arr, nested = "yes"), class = "nestimate_bad_input")
  expect_warning(mosaic_plot(arr, nested = TRUE, residuals = "asymptotic"),
                 class = "nestimate_ignored_argument")
  expect_warning(mosaic_plot(tab, nested_method = "gee", seed = 1),
                 class = "nestimate_ignored_argument")
})

test_that("report = TRUE prints and draws the test, nested or not", {
  tab <- as.table(rowSums(.na_array(), dims = 2))
  # Not nested: Pearson's chi-square on the drawn table, residuals unchanged.
  expect_output(p <- mosaic_plot(tab, residuals = "asymptotic", report = TRUE),
                "Pearson")
  E <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  x2 <- sum((tab - E)^2 / E)
  expect_match(p$labels$subtitle, sprintf("X2\\(1\\) = %.2f", x2))
  plain <- mosaic_plot(tab, residuals = "asymptotic")
  expect_identical(ggplot2::ggplot_build(p)$data[[1]]$fill,
                   ggplot2::ggplot_build(plain)$data[[1]]$fill)
  expect_null(plain$labels$subtitle)
  expect_error(mosaic_plot(tab, report = NA), class = "nestimate_bad_input")
  expect_error(mosaic_plot(tab, report = "yes"), class = "nestimate_bad_input")

  # mosaic_analysis: subtitle and printed row only with report = TRUE.
  d <- data.frame(v1 = rep(c("a", "b"), c(30, 30)),
                  v2 = rep(c("x", "y", "x", "y"), c(20, 10, 8, 22)))
  quiet <- mosaic_analysis(d, "v1", "v2", min_count = 1)
  expect_null(quiet$plot$labels$subtitle)
  expect_output(loud <- mosaic_analysis(d, "v1", "v2", min_count = 1,
                                        report = TRUE), "Chi-square")
  expect_match(loud$plot$labels$subtitle, "^Chi-square: X2\\(1\\)")
  expect_output(fis <- mosaic_analysis(d, "v1", "v2", min_count = 1,
                                       test = "fisher", report = TRUE),
                "Fisher")
  expect_match(fis$plot$labels$subtitle, "^Fisher's exact: p ")
})

test_that("verifier findings stay fixed (2026-11-16)", {
  # GEE refuses a rank-deficient covariance instead of an order-dependent Wald.
  set.seed(5)
  few <- array(stats::rpois(3 * 3 * 3, 6) + 1, c(3, 3, 3))
  expect_error(Nestimate:::.nested_association(few, "gee"),
               class = "nestimate_bad_input")
  # Rao-Scott with empty cells: estimable contrasts, no solve() error; with
  # an invertible matrix the result is unchanged (see the hand 2 x 2 test).
  sparse <- array(0, c(3, 3, 6))
  sparse[1, 1, ] <- c(3, 2, 4, 1, 2, 3); sparse[2, 2, ] <- c(2, 4, 1, 3, 2, 2)
  sparse[3, 3, ] <- c(1, 1, 3, 2, 4, 1); sparse[1, 2, ] <- c(1, 0, 1, 1, 0, 2)
  rs <- Nestimate:::.nested_association(sparse, "rao_scott")
  expect_true(is.finite(rs$test$statistic) && rs$test$statistic > 0)
  expect_false(any(is.infinite(rs$residuals)))
  # CMH with a single informative actor is a classed error.
  one <- array(0, c(2, 2, 3)); one[, , 1] <- c(3, 1, 2, 4)
  one[1, 1, 2] <- 1; one[2, 2, 3] <- 1
  expect_error(Nestimate:::.nested_association(one, "cmh"),
               class = "nestimate_bad_input")
  # A three-way table needs nested = TRUE; a plain 3-d array is accepted.
  arr <- .na_array()
  expect_error(mosaic_plot(as.table(arr)), class = "nestimate_bad_input")
  expect_true(inherits(mosaic_plot(arr, nested = TRUE), "ggplot"))
  expect_true(is.data.frame(mosaic_test(arr, method = "rao_scott")))
})

test_that("default calls keep their 1.0.12 surface", {
  tab <- as.table(rowSums(.na_array(), dims = 2))
  # Abbreviated old arguments still partially match (new ones follow ...).
  expect_true(inherits(mosaic_plot(tab, re = "asymptotic"), "ggplot"))
  expect_true(inherits(mosaic_plot(tab, n = 40, seed = 1), "ggplot"))
  # No subtitle label entry and no new legend title on default plots.
  p <- mosaic_plot(tab, seed = 1)
  expect_false("subtitle" %in% names(p$labels))
  q <- mosaic_plot(tab, residuals = "asymptotic", style = "flat",
                   legend_title = NULL)
  expect_null(ggplot2::ggplot_build(q)$plot$scales$get_scales("fill")$name)
  # print() of a default mosaic_analysis result works.
  d <- data.frame(v1 = rep(c("a", "b"), c(30, 30)),
                  v2 = rep(c("x", "y", "x", "y"), c(20, 10, 8, 22)))
  res <- mosaic_analysis(d, "v1", "v2", min_count = 1)
  expect_output(print(res), "Chi-square")
  expect_identical(names(res$stats),
                   c("test", "statistic", "df", "p_value", "cramers_v",
                     "effect_size", "n"))
})

test_that("GLMM residuals are calibrated O - E over a cluster-robust SE", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  # Hand computation of the residual of cell [1, 1].
  arr <- .na_array()
  fit <- suppressMessages(Nestimate:::.nested_association(arr, "glmm"))
  tab <- rowSums(arr, dims = 2)
  N <- sum(tab); p <- colSums(tab) / N; ni <- rowSums(tab)
  K <- dim(arr)[3]
  psi <- vapply(seq_len(K), function(k) {
    a <- arr[, , k]
    (a[1, 1] - sum(a[1, ]) * p[1]) - (ni[1] / N) * (sum(a[, 1]) - sum(a) * p[1])
  }, numeric(1))
  expected <- (tab[1, 1] - ni[1] * p[1]) /
    sqrt(sum((psi - mean(psi))^2) * K / (K - 1))
  expect_equal(fit$residuals[1, 1], unname(expected), tolerance = 1e-10)
  # A cell with no actor-to-actor variation is NA, never infinite.
  flat <- array(0, c(2, 2, 3)); flat[1, 1, ] <- 2; flat[2, 2, ] <- 2
  flat[1, 2, 1] <- 1; flat[2, 1, 2] <- 1; flat[1, 2, 3] <- 1
  res <- tryCatch(suppressMessages(suppressWarnings(
    Nestimate:::.nested_association(flat, "glmm")))$residuals,
    nestimate_bad_input = function(e) NULL)
  if (!is.null(res)) expect_false(any(is.infinite(res)))
})
