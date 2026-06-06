# Tests for mosaic_analysis() and mosaic_plot(style = "flat")

make_df <- function(n = 300, seed = 1) {
  set.seed(seed)
  data.frame(
    g = sample(c("F", "M"), n, replace = TRUE),
    lvl = sample(c("Low", "Mid", "High"), n, replace = TRUE,
                 prob = c(0.5, 0.3, 0.2)),
    stringsAsFactors = FALSE
  )
}

test_that("mosaic_analysis returns expected structure", {
  res <- mosaic_analysis(make_df(), "g", "lvl", min_count = 5)
  expect_true(inherits(res, "mosaic_analysis"))
  expect_true(inherits(res$plot, "ggplot"))
  expect_true(all(c("counts", "stats", "test", "cramers_v", "table") %in%
                    names(res)))
})

test_that("counts is tidy: one row per cell and observed sums to N", {
  res <- mosaic_analysis(make_df(), "g", "lvl", min_count = 5)
  expect_equal(nrow(res$counts), nrow(res$table) * ncol(res$table))
  expect_equal(sum(res$counts$observed), sum(res$table))
  expect_true(all(c("g", "lvl", "observed", "expected", "residual", "pct") %in%
                    names(res$counts)))
})

test_that("Cramer's V is in [0,1] and matches a manual chi-square", {
  res <- mosaic_analysis(make_df(), "g", "lvl", min_count = 5)
  expect_gte(res$cramers_v, 0)
  expect_lte(res$cramers_v, 1)
  tab <- res$table
  chi <- suppressWarnings(stats::chisq.test(tab))
  v <- sqrt(as.numeric(chi$statistic) /
              (sum(tab) * (min(dim(tab)) - 1)))
  expect_equal(res$cramers_v, v, tolerance = 1e-10)
})

test_that("fisher test path reports Fisher and still gives Cramer's V", {
  res <- mosaic_analysis(make_df(), "g", "lvl", min_count = 5, test = "fisher")
  expect_identical(res$stats$test, "Fisher's exact")
  expect_true(is.na(res$stats$statistic))
  expect_false(is.na(res$cramers_v))
})

test_that("min_count drops sparse categories", {
  df <- make_df()
  df$lvl[df$lvl == "High"] <- "High"            # keep
  # Force a rare category that must be dropped.
  df <- rbind(df, data.frame(g = "F", lvl = "Rare", stringsAsFactors = FALSE))
  res <- mosaic_analysis(df, "g", "lvl", min_count = 5)
  expect_true("Rare" %in% res$removed$var2)
  expect_false("Rare" %in% colnames(res$table))
})

test_that("input validation errors are clear", {
  expect_error(mosaic_analysis(1:10, "a", "b"), "data.frame")
  expect_error(mosaic_analysis(make_df(), "nope", "lvl"), "not found")
})

test_that("percentage_base changes the pct column base", {
  res_row <- mosaic_analysis(make_df(), "g", "lvl", min_count = 5,
                             percentage_base = "row")
  # Row-percentages within each var1 level sum to 100.
  agg <- tapply(res_row$counts$pct, res_row$counts$g, sum)
  expect_true(all(abs(agg - 100) < 1e-6))
})

test_that("plot.mosaic_analysis re-renders with overrides", {
  res <- mosaic_analysis(make_df(), "g", "lvl", min_count = 5)
  p <- plot(res, tile_label = "percent", legend_position = "bottom")
  expect_true(inherits(p, "ggplot"))
})

test_that("mosaic_plot style='flat' returns a ggplot and classic still works", {
  net <- build_network(
    data.frame(V1 = sample(LETTERS[1:4], 120, TRUE),
               V2 = sample(LETTERS[1:4], 120, TRUE),
               V3 = sample(LETTERS[1:4], 120, TRUE)),
    method = "frequency")
  p_flat    <- mosaic_plot(net, residuals = "asymptotic", style = "flat",
                           tile_label = "count")
  p_classic <- mosaic_plot(net, residuals = "asymptotic")
  expect_true(inherits(p_flat, "ggplot"))
  expect_true(inherits(p_classic, "ggplot"))
})

test_that("pct_base override does not collide with percentage_base (regression)", {
  # Was: 'formal argument "pct_base" matched by multiple actual arguments'.
  expect_no_error(
    mosaic_analysis(make_df(), "g", "lvl", min_count = 5, pct_base = "row"))
})

test_that("unknown flat styling arg gives a clear, listed error (regression)", {
  # `values` is a mosaic_plot arg, not a mosaic_analysis/flat arg -> friendly.
  expect_error(
    mosaic_analysis(make_df(), "g", "lvl", min_count = 5, values = TRUE),
    "unsupported styling argument")
  expect_error(
    mosaic_analysis(make_df(), "g", "lvl", min_count = 5, tile_lable = "count"),
    "Valid styling arguments")
})

test_that("mosaic_plot(style='flat') rejects unknown styling args clearly", {
  net <- build_network(
    data.frame(V1 = sample(LETTERS[1:4], 120, TRUE),
               V2 = sample(LETTERS[1:4], 120, TRUE)),
    method = "frequency")
  expect_error(
    mosaic_plot(net, residuals = "asymptotic", style = "flat", bogus = 1),
    "unsupported styling argument")
})

test_that("flat fill matches the reported chi-square residuals (no recompute)", {
  res <- mosaic_analysis(make_df(), "g", "lvl", min_count = 5)
  # The stored residual matrix is exactly chisq.test()$stdres on the table.
  chi <- suppressWarnings(stats::chisq.test(res$table))
  expect_equal(res$plot_parts$res, chi$stdres, tolerance = 1e-12)
})

test_that("interpret_cramers_v honours df-adjusted thresholds", {
  expect_identical(Nestimate:::.interpret_cramers_v(0.04, 1L), "negligible")
  expect_identical(Nestimate:::.interpret_cramers_v(0.20, 1L), "small")
  expect_identical(Nestimate:::.interpret_cramers_v(0.40, 1L), "medium")
  expect_identical(Nestimate:::.interpret_cramers_v(0.60, 1L), "large")
})
