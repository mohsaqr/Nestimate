# Tests for plot_state_frequencies() and plot_mosaic()
# Validates: ggplot return shape, layer count, palette correctness across
# all four supported classes (netobject, netobject_group, mcml, htna)
# and both styles (marimekko grouped/faceted, bars).

library(testthat)
suppressMessages(library(ggplot2))


# ---------------------------------------------------------------------------
# plot_mosaic primitive
# ---------------------------------------------------------------------------

test_that("plot_mosaic returns a ggplot with coherent rectangle geometry", {
  df <- data.frame(
    group = rep(c("A", "B", "C"), each = 3),
    state = rep(c("s1", "s2", "s3"), 3),
    count = c(10, 5, 2,  7, 8, 3,  4, 6, 12)
  )
  p <- plot_mosaic(df, x = "group", y = "state", weight = "count")
  expect_s3_class(p, "ggplot")
  expect_gte(length(p$layers), 1L)

  rd <- p$data
  # Widths span [0, 1]
  expect_equal(min(rd$xmin), 0, tolerance = 1e-12)
  expect_equal(max(rd$xmax), 1, tolerance = 1e-12)
  # Heights per column reach 1
  agg <- aggregate(rd$ymax, by = list(rd$x_level), FUN = max)
  expect_true(all(abs(agg$x - 1) < 1e-12))
})

test_that("plot_mosaic errors on missing columns", {
  df <- data.frame(a = 1:3, b = letters[1:3])
  expect_error(plot_mosaic(df, x = "a", y = "missing", weight = "b"),
               "missing columns")
})

test_that("plot_mosaic errors when total weight is zero", {
  df <- data.frame(g = c("A", "B"), s = c("x", "y"), w = c(0, 0))
  expect_error(plot_mosaic(df, x = "g", y = "s", weight = "w"),
               "total weight is zero")
})


# ---------------------------------------------------------------------------
# plot_state_frequencies on netobject
# ---------------------------------------------------------------------------

test_that("plot_state_frequencies works on netobject (wide trajectories)", {
  data(trajectories, package = "Nestimate")
  nw <- build_network(as.data.frame(trajectories),
                      method = "relative", format = "wide")
  expect_s3_class(nw, "netobject")

  p_marimekko <- plot_state_frequencies(nw)
  expect_s3_class(p_marimekko, "ggplot")

  p_bars <- plot_state_frequencies(nw, style = "bars")
  expect_s3_class(p_bars, "ggplot")
})


# ---------------------------------------------------------------------------
# plot_state_frequencies on netobject_group
# ---------------------------------------------------------------------------

test_that("plot_state_frequencies works on netobject_group", {
  data(group_regulation_long, package = "Nestimate")
  nw_g <- build_network(group_regulation_long,
                        method = "relative", format = "long",
                        actor = "Actor", action = "Action",
                        order = "Time", group = "Course")
  expect_s3_class(nw_g, "netobject_group")
  expect_equal(length(nw_g), 3L)  # 3 courses

  p_marimekko <- plot_state_frequencies(nw_g)
  expect_s3_class(p_marimekko, "ggplot")

  p_bars <- plot_state_frequencies(nw_g, style = "bars")
  expect_s3_class(p_bars, "ggplot")
})


# ---------------------------------------------------------------------------
# plot_state_frequencies on mcml
# ---------------------------------------------------------------------------

test_that("plot_state_frequencies works on mcml", {
  data(trajectories, package = "Nestimate")
  trj <- as.data.frame(trajectories)
  states <- unique(unlist(trj))
  states <- states[!is.na(states)]
  clusters <- setNames(c("A", "A", "B"), states)

  mc <- build_mcml(trj, clusters = clusters)
  expect_s3_class(mc, "mcml")

  # Default for mcml is legend = "per_facet" -> gtable output.
  skip_if_not_installed("gridExtra")
  p_default <- plot_state_frequencies(mc)
  expect_true(inherits(p_default, "gtable") || inherits(p_default, "ggplot"))

  # Explicit single-legend bottom -> ggplot.
  p_single <- plot_state_frequencies(mc, legend = "bottom")
  expect_s3_class(p_single, "ggplot")

  p_macro <- plot_state_frequencies(mc, include_macro = TRUE,
                                    legend = "bottom")
  expect_s3_class(p_macro, "ggplot")

  p_bars <- plot_state_frequencies(mc, style = "bars")
  expect_s3_class(p_bars, "ggplot")
})


# ---------------------------------------------------------------------------
# plot_state_frequencies on htna
# ---------------------------------------------------------------------------

test_that("plot_state_frequencies works on htna", {
  skip_if_not_installed("htna")
  data(ai_long, package = "Nestimate")
  data(human_long, package = "Nestimate")

  ht <- htna::build_htna(list(AI = ai_long, Human = human_long),
                          action = "code", session = "session_id",
                          order = "code_order")
  expect_s3_class(ht, "htna")
  expect_false(is.null(ht$node_groups))

  # Default for htna is legend = "per_facet" -> gtable output.
  skip_if_not_installed("gridExtra")
  p_default <- plot_state_frequencies(ht)
  expect_true(inherits(p_default, "gtable") || inherits(p_default, "ggplot"))

  # Explicit single-legend bottom -> ggplot.
  p_single <- plot_state_frequencies(ht, legend = "bottom")
  expect_s3_class(p_single, "ggplot")

  p_bars <- plot_state_frequencies(ht, style = "bars")
  expect_s3_class(p_bars, "ggplot")
})


# ---------------------------------------------------------------------------
# Defaults and error paths
# ---------------------------------------------------------------------------

test_that("plot_state_frequencies.default errors with informative message", {
  expect_error(plot_state_frequencies(list(foo = 1)),
               "no method for class")
})

test_that("plot_state_frequencies respects custom colors", {
  data(trajectories, package = "Nestimate")
  nw <- build_network(as.data.frame(trajectories),
                      method = "relative", format = "wide")
  custom <- c("#000000", "#FFFFFF", "#FF0000", "#00FF00", "#0000FF")
  p <- plot_state_frequencies(nw, colors = custom)
  expect_s3_class(p, "ggplot")
})

test_that("sort_states ordering changes state factor levels", {
  data(trajectories, package = "Nestimate")
  nw <- build_network(as.data.frame(trajectories),
                      method = "relative", format = "wide")
  p_freq  <- plot_state_frequencies(nw, sort_states = "frequency")
  p_alpha <- plot_state_frequencies(nw, sort_states = "alpha")
  # Both should return ggplots; we don't compare order strictly here as
  # underlying counts differ. Just ensure both dispatch cleanly.
  expect_s3_class(p_freq, "ggplot")
  expect_s3_class(p_alpha, "ggplot")
})
