# Multichannel sequence/distribution plots for mcml objects.

make_mcml_seq <- function() {
  # 8 actors, 4 states in 2 clusters, short time-ordered sequences.
  actors  <- rep(1:8, each = 5)
  states  <- c("a1", "a2", "b1", "b2")
  acts    <- states[(seq_along(actors) %% 4L) + 1L]
  times   <- as.POSIXct("2025-01-01", tz = "UTC") + seq_along(actors) * 60
  df <- data.frame(Actor = actors, Action = acts, Time = times,
                   stringsAsFactors = FALSE)
  build_mcml(df,
             clusters = list(G1 = c("a1", "a2"), G2 = c("b1", "b2")),
             actor = "Actor", action = "Action", time = "Time", type = "tna")
}

test_that("sequence_plot(mcml) returns a faceted ggplot for each type", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()

  p_index <- sequence_plot(fit)
  p_dist  <- sequence_plot(fit, type = "distribution")
  p_norm  <- sequence_plot(fit, type = "distribution", normalize = TRUE)

  expect_s3_class(p_index, "ggplot")
  expect_s3_class(p_dist,  "ggplot")
  expect_s3_class(p_norm,  "ggplot")

  # All three build without error.
  expect_silent(ggplot2::ggplot_build(p_index))
  expect_silent(ggplot2::ggplot_build(p_dist))
  expect_silent(ggplot2::ggplot_build(p_norm))
})

test_that("mcml multichannel plot carries Summary + one panel per cluster", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()
  p   <- sequence_plot(fit)
  expect_setequal(levels(p$data$channel), c("Summary", "G1", "G2"))
})

test_that("normalized distribution sums to 1 within each channel-time", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()
  p   <- sequence_plot(fit, type = "distribution", normalize = TRUE)
  totals <- tapply(p$data$prop,
                   list(p$data$channel, p$data$time), sum)
  totals <- totals[!is.na(totals)]
  # Only time points with at least one active subject reach 1.
  expect_true(all(abs(totals[totals > 0] - 1) < 1e-8))
})

test_that("prevalence distribution includes an NA band that fills to 100%", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()
  p   <- sequence_plot(fit, type = "distribution")
  expect_true("NA" %in% as.character(p$data$key))
  totals <- tapply(p$data$prop, list(p$data$channel, p$data$time), sum)
  totals <- totals[!is.na(totals)]
  expect_true(all(abs(totals - 1) < 1e-8))
})

test_that("matrix-built mcml errors cleanly (no sequence data)", {
  m <- matrix(c(0, 2, 1, 0), 2L, dimnames = list(c("A", "B"), c("A", "B")))
  fit <- build_mcml(m, clusters = list(G1 = "A", G2 = "B"))
  expect_error(sequence_plot(fit, type = "index"), "no sequence data")
})
