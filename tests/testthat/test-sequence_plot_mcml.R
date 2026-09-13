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

test_that("a cluster named after one of its own states does not break the fill scale", {
  skip_if_not_installed("ggplot2")
  # A singleton cluster named after its state is the natural way to write a
  # terminal marker; state keys and cluster keys then collide in the shared
  # `key` fill scale. Regression: factor() rejected the duplicated level.
  actors <- rep(1:6, each = 4)
  acts   <- c("a1", "a2", "b1", "End")[(seq_along(actors) %% 4L) + 1L]
  times  <- as.POSIXct("2025-01-01", tz = "UTC") + seq_along(actors) * 60
  df  <- data.frame(Actor = actors, Action = acts, Time = times,
                    stringsAsFactors = FALSE)
  fit <- build_mcml(df,
                    clusters = list(G1 = c("a1", "a2"), G2 = "b1", End = "End"),
                    actor = "Actor", action = "Action", time = "Time")

  expect_true("End" %in% names(fit$clusters))
  expect_true("End" %in% unlist(fit$cluster_members, use.names = FALSE))

  expect_silent(ggplot2::ggplot_build(sequence_plot(fit)))
  expect_silent(ggplot2::ggplot_build(sequence_plot(fit, type = "distribution")))
  expect_silent(ggplot2::ggplot_build(
    sequence_plot(fit, type = "distribution", normalize = TRUE)))

  # The colliding name appears exactly once in the legend breaks.
  ch <- .mcml_seq_channels(fit, trim = NULL)
  brks <- .mcml_grouped_breaks(ch)
  expect_equal(sum(brks == "End"), 1L)
  expect_false(anyDuplicated(brks) > 0L)
})

test_that("panel draws the macro and the channels apart, each with one palette", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()

  both     <- sequence_plot(fit, type = "distribution", panel = "both")
  summary  <- sequence_plot(fit, type = "distribution", panel = "summary")
  channels <- sequence_plot(fit, type = "distribution", panel = "channels")

  expect_s3_class(summary, "ggplot")
  expect_s3_class(channels, "ggplot")
  expect_silent(ggplot2::ggplot_build(summary))
  expect_silent(ggplot2::ggplot_build(channels))

  # The macro panel is keyed by cluster, the rest by state. Split apart, each
  # scale holds one kind of key, so a cluster and a state cannot collide.
  expect_equal(levels(summary$data$channel), "Summary")
  expect_false("Summary" %in% levels(channels$data$channel))
  # Prevalence adds an NA band to both panels; the clusters are the solid keys.
  expect_true(all(names(fit$cluster_members) %in% levels(summary$data$key)))
  expect_setequal(
    levels(sequence_plot(fit, type = "distribution", normalize = TRUE,
                         panel = "summary")$data$key),
    names(fit$cluster_members))
  # Prevalence also carries the faded "elsewhere" bands and NA; the states are
  # the solid keys. Composition drops both, leaving states alone.
  states <- unlist(fit$cluster_members, use.names = FALSE)
  expect_true(all(states %in% levels(channels$data$key)))
  norm <- sequence_plot(fit, type = "distribution", normalize = TRUE,
                        panel = "channels")
  expect_setequal(levels(norm$data$key), states)
  # the shared scale carries more keys than either panel alone
  expect_gt(nlevels(both$data$key),
            max(nlevels(summary$data$key), nlevels(norm$data$key)))
})

test_that("each panel carries a default title that main overrides", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()
  expect_match(sequence_plot(fit, type = "distribution",
                             panel = "summary")$labels$title, "Macro")
  expect_match(sequence_plot(fit, type = "distribution",
                             panel = "channels")$labels$title, "Within-cluster")
  expect_equal(sequence_plot(fit, type = "distribution", panel = "summary",
                             main = "Mine")$labels$title, "Mine")
  expect_null(sequence_plot(fit, type = "distribution",
                            panel = "both")$labels$title)
})

test_that("panel is rejected for the carpet types", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()
  expect_error(sequence_plot(fit, type = "index", panel = "summary"),
               "distribution")
})

test_that("sequence_plot(expand =) opens the Summary band only", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()

  plain <- sequence_plot(fit, type = "distribution", normalize = TRUE,
                         panel = "summary")
  open  <- sequence_plot(fit, type = "distribution", normalize = TRUE,
                         panel = "summary", expand = "G1")

  # the macro band keys by cluster, or by the states of an expanded cluster
  expect_setequal(levels(plain$data$key), c("G1", "G2"))
  expect_setequal(levels(open$data$key), c("a1", "a2", "G2"))

  # the channels are the partition either way
  expect_equal(
    levels(sequence_plot(fit, type = "distribution", normalize = TRUE,
                         panel = "channels")$data$key),
    levels(sequence_plot(fit, type = "distribution", normalize = TRUE,
                         panel = "channels", expand = "G1")$data$key))

  # and the object itself is untouched -- an expanded macro stored on the
  # mcml would be mis-drawn by cograph::plot_mcml(), which indexes the macro
  # positionally over n_clusters
  expect_equal(nrow(fit$macro$weights), length(fit$cluster_members))
})

test_that("sequence_plot(expand = all) opens every cluster in the Summary band", {
  skip_if_not_installed("ggplot2")
  fit    <- make_mcml_seq()
  states <- unlist(fit$cluster_members, use.names = FALSE)
  for (arg in list("all", TRUE)) {
    p <- sequence_plot(fit, type = "distribution", normalize = TRUE,
                       panel = "summary", expand = arg)
    expect_setequal(levels(p$data$key), states)
  }
  expect_error(sequence_plot(fit, expand = "nope"), "Unknown")
})
