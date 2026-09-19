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

# Every per-channel ggplot of a sequence_plot(mcml) figure: the figure itself
# when a single channel is drawn, else the panels stacked into it.
mcml_panels <- function(p) {
  if (inherits(p, "ggplot")) list(p) else attr(p, "panels")
}
# The long data of all panels, as one data.frame with shared factor levels.
mcml_panel_data <- function(p) {
  ds <- lapply(mcml_panels(p), function(q) q$data)
  d  <- do.call(rbind, lapply(ds, function(x) {
    x$channel <- as.character(x$channel); x$key <- as.character(x$key); x
  }))
  d$channel <- factor(d$channel, levels = unique(d$channel))
  d$key     <- factor(d$key, levels = unique(d$key))
  d
}

test_that("sequence_plot(mcml) stacks one ggplot per channel", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()

  p_index <- sequence_plot(fit)
  p_dist  <- sequence_plot(fit, type = "distribution")
  p_norm  <- sequence_plot(fit, type = "distribution", normalize = TRUE)

  lapply(list(p_index, p_dist, p_norm), function(p) {
    expect_s3_class(p, "mcml_sequence_plot")
    expect_s3_class(p, "gtable")
    lapply(mcml_panels(p), function(q) expect_silent(ggplot2::ggplot_build(q)))
  })
  # print draws without error and returns the figure
  pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(withVisible(print(p_index))$visible, FALSE)
})

test_that("each channel panel carries its own legend of its own keys", {
  skip_if_not_installed("ggplot2")
  fit    <- make_mcml_seq()
  breaks <- function(q) ggplot2::get_guide_data(q, "fill")$.label
  panels <- mcml_panels(sequence_plot(fit, type = "distribution", normalize = TRUE))
  expect_named(panels, c("Summary", "G1", "G2"))
  expect_equal(breaks(panels$Summary), c("G1", "G2"))
  expect_equal(breaks(panels$G1), c("a1", "a2"))
  expect_equal(breaks(panels$G2), c("b1", "b2"))
  # prevalence: own states first, then the faded other clusters, then NA
  prev <- mcml_panels(sequence_plot(fit, type = "distribution"))
  expect_equal(breaks(prev$G1), c("a1", "a2", "G2 (elsewhere)", "NA"))
  # carpet: no state of another cluster leaks into a channel legend
  idx <- mcml_panels(sequence_plot(fit))
  expect_false(any(c("b1", "b2") %in% breaks(idx$G1)))
  # title only on the top panel, x-axis title only on the bottom one
  expect_null(prev$G1$labels$title)
  expect_null(prev$G1$labels$x)
  expect_equal(prev$G2$labels$x, "Time")
})

test_that("mcml multichannel plot carries Summary + one panel per cluster", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()
  p   <- mcml_panel_data(sequence_plot(fit))
  expect_setequal(levels(p$channel), c("Summary", "G1", "G2"))
})

test_that("normalized distribution sums to 1 within each channel-time", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()
  p   <- mcml_panel_data(sequence_plot(fit, type = "distribution", normalize = TRUE))
  totals <- tapply(p$prop, list(p$channel, p$time), sum)
  totals <- totals[!is.na(totals)]
  # Only time points with at least one active subject reach 1.
  expect_true(all(abs(totals[totals > 0] - 1) < 1e-8))
})

test_that("prevalence distribution includes an NA band that fills to 100%", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()
  p   <- mcml_panel_data(sequence_plot(fit, type = "distribution"))
  expect_true("NA" %in% as.character(p$key))
  totals <- tapply(p$prop, list(p$channel, p$time), sum)
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

  figs <- list(sequence_plot(fit), sequence_plot(fit, type = "distribution"),
               sequence_plot(fit, type = "distribution", normalize = TRUE))
  lapply(figs, function(p) lapply(mcml_panels(p), function(q) {
    expect_silent(ggplot2::ggplot_build(q))
    # The colliding name appears at most once in any panel legend.
    brks <- ggplot2::get_guide_data(q, "fill")$.label
    expect_lte(sum(brks == "End"), 1L)
    expect_false(anyDuplicated(brks) > 0L)
  }))
})

test_that("panel draws the macro and the channels apart, each with one palette", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()

  both     <- sequence_plot(fit, type = "distribution", panel = "both")
  summary  <- sequence_plot(fit, type = "distribution", panel = "summary")
  channels <- sequence_plot(fit, type = "distribution", panel = "channels")

  expect_s3_class(summary, "ggplot")              # one channel: a plain ggplot
  expect_s3_class(channels, "mcml_sequence_plot")  # two clusters: stacked
  expect_silent(ggplot2::ggplot_build(summary))
  both     <- mcml_panel_data(both)
  channels <- mcml_panel_data(channels)

  # The macro panel is keyed by cluster, the rest by state. Split apart, each
  # scale holds one kind of key, so a cluster and a state cannot collide.
  expect_equal(levels(summary$data$channel), "Summary")
  expect_false("Summary" %in% levels(channels$channel))
  # Prevalence adds an NA band to both panels; the clusters are the solid keys.
  expect_true(all(names(fit$cluster_members) %in% levels(summary$data$key)))
  expect_setequal(
    levels(sequence_plot(fit, type = "distribution", normalize = TRUE,
                         panel = "summary")$data$key),
    names(fit$cluster_members))
  # Prevalence also carries the faded "elsewhere" bands and NA; the states are
  # the solid keys. Composition drops both, leaving states alone.
  states <- unlist(fit$cluster_members, use.names = FALSE)
  expect_true(all(states %in% levels(channels$key)))
  norm <- mcml_panel_data(sequence_plot(fit, type = "distribution",
                                       normalize = TRUE, panel = "channels"))
  expect_setequal(levels(norm$key), states)
  # the shared scale carries more keys than either panel alone
  expect_gt(nlevels(both$key),
            max(nlevels(summary$data$key), nlevels(norm$key)))
})

test_that("each panel carries a default title that main overrides", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()
  expect_match(sequence_plot(fit, type = "distribution",
                             panel = "summary")$labels$title, "Macro")
  expect_match(attr(sequence_plot(fit, type = "distribution", panel = "channels"),
                    "panels")[[1L]]$labels$title, "Within-cluster")
  expect_equal(sequence_plot(fit, type = "distribution", panel = "summary",
                             main = "Mine")$labels$title, "Mine")
  expect_null(attr(sequence_plot(fit, type = "distribution", panel = "both"),
                   "panels")[[1L]]$labels$title)
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
    levels(mcml_panel_data(sequence_plot(fit, type = "distribution", normalize = TRUE,
                                         panel = "channels"))$key),
    levels(mcml_panel_data(sequence_plot(fit, type = "distribution", normalize = TRUE,
                                         panel = "channels", expand = "G1"))$key))

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

test_that("sequence_plot(expand =) works in the default prevalence view", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq()

  # regression: the prevalence branch indexed the cluster-share table by
  # cluster name, which has no row for an expanded cluster
  p <- sequence_plot(fit, type = "distribution", expand = "G1")
  d <- mcml_panel_data(p)

  # the Summary band opens G1; the G2 panel keeps G1 as one faded band
  expect_setequal(unique(as.character(d$key[d$channel == "Summary"])),
                  c("a1", "a2", "G2", "NA"))
  expect_true("G1 (elsewhere)" %in% as.character(d$key[d$channel == "G2"]))

  # every channel still stacks to 100% at every time point
  tot <- stats::aggregate(prop ~ channel + time, d, sum)
  expect_equal(tot$prop, rep(1, nrow(tot)))

  # the faded G1 share equals the summed a1 + a2 shares of the Summary band
  faded <- d[d$channel == "G2" & d$key == "G1 (elsewhere)", ]
  own   <- stats::aggregate(prop ~ time,
                            d[d$channel == "Summary" & d$key %in% c("a1", "a2"), ],
                            sum)
  expect_equal(faded$prop[order(faded$time)], own$prop[order(own$time)])
})
