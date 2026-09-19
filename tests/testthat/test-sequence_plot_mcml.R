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
  expect_equal(breaks(prev$G1), c("a1", "a2", "G2 (Other states)", "NA"))
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
  # Prevalence also carries the faded other-cluster bands and NA; the states are
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
  expect_true("G1 (Other states)" %in% as.character(d$key[d$channel == "G2"]))

  # every channel still stacks to 100% at every time point
  tot <- stats::aggregate(prop ~ channel + time, d, sum)
  expect_equal(tot$prop, rep(1, nrow(tot)))

  # the faded G1 share equals the summed a1 + a2 shares of the Summary band
  faded <- d[d$channel == "G2" & d$key == "G1 (Other states)", ]
  own   <- stats::aggregate(prop ~ time,
                            d[d$channel == "Summary" & d$key %in% c("a1", "a2"), ],
                            sum)
  expect_equal(faded$prop[order(faded$time)], own$prop[order(own$time)])
})

make_mcml_seq3 <- function() {
  # 10 actors, 6 states in 3 clusters.
  actors <- rep(1:10, each = 6)
  states <- c("a1", "a2", "b1", "b2", "c1", "c2")
  acts   <- states[((seq_along(actors) * 7L) %% 6L) + 1L]
  times  <- as.POSIXct("2025-01-01", tz = "UTC") + seq_along(actors) * 60
  df <- data.frame(Actor = actors, Action = acts, Time = times,
                   stringsAsFactors = FALSE)
  build_mcml(df,
             clusters = list(G1 = c("a1", "a2"), G2 = c("b1", "b2"),
                             G3 = c("c1", "c2")),
             actor = "Actor", action = "Action", time = "Time", type = "tna")
}

test_that("sequence_plot(combine =) merges clusters into one channel", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()

  p <- sequence_plot(fit, type = "distribution", combine = c("G1", "G2"))
  d <- mcml_panel_data(p)

  # channels: Summary, the merged group (at G1's position), then G3
  expect_equal(levels(d$channel), c("Summary", "G1 + G2", "G3"))
  # the merged panel holds both clusters' states and one faded band for G3
  expect_setequal(unique(as.character(d$key[d$channel == "G1 + G2"])),
                  c("a1", "a2", "b1", "b2", "G3 (Other states)", "NA"))
  # the Summary band keys the group once
  expect_setequal(unique(as.character(d$key[d$channel == "Summary"])),
                  c("G1 + G2", "G3", "NA"))
  # the G3 panel fades the merged group as one band
  expect_true("G1 + G2 (Other states)" %in% as.character(d$key[d$channel == "G3"]))

  # every channel still stacks to 100% at every time point
  tot <- stats::aggregate(prop ~ channel + time, d, sum)
  expect_equal(tot$prop, rep(1, nrow(tot)))

  # invariant: the merged Summary share is the sum of the unmerged G1 and G2
  plain <- mcml_panel_data(sequence_plot(fit, type = "distribution",
                                         panel = "summary"))
  sep <- stats::aggregate(prop ~ time, plain[plain$key %in% c("G1", "G2"), ], sum)
  mer <- d[d$channel == "Summary" & d$key == "G1 + G2", ]
  expect_equal(mer$prop[order(mer$time)], sep$prop[order(sep$time)])
})

test_that("sequence_plot(combine =) takes a named list and works for every type", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()

  p <- sequence_plot(fit, type = "distribution", normalize = TRUE,
                     combine = list(Early = c("G1", "G2")))
  expect_equal(levels(mcml_panel_data(p)$channel), c("Summary", "Early", "G3"))

  idx <- sequence_plot(fit, type = "index", combine = c("G2", "G3"))
  expect_true(inherits(idx, "mcml_sequence_plot") || inherits(idx, "ggplot"))

  # expand is resolved after merging, so it can open the merged group
  open <- sequence_plot(fit, type = "distribution", panel = "summary",
                        combine = c("G1", "G2"), expand = "G1 + G2")
  expect_setequal(levels(open$data$key), c("a1", "a2", "b1", "b2", "G3", "NA"))
})

test_that("sequence_plot(combine =) rejects malformed groups", {
  fit <- make_mcml_seq3()
  expect_error(sequence_plot(fit, combine = c("G1", "nope")), "Unknown")
  expect_error(sequence_plot(fit, combine = "G1"), "at least two")
  expect_error(sequence_plot(fit, combine = list(c("G1", "G2"), c("G2", "G3"))),
               "more than one")
  expect_error(sequence_plot(fit, combine = list(G3 = c("G1", "G2"))),
               "clashes")
  expect_error(sequence_plot(fit, combine = 1:2), "character")
})

test_that("sequence_plot(rest =) sets how a panel shows the other clusters", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()

  dp <- mcml_panel_data(sequence_plot(fit, type = "distribution", rest = "pooled"))
  expect_setequal(unique(as.character(dp$key[dp$channel == "G1"])),
                  c("a1", "a2", "Other states", "NA"))
  tot <- stats::aggregate(prop ~ channel + time, dp, sum)
  expect_equal(tot$prop, rep(1, nrow(tot)))

  # the pooled band equals the sum of the per-cluster faded bands
  dc <- mcml_panel_data(sequence_plot(fit, type = "distribution"))
  per <- stats::aggregate(prop ~ time,
                          dc[dc$channel == "G1" & grepl(" (Other states)", dc$key, fixed = TRUE), ], sum)
  pool <- dp[dp$channel == "G1" & dp$key == "Other states", ]
  expect_equal(pool$prop[order(pool$time)], per$prop[order(per$time)])

  # none: own states only, and their height is the cluster's Summary share
  dn <- mcml_panel_data(sequence_plot(fit, type = "distribution", rest = "none"))
  expect_setequal(unique(as.character(dn$key[dn$channel == "G1"])), c("a1", "a2"))
  own <- stats::aggregate(prop ~ time, dn[dn$channel == "G1", ], sum)
  sm  <- dn[dn$channel == "Summary" & dn$key == "G1", ]
  expect_equal(own$prop[order(own$time)], sm$prop[order(sm$time)])

  # carpet
  ci <- mcml_panel_data(sequence_plot(fit, type = "index", rest = "pooled"))
  expect_true("Other states" %in% as.character(ci$key[ci$channel == "G2"]))
  cn <- mcml_panel_data(sequence_plot(fit, type = "index", rest = "none"))
  expect_false(any(grepl("Other states", cn$key, fixed = TRUE)))

  expect_error(sequence_plot(fit, rest = "bogus"))
})

test_that("carpet wash keys by cluster when the Summary is expanded", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()
  # regression: the wash took Summary keys (states once expanded), producing
  # keys like "a1 (Other states)" with no level, drawn as blank NA fill
  d <- mcml_panel_data(sequence_plot(fit, type = "index", expand = "G1"))
  expect_false(anyNA(d$key))
  expect_true("G1 (Other states)" %in% as.character(d$key[d$channel == "G2"]))
})

test_that("sequence_plot(mcml, na = FALSE) drops the NA band and rescales to running", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()
  d <- mcml_panel_data(sequence_plot(fit, type = "distribution",
                                     rest = "pooled", na = FALSE))
  expect_false("NA" %in% as.character(d$key))
  # every channel stacks to 100% of the running sequences at every time point
  tot <- stats::aggregate(prop ~ channel + time, d, sum)
  expect_equal(tot$prop, rep(1, nrow(tot)))
  # own-state shares are the na = TRUE shares divided by the running share
  a <- mcml_panel_data(sequence_plot(fit, type = "distribution", rest = "pooled"))
  running <- stats::aggregate(prop ~ time, a[a$channel == "Summary" & a$key != "NA", ], sum)
  g1_na  <- a[a$channel == "G1" & a$key == "a1", ]
  g1_run <- d[d$channel == "G1" & d$key == "a1", ]
  expect_equal(g1_run$prop[order(g1_run$time)],
               g1_na$prop[order(g1_na$time)] / running$prop[order(running$time)])
})

test_that("sequence_plot(rest_label =) relabels the other-cluster keys", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()
  dp <- mcml_panel_data(sequence_plot(fit, type = "distribution",
                                      rest = "pooled", rest_label = "Rest of states"))
  expect_true("Rest of states" %in% as.character(dp$key[dp$channel == "G1"]))
  dc <- mcml_panel_data(sequence_plot(fit, type = "distribution", rest_label = "Others"))
  expect_true(all(c("G2 (Others)", "G3 (Others)") %in%
                    as.character(dc$key[dc$channel == "G1"])))
  ci <- mcml_panel_data(sequence_plot(fit, type = "index", rest = "pooled",
                                      rest_label = "Others"))
  expect_true("Others" %in% as.character(ci$key))
  expect_false(anyNA(ci$key))

  expect_error(sequence_plot(fit, rest_label = "a1"), "already a state")
  expect_error(sequence_plot(fit, rest_label = ""), "non-empty")
  expect_error(sequence_plot(fit, rest_label = c("A", "B")), "single")
})

test_that("a single-channel mcml draws in every view", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()
  all3 <- c("G1", "G2", "G3")
  lapply(c("index", "heatmap", "distribution"), function(tp) {
    p <- sequence_plot(fit, type = tp, combine = all3)
    expect_true(inherits(p, "mcml_sequence_plot") || inherits(p, "ggplot"))
    expect_silent(print(p))
  })
})

test_that("expanded clusters keep distinct faded colours", {
  fit <- make_mcml_seq3()
  pal <- .mcml_seq_palettes(.mcml_seq_channels(fit, expand = "all"), NULL)
  expect_identical(anyDuplicated(unname(pal$faded)), 0L)
})

test_that("mcml-only arguments error on other input; NA label message is clear", {
  seqs <- data.frame(T1 = c("a", "b"), T2 = c("b", "a"))
  expect_error(sequence_plot(seqs, combine = c("x", "y")), "apply to an mcml")
  expect_error(sequence_plot(seqs, rest = "pooled"), "apply to an mcml")
  fit <- make_mcml_seq3()
  expect_error(sequence_plot(fit, rest_label = "NA"), "reserved")
})

# ---- named state_colors -----------------------------------------------------
# Every fill key of a panel with its resolved colour, as a named character
# vector: legend label -> fill. Reads the built guide, so it is the colour the
# figure actually draws, not the palette we handed in.
mcml_fills <- function(p, panel) {
  g <- ggplot2::get_guide_data(mcml_panels(p)[[panel]], "fill")
  stats::setNames(g$fill, g$.label)
}

test_that("named state_colors colours a combined group everywhere it appears", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()

  p <- sequence_plot(fit, type = "distribution",
                     combine = list(Merged = c("G1", "G2")),
                     state_colors = c(Merged = "#CC79A7", b1 = "#000000"))

  # the group's Summary band, and one of its member states
  expect_equal(unname(mcml_fills(p, "Summary")[["Merged"]]), "#CC79A7")
  expect_equal(unname(mcml_fills(p, "Merged")[["b1"]]), "#000000")
  # the faded band in the other panel is the 22% fade of the same colour
  m <- 0.22 * (grDevices::col2rgb("#CC79A7") / 255) + 0.78
  expect_equal(unname(mcml_fills(p, "G3")[["Merged (Other states)"]]),
               toupper(grDevices::rgb(m[1L], m[2L], m[3L])))
  # keys the user did not name keep their defaults
  expect_equal(unname(mcml_fills(p, "Merged")[["a1"]]), .okabe_ito[1L])
  expect_equal(unname(mcml_fills(p, "Summary")[["G3"]]), "#2A9D8F")
})

test_that("named state_colors reaches an expanded cluster and rest_label", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()

  # `expand` takes G1's Summary key away, but its channel and faded band remain
  p <- sequence_plot(fit, type = "distribution", expand = "G1",
                     state_colors = c(G1 = "#CC79A7"))
  m <- 0.22 * (grDevices::col2rgb("#CC79A7") / 255) + 0.78
  expect_equal(unname(mcml_fills(p, "G2")[["G1 (Other states)"]]),
               toupper(grDevices::rgb(m[1L], m[2L], m[3L])))

  pooled <- sequence_plot(fit, type = "distribution", rest = "pooled",
                          state_colors = c("Other states" = "#F0E442"))
  expect_equal(unname(mcml_fills(pooled, "G1")[["Other states"]]), "#F0E442")
})

test_that("an unnamed state_colors stays positional for an mcml", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()
  cols <- c("#111111", "#222222", "#333333", "#444444", "#555555", "#666666")

  p <- sequence_plot(fit, type = "distribution", state_colors = cols)
  # states are keyed in sorted order: a1 a2 b1 b2 c1 c2
  expect_equal(unname(mcml_fills(p, "G1")[c("a1", "a2")]), cols[1:2])
  expect_equal(unname(mcml_fills(p, "G3")[c("c1", "c2")]), cols[5:6])
  # clusters keep their own palette, untouched by a positional vector
  expect_equal(unname(mcml_fills(p, "Summary")[["G1"]]), "#264653")
})

test_that("names this figure does not draw are dropped with a message", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()

  expect_message(sequence_plot(fit, state_colors = c(nope = "red")),
                 "no name matches a key this plot draws")
  # the label of a group that was never combined is not a key either
  expect_message(sequence_plot(fit, state_colors = c("G1 + G2" = "red")),
                 "no name matches a key this plot draws")
  # a partial overlap says how many were dropped, and still draws the rest
  expect_message(sequence_plot(fit, state_colors = c(nope = "red", a1 = "#000000")),
                 "1 of 2 names are not drawn")
  # a palette that matches exactly says nothing
  expect_silent(sequence_plot(fit, combine = c("G1", "G2"),
                              state_colors = c("G1 + G2" = "red")))
})

test_that("a palette carrying keys this figure does not draw is accepted", {
  skip_if_not_installed("ggplot2")
  fit <- make_mcml_seq3()
  # One project-wide palette: states of another coding scheme, clusters this
  # figure merged away, and the states it does draw.
  master <- c(Approve = "#2CA02C", Encourage = "#B5D334", a1 = "#000000",
              G3 = "#CC79A7", "G1 + G2" = "#0072B2")

  expect_message(
    p <- sequence_plot(fit, type = "distribution", combine = c("G1", "G2"),
                       state_colors = master),
    "2 of 5 names are not drawn")

  expect_equal(unname(mcml_fills(p, "Summary")[["G1 + G2"]]), "#0072B2")
  expect_equal(unname(mcml_fills(p, "Summary")[["G3"]]), "#CC79A7")
  expect_equal(unname(mcml_fills(p, "G1 + G2")[["a1"]]), "#000000")
  # the keys it does not name are dealt the Okabe-Ito colours the palette did
  # not use: a1 is pinned, so a2 takes the first spare one
  expect_equal(unname(mcml_fills(p, "G1 + G2")[["a2"]]), .okabe_ito[1L])
})
