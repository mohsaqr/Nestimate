# Skipped on CRAN to keep the check under its time limit; runs locally and in CI.
testthat::skip_on_cran()

test_that("sequence_plot draws from a data.frame + hclust", {
  set.seed(1L)
  states <- c("A", "B", "C")
  seqs <- as.data.frame(matrix(sample(states, 30 * 10, replace = TRUE),
                               nrow = 30, ncol = 10))
  d    <- stats::dist(data.matrix(seqs))
  tree <- stats::hclust(d, method = "ward.D2")

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  res <- sequence_plot(seqs, tree = tree)
  expect_type(res, "list")
  expect_named(res, c("ord", "codes", "palette", "levels", "sort_used"))
  expect_length(res$ord, nrow(seqs))
  expect_setequal(res$ord, seq_len(nrow(seqs)))
  expect_equal(dim(res$codes), dim(seqs))
  expect_equal(res$levels, sort(states))
})

test_that("sequence_plot accepts a net_clustering directly", {
  set.seed(2L)
  seqs <- as.data.frame(matrix(sample(c("A", "B", "C", "D"), 25 * 12,
                                      replace = TRUE), 25, 12))
  cl <- build_clusters(seqs, k = 3L, dissimilarity = "hamming",
                     method = "ward.D2")

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  res <- sequence_plot(cl)
  expect_length(res$ord, nrow(seqs))
  expect_identical(res$sort_used, "net_clustering")
})

test_that("sequence_plot sort='frequency' produces a dendrogram order", {
  set.seed(3L)
  seqs <- matrix(sample(c("A", "B", "C"), 15 * 6, replace = TRUE), 15, 6)

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  res <- sequence_plot(seqs, sort = "frequency")
  expect_setequal(res$ord, seq_len(nrow(seqs)))
  expect_identical(res$sort_used, "frequency")
})

test_that("sequence_plot sort='start' sorts lexicographically forward", {
  seqs <- rbind(c("A", "B", "C"),
                c("B", "A", "A"),
                c("A", "A", "B"),
                c("B", "B", "A"))

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  res <- sequence_plot(seqs, sort = "start", legend = "none")
  # After ordering, first column should be non-decreasing.
  first_col <- res$codes[res$ord, 1]
  expect_true(all(diff(first_col) >= 0))
  expect_identical(res$sort_used, "start")
})

test_that("sequence_plot sort='end' sorts lexicographically backward", {
  seqs <- rbind(c("A", "A", "B"),
                c("B", "B", "A"),
                c("A", "B", "A"),
                c("B", "A", "B"))

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  res <- sequence_plot(seqs, sort = "end", legend = "none")
  last_col <- res$codes[res$ord, ncol(seqs)]
  expect_true(all(diff(last_col) >= 0))
})

test_that("sequence_plot routes distance sorts through build_clusters", {
  set.seed(4L)
  seqs <- matrix(sample(c("A", "B", "C"), 20 * 8, replace = TRUE), 20, 8)

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  for (metric in c("hamming", "lcs", "lv", "osa", "dl",
                   "qgram", "cosine", "jaccard", "jw")) {
    res <- sequence_plot(seqs, sort = metric, legend = "none")
    expect_identical(res$sort_used, metric, info = metric)
    expect_length(res$ord, nrow(seqs))
  }
})

test_that("sequence_plot preserves NA cells", {
  set.seed(5L)
  seqs <- matrix(sample(c("A", "B", "C"), 20 * 8, replace = TRUE), 20, 8)
  seqs[sample(length(seqs), 10)] <- NA

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  res <- sequence_plot(seqs, sort = "start", na_color = "#FF0000")
  expect_equal(sum(is.na(res$codes)), sum(is.na(seqs)))
})

test_that("sequence_plot honours all legend_position values", {
  set.seed(6L)
  seqs <- matrix(sample(c("A", "B"), 12 * 5, replace = TRUE), 12, 5)

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  for (pos in c("bottom", "right", "none")) {
    expect_silent(sequence_plot(seqs, legend = pos,
                                sort = "start"))  # no tree ⇒ exercises no-tree layout
    expect_silent(sequence_plot(seqs, legend = pos,
                                sort = "frequency"))  # tree ⇒ exercises tree layout
  }
})

test_that("sequence_plot rejects mismatched tree size", {
  seqs <- matrix(sample(c("A", "B"), 20, replace = TRUE), 10, 2)
  tree <- stats::hclust(stats::dist(matrix(stats::rnorm(16), 8)))
  expect_error(sequence_plot(seqs, tree = tree), "leaves")
})

test_that("sequence_plot rejects all-NA input", {
  seqs <- matrix(NA_character_, 5, 4)
  tree <- stats::hclust(stats::dist(matrix(stats::rnorm(10), 5)))
  expect_error(sequence_plot(seqs, tree = tree), "no non-NA values")
})

test_that("sequence_plot type='index' renders single-panel with gaps", {
  set.seed(10L)
  seqs <- matrix(sample(c("A", "B", "C"), 20 * 8, replace = TRUE), 20, 8)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  res <- sequence_plot(seqs, type = "index", row_gap = 0.2)
  expect_named(res, c("codes", "palette", "levels", "orders", "groups"))
  expect_length(res$orders, 1L)
  expect_identical(res$groups, "all")
})

test_that("sequence_plot type='index' facets by net_clustering", {
  set.seed(11L)
  seqs <- as.data.frame(matrix(sample(c("A", "B", "C"), 30 * 10,
                                      replace = TRUE), 30, 10))
  cl <- build_clusters(seqs, k = 3L, dissimilarity = "hamming",
                     method = "ward.D2")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  res <- sequence_plot(cl, type = "index")
  expect_length(res$orders, 3L)
  expect_length(res$groups, 3L)
})

test_that("sequence_plot type='distribution' dispatches to distribution_plot", {
  set.seed(12L)
  seqs <- matrix(sample(c("A", "B", "C"), 24 * 6, replace = TRUE), 24, 6)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  res <- sequence_plot(seqs, type = "distribution")
  expect_named(res, c("counts", "proportions", "levels", "palette", "groups"))
})

test_that("sequence_plot type='index' honours ncol/nrow", {
  set.seed(13L)
  seqs <- as.data.frame(matrix(sample(c("A", "B", "C", "D"), 40 * 8,
                                      replace = TRUE), 40, 8))
  cl <- build_clusters(seqs, k = 4L, dissimilarity = "hamming",
                     method = "ward.D2")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  # 4 clusters in 1x4
  expect_silent(sequence_plot(cl, type = "index", ncol = 4, nrow = 1))
  # 4 clusters in 2x2
  expect_silent(sequence_plot(cl, type = "index", ncol = 2, nrow = 2))
})

# ---- Branch-matrix coverage (task #17) ----
# Crosses type x sort for heatmap/index (sort is ignored for distribution).
# Existing tests hit each axis individually; this test guards against
# specific combinations regressing — e.g. a new sort branch that only ever
# ran with type="heatmap" silently erroring under type="index".

test_that("sequence_plot branch matrix: type x sort combinations render", {
  set.seed(17)
  seqs <- matrix(sample(c("A", "B", "C"), 20 * 6, replace = TRUE), 20, 6)
  # Mix of category-based and distance-based sorts. Distance sorts funnel
  # through build_clusters, so they cover a different backend.
  sorts <- c("lcs", "frequency", "start", "end", "hamming", "lv", "jaccard")

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  for (tp in c("heatmap", "index")) {
    for (s in sorts) {
      info <- sprintf("type=%s sort=%s", tp, s)
      res <- tryCatch(
        sequence_plot(seqs, type = tp, sort = s, legend = "none"),
        error = function(e) e
      )
      if (inherits(res, "error")) {
        fail(sprintf("%s -> %s", info, conditionMessage(res)))
        next
      }
      expect_type(res, "list")
      # Every valid sort must yield an ordering covering all rows exactly once
      ord <- if (tp == "heatmap") res$ord else unlist(res$orders)
      expect_setequal(ord, seq_len(nrow(seqs)))
    }
  }
})

test_that("sequence_plot distribution branch matrix: scale x geom combinations", {
  set.seed(17)
  seqs <- matrix(sample(c("A", "B", "C"), 24 * 6, replace = TRUE), 24, 6)
  pdf(NULL); on.exit(dev.off(), add = TRUE)

  for (scale in c("proportion", "count")) {
    for (geom in c("area", "bar")) {
      info <- sprintf("scale=%s geom=%s", scale, geom)
      res <- tryCatch(
        sequence_plot(seqs, type = "distribution", scale = scale, geom = geom),
        error = function(e) e
      )
      if (inherits(res, "error")) {
        fail(sprintf("%s -> %s", info, conditionMessage(res)))
        next
      }
      expect_type(res, "list")
      expect_true(!is.null(res$counts) || !is.null(res$proportions), info = info)
    }
  }
})

test_that("trim truncates the time axis by quantile and absolute cut", {
  # 95 short sequences (<=20 cols) + 5 long ones (200 cols).
  set.seed(1)
  n_row <- 100L; n_col <- 200L
  states <- c("A", "B", "C")
  m <- matrix(NA_character_, n_row, n_col)
  lens <- c(sample(10:20, 95, replace = TRUE), rep(200L, 5L))
  invisible(lapply(seq_len(n_row), function(i)
    m[i, seq_len(lens[i])] <<- sample(states, lens[i], replace = TRUE)))
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  colnames(df) <- paste0("T", seq_len(n_col))

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  # NULL trim plots full width.
  full <- sequence_plot(df, type = "index", trim = NULL, legend = "none")
  expect_equal(ncol(full$codes), n_col)

  # Quantile trim drops the long tail; column names are preserved.
  q <- sequence_plot(df, type = "index", trim = 0.95, legend = "none")
  expect_lt(ncol(q$codes), n_col)
  expect_identical(colnames(q$codes), paste0("T", seq_len(ncol(q$codes))))

  # Absolute cut keeps exactly the first t columns.
  abs_cut <- sequence_plot(df, type = "index", trim = 30, legend = "none")
  expect_equal(ncol(abs_cut$codes), 30L)

  # Heatmap honours trim too (z is transposed: rows = time).
  hm <- sequence_plot(df, type = "heatmap", trim = 30, legend = "none")
  expect_equal(ncol(hm$codes), 30L)

  # Cut beyond the data is a harmless no-op.
  over <- sequence_plot(df, type = "index", trim = 999, legend = "none")
  expect_equal(ncol(over$codes), n_col)

  # Invalid trim errors.
  expect_error(sequence_plot(df, type = "index", trim = -1))
})

test_that("the row half of trim ends each row earlier and removes none", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  set.seed(7)
  # 20 sequences of 30 points; rows 1-4 hold 2 states, the rest 30.
  df <- as.data.frame(matrix(sample(c("A", "B", "C"), 600, TRUE), 20, 30))
  df[1:4, 3:30] <- NA
  observed <- function(res) unname(rowSums(!is.na(res$codes)))
  idx <- sequence_plot(df, type = "index", trim = c(NA, 5), sort = "start",
                       legend = "none")
  expect_identical(dim(idx$codes), c(20L, 30L))
  expect_setequal(unlist(idx$orders), 1:20)
  expect_identical(sort(observed(idx)), c(rep(2, 4L), rep(5, 16L)))
  # A fraction: ceiling(0.5 * 2) = 1 and ceiling(0.5 * 30) = 15.
  frac <- sequence_plot(df, type = "index", trim = c(NA, 0.5), legend = "none")
  expect_identical(sort(observed(frac)), c(rep(1, 4L), rep(15, 16L)))
  # Counted over the kept time points: 10 columns, 0.5 of 10 is 5.
  both <- sequence_plot(df, type = "index", trim = c(10, 0.5), legend = "none")
  expect_identical(dim(both$codes), c(20L, 10L))
  expect_identical(sort(observed(both)), c(rep(1, 4L), rep(5, 16L)))
  # A count past a row's length changes nothing.
  expect_identical(sequence_plot(df, type = "index", trim = c(NA, 999),
                                 legend = "none")$codes,
                   sequence_plot(df, type = "index", legend = "none")$codes)
  # Heatmap and distribution take the same vector.
  hm <- sequence_plot(df, type = "heatmap", trim = c(NA, 5), legend = "none")
  expect_identical(nrow(hm$codes), 20L)
  expect_identical(sort(observed(hm)), c(rep(2, 4L), rep(5, 16L)))
  dp <- sequence_plot(df, type = "distribution", trim = c(NA, 5),
                      legend = "none")
  states <- c("A", "B", "C")
  expect_identical(sum(dp$counts[[1L]][states, 5L]), 16)
  expect_identical(sum(dp$counts[[1L]][states, 6L]), 0)
  # Grouped: every panel stays, and each counts over its own time window.
  grp <- rep(c("sparse", "dense"), c(4L, 16L))
  g <- sequence_plot(df, type = "index", group = grp, trim = c(NA, 5),
                     legend = "none")
  expect_setequal(g$groups, c("sparse", "dense"))
  expect_setequal(unlist(g$orders), 1:20)
  cw <- as.data.frame(matrix("A", 4, 100))
  cw[1:2, 3:100] <- NA
  cw[3, 11:100] <- NA
  cw_grp <- rep(c("short", "long"), each = 2L)
  cw_res <- sequence_plot(cw, type = "index", group = cw_grp,
                          trim = c(0.5, 0.5), trim_clusterwise = TRUE,
                          legend = "none")
  expect_setequal(cw_res$groups, c("short", "long"))
  expect_setequal(unlist(cw_res$orders), 1:4)
  # A tree has one leaf per row, and every row stays.
  tree <- stats::hclust(stats::dist(seq_len(20)))
  expect_identical(nrow(sequence_plot(df, type = "heatmap", tree = tree,
                                      trim = c(NA, 5), legend = "none")$codes),
                   20L)
  expect_error(sequence_plot(df, type = "index", trim = c(1, 2, 3)))
})

test_that(".trim_rows_blank keeps the first states and moves nothing", {
  codes <- rbind(c(1, NA, 2, 3), c(NA, NA, 4, NA), c(NA, NA, NA, NA))
  expect_identical(.trim_rows_blank(codes, NULL), codes)
  expect_identical(.trim_rows_blank(codes, 2),
                   rbind(c(1, NA, 2, NA), c(NA, NA, 4, NA), rep(NA_real_, 4L)))
  # ceiling(0.5 * 3) = 2, ceiling(0.5 * 1) = 1, and an empty row stays empty.
  expect_identical(.trim_rows_blank(codes, 0.5), .trim_rows_blank(codes, 2))
  one_col <- matrix(c("a", NA, "b"), ncol = 1L)
  expect_identical(.trim_rows_blank(one_col, 1), one_col)
  expect_error(.trim_rows_blank(codes, -1))
})

test_that("trim applies to the distribution and mcml paths", {
  set.seed(1)
  n_row <- 100L; n_col <- 200L
  states <- c("A", "B", "C")
  m <- matrix(NA_character_, n_row, n_col)
  lens <- c(sample(10:20, 95, replace = TRUE), rep(200L, 5L))
  invisible(lapply(seq_len(n_row), function(i)
    m[i, seq_len(lens[i])] <<- sample(states, lens[i], replace = TRUE)))
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  colnames(df) <- paste0("T", seq_len(n_col))

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  # distribution: absolute and quantile cut, via both entry points.
  d_full <- distribution_plot(df, trim = NULL, legend = "none")
  expect_equal(ncol(d_full$counts[[1L]]), n_col)
  d_abs <- distribution_plot(df, trim = 30, legend = "none")
  expect_equal(ncol(d_abs$counts[[1L]]), 30L)
  d_q <- sequence_plot(df, type = "distribution", trim = 0.95, legend = "none")
  expect_lt(ncol(d_q$counts[[1L]]), n_col)

  # mcml: every channel cut to the same width as the (trimmed) time axis.
  skip_if_not(exists("build_mcml"))
  fit <- build_mcml(
    group_regulation_long,
    clusters = list(Cognitive  = c("discuss", "synthesis", "consensus", "cohesion"),
                    Regulation = c("plan", "monitor", "adapt", "coregulate"),
                    Affective  = "emotion"),
    actor = "Actor", action = "Action", time = "Time")

  ch_full <- Nestimate:::.mcml_seq_channels(fit, NULL)
  ch_abs  <- Nestimate:::.mcml_seq_channels(fit, 15)
  expect_lt(length(ch_abs$times), length(ch_full$times))
  expect_equal(length(ch_abs$times), 15L)
  expect_true(all(vapply(ch_abs$cmats, ncol, 1L) == length(ch_abs$times)))
  expect_equal(ncol(ch_abs$summary_mat), 15L)

  expect_s3_class(sequence_plot(fit, type = "index", trim = 15), "mcml_sequence_plot")
  expect_s3_class(sequence_plot(fit, type = "distribution", trim = 0.9), "mcml_sequence_plot")

  # The row half ends each coalesced sequence earlier in every channel; rows stay.
  ch_rows <- Nestimate:::.mcml_seq_channels(fit, c(NA, 0.5))
  n_rows  <- nrow(ch_full$summary_mat)
  expect_identical(nrow(ch_rows$summary_mat), n_rows)
  expect_lt(sum(!is.na(ch_rows$summary_mat)), sum(!is.na(ch_full$summary_mat)))
  expect_identical(is.na(ch_rows$summary_cluster_mat), is.na(ch_rows$summary_mat))
  channel_cells <- Reduce(`+`, lapply(ch_rows$cmats, function(m) !is.na(m)))
  expect_identical(channel_cells > 0L, !is.na(ch_rows$summary_mat))
  expect_s3_class(sequence_plot(fit, type = "index", trim = c(15, 0.5)),
                  "mcml_sequence_plot")
})

test_that("trim_clusterwise: global keeps panels aligned, clusterwise crops per group", {
  set.seed(1)
  mk <- function(lo, hi, n) {
    m <- matrix(NA_character_, n, 100L)
    L <- sample(lo:hi, n, replace = TRUE)
    invisible(lapply(seq_len(n), function(i)
      m[i, seq_len(L[i])] <<- sample(c("A", "B", "C"), L[i], replace = TRUE)))
    m
  }
  df <- as.data.frame(rbind(mk(10, 15, 50), mk(40, 50, 50)),
                      stringsAsFactors = FALSE)
  colnames(df) <- paste0("T", seq_len(100))
  grp <- rep(c("Short", "Long"), each = 50L)

  pdf(NULL); on.exit(dev.off(), add = TRUE)

  # Default (global): both groups share one width.
  g <- distribution_plot(df, group = grp, trim = 0.95,
                         trim_clusterwise = FALSE, legend = "none")
  expect_equal(ncol(g$counts$Short), ncol(g$counts$Long))

  # Clusterwise: short group cropped much tighter than the long group.
  c1 <- distribution_plot(df, group = grp, trim = 0.95,
                          trim_clusterwise = TRUE, legend = "none")
  expect_lt(ncol(c1$counts$Short), ncol(c1$counts$Long))

  # Absolute trim is unaffected by the flag.
  a <- distribution_plot(df, group = grp, trim = 20,
                         trim_clusterwise = TRUE, legend = "none")
  expect_equal(ncol(a$counts$Short), 20L)
  expect_equal(ncol(a$counts$Long), 20L)

  # Index path honours it too, and orders still cover all rows.
  ix <- sequence_plot(df, type = "index", group = grp, trim = 0.95,
                      trim_clusterwise = TRUE, legend = "none")
  expect_setequal(unlist(ix$orders), seq_len(nrow(df)))
})

test_that("named state_colors maps by name in every non-mcml type", {
  m <- matrix(c("a", "b", "c", "a",
                "b", "c", "c", "a",
                "b", "b", "c", "a"), nrow = 3, byrow = TRUE)
  pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)

  idx  <- sequence_plot(m, type = "index", state_colors = c(b = "#000000"))
  heat <- sequence_plot(m, type = "heatmap", state_colors = c(c = "#123456"))
  dist <- sequence_plot(m, type = "distribution", state_colors = c(b = "#000000"))

  # The named key is overridden; the states left unnamed are dealt the
  # Okabe-Ito colours the palette did NOT use, one each - so "b" taking
  # #000000 (Okabe-Ito 9) leaves "a" and "c" on 1 and 2, never a repeat.
  expect_equal(idx$palette, c(.okabe_ito[1L], "#000000", .okabe_ito[2L]))
  expect_equal(heat$palette, c(.okabe_ito[1L], .okabe_ito[2L], "#123456"))
  expect_equal(dist$palette[seq_len(3L)], c(.okabe_ito[1L], "#000000",
                                            .okabe_ito[2L]))
  expect_false(anyDuplicated(idx$palette) > 0L)

  # unnamed stays positional (regression: names must not change the old path)
  pos <- sequence_plot(m, type = "index",
                       state_colors = c("#111111", "#222222", "#333333"))
  expect_equal(pos$palette, c("#111111", "#222222", "#333333"))

  # extra names are dropped with a message (one palette, many figures)
  expect_message(super <- sequence_plot(m, type = "index",
                                        state_colors = c(zz = "red", b = "#000000")),
                 "1 of 2 names are not drawn")
  expect_equal(super$palette, c(.okabe_ito[1L], "#000000", .okabe_ito[2L]))
  expect_message(sequence_plot(m, state_colors = c(zz = "red")),
                 "no name matches a key this plot draws")
})

test_that("a half-named state_colors is rejected, not half-honoured", {
  m <- matrix(c("a", "b", "c", "a", "b", "c"), nrow = 2, byrow = TRUE)
  pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(sequence_plot(m, state_colors = c("#111111", b = "#222222")),
               "mixes named and unnamed colours")
})

test_that("a fractional column trim judges each column on its own activity", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  d <- rbind(c("A", NA, "B", "C"), c("A", NA, "B", "C"), c("A", "X", "B", NA))
  colnames(d) <- paste0("T", 1:4)
  # Activity 3, 1, 3, 2: half of four columns keeps the two busiest.
  expect_identical(.trim_keep(d, 0.5), c(1L, 3L))
  expect_identical(colnames(sequence_plot(d, type = "index", trim = c(0.5, NA),
                                          legend = "none")$codes),
                   c("T1", "T3"))
  expect_identical(colnames(sequence_plot(d, type = "heatmap", trim = c(0.5, NA),
                                          legend = "none")$codes),
                   c("T1", "T3"))
  counts <- distribution_plot(d, trim = c(0.5, NA), legend = "none")$counts[[1L]]
  expect_identical(colnames(counts), c("T1", "T3"))
  expect_equal(unname(colSums(counts)), c(3, 3))

  # Ties go to the earlier column; a count still keeps the first columns.
  expect_identical(.trim_keep(d, 0.25), 1L)
  expect_identical(.trim_keep(d, 2), 1:2)
  expect_identical(.trim_keep(d, 99), 1:4)
  # A column with no data goes even when the share rounds up to all of them.
  e <- d; e[, "T2"] <- NA
  expect_identical(.trim_keep(e, 0.99), c(1L, 3L, 4L))
  # Nothing observed still leaves one column to draw.
  expect_identical(.trim_keep(matrix(NA_character_, 2, 3), 0.5), 1L)

  # Unnamed columns label by original position, not by their new place.
  expect_identical(colnames(.trim_codes(unname(d), 0.5)), c("1", "3"))
  expect_null(colnames(.trim_codes(unname(d), NULL)))
  grouped <- sequence_plot(d, type = "index", group = c("a", "a", "b"),
                           trim = c(0.5, NA), legend = "none")
  expect_identical(colnames(grouped$codes), c("T1", "T3"))
  # The input is not changed: only the drawing is.
  expect_identical(colnames(d), paste0("T", 1:4))
})

test_that("clusterwise column trim judges the columns within each group", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  d <- rbind(c("A", "B", "A", NA), c("B", "B", "A", NA),
             c("A", NA, "A", "B"), c("B", NA, "B", "B"))
  colnames(d) <- paste0("T", 1:4)
  res <- distribution_plot(d, group = c("x", "x", "y", "y"), trim = 0.6,
                           trim_clusterwise = TRUE, legend = "none")
  # Each group has three columns with data and keeps the first two of them.
  expect_identical(colnames(res$counts$x), c("T1", "T2"))
  expect_identical(colnames(res$counts$y), c("T1", "T3"))
})

test_that("the mcml axis labels the columns a fractional trim keeps", {
  label <- .mcml_time_labels(c(1L, 3L))
  expect_identical(label(c(1, 1.5, 2, NA, 3)), c("1", "", "3", "", ""))
  expect_s3_class(.mcml_time_labels(1:3), "waiver")
})

test_that("the row half of trim leaves the row order alone", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  d <- rbind(c("A", "Z"), c("B", "A"), c("A", "B"))
  order_of <- function(type, trim, group = NULL) {
    res <- sequence_plot(d, type = type, sort = "end", trim = trim, group = group,
                         legend = "none")
    unlist(res$orders %||% res$ord)
  }
  expect_identical(order_of("index", c(NA, 1)), order_of("index", NULL))
  expect_identical(order_of("heatmap", c(NA, 1)), order_of("heatmap", NULL))
  grp <- c("g", "g", "g")
  expect_identical(order_of("index", c(NA, 1), grp), order_of("index", NULL, grp))
  # Every row stays and is blank after its first state.
  codes <- sequence_plot(d, type = "index", sort = "end", trim = c(NA, 1),
                         legend = "none")$codes
  expect_identical(dim(codes), dim(d))
  expect_true(all(is.na(codes[, 2L])))
})

test_that("a fractional column trim removes a long tail", {
  # Ten rows of length 2, one of which runs on alone for eight more columns:
  # 28 states, with 20 in the first two columns.
  d <- matrix(NA_character_, 10, 10)
  d[, 1:2] <- "A"
  d[10L, ] <- "B"
  counts <- colSums(!is.na(d))
  expect_identical(unname(counts), c(10, 10, rep(1, 8)))
  # 95% of 28 is 26.6: two busy columns hold 20, then seven tail columns.
  expect_identical(.trim_keep(d, 0.95), 1:9)
  # 70% of 28 is 19.6: the first two columns are enough.
  expect_identical(.trim_keep(d, 0.7), 1:2)
  # A share of exactly what two columns hold needs no third column.
  expect_identical(.trim_keep(d, 20 / 28), 1:2)
  # The share kept never falls short of the share asked for.
  invisible(lapply(c(0.1, 0.33, 0.5, 0.8, 0.99), function(share) {
    expect_gte(sum(counts[.trim_keep(d, share)]) / sum(counts), share - 1e-8)
  }))
})
