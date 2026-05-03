testthat::skip_on_cran()

# ==============================================================================
# Tests for the unified clustering print methods:
#   - print.net_clustering           (R/cluster_data.R)
#   - print.net_mmm                  (R/mmm.R)
#   - print.net_mmm_clustering       (R/mmm.R)
#   - print.netobject_group          (R/build_network.R)
#
# Goals:
#   1. Headers / column shapes match the agreed layout.
#   2. New `digits` arg is honoured but optional (non-breaking).
#   3. Old behaviours preserved (returns invisibly, no error on minimal input).
#   4. Regression: cluster_network/cluster_mmm output flows through
#      sequence_plot()/distribution_plot() without length mismatch (the
#      .extract_seqplot_input bug fixed in this change).
# ==============================================================================

# ---- Synthetic data ----------------------------------------------------------

.make_clust_data <- function(n_per = 20, n_cols = 8, seed = 7) {
  set.seed(seed)
  g1 <- matrix(sample(c("A", "B", "C"), n_per * n_cols, replace = TRUE,
                      prob = c(0.7, 0.2, 0.1)), ncol = n_cols)
  g2 <- matrix(sample(c("A", "B", "C"), n_per * n_cols, replace = TRUE,
                      prob = c(0.1, 0.2, 0.7)), ncol = n_cols)
  out <- as.data.frame(rbind(g1, g2), stringsAsFactors = FALSE)
  colnames(out) <- paste0("T", seq_len(n_cols))
  out
}

# ==============================================================================
# print.net_clustering
# ==============================================================================

test_that("print.net_clustering: header shape and dimension line", {
  d <- .make_clust_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2", dissimilarity = "hamming")
  out <- capture.output(print(cl))

  expect_match(out[1L], "^Sequence Clustering \\[ward\\.D2\\]")
  expect_true(any(grepl("Sequences: 40", out)))
  expect_true(any(grepl("Clusters: 2", out)))
  expect_true(any(grepl("Dissimilarity: hamming", out)))
  expect_true(any(grepl("silhouette = ", out)))
})

test_that("print.net_clustering: cluster table has Cluster, N, sizes sum to total", {
  d <- .make_clust_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")
  out <- capture.output(print(cl))

  hdr <- grep("^\\s*Cluster\\b", out)
  expect_length(hdr, 1L)
  expect_true(grepl("\\bN\\b", out[hdr]))

  # extract numeric sizes from the "N (pct%)" column on the next 2 rows
  body <- out[seq.int(hdr + 1L, hdr + 2L)]
  sizes <- as.integer(sub("^\\s*\\d+\\s+(\\d+)\\s+.*", "\\1", body))
  expect_length(sizes, 2L)
  expect_equal(sum(sizes), 40L)
})

test_that("print.net_clustering: PAM emits Medoid column", {
  d <- .make_clust_data()
  cl <- build_clusters(d, k = 2, method = "pam")
  out <- capture.output(print(cl))
  hdr <- out[grep("^\\s*Cluster\\b", out)]
  expect_true(grepl("Medoid", hdr))
})

test_that("print.net_clustering: digits arg controls float decimals", {
  d <- .make_clust_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")

  o3 <- capture.output(print(cl, digits = 3))
  o6 <- capture.output(print(cl, digits = 6))

  # The silhouette line should differ between digits = 3 and digits = 6
  s3 <- o3[grep("silhouette", o3)]
  s6 <- o6[grep("silhouette", o6)]
  expect_false(identical(s3, s6))
  expect_match(s3, "silhouette = -?\\d+\\.\\d{3}")
  expect_match(s6, "silhouette = -?\\d+\\.\\d{6}")
})

test_that("print.net_clustering: returns invisibly with no extra args (non-breaking)", {
  d <- .make_clust_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")
  res <- withVisible(print(cl))
  expect_identical(res$value, cl)
  expect_false(res$visible)
})

# ==============================================================================
# print.net_mmm
# ==============================================================================

test_that("print.net_mmm: header / IC line / no leftover '%%'", {
  d <- .make_clust_data()
  mmm <- build_mmm(d, k = 2, n_starts = 2, max_iter = 30, seed = 1)
  out <- capture.output(print(mmm))

  expect_match(out[1L], "^Mixed Markov Model$")
  expect_true(any(grepl("Sequences: 40", out)))
  expect_true(any(grepl("Clusters: 2", out)))
  expect_true(any(grepl("States: 3", out)))
  expect_true(any(grepl("BIC = ", out)))
  expect_true(any(grepl("AIC = ", out)))
  expect_true(any(grepl("ICL = ", out)))

  # Regression for the literal "%%" bug in the previous print method.
  expect_false(any(grepl("%%", out)))
})

test_that("print.net_mmm: cluster table shape (Cluster, N, Mix%, AvePP)", {
  d <- .make_clust_data()
  mmm <- build_mmm(d, k = 2, n_starts = 2, max_iter = 30, seed = 1)
  out <- capture.output(print(mmm))

  hdr <- grep("^\\s*Cluster\\b", out)
  expect_length(hdr, 1L)
  expect_true(grepl("\\bN\\b",     out[hdr]))
  expect_true(grepl("Mix%",        out[hdr]))
  expect_true(grepl("AvePP",       out[hdr]))

  body <- out[seq.int(hdr + 1L, hdr + 2L)]
  sizes <- as.integer(sub("^\\s*\\d+\\s+(\\d+)\\s+.*", "\\1", body))
  expect_equal(sum(sizes), 40L)
})

test_that("print.net_mmm: digits arg honoured + invisible return", {
  d <- .make_clust_data()
  mmm <- build_mmm(d, k = 2, n_starts = 2, max_iter = 30, seed = 1)

  out6 <- capture.output(print(mmm, digits = 6))
  ll_line <- out6[grep("BIC =", out6)]
  expect_match(ll_line, "BIC = -?\\d+\\.\\d{6}")

  res <- withVisible(print(mmm))
  expect_identical(res$value, mmm)
  expect_false(res$visible)
})

# ==============================================================================
# print.net_mmm_clustering
# ==============================================================================

test_that("print.net_mmm_clustering: header + cluster table", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 2, max_iter = 30, seed = 1)

  expect_s3_class(grp, "netobject_group")
  cl <- attr(grp, "clustering")
  expect_s3_class(cl, "net_mmm_clustering")

  out <- capture.output(print(cl))
  expect_match(out[1L], "^MMM Clustering \\[k = 2\\]")
  expect_true(any(grepl("Sequences: 40", out)))
  expect_true(any(grepl("AvePP =", out)))
  expect_true(any(grepl("BIC =", out)))

  hdr <- grep("^\\s*Cluster\\b", out)
  expect_length(hdr, 1L)
  expect_true(grepl("Mix%",  out[hdr]))
  expect_true(grepl("AvePP", out[hdr]))
})

test_that("cluster_mmm stashes full sequence data on the clustering attribute", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)
  cl <- attr(grp, "clustering")

  # The invariant set by cluster_mmm + .extract_seqplot_input depends on
  # this carrying an N-row data frame, where N matches assignments.
  expect_false(is.null(cl$data))
  expect_equal(NROW(cl$data), length(cl$assignments))
  expect_equal(NROW(cl$data), nrow(d))
})

# ==============================================================================
# print.netobject_group
# ==============================================================================

test_that("print.netobject_group (group_col): legacy header preserved", {
  d <- .make_clust_data()
  d$grp <- rep(c("X", "Y"), each = 20)
  nets <- build_network(d, method = "relative", group = "grp")

  out <- capture.output(print(nets))
  expect_match(out[1L], "^Group Networks \\(2 groups, group_col: grp\\)")

  hdr <- grep("^\\s+Group\\s+Nodes", out)
  expect_length(hdr, 1L)
  expect_true(grepl("Nodes",   out[hdr]))
  expect_true(grepl("Edges",   out[hdr]))
  expect_true(grepl("Weights", out[hdr]))
  # No N column without a clustering attribute.
  expect_false(grepl("\\bN\\b", out[hdr]))
})

test_that("print.netobject_group (cluster_network): clustering source surfaced", {
  d <- .make_clust_data()
  grp <- cluster_network(d, k = 2, cluster_by = "ward.D2",
                         dissimilarity = "hamming")

  out <- capture.output(print(grp))
  expect_match(out[1L],
               "^Group Networks \\(2 clusters via ward\\.D2 / hamming\\)")

  hdr <- grep("^\\s+Group\\s+Nodes", out)
  expect_length(hdr, 1L)
  expect_true(grepl("\\bN\\b", out[hdr]))
})

test_that("print.netobject_group (cluster_mmm): MMM source surfaced", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)

  out <- capture.output(print(grp))
  expect_match(out[1L], "^Group Networks \\(2 clusters from MMM\\)")

  hdr <- grep("^\\s+Group\\s+Nodes", out)
  expect_true(grepl("\\bN\\b", out[hdr]))
})

test_that("print.netobject_group: digits arg + invisible return", {
  d <- .make_clust_data()
  grp <- cluster_network(d, k = 2, cluster_by = "ward.D2")

  out6 <- capture.output(print(grp, digits = 6))
  weights_lines <- out6[grep("\\[\\d", out6)]
  # weight ranges shown as [a, b] with `digits` decimals
  expect_true(any(grepl("\\[\\d+\\.\\d{6}, \\d+\\.\\d{6}\\]", weights_lines)))

  res <- withVisible(print(grp))
  expect_identical(res$value, grp)
  expect_false(res$visible)
})

# ==============================================================================
# Regression: .extract_seqplot_input no longer mismatches data/group length
# ==============================================================================

test_that("sequence_plot accepts cluster_network output without length mismatch", {
  d <- .make_clust_data()
  grp <- cluster_network(d, k = 2, cluster_by = "ward.D2")

  # Previous bug: x[[1]]$data was n1 rows but assignments was 40 ->
  # stopifnot(length(group) == nrow(x)) inside distribution_plot would fire,
  # and the index path would silently take only cluster 1's data.
  expect_silent(p <- sequence_plot(grp, type = "index"))
  expect_silent(p <- distribution_plot(grp))
})

test_that("sequence_plot accepts cluster_mmm output without length mismatch", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)

  expect_silent(p <- sequence_plot(grp, type = "index"))
  expect_silent(p <- distribution_plot(grp))
})

test_that("sequence_plot accepts plain netobject_group without clustering attr", {
  d <- .make_clust_data()
  d$grp <- rep(c("X", "Y"), each = 20)
  nets <- build_network(d, method = "relative", group = "grp")

  # Plain group_col split: helper rebuilds full data via rbind() and labels
  # rows by their group. Should not error.
  expect_silent(p <- sequence_plot(nets, type = "index"))
})
