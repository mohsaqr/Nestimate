testthat::skip_on_cran()

# ==============================================================================
# Tests for cluster_choice() and the cluster_choice class.
#
# Coverage:
#   1. Single-axis sweep (k only) returns expected columns and row count.
#   2. Dissimilarity sweep (fixed k) and method sweep (fixed k).
#   3. Cartesian product of k x dissimilarity gives length(k) * length(diss)
#      rows and the right columns vary across them.
#   4. "all" sentinel expands dissimilarity / method to the full canonical
#      list (9 entries each).
#   5. Per-row math: silhouette matches a manual build_clusters() call;
#      sizes / size_ratio are consistent.
#   6. Print: header adapts (single-axis -> [<method> / <dissim>]; multi-axis
#      -> "sweep: ...") and the silhouette-best row is marked.
#   7. Summary returns a data.frame with $best populated on the right row.
#   8. Plot returns a ggplot for each shape (k-only, dissim-only, k x dissim).
#   9. weighted = TRUE rejected up-front when sweeping non-hamming dissims.
# ==============================================================================

# ---- Synthetic data ----------------------------------------------------------

.make_choice_data <- function(n_per = 20, n_cols = 6, seed = 11) {
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
# Sweep shapes
# ==============================================================================

test_that("k-only sweep returns one row per k with all metric columns", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 2:4, method = "ward.D2",
                        dissimilarity = "hamming")

  expect_s3_class(ch, "cluster_choice")
  expect_s3_class(ch, "data.frame")
  expect_equal(nrow(ch), 3L)
  expect_setequal(names(ch),
                  c("k", "dissimilarity", "method",
                    "silhouette", "mean_within_dist",
                    "min_size", "max_size", "size_ratio"))
  expect_equal(ch$k, 2:4)
  expect_true(all(ch$dissimilarity == "hamming"))
  expect_true(all(ch$method == "ward.D2"))
})

test_that("dissimilarity-only sweep returns one row per dissimilarity", {
  d <- .make_choice_data()
  diss <- c("hamming", "lcs", "jaccard")
  ch <- cluster_choice(d, k = 3, dissimilarity = diss, method = "ward.D2")

  expect_equal(nrow(ch), length(diss))
  expect_setequal(ch$dissimilarity, diss)
  expect_true(all(ch$k == 3))
})

test_that("method-only sweep returns one row per method", {
  d <- .make_choice_data()
  meths <- c("ward.D2", "complete", "average")
  ch <- cluster_choice(d, k = 3, method = meths,
                        dissimilarity = "hamming")
  expect_equal(nrow(ch), length(meths))
  expect_setequal(ch$method, meths)
})

test_that("k x dissimilarity grid produces the cartesian product", {
  d <- .make_choice_data()
  diss <- c("hamming", "lcs")
  ks <- 2:4
  ch <- cluster_choice(d, k = ks, dissimilarity = diss,
                        method = "ward.D2")

  expect_equal(nrow(ch), length(ks) * length(diss))
  # Each (k, dissimilarity) appears exactly once
  combos <- paste(ch$k, ch$dissimilarity, sep = "_")
  expect_equal(length(unique(combos)), length(combos))
})

# ==============================================================================
# "all" sentinel
# ==============================================================================

test_that("dissimilarity = 'all' expands to the 9 canonical metrics", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 3, dissimilarity = "all",
                        method = "ward.D2")

  full <- c("hamming", "osa", "lv", "dl", "lcs", "qgram",
            "cosine", "jaccard", "jw")
  expect_equal(nrow(ch), length(full))
  expect_setequal(ch$dissimilarity, full)
})

test_that("method = 'all' expands to the 9 canonical methods", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 3, method = "all",
                        dissimilarity = "hamming")

  full <- c("pam", "ward.D2", "ward.D", "complete", "average",
            "single", "mcquitty", "median", "centroid")
  expect_equal(nrow(ch), length(full))
  expect_setequal(ch$method, full)
})

test_that("'all' mixed with explicit names deduplicates", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 3,
                        dissimilarity = c("hamming", "all"),
                        method = "ward.D2")
  full <- c("hamming", "osa", "lv", "dl", "lcs", "qgram",
            "cosine", "jaccard", "jw")
  expect_equal(nrow(ch), length(full))
  expect_setequal(ch$dissimilarity, full)
})

# ==============================================================================
# Per-row math
# ==============================================================================

test_that("silhouette matches a manual build_clusters() call", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 3, method = "ward.D2",
                        dissimilarity = "hamming", seed = 1)
  ref <- build_clusters(d, k = 3, method = "ward.D2",
                        dissimilarity = "hamming", seed = 1)

  expect_equal(ch$silhouette, ref$silhouette, tolerance = 1e-8)
  # Sizes also align
  expect_equal(c(ch$min_size, ch$max_size),
               c(min(ref$sizes), max(ref$sizes)))
})

test_that("size_ratio is max_size / min_size", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 2:4, method = "ward.D2")
  expect_equal(ch$size_ratio, ch$max_size / ch$min_size,
               tolerance = 1e-9)
})

# ==============================================================================
# weighted validation
# ==============================================================================

test_that("weighted = TRUE with non-hamming sweep is rejected up-front", {
  d <- .make_choice_data()
  expect_error(
    cluster_choice(d, k = 3,
                    dissimilarity = c("hamming", "lcs"),
                    weighted = TRUE),
    "weighted = TRUE requires dissimilarity"
  )
})

# ==============================================================================
# Print
# ==============================================================================

test_that("print: any swept axis triggers the (sweep: ...) header", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 2:4, method = "ward.D2",
                        dissimilarity = "hamming")
  out <- capture.output(print(ch))
  expect_match(out[1L], "^Cluster Choice \\(sweep: k\\)")
})

test_that("print: zero-axis sweep (everything fixed) reads '[<method> / <dissim>]'", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 3, method = "ward.D2",
                        dissimilarity = "hamming")
  out <- capture.output(print(ch))
  expect_match(out[1L],
               "^Cluster Choice \\[ward\\.D2 / hamming\\]")
})

test_that("print: multi-axis header reads '(sweep: ...)' and shows axis names", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 2:3, dissimilarity = c("hamming", "lcs"),
                        method = "ward.D2")
  out <- capture.output(print(ch))
  expect_match(out[1L], "^Cluster Choice \\(sweep:")
  expect_true(any(grepl("k", out[1L])))
  expect_true(any(grepl("dissimilarity", out[1L])))
})

test_that("print: silhouette-best row carries the marker", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 2:4, method = "ward.D2",
                        dissimilarity = "hamming")
  out <- capture.output(print(ch))
  expect_true(any(grepl("<-- best", out)))
})

test_that("print: constant axis columns drop out of the table", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 2:4, method = "ward.D2",
                        dissimilarity = "hamming")
  out <- capture.output(print(ch))
  # data.frame header for the table: should include `k` and metric
  # columns, but not the constant axes.
  data_hdr <- out[grep("^\\s*k\\s+silhouette", out)]
  expect_length(data_hdr, 1L)
  expect_false(grepl("dissimilarity", data_hdr))
  expect_false(grepl("method",        data_hdr))
})

# ==============================================================================
# Summary
# ==============================================================================

test_that("summary returns data.frame with $best on silhouette-max row", {
  d <- .make_choice_data()
  ch <- cluster_choice(d, k = 2:4, method = "ward.D2",
                        dissimilarity = "hamming")
  s <- summary(ch)

  expect_s3_class(s, "data.frame")
  expect_true("best" %in% names(s))
  best_idx <- which(s$best == "silhouette")
  expect_length(best_idx, 1L)
  expect_equal(best_idx, which.max(ch$silhouette))
})

# ==============================================================================
# Plot
# ==============================================================================

test_that("plot returns a ggplot for each sweep shape", {
  d <- .make_choice_data()

  ch1 <- cluster_choice(d, k = 2:4, method = "ward.D2")
  expect_true(inherits(plot(ch1), "ggplot"))

  ch2 <- cluster_choice(d, k = 3,
                        dissimilarity = c("hamming", "lcs", "jaccard"))
  expect_true(inherits(plot(ch2), "ggplot"))

  ch3 <- cluster_choice(d, k = 2:4,
                        dissimilarity = c("hamming", "lcs"))
  expect_true(inherits(plot(ch3), "ggplot"))
})
