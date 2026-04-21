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
