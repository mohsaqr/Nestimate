test_that("distribution_plot runs on a matrix with default args", {
  set.seed(1L)
  states <- c("A", "B", "C")
  seqs <- matrix(sample(states, 30 * 8, replace = TRUE), 30, 8)

  pdf(NULL); on.exit(dev.off(), add = TRUE)
  res <- distribution_plot(seqs, na = FALSE)
  expect_type(res, "list")
  expect_named(res, c("counts", "proportions", "levels", "palette", "groups"))
  expect_equal(res$levels, sort(states))
  # Column sums of proportions == 1 in every group.
  expect_true(all(abs(colSums(res$proportions[[1]]) - 1) < 1e-10))
})

test_that("distribution_plot 'count' scale preserves raw counts", {
  seqs <- matrix(c("A", "A", "B", "B", "A", "C"), nrow = 3)

  pdf(NULL); on.exit(dev.off(), add = TRUE)
  res <- distribution_plot(seqs, scale = "count", na = FALSE)
  # Column 1: all "A" except row 2 is "A" too (wait it's a 3x2 matrix)
  # Actually: matrix fills by column; seqs[,1] = c("A","A","B"), seqs[,2] = c("B","A","C")
  expect_equal(res$counts[[1]][, 1],
               c(A = 2, B = 1, C = 0), ignore_attr = TRUE)
  expect_equal(res$counts[[1]][, 2],
               c(A = 1, B = 1, C = 1), ignore_attr = TRUE)
})

test_that("distribution_plot accepts net_clustering and draws per-group panels", {
  set.seed(2L)
  seqs <- as.data.frame(matrix(sample(c("A", "B", "C", "D"), 40 * 10,
                                      replace = TRUE), 40, 10))
  cl <- cluster_data(seqs, k = 3L, dissimilarity = "hamming",
                     method = "ward.D2")

  pdf(NULL); on.exit(dev.off(), add = TRUE)
  res <- distribution_plot(cl)
  expect_length(res$groups, 3L)
  expect_length(res$counts, 3L)
  # Each per-group count matrix sums to the cluster size × n_time.
  group_sizes <- table(cl$assignments)
  for (g in seq_along(res$counts)) {
    expect_equal(sum(res$counts[[g]]), as.integer(group_sizes[g]) * ncol(seqs))
  }
})

test_that("distribution_plot include_na adds an NA band", {
  set.seed(3L)
  seqs <- matrix(sample(c("A", "B", "C"), 20 * 6, replace = TRUE), 20, 6)
  seqs[sample(length(seqs), 10)] <- NA

  pdf(NULL); on.exit(dev.off(), add = TRUE)
  res <- distribution_plot(seqs, na = TRUE)
  expect_identical(tail(res$levels, 1L), "NA")
  expect_length(res$palette, 4L)
})

test_that("distribution_plot supports geom = 'bar'", {
  set.seed(4L)
  seqs <- matrix(sample(c("A", "B"), 15 * 5, replace = TRUE), 15, 5)

  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_silent(distribution_plot(seqs, geom = "bar"))
})

test_that("distribution_plot honours legend_position values", {
  set.seed(5L)
  seqs <- matrix(sample(c("A", "B", "C"), 12 * 5, replace = TRUE), 12, 5)

  pdf(NULL); on.exit(dev.off(), add = TRUE)
  for (pos in c("right", "bottom", "none")) {
    expect_silent(distribution_plot(seqs, legend = pos))
  }
})

test_that("distribution_plot rejects group vectors of wrong length", {
  seqs <- matrix(sample(c("A", "B"), 20, replace = TRUE), 10, 2)
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(distribution_plot(seqs, group = c(1, 2, 3)), "length")
})

test_that("distribution_plot rejects all-NA input", {
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_error(distribution_plot(matrix(NA_character_, 4, 3)),
               "no non-NA values")
})
