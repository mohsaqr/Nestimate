# simplicial_features() -- tidy topology, one row per network per threshold.

dense_mat <- function() {
  matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3L, 3L,
         dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
}
sparse_mat <- function() {
  matrix(c(0, .2, 0, .2, 0, .1, 0, .1, 0), 3L, 3L,
         dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
}

test_that("a named list yields one tidy row per network", {
  f <- simplicial_features(list(dense = dense_mat(), sparse = sparse_mat()),
                           threshold = 0.3)
  expect_s3_class(f, "data.frame")
  expect_equal(nrow(f), 2L)
  expect_equal(f$network, c("dense", "sparse"))
  expect_true(all(c("b0", "b1", "euler", "max_q", "d1", "higher_order") %in%
                    names(f)))
})

test_that("a threshold vector sweeps, one row per network per threshold", {
  f <- simplicial_features(list(dense = dense_mat()),
                           threshold = c(0.1, 0.3, 0.5))
  expect_equal(nrow(f), 3L)
  expect_equal(f$threshold, c(0.1, 0.3, 0.5))
  # the triangle survives at 0.3 and is gone at 0.5 -- topology is a step
  # function of the threshold, which is why the sweep exists.
  expect_equal(f$d2[f$threshold == 0.3], 1)
  expect_equal(f$d2[f$threshold == 0.5], 0)
})

test_that("higher_order is the total of simplices above dimension 1", {
  f <- simplicial_features(list(m = dense_mat()), threshold = 0.3, max_dim = 4L)
  dcols <- paste0("d", 2:4)
  expect_equal(f$higher_order, sum(unlist(f[, dcols])))
})

test_that("normalize divides simplex counts by node count", {
  raw  <- simplicial_features(list(m = dense_mat()), threshold = 0.3)
  norm <- simplicial_features(list(m = dense_mat()), threshold = 0.3,
                              normalize = TRUE)
  expect_equal(norm$d1, raw$d1 / raw$n_nodes)
  expect_equal(norm$higher_order, raw$higher_order / raw$n_nodes)
})

test_that("an absent Betti number is reported as zero, not NA", {
  f <- simplicial_features(list(m = sparse_mat()), threshold = 0.9)
  expect_false(is.na(f$b1))
  expect_equal(f$b1, 0)
})

test_that("a single network and a netobject_group are both accepted", {
  one <- simplicial_features(dense_mat(), threshold = 0.3)
  expect_equal(nrow(one), 1L)
  expect_equal(one$network, "network_1")

  seqs <- data.frame(t1 = c("A", "B", "A", "C"), t2 = c("B", "C", "C", "A"),
                     t3 = c("C", "A", "B", "B"),
                     g  = c("x", "x", "y", "y"), stringsAsFactors = FALSE)
  grp <- build_network(seqs, method = "frequency", group = "g")
  f <- simplicial_features(grp, threshold = 0)
  expect_equal(nrow(f), length(grp))
})

test_that("input and arguments are validated", {
  expect_error(simplicial_features(list(m = dense_mat()), threshold = NA),
               "`threshold` must be")
  expect_error(simplicial_features(list(m = dense_mat()), normalize = "yes"),
               "`normalize` must be")
  expect_error(simplicial_features("not a network"), "Cannot extract networks")
})
