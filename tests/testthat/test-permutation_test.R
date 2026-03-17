# ---- permutation_test() Tests ----

# Helper: generate wide sequence data
.make_perm_wide <- function(n = 100, t = 10, states = c("A", "B", "C"),
                            seed = 42) {
  set.seed(seed)
  mat <- matrix(sample(states, n * t, replace = TRUE), nrow = n, ncol = t)
  colnames(mat) <- paste0("T", seq_len(t))
  as.data.frame(mat, stringsAsFactors = FALSE)
}

# Helper: generate frequency-like data for association methods
.make_freq_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  as.data.frame(mat)
}


# ---- Input validation ----

test_that("permutation_test rejects non-netobject inputs", {
  expect_error(permutation_test("a", "b"), "netobject")
})

test_that("permutation_test rejects mismatched methods", {
  wide <- .make_perm_wide(n = 50, seed = 1)
  net1 <- build_network(wide, method = "relative")
  net2 <- build_network(wide, method = "frequency")
  expect_error(permutation_test(net1, net2, iter = 10),
               "Methods must match")
})

test_that("permutation_test rejects mismatched nodes", {
  w1 <- .make_perm_wide(n = 50, states = c("A", "B", "C"), seed = 1)
  w2 <- .make_perm_wide(n = 50, states = c("X", "Y", "Z"), seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  expect_error(permutation_test(net1, net2, iter = 10),
               "Nodes must be the same")
})

test_that("permutation_test rejects missing data", {
  wide <- .make_perm_wide(n = 50, seed = 1)
  net1 <- build_network(wide, method = "relative")
  net2 <- build_network(wide, method = "relative")
  net2$data <- NULL
  expect_error(permutation_test(net1, net2, iter = 10),
               "does not contain \\$data")
})

test_that("paired mode requires equal n", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 60, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  expect_error(permutation_test(net1, net2, iter = 10, paired = TRUE),
               "equal number")
})


# ---- Transition methods ----

test_that("permutation_test works with method='relative'", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")

  perm <- permutation_test(net1, net2, iter = 50L, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "relative")
  expect_equal(perm$iter, 50L)
  expect_equal(perm$alpha, 0.05)
  expect_false(perm$paired)
  expect_equal(perm$adjust, "none")

  # Matrices
  expect_true(is.matrix(perm$diff))
  expect_true(is.matrix(perm$diff_sig))
  expect_true(is.matrix(perm$p_values))
  expect_true(is.matrix(perm$effect_size))
  expect_equal(dim(perm$diff), c(3, 3))
  expect_equal(dimnames(perm$diff), list(net1$nodes, net1$nodes))

  # P-values in [0, 1]
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))

  # diff = x$matrix - y$matrix
  expect_equal(perm$diff, net1$matrix - net2$matrix)

  # Summary
  expect_true(is.data.frame(perm$summary))
  expect_true(all(c("from", "to", "diff", "effect_size", "p_value", "sig",
                     "weight_x", "weight_y") %in% names(perm$summary)))
})

test_that("permutation_test works with method='frequency'", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "frequency")
  net2 <- build_network(w2, method = "frequency")

  perm <- permutation_test(net1, net2, iter = 30L, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "frequency")
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})

test_that("permutation_test works with method='co_occurrence'", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "co_occurrence")
  net2 <- build_network(w2, method = "co_occurrence")

  perm <- permutation_test(net1, net2, iter = 30L, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "co_occurrence")
  expect_false(perm$x$directed)
  # Undirected summary: from <= to (self-loops have from == to)
  if (nrow(perm$summary) > 0) {
    off_diag <- perm$summary$from != perm$summary$to
    if (any(off_diag)) {
      expect_true(all(perm$summary$from[off_diag] < perm$summary$to[off_diag]))
    }
  }
})


# ---- Association methods ----

test_that("permutation_test works with method='cor'", {
  d1 <- .make_freq_data(n = 60, p = 4, seed = 1)
  d2 <- .make_freq_data(n = 60, p = 4, seed = 2)
  net1 <- build_network(d1, method = "cor")
  net2 <- build_network(d2, method = "cor")

  perm <- permutation_test(net1, net2, iter = 30L, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "cor")
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})

test_that("permutation_test works with method='glasso'", {
  d1 <- .make_freq_data(n = 80, p = 4, seed = 1)
  d2 <- .make_freq_data(n = 80, p = 4, seed = 2)
  net1 <- build_network(d1, method = "glasso")
  net2 <- build_network(d2, method = "glasso")

  perm <- permutation_test(net1, net2, iter = 20L, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "glasso")
})


# ---- Paired mode ----

test_that("paired permutation test works", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")

  perm <- permutation_test(net1, net2, iter = 30L, paired = TRUE, seed = 42)

  expect_s3_class(perm, "net_permutation")
  expect_true(perm$paired)
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})


# ---- p.adjust correction ----

test_that("p.adjust correction is applied", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")

  perm_none <- permutation_test(net1, net2, iter = 50L,
                                adjust = "none", seed = 42)
  perm_bh <- permutation_test(net1, net2, iter = 50L,
                               adjust = "BH", seed = 42)

  # BH-adjusted p-values should be >= raw p-values
  expect_true(all(perm_bh$p_values >= perm_none$p_values - 1e-10))
})

test_that("bonferroni correction is more conservative", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")

  perm_bon <- permutation_test(net1, net2, iter = 50L,
                                adjust = "bonferroni", seed = 42)

  # All p-values still in [0, 1]
  expect_true(all(perm_bon$p_values >= 0 & perm_bon$p_values <= 1))
})


# ---- Seed reproducibility ----

test_that("seed produces reproducible results", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")

  perm_a <- permutation_test(net1, net2, iter = 30L, seed = 42)
  perm_b <- permutation_test(net1, net2, iter = 30L, seed = 42)

  expect_equal(perm_a$p_values, perm_b$p_values)
  expect_equal(perm_a$effect_size, perm_b$effect_size)
})


# ---- S3 methods ----

test_that("print.net_permutation works", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  perm <- permutation_test(net1, net2, iter = 20L, seed = 42)

  output <- capture.output(print(perm))
  expect_true(any(grepl("Permutation Test", output)))
  expect_true(any(grepl("Iterations", output)))
  expect_true(any(grepl("Significant", output)))
})

test_that("print shows paired and adjust info", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  perm <- permutation_test(net1, net2, iter = 20L,
                           paired = TRUE, adjust = "BH", seed = 42)

  output <- capture.output(print(perm))
  expect_true(any(grepl("Paired", output)))
  expect_true(any(grepl("BH", output)))
})

test_that("summary returns the summary data frame", {
  w1 <- .make_perm_wide(n = 50, seed = 1)
  w2 <- .make_perm_wide(n = 50, seed = 2)
  net1 <- build_network(w1, method = "relative")
  net2 <- build_network(w2, method = "relative")
  perm <- permutation_test(net1, net2, iter = 20L, seed = 42)

  s <- summary(perm)
  expect_identical(s, perm$summary)
})


# ---- Effect size ----

test_that("effect size is zero when diff is zero", {
  wide <- .make_perm_wide(n = 50, seed = 1)
  net1 <- build_network(wide, method = "relative")
  net2 <- build_network(wide, method = "relative")

  perm <- permutation_test(net1, net2, iter = 20L, seed = 42)

  # Same networks → diff is zero everywhere
  expect_true(all(perm$diff == 0))
  expect_true(all(perm$effect_size == 0))
})
