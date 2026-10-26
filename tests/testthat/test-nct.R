testthat::skip_on_cran()
testthat::skip_if_not_installed("Matrix")

# nct(): the network comparison test. Generators/configs live in helper-nct.R;
# fixtures/make-nct.R freezes the complete-data baseline.

nct_fit <- function(i, ...) {
  cfg <- nct_configs[i, ]
  g <- nct_groups(cfg$seed)
  set.seed(cfg$seed)
  nct(g$a, g$b, iter = nct_iter, paired = cfg$paired, abs = cfg$abs,
      weighted = cfg$weighted, ...)
}

test_that("complete-data output is frozen (input handling moves no number)", {
  base <- readRDS(test_path("fixtures", "nct-baseline.rds"))
  invisible(lapply(seq_len(nrow(nct_configs)), function(i)
    expect_identical(nct_numbers(nct_fit(i)), base[[i]],
                     label = sprintf("config %d", i))))
})

test_that("a missing cell is handled pairwise instead of failing", {
  g <- nct_groups(1)
  a_na <- g$a
  a_na[3, 2] <- NA
  set.seed(1)
  clean <- nct(g$a, g$b, iter = 5)
  set.seed(1)
  fit <- nct(a_na, g$b, iter = 5)
  expect_gt(sum(abs(fit$nw1[2, ])), 0)
  expect_lt(max(abs(fit$nw1 - clean$nw1)), 0.05)
  expect_identical(fit$params$na_method, "pairwise")
})

test_that("listwise equals running on the complete rows", {
  g <- nct_groups(2)
  a_na <- g$a
  a_na[c(3, 9), 2] <- NA
  set.seed(4)
  lw <- nct(a_na, g$b, iter = 5, na_method = "listwise")
  set.seed(4)
  cc <- nct(a_na[stats::complete.cases(a_na), ], g$b, iter = 5)
  expect_identical(nct_numbers(lw), nct_numbers(cc))
  expect_identical(lw$n, c(nrow(g$a) - 2L, nrow(g$b)))
})

test_that("a constant variable is an isolated node, with a warning when observed", {
  g <- nct_groups(1)
  g$a[, 5] <- 0
  set.seed(1)
  expect_warning(fit <- nct(g$a, g$b, iter = 5),
                 class = "nestimate_constant_variable")
  expect_identical(unname(fit$nw1[5, ]), rep(0, 5))
  expect_gt(sum(abs(fit$nw2[5, ])), 0)
})

test_that("a variable constant in some permuted splits keeps the test defined", {
  g <- nct_groups(1)
  g$a[, 5] <- 0
  g$a[1:2, 5] <- 1
  g$b[, 5] <- 0
  g$b[1, 5] <- 1
  set.seed(1)
  expect_no_warning(fit <- nct(g$a, g$b, iter = 30))
  expect_true(all(is.finite(c(fit$strength$perm, fit$structure$perm, fit$E$perm))))
})

test_that("a never co-observed pair is refused", {
  g <- nct_groups(1)
  g$a[1:60, 1] <- NA
  g$a[61:120, 2] <- NA
  expect_error(nct(g$a, g$b, iter = 3), class = "nestimate_degenerate_sample")
})

test_that("columns must match in name and order (no position-wise comparison)", {
  g <- nct_groups(1)
  expect_error(nct(g$a, g$b[, c(3, 1, 5, 2, 4)], iter = 3),
               class = "nestimate_bad_input")
  b2 <- g$b
  colnames(b2)[5] <- "W5"
  expect_error(nct(g$a, b2, iter = 3), class = "nestimate_bad_input")
  expect_error(nct(g$a, g$b[, 1:4], iter = 3), class = "nestimate_bad_input")
})

test_that("non-numeric columns are dropped", {
  g <- nct_groups(1)
  da <- data.frame(id = as.character(seq_len(nrow(g$a))), g$a)
  set.seed(1)
  with_id <- nct(da, g$b, iter = 3)
  set.seed(1)
  plain <- nct(g$a, g$b, iter = 3)
  expect_identical(nct_numbers(with_id), nct_numbers(plain))
})

test_that("malformed arguments raise nestimate_bad_input", {
  g <- nct_groups(1)
  invisible(lapply(list(0, NA, 2.7, "10", c(5, 6)), function(it)
    expect_error(nct(g$a, g$b, iter = it), class = "nestimate_bad_input")))
  expect_error(nct(g$a, g$b, paired = NA), class = "nestimate_bad_input")
  expect_error(nct(g$a, g$b, gamma = -1), class = "nestimate_bad_input")
  expect_error(nct(g$a, g$b[-1, ], paired = TRUE), class = "nestimate_bad_input")
})

test_that("summary() gives the tests, as.data.frame() the edges", {
  fit <- nct_fit(1)
  ed <- as.data.frame(fit)
  expect_named(ed, c("from", "to", "weight_1", "weight_2", "difference",
                     "p_value", "significant"))
  expect_identical(nrow(ed), 10L)
  expect_equal(ed$difference, abs(ed$weight_1 - ed$weight_2))
  expect_true(all(ed$p_value >= 1 / (nct_iter + 1) & ed$p_value <= 1))
  gl <- summary(fit)
  expect_named(gl, c("test", "statistic", "observed", "p_value", "significant"))
  expect_identical(gl$statistic, c("M", "S"))
  expect_identical(gl$p_value, c(fit$structure$p_value, fit$strength$p_value))

  g <- nct_groups(1)
  set.seed(1)
  unnamed <- as.data.frame(nct(unname(g$a), unname(g$b), iter = 3))
  expect_identical(unnamed$from[1:3], c("V1", "V1", "V2"))
  expect_identical(unnamed$to[1:3], c("V2", "V3", "V3"))
  expect_false(anyNA(unnamed$to))
})

test_that("print reports the three tests", {
  fit <- nct_fit(1)
  expect_output(print(fit), "Network Comparison Test")
  expect_output(print(fit), "S  global strength")
  expect_output(print(fit), "M  network structure")
  expect_output(print(fit), "Edges, largest difference first")
})

test_that("M and S follow the NCT notation; the old fields fail loudly", {
  fit <- nct_fit(1)
  expect_identical(fit$structure$observed, max(fit$E$observed))
  ut <- upper.tri(fit$nw1)
  expect_equal(fit$strength$by_group,
               c(sum(abs(fit$nw1[ut])), sum(abs(fit$nw2[ut]))))
  expect_error(fit$M, class = "nestimate_renamed")
  expect_error(fit$S, class = "nestimate_renamed")
  ei <- nct_fit(3)                                   # config 3 has abs = FALSE
  expect_identical(summary(ei)$test[2],
                   "global expected influence")
})

# ---- invariance / property tests ----------------------------------------------

test_that("statistics are invariant to rescaling and shifting the variables", {
  g <- nct_groups(2)
  scale_shift <- function(x) sweep(sweep(x, 2, c(2, 0.5, 10, 3, 0.1), "*"), 2,
                                   c(5, -3, 100, 0, 1), "+")
  set.seed(5)
  raw <- nct(g$a, g$b, iter = 8)
  set.seed(5)
  tx <- nct(scale_shift(g$a), scale_shift(g$b), iter = 8)
  expect_equal(tx$nw1, raw$nw1, tolerance = 1e-8)
  expect_equal(tx$strength$perm, raw$strength$perm, tolerance = 1e-8)
  expect_equal(tx$structure$perm, raw$structure$perm, tolerance = 1e-8)
  expect_identical(tx$E$p_value, raw$E$p_value)
})

test_that("identical paired groups give zero statistics and p = 1", {
  g <- nct_groups(1)
  set.seed(2)
  fit <- nct(g$a, g$a, iter = 10, paired = TRUE)
  expect_identical(fit$structure$observed, 0)
  expect_identical(fit$strength$observed, 0)
  expect_true(all(c(fit$structure$p_value, fit$strength$p_value, fit$E$p_value) == 1))
})

test_that("p-values lie on the permutation lattice and p_adjust is p.adjust", {
  raw <- nct_fit(1)
  k <- c(raw$structure$p_value, raw$strength$p_value, raw$E$p_value) * (nct_iter + 1)
  expect_true(all(abs(k - round(k)) < 1e-9))
  adj <- nct_fit(1, p_adjust = "holm")
  expect_equal(adj$E$p_value, stats::p.adjust(raw$E$p_value, method = "holm"))
  expect_identical(adj$strength$p_value, raw$strength$p_value)
})

test_that("weighted = FALSE compares binary networks", {
  fit <- nct_fit(4)
  ut <- upper.tri(fit$nw1)
  expect_true(all(c(fit$nw1, fit$nw2) %in% c(0, 1)))
  expect_equal(fit$strength$by_group, c(sum(fit$nw1[ut]), sum(fit$nw2[ut])))
})

test_that("two empty networks give M = S = 0 and p = 1 (no floating residue)", {
  # With psychnets < 0.6.0 an absent edge could come back as ~1e-17, so two
  # empty networks gave M = 5.6e-17 and a spurious p < 1. Seeds 1, 5, 6 hit it.
  invisible(lapply(c(1, 5, 6), function(seed) {
    set.seed(seed)
    a <- matrix(stats::rnorm(60 * 12), 60, 12)
    b <- matrix(stats::rnorm(60 * 12), 60, 12)
    colnames(a) <- colnames(b) <- paste0("V", 1:12)
    set.seed(1)
    fit <- nct(a, b, iter = 20)
    w <- c(fit$nw1, fit$nw2)
    expect_true(all(w == 0 | abs(w) > 1e-10), label = sprintf("seed %d weights", seed))
    if (all(w == 0)) {
      expect_identical(c(fit$structure$p_value, fit$strength$p_value), c(1, 1))
    }
  }))
})
