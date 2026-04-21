# Numerical equivalence: Nestimate::nct() vs NetworkComparisonTest::NCT()
#
# Both processes use the same RNG seed before their permutation loops, so
# the sampled label assignments and resulting networks are byte-identical.
# Implementation matches NCT defaults (abs = TRUE, weighted = TRUE,
# paired = FALSE, NCT_estimator_GGM with nearPD + EBICglasso).

set.seed(20260413)
TOL <- 1e-12

skip_if_not_installed("NetworkComparisonTest")
skip_if_not_installed("qgraph")
skip_if_not_installed("Matrix")

.gen_two_groups <- function(seed, n1, n2, p) {
  set.seed(seed)
  A <- matrix(stats::rnorm(n1 * p), n1, p)
  A <- A + 0.4 * matrix(rowSums(A), n1, p)
  colnames(A) <- paste0("V", seq_len(p))
  B <- matrix(stats::rnorm(n2 * p), n2, p)
  B <- B + 0.3 * matrix(rowSums(B), n2, p)
  colnames(B) <- paste0("V", seq_len(p))
  list(A = A, B = B)
}

cfgs <- expand.grid(seed = 1:5, n = c(150, 250), p = c(5, 8))

test_that("nct() matches NetworkComparisonTest::NCT() at machine precision", {
  invisible(lapply(seq_len(nrow(cfgs)), function(i) {
    cfg <- cfgs[i, ]
    g <- .gen_two_groups(cfg$seed, cfg$n, cfg$n, cfg$p)

    set.seed(cfg$seed * 100 + 1)
    ours <- nct(g$A, g$B, iter = 30L, gamma = 0.5,
                paired = FALSE, abs = TRUE, weighted = TRUE)

    set.seed(cfg$seed * 100 + 1)
    ref <- suppressWarnings(suppressMessages(
      NetworkComparisonTest::NCT(
        data1 = g$A, data2 = g$B, it = 30L, gamma = 0.5,
        paired = FALSE, abs = TRUE, weighted = TRUE,
        test.edges = TRUE, edges = "all", progressbar = FALSE
      )
    ))

    lbl <- sprintf("cfg(seed=%d,n=%d,p=%d)", cfg$seed, cfg$n, cfg$p)
    expect_equal(unname(ours$nw1), unname(ref$nw1), tolerance = TOL,
                 label = paste("nw1", lbl))
    expect_equal(unname(ours$nw2), unname(ref$nw2), tolerance = TOL,
                 label = paste("nw2", lbl))
    expect_equal(ours$M$observed, ref$glstrinv.real, tolerance = TOL,
                 label = paste("M_obs", lbl))
    expect_equal(unname(ours$M$perm), unname(ref$glstrinv.perm), tolerance = TOL,
                 label = paste("M_perm", lbl))
    expect_equal(ours$M$p_value, ref$glstrinv.pval, tolerance = TOL,
                 label = paste("M_pval", lbl))
    expect_equal(ours$S$observed, ref$nwinv.real, tolerance = TOL,
                 label = paste("S_obs", lbl))
    expect_equal(unname(ours$S$perm), unname(ref$nwinv.perm), tolerance = TOL,
                 label = paste("S_perm", lbl))
    expect_equal(ours$S$p_value, ref$nwinv.pval, tolerance = TOL,
                 label = paste("S_pval", lbl))
    expect_equal(unname(ours$E$p_value),
                 unname(ref$einv.pvals[, 3]),
                 tolerance = TOL,
                 label = paste("E_pval", lbl))
  }))
})

test_that("nct() print method works", {
  set.seed(1)
  g <- .gen_two_groups(1, 100, 100, 4)
  res <- nct(g$A, g$B, iter = 20L, gamma = 0.5)
  expect_s3_class(res, "net_nct")
  expect_output(print(res), "Network Comparison Test")
  expect_output(print(res), "Global strength")
  expect_output(print(res), "Network structure")
})
