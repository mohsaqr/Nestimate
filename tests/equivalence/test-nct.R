# nct(): the Network Comparison Test (van Borkulo et al. 2022) of two groups' EBIC graphical lasso networks.
# Nestimate estimates the networks with psychnets (an internal delegation); the reference is
# NetworkComparisonTest::NCT(). The two networks and the observed statistics (global strength difference S,
# largest edge difference M) are deterministic and compared; the permutation p-values are not: the permutations
# are drawn inside each function with no argument to supply them, and the two draw them differently.
# Re-expressed from local_testing_and_equivalence/test-equiv-nct.R, which compared medians of the edge
# differences and p-values within 3/31.

# Two groups with a common factor of strength 0.4 and 0.3.
psy_two_groups <- function(seed, n, p) {
  set.seed(seed)
  draw <- \(strength) {
    values <- matrix(stats::rnorm(n * p), n, p)
    values <- values + strength * matrix(rowSums(values), n, p)
    colnames(values) <- paste0("V", seq_len(p))
    values
  }
  list(a = draw(0.4), b = draw(0.3))
}

inputs <- list(list(seed = 1L, n = 150L, p = 5L), list(seed = 2L, n = 250L, p = 5L), list(seed = 3L, n = 250L, p = 8L))

test_that("the two networks and the observed statistics equal NetworkComparisonTest::NCT()", {
  skip_without_reference("NetworkComparisonTest")
  invisible(lapply(inputs, \(input) {
    groups <- psy_two_groups(input$seed, input$n, input$p)
    found <- nct(groups$a, groups$b, iter = 10L, gamma = 0.5, paired = FALSE, abs = TRUE, weighted = TRUE)
    reference <- NetworkComparisonTest::NCT(data1 = groups$a, data2 = groups$b, it = 10L, gamma = 0.5, paired = FALSE,
                                            abs = TRUE, weighted = TRUE, test.edges = FALSE, progressbar = FALSE)
    # NCT returns the networks without dimnames: they follow the columns of the data.
    named <- \(m) matrix(m, nrow(m), dimnames = list(colnames(groups$a), colnames(groups$a)))
    networks <- psy_pairs(weight_1 = named(reference$nw1), weight_2 = named(reference$nw2))
    # NCT estimates each network with qgraph::EBICglasso(), whose glasso runs to thr = 1e-4, so its partial
    # correlations are about 1e-4 from the optimum Nestimate's solver reaches; 1e-3 is the tolerance the
    # original test held. The statistics are sums and maxima of those weights.
    expect_reference(as.data.frame(found), networks, reference = "NetworkComparisonTest",
                     covers = c("nct()", "as.data.frame.net_nct()"),
                     quantity = sprintf("edge weights, seed %d, n = %d, p = %d", input$seed, input$n, input$p),
                     tolerance = 1e-3, key = c("from", "to"))
    statistics <- data.frame(test = c("network structure", "global strength"),
                             observed = c(reference$nwinv.real, reference$glstrinv.real))
    expect_reference(summary(found), statistics, reference = "NetworkComparisonTest",
                     covers = c("nct()", "summary.net_nct()"),
                     quantity = sprintf("observed statistics, seed %d, n = %d, p = %d", input$seed, input$n, input$p),
                     tolerance = 1e-3, key = "test")
  }))
})
