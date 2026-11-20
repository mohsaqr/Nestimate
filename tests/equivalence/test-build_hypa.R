# build_hypa(): HYPA, hypergeometric path anomalies (LaRock et al. 2020). Given the fitted propensity matrix Xi
# of the order-k De Bruijn graph, the count of a path is a draw of Wallenius' noncentral hypergeometric
# distribution with odds 1 (an urn of round(sum Xi) balls, round(Xi_ij) of them the path's, sum(A) drawn). The
# references: BiasedUrn::pWNCHypergeo() for both tails; the observed path counts counted from the sequences; the
# strength preservation of the fit (expected counts summed per context equal the observed ones, LaRock et al.
# eq. 3); stats::p.adjust() for the adjusted tails. Re-expressed from local_testing_and_equivalence/
# test-equiv-hypa.R (50 random inputs by their largest difference). The adjacency and Xi matrices are read from the
# result's fields `adjacency` and `xi`, the score table from `scores` (API gap: no accessor returns them).
# The fixtures of local_testing_and_equivalence/test-equiv-hypa-oracles.R (hypa_oracles.rds, said to come from the
# original tlarock/hypa code) are not used: the script that produced their JSON sources is not in the repository.

ho_hypa_sequences <- function(seed) {
  set.seed(seed)
  c(lapply(seq_len(30L), \(i) c(sample(c("X", "Y"), 1L), "A", "B", sample(c("X", "Y", "A"), 1L))),
    lapply(seq_len(30L), \(i) rep(c("X", "A", "C", "Y", "A", "D"), 2L)),
    lapply(seq_len(20L), \(i) sample(c("A", "B", "C", "D", "X", "Y"), 12L, replace = TRUE)))
}

# Both tails of every path with Xi > 0, from BiasedUrn.
ho_hypa_tails <- function(adjacency, xi) {
  cells <- which(xi > 0, arr.ind = TRUE)
  nodes <- strsplit(rownames(adjacency), "\x01", fixed = TRUE)
  pool <- round(sum(xi))
  draws <- min(sum(adjacency), pool)
  rows <- lapply(seq_len(nrow(cells)), \(r) {
    i <- cells[r, 1L]
    j <- cells[r, 2L]
    marked <- min(round(xi[i, j]), pool)
    lower <- BiasedUrn::pWNCHypergeo(adjacency[i, j], m1 = marked, m2 = pool - marked, n = draws, odds = 1,
                                     precision = 1e-12)
    data.frame(path = paste(c(nodes[[i]], utils::tail(nodes[[j]], 1L)), collapse = " -> "),
               p_under = lower, p_over = 1 - lower)
  })
  do.call(rbind, rows)
}

# Counts of every path of k + 1 states in the sequences.
ho_path_counts <- function(sequences, k) {
  paths <- unlist(lapply(sequences, \(s) {
    if (length(s) <= k) return(character())
    vapply(seq_len(length(s) - k), \(p) paste(s[p + 0:k], collapse = " -> "), character(1))
  }))
  counts <- table(paths)
  data.frame(path = names(counts), observed = as.vector(counts))
}

test_that("observed path counts equal the counts of paths in the sequences", {
  invisible(lapply(c(2L, 3L), \(k) {
    sequences <- ho_hypa_sequences(700L + k)
    found <- build_hypa(sequences, order = k, min_count = 1L)
    expect_reference(subset(found$scores, observed > 0, select = c(path, observed)), ho_path_counts(sequences, k),
                     reference = "independent: counts of paths of order + 1 states, base R", covers = "build_hypa()",
                     quantity = sprintf("observed count, order %d", k), tolerance = 0, key = "path")
  }))
})

test_that("both tail probabilities equal BiasedUrn::pWNCHypergeo() with odds 1", {
  skip_without_reference("BiasedUrn")
  invisible(lapply(c(2L, 3L), \(k) {
    found <- build_hypa(ho_hypa_sequences(700L + k), order = k, min_count = 1L)
    # BiasedUrn at precision 1e-12; 1e-6 is the original test's tolerance (BiasedUrn's default precision 1e-7
    # was seen 2.7e-6 off).
    expect_reference(subset(found$scores, select = c(path, p_under, p_over)), ho_hypa_tails(found$adjacency, found$xi),
                     reference = "BiasedUrn", covers = "build_hypa()", quantity = sprintf("tail probabilities, order %d", k),
                     tolerance = 1e-6, key = "path")
    expect_reference(subset(summary(found), select = c(path, p_tail)),
                     with(merge(subset(summary(found), select = path), ho_hypa_tails(found$adjacency, found$xi)),
                          data.frame(path = path, p_tail = pmin(p_under, p_over))),
                     reference = "BiasedUrn", covers = c("build_hypa()", "summary.net_hypa()"),
                     quantity = sprintf("anomalous paths, order %d", k), tolerance = 1e-6, measure = "relative", key = "path")
  }))
})

test_that("both tail probabilities equal stats::phyper()", {
  invisible(lapply(c(2L, 3L), \(k) {
    found <- build_hypa(ho_hypa_sequences(700L + k), order = k, min_count = 1L)
    xi <- found$xi
    adjacency <- found$adjacency
    cells <- which(xi > 0, arr.ind = TRUE)
    nodes <- strsplit(rownames(adjacency), "\x01", fixed = TRUE)
    pool <- round(sum(xi))
    marked <- pmin(round(xi[cells]), pool)
    observed <- adjacency[cells]
    draws <- min(sum(adjacency), pool)
    expected <- data.frame(
      path = vapply(seq_len(nrow(cells)), \(r) paste(c(nodes[[cells[r, 1L]]], utils::tail(nodes[[cells[r, 2L]]], 1L)),
                                                     collapse = " -> "), character(1)),
      p_under = stats::phyper(observed, marked, pool - marked, draws),
      p_over = stats::phyper(observed, marked, pool - marked, draws, lower.tail = FALSE))
    expect_reference(subset(found$scores, select = c(path, p_under, p_over)), expected, reference = "stats",
                     covers = "build_hypa()", quantity = sprintf("hypergeometric tails, order %d", k), tolerance = 1e-12,
                     measure = "relative", key = "path", note = "Wallenius with odds 1 is the hypergeometric distribution")
  }))
})

test_that("expected counts preserve the observed strength of every context", {
  sequences <- ho_hypa_sequences(702L)
  found <- build_hypa(sequences, order = 2L, min_count = 1L)
  observed <- ho_path_counts(sequences, 2L)
  observed$from <- sub(" -> [^>]*$", "", observed$path)
  expect_reference(stats::aggregate(expected ~ from, data = found$scores, FUN = sum),
                   stats::aggregate(cbind(expected = observed) ~ from, data = observed, FUN = sum),
                   reference = "independent: out-strength of each context in the counted paths", covers = "build_hypa()",
                   quantity = "expected count summed per context", tolerance = 1e-8, key = "from")
})

test_that("adjusted tails equal stats::p.adjust() of the unadjusted tails", {
  sequences <- ho_hypa_sequences(702L)
  unadjusted <- build_hypa(sequences, order = 2L, min_count = 1L, p_adjust = "none")$scores
  invisible(lapply(c("BH", "bonferroni", "holm"), \(method) {
    found <- build_hypa(sequences, order = 2L, min_count = 1L, p_adjust = method)
    expect_reference(subset(found$scores, select = c(path, p_adjusted_under, p_adjusted_over)),
                     data.frame(path = unadjusted$path, p_adjusted_under = stats::p.adjust(unadjusted$p_under, method),
                                p_adjusted_over = stats::p.adjust(unadjusted$p_over, method)),
                     reference = "stats", covers = "build_hypa()", quantity = sprintf("p_adjust = %s", method),
                     tolerance = 1e-10, key = "path")
  }))
})
