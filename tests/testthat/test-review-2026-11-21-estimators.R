testthat::skip_on_cran()

# Regression tests for the 2026-11-21 verification findings in the estimators,
# centrality and network-comparison code. Every expected value is computed by
# hand from the definition (see the comment above each expectation).

# ---- n-gram with an end state: the end column is added only when needed ----

ngram_wide <- data.frame(T1 = c("A", "B", "A"), T2 = c("B", "A", "C"),
                         T3 = c("C", NA, "A"), T4 = NA_character_)

test_that("ngram with end = TRUE puts the end state in the empty last column", {
  # With End: A B C End | B A End . | A C A End -> 4 columns, so n_gram = 3
  # gives windows 1-3 and 2-4: pair (1,2) once, (2,3) twice, (3,4) once.
  # Row 1: A->B 1, B->C 2, C->End 1. Row 2: B->A 1, A->End 2.
  # Row 3: A->C 1, C->A 2, A->End 1.
  fit <- build_network(ngram_wide, method = "ngram", end = TRUE,
                       params = list(n_gram = 3))
  states <- c("A", "B", "C", "End")
  expected <- matrix(0, 4, 4, dimnames = list(states, states))
  expected["A", c("B", "C", "End")] <- c(1, 1, 3)
  expected["B", c("A", "C")] <- c(1, 2)
  expected["C", c("A", "End")] <- c(2, 1)
  expect_equal(unclass(extract_transition_matrix(fit)), expected,
               ignore_attr = TRUE)
  # A sequence filling the last column needs a new one: without T4 the
  # first and third rows fill T3, the end state goes to a new fourth
  # column, and the windows (hence the counts) are the same.
  trimmed <- build_network(ngram_wide[c("T1", "T2", "T3")], method = "ngram",
                           end = TRUE, params = list(n_gram = 3))
  expect_equal(extract_transition_matrix(trimmed),
               extract_transition_matrix(fit))
})

test_that("an end state leaves the other methods' counts as they were", {
  # frequency with End: A->B, B->C, C->End | B->A, A->End | A->C, C->A, A->End
  fit <- build_network(ngram_wide, method = "frequency", end = TRUE)
  states <- c("A", "B", "C", "End")
  expected <- matrix(0, 4, 4, dimnames = list(states, states))
  expected["A", c("B", "C", "End")] <- c(1, 1, 2)
  expected["B", c("A", "C")] <- c(1, 1)
  expected["C", c("A", "End")] <- c(1, 1)
  expect_equal(unclass(extract_transition_matrix(fit)), expected,
               ignore_attr = TRUE)
  # n_gram = 2 counts each adjacent pair once: the frequency counts.
  bigram <- build_network(ngram_wide, method = "ngram", end = TRUE,
                          params = list(n_gram = 2))
  expect_equal(unclass(extract_transition_matrix(bigram)), expected,
               ignore_attr = TRUE)
})

test_that("the helper appends the end column only on request when no row fills the last one", {
  always <- .prepare_sequence_wide(ngram_wide, end_state = "End")
  needed <- .prepare_sequence_wide(ngram_wide, end_state = "End",
                                   end_column = "needed")
  expect_identical(ncol(always$matrix), 5L)
  expect_identical(ncol(needed$matrix), 4L)
  expect_identical(unname(needed$matrix[, 4L]), c("End", NA, "End"))
  # the appended column is the only difference
  expect_identical(needed$matrix, always$matrix[, seq_len(4L)])
})

# ---- Closeness ignoring direction with raw weights as distances ----

closeness_wide <- data.frame(T1 = c("A", "B", "C", "A"),
                             T2 = c("B", "C", "A", "C"),
                             T3 = c("A", "A", "B", "B"))

test_that("Closeness with invert = FALSE joins two nodes by their shorter edge", {
  # Counts: A->B 2, B->A 1, A->C 1, C->A 2, B->C 1, C->B 1. As distances the
  # shorter edge of each pair is 1, so every node is at 1 + 1 = 2: 1/2.
  net <- build_network(closeness_wide, method = "frequency")
  found <- suppressMessages(
    net_centrality(net, measures = "Closeness", invert = FALSE))
  expect_equal(found$Closeness, c(0.5, 0.5, 0.5))
})

test_that("Closeness with invert = TRUE is unchanged: the larger weight is the shorter edge", {
  # Distances 1/w: A-B 1/2, A-C 1/2, B-C 1. A: 1/2 + 1/2 -> 1; B: 1/2 + 1
  # (direct, equal to B-A-C) -> 2/3; C likewise 2/3.
  net <- build_network(closeness_wide, method = "frequency")
  found <- suppressMessages(
    net_centrality(net, measures = "Closeness", invert = TRUE))
  expect_equal(found$Closeness, c(1, 2 / 3, 2 / 3))
})

test_that("undirected closeness keeps a one-way edge as the only edge", {
  # A->B 3 only, B->C 1 and C->B 4 (raw distances): A-B 3, B-C min(1, 4) = 1.
  # A: 3 + 4 -> 1/7; B: 3 + 1 -> 1/4; C: 4 + 1 -> 1/5.
  W <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  W["A", "B"] <- 3
  W["B", "C"] <- 1
  W["C", "B"] <- 4
  cl <- .closeness(W, directed = TRUE, invert = FALSE)
  expect_equal(unname(cl$Closeness), c(1 / 7, 1 / 4, 1 / 5))
  # On a symmetric matrix both rules read the same edge.
  S <- W + t(W)
  expect_equal(.closeness(S, invert = FALSE)$Closeness,
               .closeness(S, invert = FALSE)$ClosenessOut)
})

# ---- min-max scaling of transition networks (the documented rule) ----

test_that("min-max on a transition network spans every cell, zeros included", {
  # Counts A->A 1, A->B 3, B->A 3, B->B 2: no zero cell, so the range is
  # [1, 3] and the smallest count maps to 0: (w - 1) / 2.
  wide <- data.frame(T1 = c("A", "A", "A", "B", "B"),
                     T2 = c("A", "B", "B", "A", "B"),
                     T3 = c(NA, "B", "A", "B", "A"))
  counts <- extract_transition_matrix(build_network(wide, method = "frequency"))
  scaled <- extract_transition_matrix(build_network(wide, method = "frequency",
                                                    scaling = "minmax"))
  expect_equal(unclass(counts), matrix(c(1, 3, 3, 2), 2), ignore_attr = TRUE)
  expect_equal(unclass(scaled), matrix(c(0, 1, 1, 0.5), 2),
               ignore_attr = TRUE)
})

# ---- compare_networks(test = "permutation") on association networks ----

cn_numeric <- function(n, seed) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * 3L), n, 3L)
  x[, 2L] <- x[, 2L] + x[, 1L]
  stats::setNames(as.data.frame(x), c("V1", "V2", "V3"))
}

test_that("compare_networks() permutes cor and pcor networks as permutation() does", {
  invisible(lapply(c("cor", "pcor"), \(method) {
    a <- build_network(cn_numeric(40L, 1L), method = method)
    b <- build_network(cn_numeric(30L, 2L), method = method)
    expect_true(is.matrix(a$data), info = method)
    cmp <- compare_networks(a, b, test = "permutation", iter = 30L, seed = 7)
    direct <- permutation(a, b, iter = 30L, seed = 7)
    run <- cmp$permutations[[1L]]
    expect_equal(run$p_values, direct$p_values, info = method)
    expect_equal(run$effect_size, direct$effect_size, info = method)
    edges <- edge_differences(cmp)
    # the edge table carries the same p-value per cell (row-major pairs)
    expect_equal(edges$perm_p[edges$from == "V1" & edges$to == "V2"],
                 direct$p_values["V1", "V2"], info = method)
  }))
})

test_that("a network without data is still refused by the permutation backend", {
  a <- build_network(cn_numeric(40L, 1L), method = "cor")
  b <- build_network(cn_numeric(30L, 2L), method = "cor")
  a$data <- NULL
  expect_error(compare_networks(a, b, test = "permutation", iter = 10L),
               class = "nestimate_compare_test_unsupported")
})

# ---- One-state networks: undefined similarities are NA, not errors ----

one_state <- build_network(data.frame(T1 = rep("A", 10L), T2 = "A", T3 = "A",
                                      T4 = "A", T5 = "A", T6 = "A"),
                           method = "relative")

test_that("the similarity of two one-node networks is NA where it is undefined", {
  # One weight (P(A -> A) = 1) in each: no spread, so correlations, the CV
  # ratio, RV and rank agreement are undefined; distances are 0 and the
  # overlap measures 1.
  sim <- .network_similarity(one_state$weights, one_state$weights)
  undefined <- c("cv_ratio", "pearson", "spearman", "kendall", "distance_cor",
                 "rv", "rank_agreement")
  expect_true(all(is.na(sim[undefined])))
  expect_false(any(is.nan(sim)))
  expect_equal(unname(sim[c("mean_abs_diff", "max_abs_diff", "euclidean")]),
               c(0, 0, 0))
  expect_equal(unname(sim[c("cosine", "jaccard", "sign_agreement")]),
               c(1, 1, 1))
})

test_that("compare_networks() and network_reliability() answer a one-state network", {
  cmp <- compare_networks(one_state, one_state)
  expect_true(inherits(cmp, "net_network_comparison"))
  global <- global_differences(cmp)
  expect_true(is.na(global$value[global$key == "pearson"]))
  expect_equal(global$value[global$key == "mean_abs_diff"], 0)
  perm <- compare_networks(one_state, one_state, test = "permutation",
                           iter = 5L, seed = 1)
  expect_true(inherits(perm, "net_network_comparison"))
  rel <- network_reliability(one_state, iter = 5L, seed = 1)
  expect_true(inherits(rel, "net_reliability"))
})

test_that("similarities of ordinary networks are unchanged by the NA guards", {
  # x = (0, 1, 2, 3), y = (0, 2, 2, 4) column-major 2 x 2: |x - y| = 0,1,0,1.
  x <- matrix(c(0, 1, 2, 3), 2)
  y <- matrix(c(0, 2, 2, 4), 2)
  sim <- .network_similarity(x, y)
  expect_equal(unname(sim["mean_abs_diff"]), 0.5)
  expect_equal(unname(sim["pearson"]), stats::cor(c(0, 1, 2, 3), c(0, 2, 2, 4)))
  # cv_ratio = sd_x * mean_y / (mean_x * sd_y)
  expect_equal(unname(sim["cv_ratio"]),
               stats::sd(c(0, 1, 2, 3)) * 2 / (1.5 * stats::sd(c(0, 2, 2, 4))))
  # rank agreement: the one row difference (1, 1) against (2, 2): signs agree
  expect_equal(unname(sim["rank_agreement"]), 1)
  # a constant y has no spread: correlations NA, as before
  flat <- .network_similarity(x, matrix(1, 2, 2))
  expect_true(is.na(flat["pearson"]))
  expect_true(is.na(flat["cv_ratio"]))
})
