# ---- Equivalence tests: cooccurrence() vs citenets::conetwork() ----
#
# Validates that Nestimate's cooccurrence() produces identical co-occurrence
# weights to citenets::conetwork() across all shared similarity measures
# and multiple random configurations.
#
# citenets differences accounted for:
#   - citenets uppercases all entity names (toupper)
#   - citenets always zeros the diagonal
#   - citenets returns an edge list; we convert to matrix for comparison
#   - citenets requires data.frame with `id` column + list-column
#   - citenets counting = "full" (binary) matches our default
#   - Shared similarities: none, jaccard, cosine, inclusion, association,
#     equivalence. (dice and relative are Nestimate-only)

skip_if_not_installed("citenets")

TOL <- 1e-12

# ---- Helper: convert citenets edge list to symmetric matrix ----
.citenets_to_matrix <- function(edges, node_names) {
  k <- length(node_names)
  M <- matrix(0, nrow = k, ncol = k, dimnames = list(node_names, node_names))
  if (nrow(edges) == 0L) return(M)
  for (r in seq_len(nrow(edges))) {
    i <- edges$from[r]
    j <- edges$to[r]
    M[i, j] <- edges$weight[r]
    M[j, i] <- edges$weight[r]
  }
  M
}

# ---- Helper: make citenets-compatible data from transactions ----
.make_citenets_data <- function(transactions) {
  d <- data.frame(id = seq_along(transactions), stringsAsFactors = FALSE)
  d$field <- transactions
  d
}

# ---- Helper: run one comparison ----
.compare_co_occurrence <- function(transactions, similarity, label) {
  # Nestimate: use uppercase to match citenets
  net <- cooccurrence(transactions, similarity = similarity, diagonal = FALSE)
  W_nest <- net$weights

  # citenets
  d <- .make_citenets_data(transactions)
  cn <- citenets::conetwork(d, "field", sep = NULL,
                            similarity = similarity, threshold = 0)

  # Get node names from Nestimate (already sorted, uppercase)
  nodes <- sort(colnames(W_nest))

  # Convert citenets edge list to matrix
  W_cite <- .citenets_to_matrix(cn, nodes)

  # Compare off-diagonal (both have zero diagonal)
  delta <- abs(W_nest - W_cite)
  max_delta <- max(delta)

  expect_true(
    max_delta < TOL,
    info = sprintf("[%s] similarity=%s: max delta = %.2e (should be < %.2e)",
                   label, similarity, max_delta, TOL)
  )

  max_delta
}

# ---- Generate random configs ----
set.seed(42)
n_configs <- 30L
all_items <- LETTERS[1:12]

configs <- lapply(seq_len(n_configs), function(i) {
  n_trans <- sample(5:20, 1)
  n_items_per <- sample(2:6, n_trans, replace = TRUE)
  lapply(seq_len(n_trans), function(j) {
    sample(all_items, n_items_per[j])
  })
})

# ---- Shared similarity methods ----
shared_sims <- c("none", "jaccard", "cosine", "inclusion",
                 "association", "equivalence")

# ========================================
# Test: all shared similarities across 30 configs
# ========================================

test_that("co_occurrence matches citenets::conetwork for similarity='none'", {
  skip_if_not_installed("citenets")
  deltas <- vapply(seq_along(configs), function(i) {
    .compare_co_occurrence(configs[[i]], "none", sprintf("config_%d", i))
  }, numeric(1))
  cat(sprintf("\n  [none] max delta across %d configs: %.2e\n",
              length(configs), max(deltas)))
})

test_that("co_occurrence matches citenets::conetwork for similarity='jaccard'", {
  skip_if_not_installed("citenets")
  deltas <- vapply(seq_along(configs), function(i) {
    .compare_co_occurrence(configs[[i]], "jaccard", sprintf("config_%d", i))
  }, numeric(1))
  cat(sprintf("\n  [jaccard] max delta across %d configs: %.2e\n",
              length(configs), max(deltas)))
})

test_that("co_occurrence matches citenets::conetwork for similarity='cosine'", {
  skip_if_not_installed("citenets")
  deltas <- vapply(seq_along(configs), function(i) {
    .compare_co_occurrence(configs[[i]], "cosine", sprintf("config_%d", i))
  }, numeric(1))
  cat(sprintf("\n  [cosine] max delta across %d configs: %.2e\n",
              length(configs), max(deltas)))
})

test_that("co_occurrence matches citenets::conetwork for similarity='inclusion'", {
  skip_if_not_installed("citenets")
  deltas <- vapply(seq_along(configs), function(i) {
    .compare_co_occurrence(configs[[i]], "inclusion", sprintf("config_%d", i))
  }, numeric(1))
  cat(sprintf("\n  [inclusion] max delta across %d configs: %.2e\n",
              length(configs), max(deltas)))
})

test_that("co_occurrence matches citenets::conetwork for similarity='association'", {
  skip_if_not_installed("citenets")
  deltas <- vapply(seq_along(configs), function(i) {
    .compare_co_occurrence(configs[[i]], "association", sprintf("config_%d", i))
  }, numeric(1))
  cat(sprintf("\n  [association] max delta across %d configs: %.2e\n",
              length(configs), max(deltas)))
})

test_that("co_occurrence matches citenets::conetwork for similarity='equivalence'", {
  skip_if_not_installed("citenets")
  deltas <- vapply(seq_along(configs), function(i) {
    .compare_co_occurrence(configs[[i]], "equivalence", sprintf("config_%d", i))
  }, numeric(1))
  cat(sprintf("\n  [equivalence] max delta across %d configs: %.2e\n",
              length(configs), max(deltas)))
})

# ========================================
# Test: min_occur equivalence
# ========================================

test_that("min_occur matches citenets min_occur filtering", {
  skip_if_not_installed("citenets")
  set.seed(99)

  deltas <- vapply(seq_len(10), function(i) {
    n_trans <- sample(8:15, 1)
    trans <- lapply(seq_len(n_trans), function(j) {
      sample(LETTERS[1:8], sample(2:5, 1))
    })
    min_occ <- 3L

    net <- cooccurrence(trans, similarity = "cosine",
                         min_occur = min_occ, diagonal = FALSE)

    d <- .make_citenets_data(trans)
    cn <- citenets::conetwork(d, "field", sep = NULL,
                              similarity = "cosine",
                              min_occur = min_occ, threshold = 0)

    nodes <- sort(colnames(net$weights))
    W_cite <- .citenets_to_matrix(cn, nodes)

    max(abs(net$weights - W_cite))
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("min_occur: max delta = %.2e", max(deltas)))
  cat(sprintf("\n  [min_occur] max delta across 10 configs: %.2e\n",
              max(deltas)))
})

# ========================================
# Test: threshold equivalence
# ========================================

test_that("threshold matches citenets threshold filtering", {
  skip_if_not_installed("citenets")
  set.seed(77)

  deltas <- vapply(seq_len(10), function(i) {
    n_trans <- sample(8:15, 1)
    trans <- lapply(seq_len(n_trans), function(j) {
      sample(LETTERS[1:8], sample(2:5, 1))
    })
    thresh <- 0.3

    net <- cooccurrence(trans, similarity = "jaccard",
                         threshold = thresh, diagonal = FALSE)

    d <- .make_citenets_data(trans)
    cn <- citenets::conetwork(d, "field", sep = NULL,
                              similarity = "jaccard",
                              threshold = thresh)

    nodes <- sort(colnames(net$weights))
    # citenets only has edges above threshold; missing pairs = 0
    W_cite <- .citenets_to_matrix(cn, nodes)

    max(abs(net$weights - W_cite))
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("threshold: max delta = %.2e", max(deltas)))
  cat(sprintf("\n  [threshold] max delta across 10 configs: %.2e\n",
              max(deltas)))
})

# ========================================
# Test: delimited input equivalence
# ========================================

test_that("delimited input matches citenets with sep splitting", {
  skip_if_not_installed("citenets")
  set.seed(55)

  deltas <- vapply(seq_len(10), function(i) {
    n_trans <- sample(5:12, 1)
    trans <- lapply(seq_len(n_trans), function(j) {
      sample(LETTERS[1:8], sample(2:5, 1))
    })

    # Nestimate: delimited data.frame
    df_nest <- data.frame(
      id = seq_along(trans),
      items = vapply(trans, paste, character(1), collapse = "; "),
      stringsAsFactors = FALSE
    )
    net <- cooccurrence(df_nest, field = "items", sep = ";",
                         similarity = "cosine", diagonal = FALSE)

    # citenets: list-column
    d <- .make_citenets_data(trans)
    cn <- citenets::conetwork(d, "field", sep = NULL,
                              similarity = "cosine", threshold = 0)

    nodes <- sort(colnames(net$weights))
    W_cite <- .citenets_to_matrix(cn, nodes)

    max(abs(net$weights - W_cite))
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("delimited: max delta = %.2e", max(deltas)))
  cat(sprintf("\n  [delimited] max delta across 10 configs: %.2e\n",
              max(deltas)))
})

# ========================================
# Test: long/bipartite input equivalence
# ========================================

test_that("long/bipartite input matches citenets", {
  skip_if_not_installed("citenets")
  set.seed(33)

  deltas <- vapply(seq_len(10), function(i) {
    n_trans <- sample(5:12, 1)
    trans <- lapply(seq_len(n_trans), function(j) {
      sample(LETTERS[1:8], sample(2:5, 1))
    })

    # Nestimate: long format
    df_long <- data.frame(
      doc = rep(seq_along(trans), lengths(trans)),
      item = unlist(trans),
      stringsAsFactors = FALSE
    )
    net <- cooccurrence(df_long, field = "item", by = "doc",
                         similarity = "jaccard", diagonal = FALSE)

    # citenets: list-column
    d <- .make_citenets_data(trans)
    cn <- citenets::conetwork(d, "field", sep = NULL,
                              similarity = "jaccard", threshold = 0)

    nodes <- sort(colnames(net$weights))
    W_cite <- .citenets_to_matrix(cn, nodes)

    max(abs(net$weights - W_cite))
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("long: max delta = %.2e", max(deltas)))
  cat(sprintf("\n  [long/bipartite] max delta across 10 configs: %.2e\n",
              max(deltas)))
})

# ========================================
# Test: binary matrix input equivalence
# ========================================

test_that("binary matrix input matches citenets", {
  skip_if_not_installed("citenets")
  set.seed(22)

  deltas <- vapply(seq_len(10), function(i) {
    n_trans <- sample(5:12, 1)
    items <- LETTERS[1:sample(4:8, 1)]
    trans <- lapply(seq_len(n_trans), function(j) {
      sample(items, sample(2:length(items), 1))
    })

    # Nestimate: binary matrix
    all_items <- sort(unique(unlist(trans)))
    bin <- matrix(0L, nrow = n_trans, ncol = length(all_items),
                  dimnames = list(NULL, all_items))
    for (j in seq_len(n_trans)) bin[j, trans[[j]]] <- 1L
    net <- cooccurrence(bin, similarity = "inclusion", diagonal = FALSE)

    # citenets: list-column
    d <- .make_citenets_data(trans)
    cn <- citenets::conetwork(d, "field", sep = NULL,
                              similarity = "inclusion", threshold = 0)

    nodes <- sort(colnames(net$weights))
    W_cite <- .citenets_to_matrix(cn, nodes)

    max(abs(net$weights - W_cite))
  }, numeric(1))

  expect_true(all(deltas < TOL),
              info = sprintf("binary: max delta = %.2e", max(deltas)))
  cat(sprintf("\n  [binary matrix] max delta across 10 configs: %.2e\n",
              max(deltas)))
})
