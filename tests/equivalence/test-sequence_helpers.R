# actor_endpoints(), net_aggregate_weights() and cooccurrence(): small numeric helpers, against their
# definitions in base R.

test_that("actor_endpoints() gives each sequence's first and last observed step and state", {
  wide <- tn_simulate(n_actors = 15L, n_states = 4L, seq_length = 8L, seed = 601L, ragged = TRUE)
  sequences <- tn_sequences(wide)
  observed <- lengths(sequences)
  expected <- data.frame(actor = seq_along(sequences),
                         first_state = vapply(sequences, \(s) s[1L], character(1L)),
                         last_state = vapply(sequences, \(s) s[length(s)], character(1L)),
                         first_step = 1L, last_step = observed, n_observed = observed,
                         dropped_out = observed < ncol(wide))
  expect_reference(actor_endpoints(wide), expected, reference = "independent: first and last non-missing cell of each row, base R",
                   covers = "actor_endpoints()", quantity = "endpoints of 15 ragged sequences", tolerance = 0, key = "actor")
})

test_that("net_aggregate_weights() equals the base R summary of the non-zero weights", {
  weights <- c(0.4, 0, 1.5, NA, 2.25, 0.1, 3)
  kept <- c(0.4, 1.5, 2.25, 0.1, 3)
  methods <- c("sum", "mean", "median", "max", "min", "prod", "geomean")
  expected <- c(sum(kept), mean(kept), stats::median(kept), max(kept), min(kept), prod(kept), exp(mean(log(kept))))
  found <- vapply(methods, \(method) net_aggregate_weights(weights, method = method), numeric(1L))
  expect_reference(found, stats::setNames(expected, methods),
                   reference = "independent: base R summaries of the non-zero, non-missing weights",
                   covers = "net_aggregate_weights()", quantity = "aggregate", tolerance = 1e-12)
  expect_reference(net_aggregate_weights(weights, method = "density", n_possible = 10), sum(kept) / 10,
                   reference = "independent: base R summaries of the non-zero, non-missing weights",
                   covers = "net_aggregate_weights()", quantity = "density over 10 possible edges", tolerance = 1e-12)
})

test_that("cooccurrence() counts and Jaccard similarities equal the incidence matrix computations", {
  set.seed(602L)
  records <- data.frame(id = rep(seq_len(20L), each = 3L), item = sample(letters[1:6], 60L, replace = TRUE))
  records <- unique(records)
  incidence <- unclass(table(records$id, records$item)) > 0
  counts <- crossprod(incidence * 1)
  dimnames(counts) <- list(colnames(incidence), colnames(incidence))
  expect_reference(extract_transition_matrix(cooccurrence(records, field = "item", by = "id")), counts,
                   reference = "independent: crossprod() of the record-by-item incidence matrix, base R",
                   covers = "cooccurrence()", quantity = "co-occurrence counts", tolerance = 0)
  occurrences <- diag(counts)
  jaccard <- counts / (outer(occurrences, occurrences, "+") - counts)
  expect_reference(extract_transition_matrix(cooccurrence(records, field = "item", by = "id", similarity = "jaccard")),
                   jaccard, reference = "independent: crossprod() of the record-by-item incidence matrix, base R",
                   covers = "cooccurrence()", quantity = "Jaccard similarity", tolerance = 1e-12)
})

# Outside references: TraMineR for sequence lengths and states, arules for co-occurrence, nnet and psych for
# indicators and the geometric mean.
test_that("actor_endpoints() lengths and end states equal TraMineR::seqlength() and the sequence object", {
  skip_without_reference("TraMineR")
  invisible(lapply(c(601L, 603L), \(seed) {
    wide <- tn_simulate(n_actors = 15L, n_states = 4L, seq_length = 8L, seed = seed, ragged = TRUE)
    sequences <- out_seqdef(wide)
    lengths <- as.integer(TraMineR::seqlength(sequences))
    cells <- as.matrix(sequences)
    expected <- data.frame(actor = seq_len(nrow(wide)), first_state = cells[, 1L],
                           last_state = cells[cbind(seq_len(nrow(wide)), lengths)], last_step = lengths,
                           n_observed = lengths, dropped_out = lengths < ncol(wide))
    expect_reference(actor_endpoints(wide), expected, reference = "TraMineR", covers = "actor_endpoints()",
                     quantity = sprintf("endpoints of 15 ragged sequences, seed %d, TraMineR", seed), tolerance = 0,
                     key = "actor")
  }))
})

test_that("cooccurrence() counts and Jaccard similarities equal arules", {
  skip_without_reference("arules")
  invisible(lapply(c(602L, 604L), \(seed) {
    set.seed(seed)
    records <- unique(data.frame(id = rep(seq_len(20L), each = 3L), item = sample(letters[1:6], 60L, replace = TRUE)))
    transactions <- arules::transactions(split(records$item, records$id))
    counts <- unclass(arules::crossTable(transactions, measure = "count", sort = FALSE))
    expect_reference(extract_transition_matrix(cooccurrence(records, field = "item", by = "id")), counts,
                     reference = "arules", covers = "cooccurrence()",
                     quantity = sprintf("co-occurrence counts, seed %d, arules crossTable()", seed), tolerance = 0)
    # Jaccard of two items: of the records holding either, the share holding both; arules' dissimilarity() of the
    # items as transactions of record ids.
    items <- arules::transactions(split(as.character(records$id), records$item))
    jaccard <- 1 - as.matrix(arules::dissimilarity(items, method = "jaccard"))
    diag(jaccard) <- 1
    expect_reference(extract_transition_matrix(cooccurrence(records, field = "item", by = "id", similarity = "jaccard")),
                     jaccard, reference = "arules", covers = "cooccurrence()",
                     quantity = sprintf("Jaccard similarity, seed %d, arules dissimilarity()", seed), tolerance = 1e-12)
  }))
})

test_that("action_to_onehot() indicators equal nnet::class.ind()", {
  skip_without_reference("nnet")
  events <- data.frame(Actor = rep(1:4, each = 3), Action = c("b", "a", "c", "a", "a", "b", "c", "b", "a", "c", "c", "a"))
  expected <- cbind(data.frame(row = seq_len(nrow(events))), as.data.frame(nnet::class.ind(events$Action)))
  found <- cbind(data.frame(row = seq_len(nrow(events))), action_to_onehot(events, sort_states = TRUE))
  expect_reference(found, expected, reference = "nnet", covers = "action_to_onehot()",
                   quantity = "indicators of 12 events, nnet class.ind()", tolerance = 0, key = "row")
})

test_that("net_aggregate_weights() equals stats and psych summaries of the non-zero weights", {
  skip_without_reference("psych")
  weights <- c(0.4, 0, 1.5, NA, 2.25, 0.1, 3)
  kept <- c(0.4, 1.5, 2.25, 0.1, 3)
  expect_reference(net_aggregate_weights(weights, method = "geomean"), psych::geometric.mean(kept), reference = "psych",
                   covers = "net_aggregate_weights()", quantity = "geometric mean, psych geometric.mean()",
                   tolerance = 1e-12)
  expect_reference(net_aggregate_weights(weights, method = "median"), stats::median(kept), reference = "stats",
                   covers = "net_aggregate_weights()", quantity = "median, stats median()", tolerance = 1e-12)
})
