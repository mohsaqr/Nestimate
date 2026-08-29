# ---- ngram / gap / reverse estimators: parity with tna::build_model() ----

.parity_seqs <- function() {
  set.seed(20260829)
  as.data.frame(
    matrix(sample(c("a", "b", "c", NA), 60, TRUE, prob = c(.3, .3, .3, .1)),
           nrow = 12, ncol = 5),
    stringsAsFactors = FALSE
  )
}

# Hand-computed reference: one sequence a b a c.
.tiny <- data.frame(V1 = "a", V2 = "b", V3 = "a", V4 = "c",
                    stringsAsFactors = FALSE)

test_that("ngram counts adjacent pairs once per containing window", {
  # n_gram = 3 on a b a c: windows (a b a), (b a c).
  # pairs: a->b (1 window), b->a (2 windows), a->c (1 window).
  net <- build_network(.tiny, method = "ngram", params = list(n_gram = 3))
  expect_equal(net$weights["a", "b"], 1)
  expect_equal(net$weights["b", "a"], 2)
  expect_equal(net$weights["a", "c"], 1)
  expect_equal(sum(net$weights), 4)
  expect_true(net$directed)
})

test_that("gap weights pairs by 1/distance up to max_gap + 1", {
  # max_gap = 1 on a b a c: distance-1 pairs a->b, b->a, a->c (w = 1);
  # distance-2 pairs a->a, b->c (w = 1/2).
  net <- build_network(.tiny, method = "gap", params = list(max_gap = 1))
  expect_equal(net$weights["a", "b"], 1)
  expect_equal(net$weights["b", "a"], 1)
  expect_equal(net$weights["a", "c"], 1)
  expect_equal(net$weights["a", "a"], 0.5)
  expect_equal(net$weights["b", "c"], 0.5)
  expect_equal(sum(net$weights), 4)
})

test_that("reverse is the transpose of frequency", {
  seqs <- .parity_seqs()
  rev <- build_network(seqs, method = "reverse")
  freq <- build_network(seqs, method = "frequency")
  expect_equal(rev$weights, t(freq$weights))
  rev_w <- build_network(seqs, method = "reverse",
                         params = list(weighted = TRUE))
  freq_w <- build_network(seqs, method = "frequency",
                          params = list(weighted = TRUE))
  expect_equal(rev_w$weights, t(freq_w$weights))
})

test_that("ngram/gap/reverse match tna::build_model() weights", {
  skip_if_not_installed("tna")
  seqs <- .parity_seqs()
  cases <- list(
    list(method = "ngram", type = "n-gram", params = list()),
    list(method = "ngram", type = "n-gram", params = list(n_gram = 3)),
    list(method = "ngram", type = "n-gram", params = list(n_gram = 4)),
    list(method = "gap", type = "gap", params = list()),
    list(method = "gap", type = "gap", params = list(max_gap = 2)),
    list(method = "gap", type = "gap", params = list(max_gap = 0)),
    list(method = "reverse", type = "reverse", params = list()),
    list(method = "reverse", type = "reverse", params = list(weighted = TRUE))
  )
  invisible(lapply(cases, function(cs) {
    ours <- build_network(seqs, method = cs$method, params = cs$params)
    ref <- tna::build_model(seqs, type = cs$type, params = cs$params)
    expect_equal(unname(ours$weights), unname(ref$weights),
                 tolerance = sqrt(.Machine$double.eps),
                 label = paste(cs$method, deparse(cs$params)))
    expect_equal(rownames(ours$weights), rownames(ref$weights))
  }))
})

test_that("ngram matches tna on group_regulation", {
  skip_if_not_installed("tna")
  seqs <- utils::head(tna::group_regulation, 200)
  ours <- build_network(seqs, method = "ngram", params = list(n_gram = 3))
  ref <- tna::build_model(seqs, type = "n-gram", params = list(n_gram = 3))
  expect_equal(unname(ours$weights), unname(ref$weights))
})

test_that("long-format input gives the same result as its wide form", {
  # NA-free on purpose: in wide form an interior NA keeps its column (so
  # `a NA b` is a distance-2 pair), whereas long form has no such slot.
  set.seed(20260830)
  seqs <- as.data.frame(
    matrix(sample(c("a", "b", "c"), 60, TRUE), nrow = 12, ncol = 5),
    stringsAsFactors = FALSE
  )
  seqs$id <- seq_len(nrow(seqs))
  long <- stats::reshape(seqs, direction = "long", varying = paste0("V", 1:5),
                         v.names = "action", timevar = "t", idvar = "id")
  long <- long[!is.na(long$action), c("id", "t", "action")]
  wide <- build_network(seqs[, paste0("V", 1:5)], method = "gap",
                        params = list(max_gap = 2))
  from_long <- build_network(long, method = "gap", actor = "id",
                             action = "action", order = "t",
                             params = list(max_gap = 2))
  expect_equal(from_long$weights, wide$weights)
})

test_that("aliases resolve to the canonical estimator names", {
  seqs <- .parity_seqs()
  expect_equal(build_network(seqs, method = "n-gram")$method, "ngram")
  expect_equal(build_network(seqs, method = "n_gram")$method, "ngram")
  expect_equal(build_network(seqs, method = "co-occurrence")$method,
               "co_occurrence")
  expect_true(all(c("ngram", "gap", "reverse") %in% list_estimators()$name))
})

test_that("normalize scaling turns ngram counts into row probabilities", {
  seqs <- .parity_seqs()
  net <- build_network(seqs, method = "ngram", scaling = "normalize")
  rs <- rowSums(net$weights)
  expect_true(all(abs(rs[rs > 0] - 1) < 1e-10))
})

test_that("invalid params fail loudly", {
  seqs <- .parity_seqs()
  expect_error(build_network(seqs, method = "ngram", params = list(n_gram = 1)),
               "n_gram")
  expect_error(build_network(seqs, method = "gap", params = list(max_gap = -1)),
               "max_gap")
  expect_error(build_network(seqs, method = "reverse",
                             params = list(weighted = "yes")), "weighted")
})

test_that("sequences shorter than n_gram give an all-zero matrix", {
  short <- data.frame(V1 = c("a", "b"), V2 = c("b", "a"),
                      stringsAsFactors = FALSE)
  net <- build_network(short, method = "ngram", params = list(n_gram = 5))
  expect_equal(sum(net$weights), 0)
  expect_equal(dim(net$weights), c(2L, 2L))
})

test_that("column permutation of states leaves the weight matrix invariant", {
  seqs <- .parity_seqs()
  relabel <- c(a = "z", b = "y", c = "x")
  seqs_perm <- as.data.frame(lapply(seqs, function(v) unname(relabel[v])),
                             stringsAsFactors = FALSE)
  w1 <- build_network(seqs, method = "gap")$weights
  w2 <- build_network(seqs_perm, method = "gap")$weights
  # Same counts, states reordered: map back and compare.
  w2 <- w2[relabel[rownames(w1)], relabel[colnames(w1)]]
  expect_equal(unname(w1), unname(w2))
})
