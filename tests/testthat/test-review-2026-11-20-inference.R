# Regression tests for the 2026-11-20 review of construction, latency, tna
# input, bootstrap and permutation (findings C1-C14, I1-I11). Expected values
# come from hand computations or from the definition of the quantity.

# ---- Latency: not resampled (C1, C2/I1, C5, C6) ----

test_that("latency networks built with different settings are not permuted", {
  # C2/I1: identical logs, different floors used to give p = 0.01
  log <- data.frame(chat = rep(c("c1", "c2", "c3"), each = 6),
                    sender = rep(c("A", "B"), 9),
                    at = rep(10 * (0:5), 3))
  fit <- function(floor) {
    build_network(log, method = "latency", actor = "chat", action = "sender",
                  time = "at", latency_floor = floor, latency_cap = "20s")
  }
  expect_error(permutation(fit(0.5), fit(1), iter = 10L, seed = 1),
               class = "nestimate_latency_resampling")
  # a latency network beside a network of another method is refused too
  rel <- build_network(log, method = "relative", actor = "chat",
                       action = "sender")
  expect_error(network_reliability(rel, fit(0.5), iter = 10L, seed = 1),
               class = "nestimate_latency_resampling")
  # the refusal is a bad input, as the other refusals of these verbs are
  expect_error(centrality_stability(fit(0.5), iter = 10L),
               class = "nestimate_bad_input")
})

# ---- C8: actor-scoped caps of hyphenated compound actors ----

test_that("compound actors that join to the same label keep their own caps", {
  # ("a-b", "c") replies every 5 s, ("a", "b-c") every 50 s. Joined with "-"
  # both read "a-b-c"; each must get the cap of its own latencies.
  d <- data.frame(a = rep(c("a-b", "a"), each = 6),
                  b = rep(c("c", "b-c"), each = 6),
                  sender = rep(c("A", "B"), 6),
                  t = c(5 * (0:5), 50 * (0:5)))
  net <- build_network(d, "latency", latency_aggregate = "sum", actor = c("a", "b"), action = "sender",
                       time = "t", latency_scope = "actor", latency_cap = 1)
  expect_equal(net$latency_unit_caps[c("a\\-b-c", "a-b\\-c")],
               c("a\\-b-c" = 5, "a-b\\-c" = 50))
  # every reply then sits at its own unit's cap: weight 1 - 0.5 = 0.5.
  # Each actor has 3 replies B -> A and 2 replies A -> B.
  expect_equal(net$weights["B", "A"], 6 * 0.5, tolerance = 1e-12)
  expect_equal(net$weights["A", "B"], 4 * 0.5, tolerance = 1e-12)

  # the same with the two actors in different groups: the caps are pooled
  # over the groups by unit, and the group networks add up to the full one
  d$team <- rep(c("x", "y"), each = 6)
  grouped <- build_network(d, "latency", latency_aggregate = "sum", actor = c("a", "b"),
                           action = "sender", time = "t",
                           latency_scope = "actor", latency_cap = 1,
                           group = "team")
  invisible(lapply(grouped, function(n) {
    expect_equal(n$latency_unit_caps[names(net$latency_unit_caps)],
                 net$latency_unit_caps)
  }))
  total <- Reduce(`+`, lapply(grouped, function(n) n$weights))
  expect_equal(total, net$weights, tolerance = 1e-12)

  # without a "-" in any value the label is the plain join
  plain <- transform(d, a = rep(c("ab", "a"), each = 6),
                     b = rep(c("c", "bc"), each = 6))
  net2 <- build_network(plain, "latency", latency_aggregate = "sum", actor = c("a", "b"),
                        action = "sender", time = "t",
                        latency_scope = "actor", latency_cap = 1)
  expect_setequal(names(net2$latency_unit_caps), c("ab-c", "a-bc"))
})

# ---- C3/I2: windowed bootstrap with concat ----

test_that("a windowed bootstrap resamples the sequences concat merges", {
  # rows A, B, A, B merged in pairs: two sequences (A, B), each with A -> B
  n <- build_network(data.frame(T1 = c("A", "B", "A", "B")), "frequency",
                     params = list(window_size = 1L, concat = 2L))
  expect_equal(n$weights["A", "B"], 2)
  # both sequences are alike, so every resample is the data
  b <- bootstrap_network(n, iter = 10L, seed = 1)
  expect_equal(unname(b$mean), unname(n$weights))

  d <- data.frame(T1 = rep(c("A", "B"), 2), T2 = rep(c("A", "B"), 2))
  x <- build_network(d, method = "frequency",
                     params = list(window_size = 2, concat = 2))
  expect_equal(x$weights["A", "B"], 8)
  b2 <- bootstrap_network(x, iter = 10L, seed = 1)
  expect_equal(unname(b2$mean), unname(x$weights))

  # the per-sequence counts add up to the network's counts, one row per
  # merged sequence (ceiling(5 / 2) = 3)
  d3 <- data.frame(T1 = c("A", "B", "B", "A", "C"))
  n3 <- build_network(d3, "frequency",
                      params = list(window_size = 1L, concat = 2L))
  states <- n3$nodes$label
  per_seq <- .precompute_per_sequence(n3$data, "frequency", n3$params, states)
  expect_identical(nrow(per_seq), 3L)
  expect_equal(matrix(colSums(per_seq), length(states), byrow = TRUE,
                      dimnames = list(states, states)),
               n3$weights)
})

# ---- C4/I5: wtna scaling domain under resampling ----

test_that("a scaled wtna network is resampled on its own scaling domain", {
  # two identical actors: every resample is the data, so the bootstrap mean
  # is the fitted network, cell for cell
  d <- data.frame(id = rep(1:2, each = 4), A = rep(c(1, 0), 4),
                  B = rep(c(0, 1), 4))
  n <- suppressWarnings(
    build_network(d, "wtna", actor = "id", codes = c("A", "B"),
                  window_size = 1, scaling = "minmax"),
    classes = "nestimate_sparse_onehot")
  b <- bootstrap_network(n, iter = 10L, seed = 1)
  expect_equal(unname(b$mean), unname(n$weights), tolerance = 1e-12)

  d2 <- data.frame(id = rep(1:2, each = 2), A = rep(c(1, 0), 2),
                   B = rep(c(0, 1), 2))
  x <- suppressWarnings(
    build_network(d2, method = "wtna", format = "onehot", actor = "id",
                  codes = c("A", "B"), window_size = 1, scaling = "minmax"),
    classes = "nestimate_sparse_onehot")
  b2 <- bootstrap_network(x, iter = 10L, seed = 1)
  expect_equal(unname(b2$mean), unname(x$weights), tolerance = 1e-12)
})

# ---- C7, C11, C12, C14: the columns a network is estimated from ----

test_that("$data keeps an all-missing position the estimator read", {
  d <- data.frame(T1 = c("A", "A"), T2 = NA_character_, T3 = c("B", "B"),
                  T4 = c("C", "C"))
  n <- build_ftna(d)
  # A, NA, B has no adjacent pair: the only transition is B -> C
  expect_equal(n$weights["A", "B"], 0)
  expect_equal(n$weights["B", "C"], 2)
  expect_named(n$data, c("T1", "T2", "T3", "T4"))
  expect_named(extract_sequences(n), c("T1", "T2", "T3", "T4"))
  b <- bootstrap_network(n, iter = 10L, seed = 1)
  expect_equal(unname(b$mean), unname(n$weights))
})

test_that("a time column of wide data is never a state or a variable", {
  d <- data.frame(T1 = "A", T2 = "B", T3 = "2024-01-01")
  n <- build_ftna(d, time = "T3")
  expect_setequal(n$nodes$label, c("A", "B"))
  expect_equal(n$weights["A", "B"], 1)
  expect_equal(sum(n$weights), 1)
  expect_named(n$data, c("T1", "T2"))
  expect_true("T3" %in% names(n$metadata))

  .local_seed(1)
  num <- data.frame(a = stats::rnorm(20), b = stats::rnorm(20), tt = 1:20)
  expect_setequal(colnames(build_cor(num, time = "tt")$weights), c("a", "b"))
})

test_that("state_cols and params$cols that disagree are refused", {
  d <- data.frame(T1 = "A", T2 = "B", T3 = "C")
  expect_error(build_ftna(d, state_cols = c("T1", "T2"),
                          params = list(cols = c("T2", "T3"))),
               class = "nestimate_bad_input")
  # the same positions given twice agree
  n <- build_ftna(d, state_cols = c("T1", "T2"),
                  params = list(cols = c("T1", "T2")))
  expect_equal(n$weights["A", "B"], 1)
  expect_named(n$data, c("T1", "T2"))
})

test_that("build_cor() converts the declared state_cols to frequencies", {
  d <- data.frame(T1 = c("A", "A", "B", "B"), T2 = c("A", "B", "A", "B"),
                  score = c(1, 2, 4, 3), age = c(11, 13, 12, 14))
  r <- build_cor(d, state_cols = c("T1", "T2"))
  expect_setequal(r$nodes$label, c("A", "B"))
  # per row the A and B counts are (2, 0), (1, 1), (1, 1), (0, 2): B = 2 - A
  expect_equal(r$weights["A", "B"], -1, tolerance = 1e-12)
})

# ---- C9, C13: tna models converted to netobjects ----

test_that("a co-occurrence tna model converts to an undirected network", {
  s <- c("A", "B")
  m <- structure(list(weights = matrix(c(0, 1, 1, 0), 2,
                                       dimnames = list(s, s)),
                      labels = s),
                 type = "co-occurrence", class = "tna")
  cv <- as_netobject(m)
  expect_false(cv$directed)
  expect_identical(nrow(extract_edges(cv)), nrow(extract_edges(m)))
  expect_identical(nrow(extract_edges(cv)), 1L)
})

test_that("a tna model keeps its own initial probabilities", {
  skip_if_not_installed("tna")
  m <- tna::tna(data.frame(T1 = c("A", NA_character_), T2 = c("B", "B"),
                           T3 = c("A", "A")))
  cv <- as_netobject(m)
  # tna counts the first position: A once, the NA not at all
  expect_equal(cv$initial[c("A", "B")], c(A = 1, B = 0))
  expect_equal(extract_initial_probs(cv), extract_initial_probs(m))
  expect_equal(unname(cv$weights), unname(m$weights))
})

# ---- C10/I4: relative WTNA co-occurrence ----

test_that("relative WTNA co-occurrence is directed, whatever the code order", {
  d <- data.frame(A = c(1, 1), B = c(1, 0))
  w <- function(codes) {
    wtna(d, method = "cooccurrence", type = "relative", codes = codes,
         window_size = 1)
  }
  # counts AA 2, AB 1, BA 1, BB 1; rows divided by 3 and 2
  expected <- matrix(c(2 / 3, 1 / 2, 1 / 3, 1 / 2), 2,
                     dimnames = list(c("A", "B"), c("A", "B")))
  expect_equal(w(c("A", "B"))$weights, expected, tolerance = 1e-12)
  expect_true(w(c("A", "B"))$directed)
  edges_ab <- extract_edges(w(c("A", "B")), threshold = 0.4)
  edges_ba <- extract_edges(w(c("B", "A")), threshold = 0.4)
  expect_identical(nrow(edges_ab), 1L)
  expect_identical(c(edges_ab$from, edges_ab$to), c("B", "A"))
  expect_equal(edges_ab$weight, 0.5)
  expect_equal(edges_ba$weight, edges_ab$weight)
  # the frequency counts stay symmetric and undirected
  expect_false(wtna(d, method = "cooccurrence", codes = c("A", "B"),
                    window_size = 1)$directed)
})

test_that("BH correction of relative co-occurrence tests every cell", {
  fit <- function(last) {
    d <- data.frame(id = rep(1:2, each = 3), A = rep(c(1, 1, 0), 2),
                    B = rep(c(1, 0, last), 2))
    wtna(d, method = "cooccurrence", type = "relative", actor = "id",
         window_size = 1)
  }
  x <- fit(0)
  y <- fit(1)
  raw <- permutation(x, y, iter = 199L, seed = 1)
  bh <- permutation(x, y, iter = 199L, seed = 1, adjust = "BH")
  # all four cells are non-zero, so all four are one family
  expected <- matrix(stats::p.adjust(as.vector(raw$p_values), "BH"), 2,
                     dimnames = dimnames(raw$p_values))
  expect_equal(bh$p_values, expected)
  # an undirected flag on asymmetric weights does not mirror p-values
  p <- matrix(c(0.01, 0.5, 0.02, 0.6), 2)
  w <- matrix(c(1, 2, 3, 4), 2)
  expect_equal(.adjust_edge_p(p, w, w, directed = FALSE, adjust = "BH"),
               matrix(stats::p.adjust(as.vector(p), "BH"), 2))
})

# ---- I3: paired permutation of one-hot networks ----

test_that("paired one-hot units pair by first appearance, not factor levels", {
  d <- data.frame(id = factor(rep(c("a", "b"), each = 2)),
                  A = c(1, 1, 0, 0), B = c(0, 0, 1, 1))
  fit <- function(z) suppressWarnings(wtna(z, actor = "id", window_size = 1),
                                      classes = "nestimate_sparse_onehot")
  y <- transform(d, id = factor(id, levels = c("b", "a")))
  # the same rows: every pair holds two identical units, so every swap
  # leaves the difference at 0
  p <- permutation(fit(d), fit(y), paired = TRUE, iter = 99L, seed = 1)
  expect_equal(unname(p$null_sd), matrix(0, 2, 2))
  expect_true(all(p$p_values == 1))
})

# ---- I6: permutation of a registered estimator with an auxiliary column ----

test_that("permutation gives a registered estimator every column it reads", {
  register_estimator("review_aux", function(data, ...) {
    stopifnot("aux" %in% names(data))
    w <- sum(data$aux == "A")
    m <- matrix(c(0, w, w, 0), 2, dimnames = list(c("A", "B"), c("A", "B")))
    list(matrix = m, nodes = colnames(m), directed = FALSE,
         cleaned_data = data)
  }, "Estimator with an auxiliary input", directed = FALSE)
  on.exit(remove_estimator("review_aux"), add = TRUE)
  d <- data.frame(A = c("A", "B"), B = c("B", "A"), aux = c("A", "B"))
  x <- build_network(d, method = "review_aux")
  expect_named(x$data, c("A", "B", "aux"))
  r <- permutation(x, x, iter = 5L, seed = 1)
  # x against itself: the observed difference is 0 and every shuffle is at
  # least as large, so p = (iter + 1) / (iter + 1)
  expect_equal(unname(r$p_values), matrix(1, 2, 2))
})

# ---- I7: boot_glasso() reads the network's threshold ----

test_that("boot_glasso() thresholds as the glasso network was", {
  d <- data.frame(A = 1:20, B = 1:20 + sin(1:20) / 10)
  x <- build_network(d, method = "glasso", threshold = 1)
  b <- boot_glasso(x, iter = 4L, cs_iter = 2L, seed = 1)
  expect_equal(b$original_pcor["A", "B"], 0)
  expect_equal(unname(b$original_pcor), unname(x$weights))

  .local_seed(1)
  z <- matrix(stats::rnorm(800), 200)
  z[, 2] <- z[, 1] + stats::rnorm(200)
  z[, 3] <- z[, 2] * 0.3 + stats::rnorm(200)
  d3 <- stats::setNames(as.data.frame(z), c("V1", "V2", "V3", "V4"))
  x3 <- build_network(d3, method = "glasso", threshold = 0.2)
  b3 <- boot_glasso(x3, iter = 10L, cs_iter = 2L, seed = 1)
  expect_equal(unname(b3$original_pcor), unname(x3$weights),
               tolerance = 1e-10)
  # every resample is thresholded alike: no bootstrap edge in (0, 0.2)
  boot_abs <- abs(b3$boot_edges)
  expect_false(any(boot_abs > 0 & boot_abs < 0.2, na.rm = TRUE))
})

# ---- I8: the single-sequence guard counts what the resampler draws ----

test_that("wide rows sharing an id are each a resampled sequence", {
  d <- data.frame(id = "same", T1 = c("A", "B", "A"), T2 = c("B", "A", "A"))
  x <- build_network(d, method = "frequency", format = "wide",
                     params = list(id = "id"))
  # every row counts as a sequence: three transitions
  expect_equal(sum(x$weights), 3)
  per_seq <- .precompute_per_sequence(x$data, "frequency", x$params,
                                      x$nodes$label)
  expect_identical(nrow(per_seq), 3L)
  expect_identical(.transition_resampling_n_sequences(x$data, x$params),
                   nrow(per_seq))
  expect_s3_class(bootstrap_network(x, iter = 4L, seed = 1), "net_bootstrap")
  expect_s3_class(network_reliability(x, iter = 4L, seed = 1),
                  "net_reliability")
  # one-hot rows without an actor remain one sequence
  oh <- data.frame(A = c(1, 0, 1), B = c(0, 1, 0))
  expect_identical(.transition_resampling_n_sequences(
    oh, list(codes = c("A", "B"))), 1L)
})

# ---- I9: the constant-measure tolerance is relative ----

test_that("centrality stability does not depend on the measure's units", {
  .local_seed(2)
  s <- as.data.frame(matrix(sample(c("A", "B", "C", "D"), 60 * 8, TRUE,
                                   prob = c(0.4, 0.3, 0.2, 0.1)), 60))
  x <- build_network(s, method = "frequency")
  cs <- function(k) {
    centrality_stability(x, measures = "custom", iter = 4L, drop_prop = 0.5,
                         seed = 1, centrality_fn = function(w) {
                           list(custom = k * rowSums(w))
                         })$cs
  }
  # a scale change leaves every correlation, so the CS coefficient, unchanged
  expect_equal(cs(1e-10), cs(1))
  expect_equal(cs(1e-12), cs(1))
  expect_false(.is_constant_vector(c(1e-10, 2e-10)))
  expect_true(.is_constant_vector(c(1, 1 + 1e-12)))
  expect_true(.is_constant_vector(c(0, 0, 0)))
  expect_true(.is_constant_vector(c(1e-10, 1e-10 * (1 + 1e-12))))
})

# ---- I11: what min_p means under Monte Carlo draws ----

test_that("the few-actors note speaks of the exact test", {
  d <- data.frame(team = rep(1:3, each = 2), person = 1:6,
                  T1 = rep(c("B", "A", "A"), each = 2),
                  T2 = rep(c("A", "B", "B"), each = 2))
  fit <- function(z) {
    build_network(z, method = "frequency", metadata_cols = c("team", "person"))
  }
  x <- fit(d[1:2, ])
  y <- fit(d[3:6, ])
  # 3 actors, 1 in x: 3 arrangements, so the exact test stops at 1/3
  expect_warning(
    p <- permutation(x, y, actor = "team", iter = 50L, seed = 1),
    "exact test cannot give a p-value below 0.333",
    class = "nestimate_few_actors")
  expect_equal(p$clustering$min_p, 1 / 3)
  expect_output(print(p), "with 3 actors the exact test gives no p-value")
})

test_that("bootstrap_network() on a one-state network returns 1 x 1 results", {
  d <- data.frame(T1 = c("A", "A", "A"), T2 = c("A", "A", "A"),
                  T3 = c("A", NA, "A"))
  net <- suppressMessages(build_network(d, method = "frequency"))
  b <- expect_no_warning(bootstrap_network(net, iter = 20, seed = 1))
  expect_identical(dim(b$mean), c(1L, 1L))
  expect_identical(dim(b$ci_lower), c(1L, 1L))
  # Every resample counts between 3 and 6 A -> A transitions (rows hold 2, 1, 2).
  expect_true(b$ci_lower[1, 1] >= 3 && b$ci_upper[1, 1] <= 6)
  rel <- suppressMessages(build_network(d, method = "relative"))
  expect_equal(bootstrap_network(rel, iter = 20, seed = 1)$mean[1, 1], 1)
})
