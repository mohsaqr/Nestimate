testthat::skip_on_cran()

# Review 2026-11-20, wave 2: one-hot cluster rebuilds, undeclared direction,
# redundant column declarations, one-step windows in the landscapes.

# ---- C5: clustering a one-hot relative co-occurrence network ----

.w2_onehot <- function(s, n_actors) {
  data.frame(id = rep(sprintf("u%02d", seq_len(n_actors)),
                      each = length(s) / n_actors),
             A = as.integer(s == "A"), B = as.integer(s == "B"))
}

.w2_direct <- function(d, ids) {
  build_network(subset(d, id %in% ids), "co_occurrence", actor = "id",
                codes = c("A", "B"), window_size = 3,
                params = list(type = "relative"))
}

test_that("cluster_network() rebuilds a one-hot relative co-occurrence network with its settings", {
  s <- c("A", "A", "B", "A", "A", "B", "A", "B", "B", "A", "B", "B")
  d <- .w2_onehot(s, 4L)
  net <- build_network(d, "co_occurrence", actor = "id", codes = c("A", "B"),
                       window_size = 3, params = list(type = "relative"))
  expect_true(net$directed)
  cl <- cluster_network(net, k = 2)
  asg <- attr(cl, "clustering")$assignments
  ids <- sprintf("u%02d", seq_len(4L))
  invisible(Map(function(z, k) {
    direct <- .w2_direct(d, ids[asg == k])
    expect_true(z$directed)
    expect_equal(unname(z$weights), unname(direct$weights), tolerance = 1e-12)
    expect_equal(unname(rowSums(z$weights)), c(1, 1), tolerance = 1e-12)
  }, cl, seq_along(cl)))
  # the AAB actors: both rows (2/3, 1/3)
  aab <- cl[[asg[1L]]]
  expect_equal(unname(aab$weights),
               rbind(c(2, 1), c(2, 1)) / 3, tolerance = 1e-12)
  # build_network() on the clustering gives the same networks
  again <- build_network(build_clusters(net, k = 2))
  invisible(Map(function(a, b) {
    expect_equal(a$weights, b$weights, tolerance = 1e-12)
    expect_identical(a$directed, b$directed)
  }, cl, again))
})

test_that("cluster_network(cluster_by = \"mmm\") of a one-hot relative co-occurrence network keeps its settings", {
  a_type <- rep(c("A", "A", "A", "B"), 3L)
  b_type <- rep(c("B", "B", "B", "A"), 3L)
  d <- .w2_onehot(c(rep(a_type, 5L), rep(b_type, 5L)), 10L)
  net <- build_network(d, "co_occurrence", actor = "id", codes = c("A", "B"),
                       window_size = 3, params = list(type = "relative"))
  cl <- cluster_network(net, k = 2, cluster_by = "mmm", seed = 1)
  asg <- attr(cl, "clustering")$assignments
  ids <- sprintf("u%02d", seq_len(10L))
  invisible(Map(function(z, k) {
    expect_true(z$directed)
    expect_equal(unname(z$weights), unname(.w2_direct(d, ids[asg == k])$weights),
                 tolerance = 1e-12)
  }, cl, seq_along(cl)))
})

test_that("a caller's params or another method replace the one-hot settings", {
  s <- c("A", "A", "B", "A", "A", "B", "A", "B", "B", "A", "B", "B")
  d <- .w2_onehot(s, 4L)
  net <- build_network(d, "co_occurrence", actor = "id", codes = c("A", "B"),
                       window_size = 3, params = list(type = "relative"))
  counts <- cluster_network(net, k = 2, params = list(type = "frequency"))
  invisible(lapply(counts, function(z) {
    expect_false(z$directed)
    expect_true(isSymmetric(unname(z$weights)))
  }))
  # the counts of the clusters add up to the whole network's counts
  whole <- build_network(d, "co_occurrence", actor = "id",
                         codes = c("A", "B"), window_size = 3)
  expect_equal(Reduce(`+`, lapply(counts, `[[`, "weights")), whole$weights)
  freq <- build_network(build_clusters(net, k = 2), method = "frequency")
  expect_true(all(vapply(freq, function(z) z$method == "frequency",
                         logical(1L))))
})

test_that("co-occurrence of sequences honours type = \"relative\"", {
  seqs <- data.frame(T1 = c("A", "B", "A"), T2 = c("B", "B", "A"),
                     T3 = c("A", "A", "B"))
  raw <- build_network(seqs, "co_occurrence")
  rel <- build_network(seqs, "co_occurrence", params = list(type = "relative"))
  expect_false(raw$directed)
  expect_true(rel$directed)
  expect_equal(rel$weights, raw$weights / rowSums(raw$weights),
               tolerance = 1e-12)
  expect_error(build_network(seqs, "co_occurrence",
                             params = list(type = "share")),
               class = "nestimate_bad_input")
})

# ---- C6: a bare cograph_network without a direction flag ----

test_that("subtract_networks() infers an undeclared direction from symmetry", {
  m <- matrix(c(0, 1, 1, 0), 2L, dimnames = list(c("A", "B"), c("A", "B")))
  nodes <- data.frame(id = 1:2, name = c("A", "B"), label = c("A", "B"))
  bare <- structure(list(weights = m, nodes = nodes), class = "cograph_network")
  diff <- subtract_networks(bare, 0 * m)
  expect_false(diff$directed)
  expect_identical(diff$n_edges, 1L)
  # an asymmetric one is read as directed
  a <- matrix(c(0, 1, 0, 0), 2L, dimnames = list(c("A", "B"), c("A", "B")))
  asym <- structure(list(weights = a, nodes = nodes), class = "cograph_network")
  expect_true(subtract_networks(asym, 0 * a)$directed)
  # a declared flag still decides, also on symmetric weights
  declared <- structure(list(weights = m, nodes = nodes, directed = TRUE),
                        class = "cograph_network")
  expect_true(subtract_networks(declared, 0 * m)$directed)
})

test_that("compare_model() infers an undeclared direction from symmetry", {
  m <- matrix(c(0, 1, 1, 0), 2L, dimnames = list(c("A", "B"), c("A", "B")))
  nodes <- data.frame(id = 1:2, name = c("A", "B"), label = c("A", "B"))
  bare <- structure(list(weights = m, nodes = nodes), class = "cograph_network")
  undeclared <- compare_model(bare, m / 2)
  as_matrix <- compare_model(m, m / 2)
  expect_equal(undeclared$network_metrics, as_matrix$network_metrics)
})

# ---- C8: redundant column declarations ----

test_that("build_network() compares the effective column declarations", {
  d <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"))
  p <- c("T1", "T2")
  # `cols` in the call replaces params$cols, and then agrees with state_cols
  net <- build_network(d, "frequency", state_cols = p,
                       params = list(cols = rev(p)), cols = p)
  expect_identical(net$n_edges, 2L)
  plain <- build_network(d, "frequency", state_cols = p)
  expect_equal(net$weights, plain$weights)
  # a real conflict still stops, by params or by the effective `cols`
  expect_error(build_network(d, "frequency", state_cols = p,
                             params = list(cols = rev(p))),
               class = "nestimate_bad_input")
  expect_error(build_network(d, "frequency", state_cols = p,
                             params = list(cols = p), cols = rev(p)),
               class = "nestimate_bad_input")
})

# ---- C9: one-step windows in the landscapes ----

.w2_landscape_data <- function() {
  seqs <- c("AABBA", "ABBBA", "BBAAB", "AAABB", "BABBA", "ABABB", "BBBAA",
            "AABAB", "AAAAB", "BBBBA", "ABBBB", "BAAAA")
  d <- as.data.frame(do.call(rbind, strsplit(seqs, "")),
                     stringsAsFactors = FALSE)
  names(d) <- paste0("T", seq_len(5L))
  d
}

test_that("stability_landscape() accepts a width-1 window, equal to the unwindowed fit", {
  d <- .w2_landscape_data()
  one <- build_network(d, "relative", params = list(window_size = 1L))
  one_span <- build_network(d, "relative",
                            params = list(window_size = 1L, window_span = 1L))
  plain <- build_network(d, "relative")
  ref <- stability_landscape(plain, iter = 20, seed = 3)
  expect_equal(stability_landscape(one, iter = 20, seed = 3), ref)
  expect_equal(stability_landscape(one_span, iter = 20, seed = 3), ref)
  expect_equal(summary(landscape_trend(one, window = 3L, iter = 10, seed = 2)),
               summary(landscape_trend(plain, window = 3L, iter = 10, seed = 2)))
  expect_warning(surface_one <- landscape_surface(one, k = 5L),
                 class = "nestimate_thin_slice")
  expect_warning(surface_plain <- landscape_surface(plain, k = 5L),
                 class = "nestimate_thin_slice")
  expect_equal(surface_one, surface_plain)
})

test_that("the landscapes still refuse windows wider than one step", {
  d <- .w2_landscape_data()
  wide <- build_network(d, "relative", params = list(window_size = 2L))
  spanned <- build_network(d, "relative",
                           params = list(window_size = 1L, window_span = 2L))
  windowed <- build_network(d, "relative", params = list(windowed = TRUE))
  invisible(lapply(list(wide, spanned, windowed), function(net) {
    expect_error(stability_landscape(net, iter = 0),
                 class = "nestimate_bad_input")
    expect_error(landscape_trend(net, window = 3L, iter = 10),
                 class = "nestimate_bad_input")
    expect_error(landscape_surface(net, k = 5L),
                 class = "nestimate_bad_input")
  }))
})

test_that("compare_networks() infers an undeclared direction from symmetry", {
  weights <- matrix(c(0, 1, 1, 0), 2,
                    dimnames = list(c("A", "B"), c("A", "B")))
  bare <- structure(list(weights = weights,
                         nodes = data.frame(id = 1:2, name = c("A", "B"),
                                            label = c("A", "B"))),
                    class = "cograph_network")
  comparison <- compare_networks(bare = bare, half = weights / 2)
  expect_false(any(comparison$directed))
  edges <- subset(network_metrics(comparison), metric == "Edge Count")
  expect_identical(as.numeric(edges$value), c(1, 1))
})

test_that("build_network() of an mmm fit keeps a one-hot network's settings", {
  set.seed(3)
  pattern <- function(p) {
    unlist(lapply(1:6, \(i) sample(c("A", "B"), 3, replace = TRUE,
                                   prob = c(p, 1 - p))))
  }
  states <- c(unlist(lapply(1:10, \(i) pattern(0.85))),
              unlist(lapply(1:10, \(i) pattern(0.15))))
  onehot <- data.frame(id = rep(sprintf("u%02d", 1:20), each = 18),
                       A = as.integer(states == "A"),
                       B = as.integer(states == "B"))
  net <- build_network(onehot, "co_occurrence", actor = "id",
                       codes = c("A", "B"), window_size = 3,
                       params = list(type = "relative"))
  fit <- suppressWarnings(build_mmm(net, k = 2, n_starts = 5, seed = 1))
  rebuilt <- build_network(fit)
  expect_true(all(vapply(rebuilt, \(z) isTRUE(z$directed), logical(1))))
  row_sums <- unlist(lapply(rebuilt, \(z) rowSums(z$weights)))
  expect_equal(unname(row_sums), rep(1, length(row_sums)))
})
