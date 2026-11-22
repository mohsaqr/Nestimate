testthat::skip_on_cran()

# Refusals declared as defects by the layer-7 refusal-parity suite
# (tests/equivalence/test-layer7-refusal-*.R): each input below gave garbage or
# an internal error, and now raises the class its siblings raise.

.r21_states9 <- c("a1", "a2", "a3", "b1", "b2", "b3", "c1", "c2", "c3")
.r21_clusters9 <- list(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"),
                       C = c("c1", "c2", "c3"))
.r21_wide9 <- function(seed = 21L) {
  set.seed(seed)
  probability <- sort(stats::runif(9L))
  wide <- as.data.frame(matrix(sample(.r21_states9, 40L * 12L, replace = TRUE,
                                      prob = probability / sum(probability)),
                               40L, 12L), stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(12L))
  wide
}
.r21_edges9 <- function() {
  data.frame(from = c("a1", "a2", "a3", "b1", "b2", "b3", "c1", "c2", "c3",
                      "a1", "b1", "c1", "a2"),
             to = c("a2", "a3", "a1", "b2", "b3", "b1", "c2", "c3", "c1",
                    "b1", "c1", "a1", "a3"),
             w = c(rep(1, 12L), 50))
}
.r21_wide <- function(n = 30L, length = 10L, seed = 11L, group = FALSE) {
  set.seed(seed)
  wide <- as.data.frame(matrix(sample(c("A", "B", "C", "D"), n * length,
                                      replace = TRUE), n, length),
                        stringsAsFactors = FALSE)
  names(wide) <- paste0("T", seq_len(length))
  if (group) wide$g <- rep(c("x", "y"), each = n %/% 2L, length.out = n)
  wide
}
.r21_event_log <- function(seed = 7L) {
  set.seed(seed)
  t0 <- as.POSIXct("2024-03-01", tz = "UTC")
  do.call(rbind, lapply(seq_len(12L), \(i) {
    data.frame(team = if (i <= 6L) "red" else "blue", chat = paste0("chat", i),
               sender = sample(c("A", "B", "C", "D"), 25L, replace = TRUE),
               at = t0 + cumsum(stats::rexp(25L, 1 / 60)),
               stringsAsFactors = FALSE)
  }))
}
.r21_latency <- function(group = NULL) {
  build_network(.r21_event_log(), method = "latency", actor = "chat",
                action = "sender", time = "at", group = group)
}
.r21_edge_fit <- function() {
  build_mcml(.r21_edges9(), clusters = .r21_clusters9, from = "from",
             to = "to", weight = "w")
}
.r21_one_sequence <- function() {
  suppressMessages(
    build_tna(data.frame(T1 = "A", T2 = "B", T3 = "C", T4 = "A", T5 = "B",
                         T6 = "C", T7 = "A", T8 = "A")),
    classes = "nestimate_single_sequence")
}

test_that("markov_order_test() refuses an mcml instead of reading its fields as sequences", {
  sum_fit <- build_mcml(.r21_wide9(), clusters = .r21_clusters9)
  expect_error(markov_order_test(sum_fit, max_order = 1, n_perm = 5),
               class = "nestimate_bad_input")
  expect_error(markov_order_test(.r21_edge_fit(), max_order = 1, n_perm = 5),
               class = "nestimate_bad_input")
  # a classed list of any kind is not a list of sequences
  expect_error(markov_order_test(stats::t.test(1:10), max_order = 1, n_perm = 5),
               class = "nestimate_bad_input")
})

test_that("markov_order_test() reads a plain list of sequences as the same rows of a frame", {
  wide <- .r21_wide()
  rows <- unname(split(as.matrix(wide), row(as.matrix(wide))))
  from_list <- markov_order_test(rows, max_order = 2, n_perm = 20, seed = 4)
  from_frame <- markov_order_test(wide, max_order = 2, n_perm = 20, seed = 4)
  expect_s3_class(from_list, "net_markov_order")
  expect_equal(from_list$test_table, from_frame$test_table)
})

test_that("an edge-list mcml layer is not read as sequences", {
  layers <- as_tna(.r21_edge_fit())
  layer <- layers$macro
  expect_identical(attr(layer$data, "source"), "edgelist")
  readers <- list(
    state_frequencies = \(x) state_frequencies(x),
    sequence_metrics = \(x) sequence_metrics(x),
    sequence_patterns = \(x) sequence_patterns(x),
    extract_sequences = \(x) extract_sequences(x),
    build_hon = \(x) build_hon(x, max_order = 2L),
    build_hypa = \(x) build_hypa(x, order = 2L),
    build_mogen = \(x) build_mogen(x, max_order = 2L),
    markov_order_test = \(x) markov_order_test(x, max_order = 1L, n_perm = 5L),
    path_counts = \(x) path_counts(x),
    stability_landscape = \(x) stability_landscape(x, iter = 5L, seed = 1)
  )
  invisible(Map(\(reader, name) {
    expect_error(reader(layer), class = "nestimate_bad_input", info = name)
  }, readers, names(readers)))
  # the whole set of layers too, through each of its edge-list layers
  expect_error(sequence_metrics(layers), class = "nestimate_bad_input")
})

test_that("the sequence readers still read a sequence-built mcml layer", {
  layer <- as_tna(build_mcml(.r21_wide9(), clusters = .r21_clusters9))$macro
  expect_s3_class(state_frequencies(layer), "data.frame")
  expect_s3_class(extract_sequences(layer), "data.frame")
})

test_that("the vertex bootstrap refuses a latency network", {
  latency <- .r21_latency()
  expect_error(vertex_bootstrap(latency, iter = 5L, seed = 1),
               class = "nestimate_latency_resampling")
  expect_error(vertex_compare(latency, latency, iter = 5L, seed = 1),
               class = "nestimate_latency_resampling")
  net <- build_network(.r21_wide(), method = "relative")
  expect_error(vertex_compare(net, latency, iter = 5L, seed = 1),
               class = "nestimate_latency_resampling")
  expect_s3_class(vertex_bootstrap(net, iter = 5L, seed = 1),
                  "net_vertex_bootstrap")
})

test_that("the landscape verbs refuse latency networks with the latency class", {
  single <- .r21_latency()
  grouped <- .r21_latency(group = "team")
  invisible(lapply(list(single, grouped), \(net) {
    expect_error(stability_landscape(net, iter = 5L, seed = 1),
                 class = "nestimate_latency_resampling")
    expect_error(landscape_trend(net, iter = 5L, seed = 1, window = 3L),
                 class = "nestimate_latency_resampling")
  }))
})

test_that("the landscape verbs refuse to resample one sequence", {
  net <- .r21_one_sequence()
  expect_error(stability_landscape(net, iter = 5L, seed = 1),
               class = "nestimate_single_sequence")
  expect_error(landscape_trend(net, iter = 5L, seed = 1, window = 3L),
               class = "nestimate_single_sequence")
  one_row <- data.frame(T1 = "A", T2 = "B", T3 = "C", T4 = "A")
  expect_error(stability_landscape(one_row, iter = 5L, seed = 1),
               class = "nestimate_single_sequence")
  # a group with one one-sequence member is refused, naming the group
  wide <- .r21_wide(n = 11L)
  wide$g <- c(rep("x", 10L), "y")
  grouped <- suppressMessages(
    build_network(wide, method = "relative", group = "g"),
    classes = "nestimate_single_sequence")
  expect_error(stability_landscape(grouped, iter = 5L, seed = 1),
               regexp = "^y: ", class = "nestimate_single_sequence")
  # without resamples there is no interval, so the estimates are still given
  expect_s3_class(stability_landscape(net, iter = 0L), "net_landscape")
  expect_s3_class(landscape_trend(net, iter = 0L, window = 3L),
                  "net_landscape_trend")
})

test_that("stability_landscape() refuses an mcml layer combined by a mean", {
  mean_layers <- as_tna(build_mcml(.r21_wide9(), clusters = .r21_clusters9,
                                   method = "mean", type = "frequency"))
  expect_error(stability_landscape(mean_layers$A, iter = 5L, seed = 1),
               class = "nestimate_bad_input")
  expect_error(landscape_trend(mean_layers$A, iter = 5L, seed = 1, window = 3L),
               class = "nestimate_bad_input")
  # the summed layer is read
  sum_layers <- as_tna(build_mcml(.r21_wide9(), clusters = .r21_clusters9,
                                  type = "frequency"))
  expect_s3_class(stability_landscape(sum_layers$A, iter = 5L, seed = 1),
                  "net_landscape")
})

test_that("casedrop_reliability() refuses a mean-combined layer before its constant-edge answer", {
  mean_fit <- build_mcml(.r21_wide9(), clusters = .r21_clusters9,
                         method = "mean", type = "frequency")
  mean_layers <- as_tna(mean_fit)
  expect_error(casedrop_reliability(mean_layers$A, iter = 5L, seed = 1),
               class = "nestimate_bad_input")
  expect_error(casedrop_reliability(mean_fit, iter = 5L, seed = 1),
               class = "nestimate_bad_input")
  # a summed network with constant edges is still answered "CS undefined"
  one_state <- build_network(
    as.data.frame(matrix("A", 10L, 6L, dimnames = list(NULL, paste0("T", 1:6)))),
    method = "relative")
  expect_warning(res <- casedrop_reliability(one_state, iter = 5L, seed = 1),
                 class = "nestimate_constant_measure")
  expect_true(is.na(res$cs))
})

test_that("entropy_bayes() refuses a network with no observed transition", {
  no_transition <- build_network(data.frame(T1 = c("A", "B", "C", "A"),
                                            T2 = NA_character_),
                                 method = "relative")
  expect_error(entropy_bayes(no_transition, draws = 200, seed = 1),
               class = "nestimate_bad_input")
  empty <- .r21_wide(group = TRUE)
  empty[empty$g == "y", paste0("T", seq_len(10L))] <- NA_character_
  empty_group <- build_network(empty, method = "relative", group = "g")
  expect_error(entropy_bayes(empty_group, draws = 200, seed = 1),
               class = "nestimate_bad_input")
  # a zero count matrix the same way
  zero <- matrix(0, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(entropy_bayes(zero, draws = 200, seed = 1),
               class = "nestimate_bad_input")
})

test_that("entropy_bayes() reads a group_tna as the netobject_group it converts to", {
  skip_if_not_installed("tna")
  wide <- .r21_wide(group = TRUE)
  from_tna <- entropy_bayes(tna::group_tna(wide, group = "g",
                                           cols = paste0("T", seq_len(10L))),
                            draws = 300, seed = 9)
  from_group <- entropy_bayes(build_network(wide, method = "relative",
                                            group = "g"),
                              draws = 300, seed = 9)
  expect_s3_class(from_tna, "net_entropy_bayes_group")
  expect_identical(names(from_tna), c("x", "y"))
  expect_equal(from_tna, from_group)
})
