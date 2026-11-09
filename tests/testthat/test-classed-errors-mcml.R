# Input errors in mcml.R, sequence_plot_mcml.R and macro_network.R carry nestimate_bad_input.
testthat::skip_on_cran()

.ce_seqs <- data.frame(
  T1 = c("A", "C", "A", "B"), T2 = c("B", "D", "C", "A"),
  T3 = c("C", "A", "D", "C"), stringsAsFactors = FALSE
)
.ce_clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))
.ce_mat <- matrix(c(0, 1, 2, 0), 2, 2, dimnames = list(c("A", "B"), c("A", "B")))

test_that("net_aggregate_weights() classes its input errors", {
  expect_error(net_aggregate_weights("a"), class = "nestimate_bad_input")
  expect_error(net_aggregate_weights(1:3, method = "nope"),
               class = "nestimate_bad_input")
})

test_that("cluster_summary() classes its input errors", {
  expect_error(cluster_summary(.ce_mat), class = "nestimate_bad_input")
  expect_error(cluster_summary(.ce_mat, list(G = "A", H = "B"), directed = NA),
               class = "nestimate_bad_input")
  bad <- .ce_mat
  bad[1, 2] <- NA
  expect_error(cluster_summary(bad, list(G = "A", H = "B")),
               class = "nestimate_bad_input")
})

test_that("build_mcml() classes its input errors", {
  expect_error(build_mcml(.ce_seqs, .ce_clusters, compute_within = NA),
               class = "nestimate_bad_input")
  expect_error(build_mcml(.ce_seqs, .ce_clusters, trim = -1),
               class = "nestimate_bad_input")
  expect_error(build_mcml(.ce_seqs, list(G1 = c("A", "B"), G2 = "C")),
               class = "nestimate_bad_input")
  expect_error(build_mcml(.ce_seqs, .ce_clusters, type = "semi_markov"),
               class = "nestimate_bad_input")
  edges <- data.frame(src = c("A", "B"), dst = c("B", "C"), w = c(1, -1))
  expect_error(build_mcml(edges, list(G1 = c("A", "B"), G2 = "C"),
                          from = "src", to = "dst", weight = "w"),
               class = "nestimate_bad_input")
  mc <- build_mcml(.ce_seqs, .ce_clusters)
  expect_error(build_mcml(mc, clusters = .ce_clusters, trim = 2),
               class = "nestimate_bad_input")
})

test_that("print.mcml() classes unsupported arguments", {
  mc <- build_mcml(.ce_seqs, .ce_clusters)
  expect_error(print(mc, foo = 1), class = "nestimate_bad_input")
})

test_that("macro_network() classes its input errors", {
  mc <- build_mcml(.ce_seqs, .ce_clusters)
  expect_error(macro_network(.ce_mat), class = "nestimate_bad_input")
  expect_error(macro_network(mc, expand = "nope"),
               class = "nestimate_bad_input")
})

test_that("as_tna() classes an unconvertible input", {
  expect_error(as_tna(1:3), class = "nestimate_bad_input")
})

test_that("as_htna() classes its input errors", {
  expect_error(as_htna(.ce_seqs, clusters = list(G1 = c("A", "B"), G2 = "Z")),
               class = "nestimate_bad_input")
  mc <- build_mcml(.ce_mat, list(G = "A", H = "B"))
  expect_error(as_htna(mc), class = "nestimate_bad_input")
})

test_that("sequence_plot() on an mcml classes its input errors", {
  mc <- build_mcml(.ce_seqs, .ce_clusters)
  expect_error(sequence_plot(mc, type = "distribution", expand = "nope"),
               class = "nestimate_bad_input")
  expect_error(sequence_plot(mc, type = "distribution", rest_label = ""),
               class = "nestimate_bad_input")
})
