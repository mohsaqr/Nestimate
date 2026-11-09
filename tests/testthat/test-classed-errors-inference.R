# Input errors of the inference / comparison verbs carry class nestimate_bad_input.
testthat::skip_on_cran()

.cei_seqs <- data.frame(
  V1 = c("A", "B", "A", "C", "B"), V2 = c("B", "C", "B", "A", "C"),
  V3 = c("C", "A", "C", "B", "A")
)
.cei_net <- function(method = "relative") {
  build_network(.cei_seqs, method = method)
}

test_that("permutation() classes its input errors", {
  net <- .cei_net()
  expect_error(permutation(net, net, iter = 1), class = "nestimate_bad_input")
  expect_error(permutation(net, .cei_net("frequency"), iter = 10),
               class = "nestimate_bad_input")
})

test_that("bayes_compare() classes its input errors", {
  net <- .cei_net()
  expect_error(bayes_compare(net, net, prior = -1),
               class = "nestimate_bad_input")
})

test_that("compare_model() classes its input errors", {
  expect_error(compare_model(matrix(1:6, 2), matrix(1:6, 2)),
               class = "nestimate_bad_input")
  expect_error(compare_model(diag(3), diag(4)), class = "nestimate_bad_input")
})

test_that("boot_glasso() classes its input errors", {
  expect_error(boot_glasso(data.frame(a = 1:5, b = 5:1), iter = 1),
               class = "nestimate_bad_input")
  expect_error(boot_glasso(letters), class = "nestimate_bad_input")
})

test_that("net_prune() and its companions class their input errors", {
  net <- .cei_net()
  expect_error(net_prune(net, threshold = "a"), class = "nestimate_bad_input")
  expect_error(net_prune(1), class = "nestimate_bad_input")
  expect_error(net_deprune(net), class = "nestimate_bad_input")
  expect_error(net_reprune(net), class = "nestimate_bad_input")
  expect_error(net_pruning_details(net), class = "nestimate_bad_input")
})

test_that("vertex_bootstrap() and vertex_compare() class their input errors", {
  expect_error(vertex_bootstrap("a"), class = "nestimate_bad_input")
  expect_error(vertex_bootstrap(diag(2)), class = "nestimate_bad_input")
  net <- .cei_net()
  expect_error(vertex_compare(net, net, labels = "x"),
               class = "nestimate_bad_input")
})

test_that("compare_networks() and its verbs class their input errors", {
  net1 <- .cei_net()
  net2 <- .cei_net()
  expect_error(compare_networks(net1, net2, iter = 1),
               class = "nestimate_bad_input")
  expect_error(edge_differences(1), class = "nestimate_bad_input")
  cmp <- compare_networks(net1, net2)
  expect_error(plot(cmp, combined = NA), class = "nestimate_bad_input")
})

test_that("permutation_diagnostics() classes its input errors", {
  net <- .cei_net()
  expect_error(permutation_diagnostics(net, net, actor = 1),
               class = "nestimate_bad_input")
})

test_that("subtract_networks() and as_netdifference() class their input errors", {
  expect_error(subtract_networks(diag(3), diag(4)),
               class = "nestimate_bad_input")
  expect_error(as_netdifference(1), class = "nestimate_bad_input")
})

test_that("magnitude_difference() classes its input errors", {
  expect_error(magnitude_difference(1), class = "nestimate_bad_input")
})

test_that("bootstrap_network() classes its input errors", {
  expect_error(bootstrap_network(1), class = "nestimate_bad_input")
  expect_error(bootstrap_network(.cei_net(), iter = 1),
               class = "nestimate_bad_input")
})
